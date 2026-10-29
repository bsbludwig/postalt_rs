// SPDX-License-Identifier: GPL-3.0-only
// Part of the BWA/bwakit postalt Rust port; see NOTICE.md for provenance.

//! Batch-parallel streaming: read groups, process them, write in input order.
//!
//! Deviation from the single-threaded script: input is cut into [`Batch`]es of
//! whole read-end groups. Groups are independent in r985, so each batch can be
//! processed by any worker and the results written in input order; output is
//! byte-identical to sequential processing. With `jobs == 0` everything runs on
//! the calling thread; otherwise one thread reads, `jobs` threads process and
//! the calling thread writes, with a bounded number of batches in flight.
//!
//! Errors keep sequential semantics: everything emitted before the failing
//! record is written, and a group is written only if sequential processing
//! would have flushed it (see [`Processed::tail`]).
use crate::io::{Batch, Format, Input, Output, Record, Serializer, is_header_line};
use crate::postalt::{AltIndex, Processor};
use anyhow::{Context, Result, anyhow};
use seqair::bam::{BamHeader, RawBamRecord};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufWriter, Write},
    panic::{AssertUnwindSafe, catch_unwind},
    sync::{
        Arc, Mutex,
        atomic::{AtomicBool, Ordering},
        mpsc,
    },
    thread,
};

/// Target input bytes per batch. Larger batches amortize hand-off costs;
/// smaller ones bound memory, which is a few times `BATCH_BYTES` per batch in
/// flight (input, output, and generated ALT records).
pub const BATCH_BYTES: usize = 512 << 10;

/// Byte buffers returned for reuse between pipeline stages.
///
/// Batches and their output are a few MiB each; reusing them avoids allocator
/// churn and keeps resident memory near what is actually in flight.
#[derive(Debug, Default)]
pub struct Buffers(Mutex<Vec<Vec<u8>>>);
impl Buffers {
    /// Spare buffers kept at most; more are freed.
    const LIMIT: usize = 64;
    /// Capacity reserved up front, so that a typical batch or its output never
    /// regrows (regrowth would copy and churn the allocator).
    const CAPACITY: usize = 2 * BATCH_BYTES + (64 << 10);
    /// An empty buffer, reusing a returned allocation if available.
    pub fn take(&self) -> Vec<u8> {
        let mut buffer = self
            .0
            .lock()
            .ok()
            .and_then(|mut v| v.pop())
            .unwrap_or_default();
        buffer.reserve(Self::CAPACITY);
        buffer
    }
    /// Return a buffer for reuse.
    pub fn give(&self, mut buffer: Vec<u8>) {
        buffer.clear();
        // One huge group must not keep its allocation resident for the rest of
        // the run.
        if buffer.capacity() > 4 * Self::CAPACITY {
            buffer.shrink_to(Self::CAPACITY);
        }
        if let Ok(mut spare) = self.0.lock()
            && spare.len() < Self::LIMIT
        {
            spare.push(buffer);
        }
    }
}

/// Processing options shared by every worker.
#[derive(Clone, Debug)]
pub struct Settings {
    /// Minimum `pa` ratio, as for [`Processor::new`].
    pub ratio: f64,
    /// Installed-k8 reversal behavior, as for [`Processor::set_k8_compat`].
    pub k8_compat: bool,
    /// Whether HLA FASTQ text is produced.
    pub fastq: bool,
}

/// Output produced for some consecutive groups.
#[derive(Debug, Default)]
pub struct Chunk {
    /// Serialized alignments in the output format.
    pub records: Vec<u8>,
    /// FASTQ text per HLA gene index.
    pub fastq: Vec<Vec<u8>>,
}
impl Chunk {
    fn new(genes: usize, records: Vec<u8>) -> Self {
        Self {
            records,
            fastq: vec![Vec::new(); genes],
        }
    }
}

/// A processed batch.
///
/// Sequentially, a group is flushed when the first record of the next group
/// parses, or at a clean end of input. The final group of a batch is therefore
/// kept apart in `tail`, to be written only if the next batch's first record
/// parsed (`first_ok`) or input ended cleanly.
#[derive(Debug)]
pub struct Processed {
    /// The batch's input position.
    pub index: u64,
    /// Whether the batch's first record parsed successfully.
    pub first_ok: bool,
    /// Groups flushed within the batch.
    pub body: Chunk,
    /// The batch's final group, and an error raised while flushing it.
    pub tail: (Chunk, Option<anyhow::Error>),
    /// The error that stopped processing within the batch, if any.
    pub error: Option<anyhow::Error>,
}

/// Parses, processes and serializes batches; one per processing thread.
pub struct Worker {
    processor: Processor,
    serializer: Serializer,
    header: Option<Arc<BamHeader>>,
    genes: usize,
}
impl Worker {
    /// Create a worker for this BAM input dictionary (`None` for SAM input) and
    /// BAM output dictionary (`None` for SAM output).
    pub fn new(
        index: Arc<AltIndex>,
        settings: &Settings,
        input: Option<Arc<BamHeader>>,
        output: Option<Arc<BamHeader>>,
    ) -> Self {
        let genes = if settings.fastq {
            index.hla_genes().len()
        } else {
            0
        };
        let mut processor = Processor::new(index, settings.ratio);
        processor.set_k8_compat(settings.k8_compat);
        processor.set_hla_fastq(settings.fastq);
        processor.set_bam_header(input.clone());
        Self {
            processor,
            serializer: Serializer::new(input.clone(), output),
            header: input,
            genes,
        }
    }

    /// Process every group in `batch`, taking output buffers from `buffers`.
    pub fn process(&mut self, batch: &Batch, buffers: &Buffers) -> Processed {
        let mut body = Chunk::new(self.genes, buffers.take());
        // A tail holds a single group, so it does not use a pooled batch buffer.
        let mut tail = Chunk::new(self.genes, Vec::new());
        let mut first_ok = false;
        let mut error = None;
        for (i, bytes) in batch.records().enumerate() {
            let number = batch.first + i as u64;
            // r985 prints a header line wherever it appears, at once, while the
            // current group is still buffered.
            if is_header_line(batch.format, bytes) {
                first_ok |= i == 0;
                if let Err(e) = self.serializer.write_header_line(bytes, &mut body.records) {
                    error = Some(e.context(format!("SAM line {number}")));
                    break;
                }
                continue;
            }
            let record = match self.parse(batch.format, bytes) {
                Ok(record) => record,
                Err(e) => {
                    error = Some(e.context(match batch.format {
                        Format::Sam => format!("SAM line {number}"),
                        Format::Bam => format!("BAM record {number}"),
                    }));
                    break;
                }
            };
            first_ok |= i == 0;
            // A later ALT alignment can change earlier MAPQs within this read end;
            // emission must wait until the group boundary, not each input record.
            if self.processor.starts_new_group(&record) {
                let serializer = &mut self.serializer;
                let flushed = self
                    .processor
                    .flush(|r| serializer.write(r, &mut body.records), &mut body.fastq);
                if let Err(e) = flushed {
                    error = Some(e);
                    break;
                }
            }
            if let Err(e) = self.processor.push(record) {
                error = Some(e);
                break;
            }
        }
        let mut tail_error = None;
        if error.is_none() {
            let serializer = &mut self.serializer;
            tail_error = self
                .processor
                .flush(|r| serializer.write(r, &mut tail.records), &mut tail.fastq)
                .err();
        }
        // Leave no state behind for the next batch, even after an error.
        self.processor.clear();
        Processed {
            index: batch.index,
            first_ok,
            body,
            tail: (tail, tail_error),
            error,
        }
    }

    /// Parse one framed record into a pooled record.
    fn parse(&mut self, format: Format, bytes: &[u8]) -> Result<Record> {
        match format {
            Format::Sam => {
                let mut record = self.processor.spare_sam();
                record.parse_into(bytes)?;
                Ok(Record::Sam(record))
            }
            Format::Bam => {
                let header = self.header.as_deref().context("BAM batch without header")?;
                // Reuse a spare record's allocation for this payload.
                let mut buffer = self.processor.spare_bam().into_payload();
                buffer.clear();
                buffer.extend_from_slice(bytes);
                Ok(Record::Bam(RawBamRecord::from_payload(buffer, header)?))
            }
        }
    }
}

/// Writes processed batches in order, deciding when each batch's tail is due.
struct Sink<'a> {
    output: &'a mut Output,
    fastq: &'a mut [BufWriter<File>],
    tail: Option<(Chunk, Option<anyhow::Error>)>,
    buffers: &'a Buffers,
}
impl Sink<'_> {
    /// Write a chunk, returning a batch-sized alignment buffer for reuse.
    fn write(&mut self, chunk: Chunk) -> Result<()> {
        self.output.write_serialized(&chunk.records)?;
        for (file, text) in self.fastq.iter_mut().zip(&chunk.fastq) {
            file.write_all(text)?;
        }
        self.recycle(chunk);
        Ok(())
    }
    fn recycle(&self, chunk: Chunk) {
        if chunk.records.capacity() >= Buffers::CAPACITY {
            self.buffers.give(chunk.records);
        }
    }
    /// Write the pending tail, if any, and report its flush error.
    fn write_tail(&mut self) -> Result<()> {
        if let Some((tail, error)) = self.tail.take() {
            self.write(tail)?;
            if let Some(e) = error {
                return Err(e);
            }
        }
        Ok(())
    }
    /// Accept the next batch in input order.
    fn accept(&mut self, processed: Processed) -> Result<()> {
        if processed.first_ok {
            self.write_tail()?;
        }
        if let Some((tail, _)) = self.tail.take() {
            self.recycle(tail);
        }
        self.write(processed.body)?;
        if let Some(e) = processed.error {
            return Err(e);
        }
        self.tail = Some(processed.tail);
        Ok(())
    }
}

/// Messages from the reader and workers to the writing thread.
enum Message {
    Done(Box<Processed>),
    /// Input ended cleanly after this many batches.
    End(u64),
    /// Reading failed after this many batches.
    Failed(u64, anyhow::Error),
}

/// Stream `input` to `output`, processing groups with `jobs` worker threads.
///
/// `fastq[i]` receives FASTQ text for HLA gene `i`. `make_worker` is called
/// once per worker (once in total when `jobs == 0`).
///
/// # Errors
/// Returns the first error in input order: read, parse, processing,
/// serialization or write failures. Output may already contain earlier groups.
pub fn run(
    mut input: Input,
    output: &mut Output,
    fastq: &mut [BufWriter<File>],
    jobs: usize,
    mut make_worker: impl FnMut() -> Worker,
) -> Result<()> {
    let buffers = Buffers::default();
    let mut sink = Sink {
        output,
        fastq,
        tail: None,
        buffers: &buffers,
    };
    if jobs == 0 {
        let mut worker = make_worker();
        while let Some(batch) = input.next_batch_into(BATCH_BYTES, buffers.take())? {
            let processed = worker.process(&batch, &buffers);
            buffers.give(batch.data);
            sink.accept(processed)?;
        }
        // EOF does not produce a boundary record, so explicitly drain the
        // last group.
        return sink.write_tail();
    }
    let in_flight = 2 * jobs + 2;
    // Shared by the scoped threads, so declared outside the scope.
    let (batch_tx, batch_rx) = mpsc::sync_channel::<Batch>(jobs);
    let batch_rx = Mutex::new(batch_rx);
    let stop = AtomicBool::new(false);
    thread::scope(|scope| {
        let (result_tx, result_rx) = mpsc::channel::<Message>();
        // Tokens bound the batches between reading and writing, and therefore
        // memory, regardless of how far workers run ahead of a slow batch.
        let (token_tx, token_rx) = mpsc::sync_channel::<()>(in_flight);
        for _ in 0..in_flight {
            token_tx.send(()).expect("receiver is alive");
        }
        for i in 0..jobs {
            let (batch_rx, stop, buffers) = (&batch_rx, &stop, &buffers);
            let result_tx = result_tx.clone();
            let mut worker = make_worker();
            thread::Builder::new()
                .name(format!("postalt-worker-{i}"))
                .spawn_scoped(scope, move || {
                    loop {
                        let batch = batch_rx.lock().map(|rx| rx.recv());
                        let Ok(Ok(batch)) = batch else { break };
                        if stop.load(Ordering::Relaxed) {
                            buffers.give(batch.data);
                            continue;
                        }
                        // A panic becomes this batch's error, so the writer
                        // never waits for a batch that will not arrive.
                        let processed =
                            catch_unwind(AssertUnwindSafe(|| worker.process(&batch, buffers)))
                                .unwrap_or_else(|_| Processed {
                                    index: batch.index,
                                    first_ok: false,
                                    body: Chunk::default(),
                                    tail: (Chunk::default(), None),
                                    error: Some(anyhow!("worker panicked")),
                                });
                        buffers.give(batch.data);
                        // After the writer stops, keep draining so the reader
                        // is never blocked on a full batch queue.
                        let _ = result_tx.send(Message::Done(Box::new(processed)));
                    }
                })?;
        }
        // Spawn workers first: if a spawn fails, dropping the unsent `batch_tx`
        // (owned by this closure) releases the workers already started.
        let reader_tx = result_tx.clone();
        let buffers = &buffers;
        thread::Builder::new()
            .name("postalt-reader".into())
            .spawn_scoped(scope, move || {
                let mut count = 0;
                while token_rx.recv().is_ok() {
                    let message = match input.next_batch_into(BATCH_BYTES, buffers.take()) {
                        Ok(Some(batch)) => {
                            count += 1;
                            if batch_tx.send(batch).is_err() {
                                return;
                            }
                            continue;
                        }
                        Ok(None) => Message::End(count),
                        Err(e) => Message::Failed(count, e),
                    };
                    let _ = reader_tx.send(message);
                    return;
                }
            })?;
        drop(result_tx);
        let mut pending = BTreeMap::new();
        let mut next = 0;
        let mut end: Option<(u64, Option<anyhow::Error>)> = None;
        let result = (|| loop {
            while let Some(processed) = pending.remove(&next) {
                sink.accept(processed)?;
                next += 1;
                let _ = token_tx.send(());
            }
            if end.as_ref().is_some_and(|(n, _)| *n == next) {
                return match end.take() {
                    Some((_, Some(e))) => Err(e),
                    _ => sink.write_tail(),
                };
            }
            match result_rx.recv() {
                Ok(Message::Done(processed)) => {
                    pending.insert(processed.index, *processed);
                }
                Ok(Message::End(n)) => end = Some((n, None)),
                Ok(Message::Failed(n, e)) => end = Some((n, Some(e))),
                Err(_) => return Err(anyhow!("processing threads stopped unexpectedly")),
            }
        })();
        // Stop the other threads before the scope joins them: workers drain
        // remaining batches without processing them, and the reader stops at
        // its next token.
        stop.store(true, Ordering::Relaxed);
        drop(result_rx);
        drop(token_tx);
        result
    })
}
