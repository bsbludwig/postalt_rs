//! Text-preserving sequential SAM/BAM adapter for postalt_rs.
//!
//! Keeps original SAM fields (including duplicate tags and decimal spelling).
//! Uses seqair's CIGAR/aux parsers, BAM wire validation and sequence codecs.
//!
//! This module is a local extension, not an upstream seqair streaming API. The
//! indexed/pileup APIs cannot preserve all required SAM text, unplaced records
//! and IUPAC bases. BAM conversion is necessarily semantic rather than lexical:
//! integer widths, float spelling, sequence case and RNEXT aliases can normalize.
//! Parsing is deliberately layered; `SamRecord::parse` does not validate every
//! optional field or eagerly decode the sequence/CIGAR of pass-through records.
use crate::bam::aux::{self, AuxValue};
use crate::bam::{self, AuxData, BamHeader, CigarOp, MateFields, OwnedBamRecord};
use seqair_types::{BamFlags, Base, BaseQuality, MappingQuality, Pos0};
use std::{
    io::{self, Write},
    num::NonZeroI32,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

/// Owned SAM fields with separately mutable numeric values used by postalt.
///
/// # Invariants
///
/// `fields` must contain at least the eleven mandatory SAM columns. `flag`,
/// `pos` and `mapq` are authoritative: serialization uses them instead of the
/// potentially stale strings at indices 1, 3 and 4. All other fields, including
/// duplicate tags, retain their order and spelling until BAM conversion.
/// Construct through the parsers or maintain these invariants when editing.
#[derive(Clone, Debug)]
pub struct SamRecord {
    /// Eleven core SAM columns followed by optional fields, in original order.
    pub fields: Vec<String>,
    /// SAM FLAG bits, overriding `fields[1]` on output.
    pub flag: u16,
    /// One-based SAM POS, or zero for unavailable; overrides `fields[3]`.
    pub pos: i64,
    /// SAM MAPQ including the 255 unavailable sentinel; overrides `fields[4]`.
    pub mapq: u8,
}
impl SamRecord {
    /// Split one alignment line, stripping CR/LF but retaining optional-field text.
    ///
    /// # Errors
    /// Rejects fewer than eleven fields, invalid FLAG/MAPQ integers, and POS
    /// outside 0..=i32::MAX. CIGAR, sequence, mate fields and aux values are decoded
    /// only when needed later. This is stricter than JS numeric coercion but is
    /// not a comprehensive SAM validator; header lines belong to the caller.
    pub fn parse(line: &str) -> Result<Self> {
        let fields: Vec<String> =
            line.trim_end_matches(['\r', '\n']).split('\t').map(str::to_owned).collect();
        if fields.len() < 11 {
            return Err("expected at least 11 SAM fields".into());
        }
        let flag = fields[1].parse()?;
        let pos = fields[3].parse()?;
        let mapq = fields[4].parse()?;
        if !(0..=i32::MAX as i64).contains(&pos) {
            return Err("SAM POS out of range".into());
        }
        Ok(Self { fields, flag, pos, mapq })
    }
    /// Write TAB-separated fields and a newline, using the numeric overrides.
    ///
    /// Duplicate optional tags and their decimal spelling are retained. Returns
    /// an I/O error on failure; some fields may already have been written.
    pub fn write_sam(&self, out: &mut impl Write) -> io::Result<()> {
        for (i, field) in self.fields.iter().enumerate() {
            if i > 0 {
                out.write_all(b"\t")?;
            }
            match i {
                1 => write!(out, "{}", self.flag)?,
                3 => write!(out, "{}", self.pos)?,
                4 => write!(out, "{}", self.mapq)?,
                _ => out.write_all(field.as_bytes())?,
            }
        }
        out.write_all(b"\n")
    }
    /// Return the first matching optional-field suffix, e.g. for `"XA:Z:"`.
    /// This mirrors r985's first regex match for NM/XA/RG, not its last-match scan
    /// of buffered om/pa tags. Duplicate tags are neither collapsed nor removed.
    pub fn tag(&self, prefix: &str) -> Option<&str> {
        self.fields[11..].iter().find_map(|s| s.strip_prefix(prefix))
    }
    /// Decode field 6 through seqair's CIGAR parser; `*` yields no operations.
    /// Returns recognized parse errors without applying postalt scoring rules.
    pub fn cigar(&self) -> Result<Vec<CigarOp>> {
        let mut out = Vec::new();
        super::reader::parse_cigar(self.fields[5].as_bytes(), &mut out)?;
        Ok(out)
    }
    /// Resolve SAM references and encode typed fields for seqair's owned writer.
    ///
    /// # Sequence limitation
    /// `OwnedBamRecord` uses the five-base pileup alphabet, so this intermediate
    /// value loses non-ACGT ambiguity codes. Use `BamWriter::write_sam_record` to
    /// restore the original IUPAC nibbles before emitting a BAM record. Do not
    /// directly serialize this intermediate when full sequence fidelity matters.
    ///
    /// Aux fields are encoded in order, preserving duplicates. Missing sequence,
    /// quality and positions use the owned representation's empty/None sentinels.
    ///
    /// # Errors
    /// Returns unknown-reference, mate-coordinate, CIGAR, quality, auxiliary-value
    /// or owned-record construction errors. Requires the `SamRecord` invariants.
    pub fn to_owned_bam(&self, header: &BamHeader) -> Result<OwnedBamRecord> {
        let resolve = |name: &str| -> Result<_> {
            if name == "*" {
                Ok(None)
            } else {
                Ok(Some(
                    header
                        .tid(name)
                        .ok_or_else(|| format!("reference {name} missing from SAM @SQ header"))?,
                ))
            }
        };
        let tid = resolve(&self.fields[2])?;
        let pos = if self.pos == 0 { None } else { Some(Pos0::try_from(self.pos - 1)?) };
        let mate_tid = if self.fields[6] == "=" { tid } else { resolve(&self.fields[6])? };
        let pn: i64 = self.fields[7].parse()?;
        let mate_pos = if pn == 0 { None } else { Some(Pos0::try_from(pn - 1)?) };
        let seq = if self.fields[9] == "*" {
            Vec::new()
        } else {
            Base::from_ascii_vec(self.fields[9].as_bytes().to_vec())
        };
        let qual = if self.fields[10] == "*" {
            Vec::new()
        } else {
            self.fields[10]
                .bytes()
                .map(|b| b.checked_sub(33).map(BaseQuality::from_byte).ok_or("invalid SAM quality"))
                .collect::<std::result::Result<Vec<_>, _>>()?
        };
        // Encode each field independently. Setting all tags on one AuxData map
        // would replace duplicates, unlike the reference's append-only behavior.
        let mut aux = Vec::new();
        for field in &self.fields[11..] {
            // BAM's double-precision extension is not handled by the indexed
            // SAM parser. Retain it through the existing AuxData encoder.
            if field.as_bytes().get(2..5) == Some(b":d:") {
                let mut tag = AuxData::new();
                tag.set_double(field.as_bytes()[..2].try_into()?, field[5..].parse()?);
                aux.extend_from_slice(tag.as_bytes());
            } else if field.len() == 6 && field.as_bytes().get(2..5) == Some(b":B:") {
                // Empty BAM arrays have a subtype and zero elements. The indexed
                // parser treats the absent element text as one invalid element.
                let mut tag = AuxData::new();
                let key = field.as_bytes()[..2].try_into()?;
                match field.as_bytes()[5] {
                    b'c' => tag.set_array_i8(key, &[])?,
                    b'C' => tag.set_array_u8(key, &[])?,
                    b's' => tag.set_array_i16(key, &[])?,
                    b'S' => tag.set_array_u16(key, &[])?,
                    b'i' => tag.set_array_i32(key, &[])?,
                    b'I' => tag.set_array_u32(key, &[])?,
                    b'f' => tag.set_array_f32(key, &[])?,
                    _ => return Err("invalid empty BAM array subtype".into()),
                }
                aux.extend_from_slice(tag.as_bytes());
            } else {
                super::reader::parse_aux_tags(field.as_bytes(), &mut aux)?;
            }
        }
        Ok(OwnedBamRecord::builder(tid, pos, self.fields[0].as_bytes().to_vec())
            .flags(BamFlags::from(self.flag))
            .mapq(MappingQuality::from_byte(self.mapq))
            .mate_fields(MateFields {
                tid: mate_tid,
                pos: mate_pos,
                template_len: NonZeroI32::new(self.fields[8].parse()?),
            })
            .cigar(self.cigar()?)
            .seq(seq)
            .qual(qual)
            .aux(AuxData::from_bytes(aux))
            .build()?)
    }
    /// Decode one BAM payload (without its block_size prefix) into SAM fields.
    ///
    /// Unplaced records are retained, unlike seqair's indexed read path. Temporary
    /// POS normalization is restored before propagating any parser error, leaving
    /// the caller's payload unchanged. Sequence decoding retains full IUPAC codes.
    ///
    /// # Errors
    /// Returns truncated-layout, reference/position, read-name or text-conversion
    /// errors. Auxiliary iteration uses upstream seqair's permissive iterator;
    /// this is not a complete malformed-aux validator.
    pub fn from_bam(raw: &mut [u8], header: &BamHeader) -> Result<Self> {
        // Indexed seqair's wire parser requires a placed position. Normalize only
        // for validation, then restore BAM's -1 sentinel for streaming unmapped reads.
        if raw.len() < 32 {
            return Err("truncated BAM record".into());
        }
        let pos = i32::from_le_bytes(raw[4..8].try_into()?);
        if pos == -1 {
            raw[4..8].copy_from_slice(&0i32.to_le_bytes());
        }
        let parsed = bam::record::parse_header(raw);
        raw[4..8].copy_from_slice(&pos.to_le_bytes());
        let h = parsed?;
        let refs = h.refs.resolve(header)?;
        if h.var_start <= 32 || raw[h.var_start - 1] != 0 {
            return Err("invalid BAM read name".into());
        }
        let name = |tid| header.target_name(tid).unwrap_or("*").to_owned();
        let mut cigar = Vec::new();
        CigarOp::extend_from_bam_bytes(&mut cigar, &raw[h.var_start..h.cigar_end]);
        // Decode raw nibbles directly: RecordStore's Base conversion would turn
        // ambiguity codes into N before postalt could emit or reverse them.
        let seq = bam::seq::decode_seq(&raw[h.cigar_end..h.seq_end], h.seq_len as usize);
        let qual = &raw[h.seq_end..h.qual_end];
        let mut fields = vec![
            String::from_utf8(raw[32..h.var_start - 1].to_vec())?,
            h.flags.to_string(),
            refs.tid.map(name).unwrap_or_else(|| "*".into()),
            (pos as i64 + 1).to_string(),
            h.mapq.as_byte().to_string(),
            if cigar.is_empty() { "*".into() } else { bam::CigarStr(&cigar).to_string() },
            if refs.mate_fields.tid.is_some() && refs.mate_fields.tid == refs.tid {
                "=".into()
            } else {
                refs.mate_fields.tid.map(name).unwrap_or_else(|| "*".into())
            },
            refs.mate_fields.pos.map(|p| p.as_u64() + 1).unwrap_or(0).to_string(),
            refs.mate_fields.template_len.map(|n| n.get()).unwrap_or(0).to_string(),
            if seq.is_empty() { "*".into() } else { String::from_utf8(seq)? },
            if qual.is_empty() || qual[0] == 255 {
                "*".into()
            } else {
                String::from_utf8(qual.iter().map(|q| q.saturating_add(33)).collect())?
            },
        ];
        for tag in aux::iter_tags(&raw[h.qual_end..]) {
            let (tag, value) = tag;
            let prefix = std::str::from_utf8(&tag)?;
            fields.push(format_aux(prefix, value)?);
        }
        Ok(Self { flag: u16::from(h.flags), pos: pos as i64 + 1, mapq: h.mapq.as_byte(), fields })
    }
}
/// Format a decoded BAM value using its SAM representation.
///
/// All scalar integer widths become `i`; arrays retain their subtype and order.
/// Rust's float display may differ lexically from k8/samtools while representing
/// the same stored value. Empty arrays emit only their subtype, e.g. `XX:B:i`.
fn format_aux(tag: &str, value: AuxValue<'_>) -> Result<String> {
    use AuxValue::*;
    if let Some(n) = value.as_i64() {
        return Ok(format!("{tag}:i:{n}"));
    }
    let mut s;
    macro_rules! array {
        ($bytes:expr, $typ:ty, $sub:literal) => {{
            s = format!("{tag}:B:{}", $sub);
            for chunk in $bytes.chunks_exact(size_of::<$typ>()) {
                use std::fmt::Write;
                write!(s, ",{}", <$typ>::from_le_bytes(chunk.try_into()?))?;
            }
            s
        }};
    }
    Ok(match value {
        Char(c) => format!("{tag}:A:{}", c as char),
        Float(f) => format!("{tag}:f:{f}"),
        Double(f) => format!("{tag}:d:{f}"),
        String(v) => format!("{tag}:Z:{}", std::str::from_utf8(v)?),
        Hex(v) => format!("{tag}:H:{}", std::str::from_utf8(v)?),
        ArrayI8(b) => array!(b, i8, 'c'),
        ArrayU8(b) => array!(b, u8, 'C'),
        ArrayI16(b) => array!(b, i16, 's'),
        ArrayU16(b) => array!(b, u16, 'S'),
        ArrayI32(b) => array!(b, i32, 'i'),
        ArrayU32(b) => array!(b, u32, 'I'),
        ArrayFloat(b) => array!(b, f32, 'f'),
        _ => unreachable!(),
    })
}

/// Decode SAM CIGAR operations using the indexed reader's parser.
///
/// Returns an empty vector for `*` and propagates seqair's parse errors.
/// Recognizes `=`, `X` and `P`; ignoring them for r985 scoring is the caller's job.
pub fn parse_cigar(text: &str) -> Result<Vec<CigarOp>> {
    let mut out = Vec::new();
    super::reader::parse_cigar(text.as_bytes(), &mut out)?;
    Ok(out)
}

/// Read a sequential BAM header through the existing seqair header decoder.
///
/// Added to expose header parsing without enabling upstream's fuzz feature or
/// opening an indexed reader. Returns BGZF/header decoding errors and leaves the
/// stream positioned at the first record. Sequential use does not seek backward.
pub fn read_bam_header<R: std::io::Read + std::io::Seek>(
    reader: &mut bam::bgzf::BgzfReader<R>,
) -> Result<BamHeader> {
    Ok(BamHeader::parse(reader)?)
}
