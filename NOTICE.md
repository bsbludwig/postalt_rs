# Source attribution and licensing

`postalt_rs` is a Rust port of `bwa-postalt.js` from the **BWA/bwakit** project.
The original postalt implementation is by **Heng Li**. The upstream change to
the SAM input loop's empty-line check is by **Imran Haque**.

## Pinned upstream source

The relevant BWA revision is
`27dd1da7020151be9358fce67cbe2c10801c348d` (2018-04-07):

- [bwakit/bwa-postalt.js](https://github.com/lh3/bwa/blob/27dd1da7020151be9358fce67cbe2c10801c348d/bwakit/bwa-postalt.js)
- [README license declaration](https://github.com/lh3/bwa/blob/27dd1da7020151be9358fce67cbe2c10801c348d/README.md#availability)
- [COPYING: GNU GPL version 3](https://github.com/lh3/bwa/blob/27dd1da7020151be9358fce67cbe2c10801c348d/COPYING)
- [Empty-line handling change and authorship](https://github.com/lh3/bwa/commit/27dd1da7020151be9358fce67cbe2c10801c348d)
- [Earlier postalt change by Heng Li](https://github.com/lh3/bwa/commit/6c668d3e1a37ba65385b5dba76f94986933f8d86)

The reference used for compatibility testing is installed at
`/opt/apps/manual/bwa-mem2/2.3.0/noarch/bin/bwa-postalt.js`, with SHA-256
`50660424f540ad38e08c181f628289754a39b38e2092f3010cf7478e1d0d83b7`.
Comparison confirms that it matches the pinned upstream script after excluding
the installed `#!/bin/env k8` shebang and trailing blank lines. The installation
directory identifies how the script was supplied, not its source license.

## License and modifications

The port uses **GPL-3.0-only** metadata. BWA's README declares GPLv3, its COPYING
contains GPL version 3, and this script has no separate permissive or
later-version license notice. [LICENSE](LICENSE) contains the upstream GPLv3
text. The previous Intel/BWA-MEM2 MIT notice was incorrectly applied and has
been removed.

This is a modified implementation, not an upstream BWA release. Changes made on
2026-10-27 include translation to Rust, precomputed projection/index structures,
seqair-based streaming I/O, direct BAM support, corrected reversal by default,
and an explicit mode reproducing the supplied k8 runtime's reversal behavior.
Preserved quirks and intentional differences are described in [README.md](README.md)
and the source comments. Test-only reference repairs are documented in
`scripts/compatibility.py`; they do not modify the installed script.

## Vendored dependency

The upstream seqair code is by Pascal Hertleif and retains its original
[MIT](vendor/seqair/LICENSE-MIT) and
[Apache-2.0](vendor/seqair/LICENSE-APACHE) license notices. Its provenance and
local streaming additions are documented in
[POSTALT_PATCHES.md](vendor/seqair/POSTALT_PATCHES.md). The project's GPL license
does not replace the upstream dependency's notices.
