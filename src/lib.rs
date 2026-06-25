//! Annotate ChIP/ATAC peaks with nearest gene, TSS distance, and a
//! genomic-feature category — Rust port of Bioconductor ChIPseeker
//! `annotatePeak` (level = "transcript", default tssRegion = (-3000, 3000)).
//!
//! Input peaks (BED / narrowPeak) are 0-based half-open and pass through
//! unmodified; nine annotation columns are appended.

pub mod annotate;
pub mod gtf;
pub mod model;

use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Context, Result, RsomicsError};

use annotate::{Annotator, Peak};

pub const APPENDED_COLUMNS: &[&str] = &[
    "annotation",
    "geneChr",
    "geneStart",
    "geneEnd",
    "geneLength",
    "geneStrand",
    "geneId",
    "transcriptId",
    "distanceToTSS",
];

/// Read peaks, build the gene model, annotate, and write the tab-separated
/// table (input columns passed through, nine annotation columns appended).
pub fn annotate_peaks(
    peaks_path: &Path,
    gtf_path: &Path,
    tss_region: (i64, i64),
    out: &mut dyn Write,
) -> Result<usize> {
    let peaks = read_peaks(peaks_path)?;
    let model = gtf::read(gtf_path)?;
    let annotator = Annotator::new(&model, tss_region);

    let mut w = BufWriter::new(out);
    let mut written = 0usize;
    for peak in &peaks {
        let Some(a) = annotator.annotate(peak) else {
            return Err(RsomicsError::InvalidInput(format!(
                "no gene on chromosome {:?} to annotate peak at {}:{}-{}",
                peak.chrom, peak.chrom, peak.start, peak.end
            )));
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            peak.raw,
            a.annotation,
            a.gene_chr,
            a.gene_start,
            a.gene_end,
            a.gene_length,
            a.gene_strand,
            a.gene_id,
            a.transcript_id,
            a.distance_to_tss,
        )
        .rs_context("writing annotated peak")?;
        written += 1;
    }
    w.flush().rs_context("flushing output")?;
    Ok(written)
}

/// Parse a header line of column names, used only to echo a header when the
/// caller wants one. The appended column names are [`APPENDED_COLUMNS`].
#[must_use]
pub fn header_line(peak_columns: &[&str]) -> String {
    let mut cols: Vec<&str> = peak_columns.to_vec();
    cols.extend_from_slice(APPENDED_COLUMNS);
    cols.join("\t")
}

fn read_peaks(path: &Path) -> Result<Vec<Peak>> {
    let data =
        std::fs::read(path).rs_with_context(|| format!("reading peaks {}", path.display()))?;
    let text = std::str::from_utf8(&data)
        .map_err(|e| RsomicsError::InvalidInput(format!("non-UTF8 peak file: {e}")))?;

    let mut peaks = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        let trimmed = line.trim_end();
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("track") {
            continue;
        }
        let mut f = trimmed.split('\t');
        let chrom = f.next().ok_or_else(|| bad_peak(lineno, "missing chrom"))?;
        let start = f.next().ok_or_else(|| bad_peak(lineno, "missing start"))?;
        let end = f.next().ok_or_else(|| bad_peak(lineno, "missing end"))?;
        let start: u64 = start
            .parse()
            .map_err(|_| bad_peak(lineno, "non-numeric start"))?;
        let end: u64 = end
            .parse()
            .map_err(|_| bad_peak(lineno, "non-numeric end"))?;
        if end <= start {
            return Err(bad_peak(lineno, "end <= start"));
        }
        peaks.push(Peak {
            chrom: chrom.to_string(),
            start,
            end,
            raw: trimmed.to_string(),
        });
    }
    Ok(peaks)
}

fn bad_peak(lineno: usize, msg: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!("peak line {}: {msg}", lineno + 1))
}

/// Parse a `up,down` tssRegion string into a signed `(up, down)` pair.
#[allow(clippy::missing_errors_doc)]
pub fn parse_tss_region(s: &str) -> Result<(i64, i64)> {
    let (a, b) = s
        .split_once(',')
        .ok_or_else(|| RsomicsError::InvalidInput(format!("tss-region {s:?} must be 'up,down'")))?;
    let up: i64 = a
        .trim()
        .parse()
        .map_err(|_| RsomicsError::InvalidInput(format!("tss-region upstream {a:?} not an int")))?;
    let down: i64 = b.trim().parse().map_err(|_| {
        RsomicsError::InvalidInput(format!("tss-region downstream {b:?} not an int"))
    })?;
    if up > down {
        return Err(RsomicsError::InvalidInput(format!(
            "tss-region upstream {up} must be <= downstream {down}"
        )));
    }
    Ok((up, down))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tss_region_parses_signed_pair() {
        assert_eq!(parse_tss_region("-3000,3000").unwrap(), (-3000, 3000));
        assert_eq!(parse_tss_region("-1000, 500").unwrap(), (-1000, 500));
    }

    #[test]
    fn tss_region_rejects_inverted() {
        assert!(parse_tss_region("3000,-3000").is_err());
        assert!(parse_tss_region("garbage").is_err());
    }

    #[test]
    fn header_appends_nine_columns() {
        let h = header_line(&["chrom", "start", "end", "name"]);
        let cols: Vec<&str> = h.split('\t').collect();
        assert_eq!(cols.len(), 4 + 9);
        assert_eq!(cols[4], "annotation");
        assert_eq!(*cols.last().unwrap(), "distanceToTSS");
    }
}
