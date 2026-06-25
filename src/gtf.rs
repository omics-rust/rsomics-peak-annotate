//! GTF parser → [`GeneModel`]. Fail-loud on any malformed data line; only
//! `#`-comment and blank lines are skipped. Both GTF (`gene_id "X";`) and GFF3
//! (`gene_id=X`) attribute styles are accepted.

use std::collections::BTreeMap;
use std::path::Path;

use rsomics_common::{Context, Result, RsomicsError};
use rsomics_intervals::Strand;

use crate::model::{GeneModel, Region, Transcript, tss_of};

/// A list of `(start, end)` intervals on one transcript.
type Spans = Vec<(u64, u64)>;

struct GtfRecord {
    chrom: String,
    feature: String,
    /// 0-based half-open, already converted from the GTF's 1-based inclusive.
    start: u64,
    end: u64,
    strand: Strand,
    gene_id: String,
    transcript_id: String,
    exon_number: Option<u32>,
}

/// Per-transcript accumulator while scanning records in file order.
struct TxBuilder {
    gene_id: String,
    chrom: String,
    strand: Strand,
    /// `(start, end, exon_number_attr)` in file order.
    exons: Vec<(u64, u64, Option<u32>)>,
    cds: Vec<(u64, u64)>,
    utr5: Vec<(u64, u64)>,
    utr3: Vec<(u64, u64)>,
}

pub fn read(path: &Path) -> Result<GeneModel> {
    let data = std::fs::read(path).rs_with_context(|| format!("reading GTF {}", path.display()))?;
    parse(&data)
}

fn parse(data: &[u8]) -> Result<GeneModel> {
    let mut txs: BTreeMap<String, TxBuilder> = BTreeMap::new();

    for (lineno, raw) in data.split(|&b| b == b'\n').enumerate() {
        let line = trim_ascii_end(raw);
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        let rec = parse_line(line)
            .map_err(|e| RsomicsError::InvalidInput(format!("GTF line {}: {e}", lineno + 1)))?;
        ingest(rec, &mut txs);
    }

    Ok(assemble(txs))
}

fn ingest(rec: GtfRecord, txs: &mut BTreeMap<String, TxBuilder>) {
    if rec.transcript_id.is_empty() {
        return;
    }
    let tx = txs
        .entry(rec.transcript_id.clone())
        .or_insert_with(|| TxBuilder {
            gene_id: rec.gene_id.clone(),
            chrom: rec.chrom.clone(),
            strand: rec.strand,
            exons: Vec::new(),
            cds: Vec::new(),
            utr5: Vec::new(),
            utr3: Vec::new(),
        });
    match rec.feature.as_str() {
        "exon" => tx.exons.push((rec.start, rec.end, rec.exon_number)),
        "CDS" => tx.cds.push((rec.start, rec.end)),
        "five_prime_utr" | "5UTR" | "five_prime_UTR" => tx.utr5.push((rec.start, rec.end)),
        "three_prime_utr" | "3UTR" | "three_prime_UTR" => tx.utr3.push((rec.start, rec.end)),
        _ => {}
    }
}

fn assemble(txs: BTreeMap<String, TxBuilder>) -> GeneModel {
    let mut model = GeneModel::default();

    for (transcript_id, tx) in txs {
        let (Some(start), Some(end)) = (
            tx.exons.iter().map(|&(s, _, _)| s).min(),
            tx.exons.iter().map(|&(_, e, _)| e).max(),
        ) else {
            continue; // a transcript with no exons cannot anchor or be annotated
        };
        model.transcripts.push(Transcript {
            chrom: tx.chrom.clone(),
            start,
            end,
            strand: tx.strand,
            gene_id: tx.gene_id.clone(),
            transcript_id: transcript_id.clone(),
            tss: tss_of(start, end, tx.strand),
        });
        for r in rank_exons(&transcript_id, &tx) {
            model.exons.push(r);
        }
        for r in rank_introns(&transcript_id, &tx) {
            model.introns.push(r);
        }
        let (utr5, utr3) = transcript_utrs(&tx);
        for (s, e) in utr5 {
            model.utr5.push(plain_region(&transcript_id, &tx, s, e));
        }
        for (s, e) in utr3 {
            model.utr3.push(plain_region(&transcript_id, &tx, s, e));
        }
    }
    model.transcripts.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.tss.cmp(&b.tss))
            .then(a.transcript_id.cmp(&b.transcript_id))
    });
    model
}

/// Introns are the gaps between consecutive (genomically sorted) exons, ranked
/// like exons (1-based, minus-strand counts from the 3' genomic end), with a
/// total of `exons - 1`.
fn rank_introns(transcript_id: &str, tx: &TxBuilder) -> Vec<Region> {
    let mut sorted: Vec<(u64, u64)> = tx.exons.iter().map(|&(s, e, _)| (s, e)).collect();
    sorted.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    if sorted.len() < 2 {
        return Vec::new();
    }
    let total = u32::try_from(sorted.len() - 1).unwrap_or(u32::MAX);
    let mut out = Vec::new();
    for i in 0..sorted.len() - 1 {
        let (s, e) = (sorted[i].1, sorted[i + 1].0);
        if s >= e {
            continue;
        }
        let genomic_rank = u32::try_from(i).unwrap_or(u32::MAX) + 1;
        let rank = if tx.strand == Strand::Forward {
            genomic_rank
        } else {
            total - genomic_rank + 1
        };
        out.push(Region {
            chrom: tx.chrom.clone(),
            start: s,
            end: e,
            gene_id: tx.gene_id.clone(),
            transcript_id: transcript_id.to_string(),
            exon_rank: Some(rank),
            exon_total: Some(total),
        });
    }
    out
}

/// 5'/3' UTR intervals. With CDS present, UTRs are the exonic regions outside the
/// coding span (strand-aware), matching TxDb's derivation; without CDS, explicit
/// UTR features are used.
fn transcript_utrs(tx: &TxBuilder) -> (Spans, Spans) {
    if tx.cds.is_empty() {
        return (tx.utr5.clone(), tx.utr3.clone());
    }
    let cds_min = tx.cds.iter().map(|&(s, _)| s).min().unwrap();
    let cds_max = tx.cds.iter().map(|&(_, e)| e).max().unwrap();
    let mut low = Vec::new();
    let mut high = Vec::new();
    for &(s, e, _) in &tx.exons {
        if s < cds_min {
            low.push((s, e.min(cds_min)));
        }
        if e > cds_max {
            high.push((s.max(cds_max), e));
        }
    }
    if tx.strand == Strand::Forward {
        (low, high)
    } else {
        (high, low)
    }
}

/// Exon rank: order exons by genomic position, number 1..N, then reverse the
/// numbering for minus-strand transcripts so rank 1 is the biological 5' exon.
fn rank_exons(transcript_id: &str, tx: &TxBuilder) -> Vec<Region> {
    let mut sorted: Vec<(u64, u64)> = tx.exons.iter().map(|&(s, e, _)| (s, e)).collect();
    sorted.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let total = u32::try_from(sorted.len()).unwrap_or(u32::MAX);

    sorted
        .into_iter()
        .enumerate()
        .map(|(i, (s, e))| {
            let genomic_rank = u32::try_from(i).unwrap_or(u32::MAX) + 1;
            let rank = if tx.strand == Strand::Forward {
                genomic_rank
            } else {
                total - genomic_rank + 1
            };
            Region {
                chrom: tx.chrom.clone(),
                start: s,
                end: e,
                gene_id: tx.gene_id.clone(),
                transcript_id: transcript_id.to_string(),
                exon_rank: Some(rank),
                exon_total: Some(total),
            }
        })
        .collect()
}

fn plain_region(transcript_id: &str, tx: &TxBuilder, start: u64, end: u64) -> Region {
    Region {
        chrom: tx.chrom.clone(),
        start,
        end,
        gene_id: tx.gene_id.clone(),
        transcript_id: transcript_id.to_string(),
        exon_rank: None,
        exon_total: None,
    }
}

fn parse_line(line: &[u8]) -> std::result::Result<GtfRecord, String> {
    let s = std::str::from_utf8(line).map_err(|e| format!("non-UTF8: {e}"))?;
    let mut f = s.split('\t');
    let chrom = f.next().ok_or("missing seqname")?;
    let _source = f.next().ok_or("missing source")?;
    let feature = f.next().ok_or("missing feature")?;
    let start_1 = f.next().ok_or("missing start")?;
    let end_1 = f.next().ok_or("missing end")?;
    let _score = f.next().ok_or("missing score")?;
    let strand_s = f.next().ok_or("missing strand")?;
    let _frame = f.next().ok_or("missing frame")?;
    let attrs = f.next().ok_or("missing attributes")?;

    let start_1: u64 = start_1
        .parse()
        .map_err(|_| format!("bad start {start_1:?}"))?;
    let end_1: u64 = end_1.parse().map_err(|_| format!("bad end {end_1:?}"))?;
    if start_1 == 0 || end_1 < start_1 {
        return Err(format!("inverted/zero coordinates {start_1}..{end_1}"));
    }
    let strand = Strand::from_byte(strand_s.bytes().next().unwrap_or(b'.'))
        .ok_or_else(|| format!("unstranded feature (strand {strand_s:?})"))?;

    let gene_id = attr(attrs, "gene_id").ok_or("missing gene_id attribute")?;
    let transcript_id = attr(attrs, "transcript_id").unwrap_or_default();
    let exon_number = attr(attrs, "exon_number").and_then(|v| v.parse().ok());

    Ok(GtfRecord {
        chrom: chrom.to_string(),
        feature: feature.to_string(),
        start: start_1 - 1,
        end: end_1,
        strand,
        gene_id,
        transcript_id,
        exon_number,
    })
}

/// Pull `key` from a GTF (`key "value";`) or GFF3 (`key=value;`) attribute field.
fn attr(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let part = part.trim();
        let Some(rest) = part.strip_prefix(key).filter(|r| r.starts_with([' ', '='])) else {
            continue;
        };
        let value = rest.trim_start_matches([' ', '=']).trim();
        return Some(value.trim_matches('"').to_string());
    }
    None
}

fn trim_ascii_end(raw: &[u8]) -> &[u8] {
    let mut len = raw.len();
    while len > 0 && raw[len - 1].is_ascii_whitespace() {
        len -= 1;
    }
    &raw[..len]
}

#[cfg(test)]
mod tests {
    use super::*;

    const GTF: &str = "\
# comment line
chr1\thavana\tgene\t1001\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\thavana\ttranscript\t1001\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\thavana\texon\t1001\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_number 1;
chr1\thavana\texon\t2001\t2300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_number 2;
chr1\thavana\tfive_prime_utr\t1001\t1050\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr2\thavana\tgene\t8000\t9000\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\";
chr2\thavana\texon\t8000\t8400\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\";
chr2\thavana\texon\t8600\t9000\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\";
";

    #[test]
    fn parses_transcripts_with_zero_based_coords() {
        let m = parse(GTF.as_bytes()).unwrap();
        assert_eq!(m.transcripts.len(), 2);
        let t1 = m
            .transcripts
            .iter()
            .find(|t| t.transcript_id == "T1")
            .unwrap();
        // span is the exon extent: 1-based 1001..2300 → 0-based half-open 1000..2300
        assert_eq!((t1.start, t1.end), (1000, 2300));
        assert_eq!(t1.strand, Strand::Forward);
        assert_eq!(t1.tss, 1000, "plus-strand TSS is genomic start");

        let t2 = m
            .transcripts
            .iter()
            .find(|t| t.transcript_id == "T2")
            .unwrap();
        assert_eq!((t2.start, t2.end), (7999, 9000));
        assert_eq!(t2.tss, 8999, "minus-strand TSS is genomic end-1");
    }

    #[test]
    fn plus_strand_exon_rank_is_genomic_order() {
        let m = parse(GTF.as_bytes()).unwrap();
        let mut ex: Vec<_> = m.exons.iter().filter(|r| r.gene_id == "G1").collect();
        ex.sort_by_key(|r| r.start);
        assert_eq!(ex[0].exon_rank, Some(1));
        assert_eq!(ex[1].exon_rank, Some(2));
        assert_eq!(ex[0].exon_total, Some(2));
    }

    #[test]
    fn minus_strand_exon_rank_counts_from_genomic_3prime() {
        let m = parse(GTF.as_bytes()).unwrap();
        let mut ex: Vec<_> = m.exons.iter().filter(|r| r.gene_id == "G2").collect();
        ex.sort_by_key(|r| r.start);
        // genomic-first exon is the biological last → rank 2; last → rank 1
        assert_eq!(ex[0].exon_rank, Some(2));
        assert_eq!(ex[1].exon_rank, Some(1));
    }

    #[test]
    fn utr_regions_captured() {
        let m = parse(GTF.as_bytes()).unwrap();
        assert_eq!(m.utr5.len(), 1);
        assert_eq!((m.utr5[0].start, m.utr5[0].end), (1000, 1050));
    }

    #[test]
    fn gff3_attribute_style() {
        let g = "chr1\tx\texon\t10\t20\t.\t+\t.\tgene_id=GX;transcript_id=TX";
        let m = parse(g.as_bytes()).unwrap();
        assert_eq!(m.transcripts[0].gene_id, "GX");
        assert_eq!(m.transcripts[0].transcript_id, "TX");
    }

    #[test]
    fn malformed_data_line_fails_loud() {
        let bad = "chr1\thavana\tgene\tNOTANUMBER\t20\t.\t+\t.\tgene_id \"G\";";
        assert!(parse(bad.as_bytes()).is_err());
    }

    #[test]
    fn missing_gene_id_fails_loud() {
        let bad = "chr1\thavana\tgene\t10\t20\t.\t+\t.\tfoo \"bar\";";
        assert!(parse(bad.as_bytes()).is_err());
    }
}
