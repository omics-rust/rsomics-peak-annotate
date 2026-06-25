//! Peak annotation: nearest strand-aware TSS, signed distance, and a single
//! genomic-feature category per peak, following ChIPseeker `annotatePeak`
//! (level="transcript", default tssRegion = (-3000, 3000)).

use std::collections::BTreeMap;

use rsomics_intervals::{Interval, IntervalIndex, IntervalSet, Strand};

use crate::model::{GeneModel, Region, Transcript};

pub const DOWNSTREAM_MAX: i64 = 300;

/// A peak to annotate. `start`/`end` are 0-based half-open. `raw` is the input
/// line verbatim, passed through unmodified into the output.
#[derive(Debug, Clone)]
pub struct Peak {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub raw: String,
}

/// One peak's full annotation. `gene_*` are reported in ChIPseeker's 1-based
/// TxDb convention; `distance_to_tss` is signed.
#[derive(Debug, Clone)]
pub struct Annotation {
    pub annotation: String,
    pub gene_chr: String,
    pub gene_start: u64,
    pub gene_end: u64,
    pub gene_length: u64,
    pub gene_strand: char,
    pub gene_id: String,
    pub transcript_id: String,
    pub distance_to_tss: i64,
}

/// Strand-aware signed distance from a peak to a TSS, measured from the peak
/// edge nearest the TSS. Negative = upstream (5') of the TSS, positive =
/// downstream (3'), 0 when the peak spans the TSS. Inputs are 0-based half-open
/// (peak) and 0-based (tss); the comparison is done in 1-based coordinates to
/// match ChIPseeker.
#[must_use]
pub fn signed_distance(peak_start: u64, peak_end: u64, tss: u64, strand: Strand) -> i64 {
    let pstart = peak_start as i64 + 1;
    let pend = peak_end as i64;
    let tss = tss as i64 + 1;
    if pstart <= tss && tss <= pend {
        return 0;
    }
    if strand == Strand::Forward {
        if pend < tss { pend - tss } else { pstart - tss }
    } else if pstart > tss {
        tss - pstart
    } else {
        tss - pend
    }
}

/// Indexes the gene model once for repeated peak queries: per-chromosome
/// TSS-sorted gene order plus COITree overlap indices for feature regions.
pub struct Annotator<'m> {
    model: &'m GeneModel,
    tss_by_chrom: BTreeMap<&'m str, Vec<usize>>,
    exon_index: IntervalIndex,
    intron_index: IntervalIndex,
    utr5_index: IntervalIndex,
    utr3_index: IntervalIndex,
    promoter_up: i64,
    promoter_down: i64,
}

impl<'m> Annotator<'m> {
    #[must_use]
    pub fn new(model: &'m GeneModel, tss_region: (i64, i64)) -> Self {
        let mut tss_by_chrom: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
        for (i, t) in model.transcripts.iter().enumerate() {
            tss_by_chrom.entry(t.chrom.as_str()).or_default().push(i);
        }
        for ids in tss_by_chrom.values_mut() {
            ids.sort_by_key(|&i| model.transcripts[i].tss);
        }

        Self {
            model,
            tss_by_chrom,
            exon_index: build_index(&model.exons),
            intron_index: build_index(&model.introns),
            utr5_index: build_index(&model.utr5),
            utr3_index: build_index(&model.utr3),
            promoter_up: tss_region.0,
            promoter_down: tss_region.1,
        }
    }

    /// Annotate one peak. `None` only when the model has no gene on the peak's
    /// chromosome (no candidate to anchor against).
    #[must_use]
    pub fn annotate(&self, peak: &Peak) -> Option<Annotation> {
        let nearest = self.nearest_transcript(peak)?;
        let tx = &self.model.transcripts[nearest.tx_idx];
        let category = self.category(peak, tx, nearest.distance);

        Some(Annotation {
            annotation: category,
            gene_chr: tx.chrom.clone(),
            // 0-based half-open → 1-based inclusive, as ChIPseeker reports.
            gene_start: tx.start + 1,
            gene_end: tx.end,
            gene_length: tx.length_1based(),
            gene_strand: tx.strand.as_byte() as char,
            gene_id: tx.gene_id.clone(),
            transcript_id: tx.transcript_id.clone(),
            distance_to_tss: nearest.distance,
        })
    }

    /// Among the transcript whose TSS precedes the peak and the one whose TSS
    /// follows, pick the smaller |distance|; ties go upstream (first index).
    fn nearest_transcript(&self, peak: &Peak) -> Option<Nearest> {
        let order = self.tss_by_chrom.get(peak.chrom.as_str())?;
        let pos = peak.start;

        // Partition by genomic TSS position; the candidates are the transcripts
        // straddling the peak start.
        let split = order.partition_point(|&i| self.model.transcripts[i].tss <= pos);

        let mut candidates: Vec<usize> = Vec::with_capacity(2);
        if split > 0 {
            candidates.push(order[split - 1]);
        }
        if split < order.len() {
            candidates.push(order[split]);
        }

        candidates
            .into_iter()
            .map(|ti| {
                let t = &self.model.transcripts[ti];
                Nearest {
                    tx_idx: ti,
                    distance: signed_distance(peak.start, peak.end, t.tss, t.strand),
                }
            })
            .min_by_key(|n| n.distance.abs())
    }

    /// The single highest-priority category for this peak:
    /// Promoter > 5'UTR > 3'UTR > Exon > Intron > Downstream > Distal Intergenic.
    fn category(&self, peak: &Peak, tx: &Transcript, distance: i64) -> String {
        if distance >= self.promoter_up && distance <= self.promoter_down {
            return promoter_label(distance);
        }
        if first_overlap(&self.utr5_index, &self.model.utr5, peak).is_some() {
            return "5' UTR".to_string();
        }
        if first_overlap(&self.utr3_index, &self.model.utr3, peak).is_some() {
            return "3' UTR".to_string();
        }
        if let Some(r) = first_overlap(&self.exon_index, &self.model.exons, peak) {
            return feature_label("Exon", "exon", r);
        }
        if let Some(r) = first_overlap(&self.intron_index, &self.model.introns, peak) {
            return feature_label("Intron", "intron", r);
        }
        if let Some(label) = downstream_label(peak, tx) {
            return label;
        }
        "Distal Intergenic".to_string()
    }
}

struct Nearest {
    tx_idx: usize,
    distance: i64,
}

fn build_index(regions: &[Region]) -> IntervalIndex {
    let mut set = IntervalSet::new();
    for r in regions {
        if r.start < r.end {
            let mut iv = Interval::new(r.chrom.clone(), r.start, r.end)
                .expect("region start < end checked above");
            iv.strand = None;
            set.push(iv);
        }
    }
    IntervalIndex::build(&set)
}

/// First region (lowest start) overlapping the peak. The index lookup uses the
/// peak coordinates; the returned region is resolved from `regions` by position.
fn first_overlap<'a>(
    index: &IntervalIndex,
    regions: &'a [Region],
    peak: &Peak,
) -> Option<&'a Region> {
    if peak.start >= peak.end {
        return None;
    }
    let hits = index.query(&peak.chrom, peak.start, peak.end);
    let lead = hits
        .into_iter()
        .min_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)))?;
    regions
        .iter()
        .find(|r| r.chrom == lead.chrom && r.start == lead.start && r.end == lead.end)
}

/// `Downstream (<=300bp)`. ChIPseeker's downstream is genomic, not strand-
/// corrected: the peak must sit within `DOWNSTREAM_MAX` bp to the genomic right
/// of the transcript's higher coordinate. A minus-strand gene's 3' (genomic-
/// left) side is therefore never Downstream.
fn downstream_label(peak: &Peak, tx: &Transcript) -> Option<String> {
    let gap = (peak.start as i64 + 1) - tx.end as i64;
    (gap > 0 && gap <= DOWNSTREAM_MAX).then(|| "Downstream (<=300bp)".to_string())
}

/// Promoter sub-bin by |distanceToTSS|: `<=1kb`, `1-2kb`, `2-3kb`, … (1kb steps).
#[must_use]
pub fn promoter_label(distance: i64) -> String {
    let kb = (distance.abs() as f64 / 1000.0).ceil().max(1.0) as i64;
    if kb <= 1 {
        "Promoter (<=1kb)".to_string()
    } else {
        format!("Promoter ({}-{}kb)", kb - 1, kb)
    }
}

/// `Exon (txId/geneId, exon N of M)` / `Intron (txId/geneId, intron N of M)` —
/// ChIPseeker's detail string anchored on the overlapping exon or intron.
fn feature_label(kind: &str, noun: &str, r: &Region) -> String {
    match (r.exon_rank, r.exon_total) {
        (Some(rank), Some(total)) => {
            format!(
                "{kind} ({}/{}, {noun} {rank} of {total})",
                r.transcript_id, r.gene_id
            )
        }
        _ => format!("{kind} ({}/{})", r.transcript_id, r.gene_id),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plus_strand_distance_uses_nearest_edge() {
        // + gene TSS at 0-based 999 (1-based 1000). Upstream peak [800,900):
        // nearest edge = end 900 → 900 - 1000 = -100.
        assert_eq!(signed_distance(800, 900, 999, Strand::Forward), -100);
        // Downstream peak [1100,1200): nearest edge = start 1101 → 1101-1000 = 101.
        assert_eq!(signed_distance(1100, 1200, 999, Strand::Forward), 101);
        // Peak spanning the TSS → 0.
        assert_eq!(signed_distance(950, 1050, 999, Strand::Forward), 0);
    }

    #[test]
    fn minus_strand_distance_is_flipped() {
        // - gene TSS at 0-based 1999 (1-based 2000). Genomic-right peak [2100,2200)
        // is upstream in transcription → 2000 - 2101 = -101.
        assert_eq!(signed_distance(2100, 2200, 1999, Strand::Reverse), -101);
        // Genomic-left peak [1800,1900) is downstream → 2000 - 1900 = 100.
        assert_eq!(signed_distance(1800, 1900, 1999, Strand::Reverse), 100);
    }

    #[test]
    fn promoter_subbins() {
        assert_eq!(promoter_label(0), "Promoter (<=1kb)");
        assert_eq!(promoter_label(-500), "Promoter (<=1kb)");
        assert_eq!(promoter_label(1000), "Promoter (<=1kb)");
        assert_eq!(promoter_label(-1500), "Promoter (1-2kb)");
        assert_eq!(promoter_label(2500), "Promoter (2-3kb)");
    }

    #[test]
    fn downstream_is_genomic_right_within_300_strand_independent() {
        let tx = |strand| Transcript {
            chrom: "c".to_string(),
            start: 1000,
            end: 2000, // 0-based half-open → 1-based last base 2000
            strand,
            gene_id: "g".to_string(),
            transcript_id: "t".to_string(),
            tss: if strand == Strand::Forward {
                1000
            } else {
                1999
            },
        };
        let peak = |s, e| Peak {
            chrom: "c".to_string(),
            start: s,
            end: e,
            raw: String::new(),
        };
        // 51 bp to the genomic right of the end → Downstream for either strand.
        assert!(downstream_label(&peak(2050, 2100), &tx(Strand::Forward)).is_some());
        assert!(downstream_label(&peak(2050, 2100), &tx(Strand::Reverse)).is_some());
        // boundary: gap 300 in, 301 out.
        assert!(downstream_label(&peak(2299, 2400), &tx(Strand::Forward)).is_some());
        assert!(downstream_label(&peak(2300, 2400), &tx(Strand::Forward)).is_none());
        // genomic-left of the gene (a minus gene's 3' side) is never Downstream.
        assert!(downstream_label(&peak(700, 900), &tx(Strand::Reverse)).is_none());
    }
}
