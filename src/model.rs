//! Gene model assembled from a GTF: genes carry a strand-aware TSS, transcripts
//! carry rank-ordered exons and UTR intervals. All coordinates are 0-based
//! half-open internally; the GTF's 1-based start is converted on ingest.

use rsomics_intervals::Strand;

/// A genomic feature region used for category overlap tests. `[start, end)`.
#[derive(Debug, Clone)]
pub struct Region {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub gene_id: String,
    pub transcript_id: String,
    /// 1-based rank along the transcript (strand-aware: minus-strand exons
    /// count from the 3' genomic end). `None` for non-exon regions.
    pub exon_rank: Option<u32>,
    /// Total exon count of the owning transcript. `None` for non-exon regions.
    pub exon_total: Option<u32>,
}

/// One transcript: its span, strand, and the strand-aware transcription start
/// site. ChIPseeker `level="transcript"` anchors each peak on the nearest
/// transcript's TSS and reports that transcript's span as `geneStart/geneEnd`.
#[derive(Debug, Clone)]
pub struct Transcript {
    pub chrom: String,
    /// 0-based half-open transcript span (exon extent).
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub gene_id: String,
    pub transcript_id: String,
    /// Strand-aware 5' end: `+` → genomic start, `-` → genomic end-1 (the last
    /// base of the span). Mirrors GRanges `resize(width = 1)`.
    pub tss: u64,
}

impl Transcript {
    /// 1-based inclusive span length, matching the value ChIPseeker reports
    /// (`geneEnd - geneStart + 1` over TxDb 1-based coordinates).
    #[must_use]
    pub fn length_1based(&self) -> u64 {
        self.end - self.start
    }
}

/// The fully assembled model: transcripts for nearest-TSS search, plus typed
/// feature regions for category overlap.
#[derive(Debug, Default)]
pub struct GeneModel {
    pub transcripts: Vec<Transcript>,
    pub exons: Vec<Region>,
    pub introns: Vec<Region>,
    pub utr5: Vec<Region>,
    pub utr3: Vec<Region>,
}

/// Strand-aware TSS from a 0-based half-open span.
#[must_use]
pub fn tss_of(start: u64, end: u64, strand: Strand) -> u64 {
    if strand == Strand::Forward {
        start
    } else {
        end - 1
    }
}
