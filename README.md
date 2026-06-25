# rsomics-peak-annotate

Annotate ChIP-seq / ATAC-seq peaks with their nearest gene, signed distance to
the TSS, and a single genomic-feature category — a Rust port of Bioconductor
**ChIPseeker** `annotatePeak`. Output is **value-exact** to ChIPseeker 1.46.1,
at roughly **160× the throughput** and **a fraction of the memory**.

```sh
cargo install rsomics-peak-annotate
```

## Usage

```sh
rsomics-peak-annotate --peaks peaks.narrowPeak --gtf genes.gtf > annotated.tsv
rsomics-peak-annotate --peaks peaks.bed --gtf genes.gtf --tss-region=-3000,3000 --output out.tsv
```

Peaks are read as BED / narrowPeak (0-based half-open; all input columns pass
through unchanged). The gene model is a GTF (1-based). For each peak the tool
appends: `annotation`, `geneChr`, `geneStart`, `geneEnd`, `geneLength`,
`geneStrand`, `geneId`, `transcriptId`, `distanceToTSS`.

`--tss-region up,down` sets the promoter window (default `-3000,3000`).

## Semantics (matched to ChIPseeker 1.46.1)

- Each peak is anchored on the nearest **transcript** TSS (strand-aware 5' end);
  the reported `geneStart`/`geneEnd` are that transcript's span (`level="transcript"`).
- `distanceToTSS` is measured from the peak edge nearest the TSS, signed
  (negative = upstream / 5', positive = downstream / 3', 0 when the peak spans it).
- Category priority (one label per peak): Promoter (within the TSS region,
  sub-binned `<=1kb` / `1-2kb` / `2-3kb`) → 5'/3' UTR → Exon `(tx/gene, exon N of M)`
  → Intron `(tx/gene, intron N of M)` → Downstream `(<=300bp)` → Distal Intergenic.
  UTRs are derived from the CDS/exon boundaries. Downstream is genomic (within
  300 bp to the right of a transcript's genomic end), matching ChIPseeker.

## Performance

On a 2000-gene GTF + 20,000 peaks (mini_m2, Apple Silicon), versus ChIPseeker
1.46.1 (the full `makeTxDbFromGFF` + `annotatePeak` pipeline):

| | ours | ChIPseeker 1.46.1 | ratio |
|---|---|---|---|
| wall | 0.07 s | 11.25 s | ~160× faster |
| peak RSS | 11.5 MB | 1015 MB | ~88× smaller |

Every one of the 20,000 peaks' annotation, nearest transcript, gene bounds,
strand, and signed TSS distance is identical to ChIPseeker. The committed golden
(`tests/golden/`, covering all categories on both strands plus a multi-transcript
gene) is checked in CI via `tests/compat.rs` — no R needed at test time.

## Origin

This crate is an independent Rust reimplementation of ChIPseeker `annotatePeak`,
informed by:

- The published method: Yu, Wang & He, *ChIPseeker: an R/Bioconductor package for
  ChIP peak annotation, comparison and visualization*, Bioinformatics 2015,
  [doi:10.1093/bioinformatics/btv145](https://doi.org/10.1093/bioinformatics/btv145).
- The ChIPseeker source (Artistic-2.0, a permissive license that allows reading
  and citing): `getNearestFeatureIndicesAndDistances`, `getGenomicAnnotation`,
  `annotatePeak`.
- Black-box behaviour testing against ChIPseeker 1.46.1 (`tests/golden/run_chipseeker.R`).

License: MIT OR Apache-2.0.
Upstream credit: [ChIPseeker](https://github.com/YuLab-SMU/ChIPseeker) (Artistic-2.0).
