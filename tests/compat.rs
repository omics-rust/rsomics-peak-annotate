//! Value-exact compat against ChIPseeker 1.46.1. Runs the binary on the
//! committed fixture (`tests/golden/{peaks.bed,anno.gtf}`) and checks every
//! peak's category, nearest gene, gene bounds, strand, and signed TSS distance
//! against the checked-in `golden.tsv` (ChIPseeker `annotatePeak`,
//! tssRegion=(-3000,3000), level="transcript"). Always-run; no R needed.
//!
//! ChIPseeker's `as.data.frame` encodes geneChr as a factor int and geneStrand
//! as 1(+)/2(-); we emit chrom names and +/-, so those two columns are mapped.

use std::collections::HashMap;
use std::path::Path;
use std::process::Command;

const GOLDEN_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden");

/// The value-exact fields, keyed by peak name.
#[derive(PartialEq, Eq, Debug)]
struct Fields {
    annotation: String,
    gene_start: String,
    gene_end: String,
    gene_length: String,
    gene_strand: String, // normalized to "+"/"-"
    gene_id: String,
    transcript_id: String,
    distance: String,
}

fn parse_golden(text: &str) -> HashMap<String, Fields> {
    let mut out = HashMap::new();
    for line in text.lines().skip(1) {
        let c: Vec<&str> = line.split('\t').collect();
        // seqnames start end V4 annotation geneChr geneStart geneEnd geneLength
        // geneStrand geneId transcriptId distanceToTSS
        out.insert(
            c[3].to_string(),
            Fields {
                annotation: c[4].to_string(),
                gene_start: c[6].to_string(),
                gene_end: c[7].to_string(),
                gene_length: c[8].to_string(),
                gene_strand: if c[9] == "1" { "+" } else { "-" }.to_string(),
                gene_id: c[10].to_string(),
                transcript_id: c[11].to_string(),
                distance: c[12].to_string(),
            },
        );
    }
    out
}

fn parse_ours(text: &str) -> HashMap<String, Fields> {
    let mut out = HashMap::new();
    for line in text.lines() {
        let c: Vec<&str> = line.split('\t').collect();
        // BED6 passthrough (chrom start end name score strand) + annotation
        // geneChr geneStart geneEnd geneLength geneStrand geneId transcriptId distanceToTSS
        out.insert(
            c[3].to_string(),
            Fields {
                annotation: c[6].to_string(),
                gene_start: c[8].to_string(),
                gene_end: c[9].to_string(),
                gene_length: c[10].to_string(),
                gene_strand: c[11].to_string(),
                gene_id: c[12].to_string(),
                transcript_id: c[13].to_string(),
                distance: c[14].to_string(),
            },
        );
    }
    out
}

#[test]
fn matches_chipseeker_golden() {
    let dir = Path::new(GOLDEN_DIR);
    let out_dir = tempfile::tempdir().expect("tempdir");
    let out_path = out_dir.path().join("ours.tsv");

    let status = Command::new(env!("CARGO_BIN_EXE_rsomics-peak-annotate"))
        .arg("--peaks")
        .arg(dir.join("peaks.bed"))
        .arg("--gtf")
        .arg(dir.join("anno.gtf"))
        .arg("--tss-region=-3000,3000")
        .arg("--output")
        .arg(&out_path)
        .status()
        .expect("run rsomics-peak-annotate");
    assert!(status.success(), "binary exited non-zero");

    let golden =
        parse_golden(&std::fs::read_to_string(dir.join("golden.tsv")).expect("read golden.tsv"));
    let ours = parse_ours(&std::fs::read_to_string(&out_path).expect("read ours.tsv"));

    assert_eq!(ours.len(), golden.len(), "peak count");
    for (name, want) in &golden {
        let got = ours
            .get(name)
            .unwrap_or_else(|| panic!("peak {name} missing from output"));
        assert_eq!(got, want, "peak {name} diverges from ChIPseeker golden");
    }
}
