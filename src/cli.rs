use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Context, Result, Tool, ToolMeta};

use rsomics_peak_annotate::{annotate_peaks, parse_tss_region};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-peak-annotate",
    version,
    about = "Annotate ChIP/ATAC peaks with nearest gene, TSS distance, and genomic-feature category — port of ChIPseeker annotatePeak",
    long_about = None
)]
pub struct Cli {
    /// Peak file (BED / narrowPeak, 0-based half-open). Columns pass through.
    #[arg(long = "peaks")]
    pub peaks: PathBuf,

    /// Gene model GTF (1-based inclusive).
    #[arg(long = "gtf")]
    pub gtf: PathBuf,

    /// TSS region window as `upstream,downstream` (signed bp).
    #[arg(long = "tss-region", default_value = "-3000,3000")]
    pub tss_region: String,

    /// Output path (default: stdout).
    #[arg(long = "output")]
    pub output: Option<PathBuf>,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let tss_region = parse_tss_region(&self.tss_region)?;

        let mut out: Box<dyn Write> = match &self.output {
            Some(p) => {
                Box::new(BufWriter::new(File::create(p).rs_with_context(|| {
                    format!("creating output {}", p.display())
                })?))
            }
            None => Box::new(io::stdout().lock()),
        };

        annotate_peaks(&self.peaks, &self.gtf, tss_region, &mut out)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn parses_required_flags() {
        let cli = Cli::parse_from([
            "rsomics-peak-annotate",
            "--peaks",
            "p.bed",
            "--gtf",
            "g.gtf",
        ]);
        assert_eq!(cli.peaks, PathBuf::from("p.bed"));
        assert_eq!(cli.gtf, PathBuf::from("g.gtf"));
        assert_eq!(cli.tss_region, "-3000,3000");
        assert!(cli.output.is_none());
    }
}
