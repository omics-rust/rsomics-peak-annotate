mod cli;

use std::process::ExitCode;

use clap::Parser;
use rsomics_common::Tool;

fn main() -> ExitCode {
    cli::Cli::parse().run()
}
