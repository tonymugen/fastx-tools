# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A Rust library and CLI toolset for manipulating FASTA and FASTQ files (common bioinformatics sequence formats). The library is in `src/lib.rs`; binary tools live in `src/bin/`.

## Commands

```bash
# Build
cargo build

# Build release
cargo build --release

# Run tests
cargo test

# Run a single test
cargo test <test_name>

# Run a specific binary
cargo run --bin subsetfa -- <args>

# Lint
cargo clippy

# Format
cargo fmt
```

## Architecture

- `src/lib.rs` — shared parsing and manipulation logic for FASTA/FASTQ records; exposed as a library crate
- `src/bin/subsetfa.rs` — CLI binary for subsetting FASTA files; links against the library crate

New binaries for other tools should be added as `src/bin/<toolname>.rs` and share common logic through `src/lib.rs`.
