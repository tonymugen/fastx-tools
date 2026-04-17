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

## Library structure (`src/lib.rs`)

All types live in the `fastx` module (`pub mod fastx`).

### Types

- `NameWithRange` — struct with `name: String`, `start: usize`, `end: usize`; used to pass per-record ranges to subsetting methods
- `IndexedSequence` — FASTA record: sequence string plus its base-0 index in the source file; derives `Clone`, `Debug`
- `IndexedSequenceWithQuality` — FASTQ record: sequence, quality scores, and base-0 index; derives `Clone`, `Debug`
- `FastaRecords` — collection of `IndexedSequence` records keyed by header name; derives `Debug`
- `FastqRecords` — collection of `IndexedSequenceWithQuality` records keyed by header name; derives `Debug`

### Key design decisions

- FASTA headers are stored including the leading `>` as the HashMap key.
- FASTQ headers are stored including the leading `@` as the HashMap key.
- `subsequence()` on individual records returns a new object of the same type (not a `&str`); use `.get_sequence()` / `.get_quality_scores()` to extract strings.
- `records_by_name` and `subsequences_by_name` return `(RecordCollection, String)` where the second element is absent names joined by `\n`.
- `FastaRecords::new` returns an error for empty files or files with no valid records (header + sequence pairs).
- `FastqRecords::new` parses strict 4-line-per-record FASTQ; returns errors for missing lines, non-`@` headers, and quality/sequence length mismatches.

### Test data

- `tests/data/correct.fasta` — 18 HIV sequence records
- `tests/data/correct.fastq` — 3 synthetic records for unit testing

### Dependencies

- `[dev-dependencies]`: `tempfile = "3"` — used in error-path tests to create and auto-delete temporary files
