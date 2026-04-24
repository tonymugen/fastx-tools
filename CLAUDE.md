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
cargo run --bin subsetfx -- <args>

# Lint
cargo clippy

# Format
cargo fmt
```

## Architecture

- `src/lib.rs` — shared parsing and manipulation logic for FASTA/FASTQ records; exposed as a library crate
- `src/bin/subsetfx.rs` — CLI binary for subsetting FASTA and FASTQ files; links against the library crate

New binaries for other tools should be added as `src/bin/<toolname>.rs` and share common logic through `src/lib.rs`.

## Library structure (`src/lib.rs`)

All types live in the `fastx` module (`pub mod fastx`).

### Types

- `NameWithRange` — struct with `name: String`, `start: usize`, `end: usize`; used to pass per-record ranges to subsetting methods; `start` is inclusive, `end` is exclusive
- `IndexedSequence` — FASTA record: sequence string plus its base-0 index in the source file; implements `FastxRecord`; derives `Clone`, `Debug`
- `IndexedSequenceWithQuality` — FASTQ record: sequence, quality scores, and base-0 index; implements `FastxRecord`; derives `Clone`, `Debug`
- `FastxRecords<T: FastxRecord>` — generic collection of records keyed by header name; replaces the former `FastaRecords` and `FastqRecords` types; derives `Debug`

### Trait

- `FastxRecord` — required bound for `T` in `FastxRecords<T>`; requires `Clone` and the methods `get_index() -> u32`, `get_sequence() -> &str`, `subsequence(usize, usize) -> Self`, `format_output(&str) -> String`

### Free constructor functions

- `read_fasta(path: &str) -> Result<FastxRecords<IndexedSequence>, String>`
- `read_fastq(path: &str) -> Result<FastxRecords<IndexedSequenceWithQuality>, String>`

### Key design decisions

- FASTA headers are stored including the leading `>` as the HashMap key.
- FASTQ headers are stored including the leading `@` as the HashMap key.
- `subsequence()` on individual records returns a new object of the same type (not a `&str`); use `.get_sequence()` / `.get_quality_scores()` to extract strings.
- `records_by_name` and `subsequences_by_name` return `(FastxRecords<T>, String)` where the second element is absent names joined by `\n`. Duplicate present names are deduplicated by the HashMap; duplicate absent names appear multiple times in the string.
- `read_fasta` returns an error for empty files or files with no valid records (header + sequence pairs). Multi-line sequences are concatenated.
- `read_fastq` parses strict 4-line-per-record FASTQ; returns errors for missing lines, non-`@` headers, invalid `+` separator, and quality/sequence length mismatches. Blank lines between records are skipped.
- `FastxRecords` tracks `max_sequence_length_`; accessible via `get_max_length() -> usize`. Updated in all methods that produce a new `FastxRecords`.
- `clamp_slice` is a private helper that safely slices strings with out-of-bounds and inverted range inputs.

### Test data

- `tests/data/correct.fasta` — 18 HIV sequence records
- `tests/data/correct.fastq` — 3 synthetic records for unit testing

### Dependencies

- `[dev-dependencies]`: `tempfile = "3"` — used in error-path tests to create and auto-delete temporary files

## Binary structure (`src/bin/subsetfx.rs`)

### Types

- `NamesWithRanges` — struct with `name: String`, `start: usize`, `end: usize`; holds a sequence name and optional subset range parsed from the name file; `start` is inclusive, `end` is exclusive
- `ParsedCLIFlags` — holds parsed and type-classified command line flags in three internal `HashMap`s: string, integer (`usize`), and bool values

### Functions

- `parse_name_file(path: &str) -> Result<Vec<NamesWithRanges>, String>` — parses a whitespace-delimited name file; lines starting with `#` and blank lines are skipped; one field = name only (start/end = 0); two fields = error; three or more fields = name, start, end (user-provided end is 1-based inclusive, stored as `end + 1`)

### `ParsedCLIFlags` methods

- `new(command_line: Vec<String>) -> Result<ParsedCLIFlags, String>` — parses and validates `std::env::args()`
- `get_string_value(flag: &str) -> Result<&str, String>`
- `get_int_value(flag: &str) -> Result<usize, String>`
- `get_bool_value(flag: &str) -> Result<bool, String>`

### Supported CLI flags

| Flag | Type | Required | Description |
|------|------|----------|-------------|
| `--input-fx-file` | string | yes | Input FASTA/FASTQ file path |
| `--output-fx-file` | string | yes | Output file path |
| `--file-type` | string | yes | `fasta` or `fastq` |
| `--subset-file` | string | no | File with sequence names and optional ranges |
| `--start` | usize | no | Global start position (0-based inclusive) |
| `--end` | usize | no | Global end position (0-based exclusive) |
| `--keep-record-order` | bool | no | Preserve input file record order in output |

### CLI parsing design decisions

- `--keep-record-order` is a presence flag; providing a value is an error.
- Duplicate flags silently overwrite (last value wins); documented in usage instructions.
- Unknown flags, values before any flag name, and multiple values for one flag all return errors.
- If none of `--subset-file`, `--start`, or `--end` are provided, the caller copies the full input file with a warning.
