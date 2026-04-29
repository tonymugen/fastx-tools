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
- `FastxRecords<T: FastxRecord>` — generic collection of records keyed by header name; derives `Debug`

### Trait

- `FastxRecord` — required bound for `T` in `FastxRecords<T>`; requires `Clone` and the methods `get_index() -> u32`, `get_sequence() -> &str`, `subsequence(usize, usize) -> Self`, `format_output(&str) -> String`

### `FastxRecords<T>` methods

- `num_records() -> usize`
- `get_max_length() -> usize`
- `records_by_name(names: Vec<String>) -> (FastxRecords<T>, String)` — returns matching records plus absent names joined by `\n`
- `subsequences(start: usize, end: usize) -> FastxRecords<T>` — returns all records trimmed to `[start, end)`
- `subsequences_by_name(names_ranges: Vec<NameWithRange>) -> (FastxRecords<T>, String)` — returns per-record subsequences plus absent names joined by `\n`
- `save_records(output_path: &str) -> Result<(), String>` — writes records to file in arbitrary (HashMap) order
- `save_sorted_records(output_path: &str) -> Result<(), String>` — writes records in original source-file order (sorted by `get_index()`)

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

- `tests/data/small.fasta` — 4 synthetic FASTA records for integration testing
- `tests/data/small.fastq` — 3 synthetic FASTQ records for integration testing

### Dependencies

- `[dev-dependencies]`: `tempfile = "3"` — used in error-path tests to create and auto-delete temporary files

## Binary structure (`src/bin/subsetfx.rs`)

### Types

- `ParsedCLIFlags` — holds parsed and type-classified command line flags in three internal `HashMap`s: string, integer (`usize`), and bool values

The binary does not define its own record struct; it uses `fastx::NameWithRange` from the library directly.

### Functions

- `parse_name_file(path: &str) -> Result<Vec<fastx::NameWithRange>, String>` — parses a whitespace-delimited name file; lines starting with `#` and blank/whitespace-only lines are skipped; one field = name only (start/end = 0); two fields = error; three or more fields = name, start, end stored as raw 1-based values (zero is rejected for start and end); extra fields beyond three are ignored
- `split_subset_list(records: Vec<fastx::NameWithRange>) -> (Vec<String>, Vec<fastx::NameWithRange>)` — splits records into (names to copy in full, records with ranges); converts the 1-based inclusive start from `parse_name_file` to 0-based by subtracting 1; end is left unchanged (1-based inclusive == 0-based exclusive)
- `run<T: fastx::FastxRecord>(input_records: FastxRecords<T>, cli_args: &ParsedCLIFlags) -> Result<(), String>` — generic subsetting driver; dispatched separately per file type in `main()` because `FastxRecord` is not object-safe

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
| `--start` | usize | no | Global start position (1-based inclusive) |
| `--end` | usize | no | Global end position (1-based inclusive) |
| `--keep-record-order` | bool | no | Preserve input file record order in output |

### CLI parsing design decisions

- `--keep-record-order` is a presence flag; providing a value is an error.
- Duplicate flags silently overwrite (last value wins).
- Unknown flags, values before any flag name, and multiple values for one flag all return errors.
- `--start` and `--end` are 1-based inclusive on the command line, matching the subset file convention.
- Missing `--end` defaults to `usize::MAX`, clamped to `get_max_length()` before slicing, so omitting `--end` runs to the end of each sequence.
- If neither `--subset-file`, `--start`, nor `--end` is provided, all records are copied unchanged and a warning is printed to stderr.
- If `--start` or `--end` is present alongside `--subset-file`, the global range overrides any per-record ranges in the subset file.

### `run()` flag-combination behaviour

| `--subset-file` | `--start`/`--end` | Behaviour |
|---|---|---|
| absent | absent | copy all records unchanged (warning to stderr) |
| absent | present | apply global subsequence to all records |
| present (name-only lines) | absent | copy named records in full |
| present (range lines) | absent | apply per-record subsequences |
| present (mixed) | absent | full-copy name-only lines; per-record subseq for range lines |
| present | present | filter to named records, then apply global subsequence (per-record ranges ignored) |

### Integration tests

Integration tests live in `tests/subsetfx_integration.rs` and use `assert_cmd::Command::cargo_bin("subsetfx")`. They cover all flag combinations above plus three error-exit cases (missing required flag, bad input file, bad subset file).
