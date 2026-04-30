# fastx-tools

Command-line tools for simple manipulation of FASTA and FASTQ sequence files. Currently only a subsetting tool is available.

## Requirements

- [Rust](https://www.rust-lang.org/tools/install) 1.85 or later (edition 2024)

## Installation

Clone the repository and build a release binary:

```sh
git clone <repository-url>
cd fastx-tools
cargo install --path .
```

This compiles all binaries and installs them into `~/.cargo/bin/`, which is on your `PATH` if you installed Rust with `rustup`.

To install to a custom location, use the `--root` flag:

```sh
cargo install --path . --root /usr/local
```

To build without installing (output in `target/release/`):

```sh
cargo build --release
```

## Tools

### `subsetfx` — subset FASTA/FASTQ files

Extract a subset of records from a FASTA or FASTQ file, with optional subsetting of each sequence.

```
subsetfx --input-fx-file <input> --output-fx-file <output> --file-type <fasta|fastq> [options]
```

**Required flags**

| Flag               | Description               |
|--------------------|---------------------------|
| `--input-fx-file`  | Input FASTA or FASTQ file |
| `--output-fx-file` | Output file path          |
| `--file-type`      | `fasta` or `fastq`        |

**Optional flags**

| Flag                  | Description                                              |
|-----------------------|----------------------------------------------------------|
| `--subset-file`       | File listing record names to keep (see format below)     |
| `--start`             | Subsequence start position (1-based inclusive)           |
| `--end`               | Subsequence end position (1-based inclusive)             |
| `--keep-record-order` | Write output records in the same order as the input file |

**Subset file format**

Each line names one record to keep. Optionally include start and end positions to extract a subsequence of that record.

Names must include the record-type prefix: `>` for FASTA records, `@` for FASTQ records.

```
>seq1
>seq2 4 10
>seq3
```

Lines beginning with `#` and blank lines are ignored. Positions are 1-based inclusive. A line with only a name copies the full sequence. A line with name, start, and end extracts that subsequence. Two fields (name + one position) are an error. If more than three fields are provided, the additional fields are silently ignored.

**Flag interaction rules**

- If `--start` or `--end` is given alongside `--subset-file`, the global range is applied to all selected records and any per-record ranges in the subset file are ignored. Note that in this case a subset file still cannot have two fields on a line as explained above.
- If neither `--subset-file`, `--start`, nor `--end` is provided, all records are copied and a warning is printed to `stderr`.
- Missing `--end` runs to the end of each sequence; missing `--start` runs from the beginning.

**Examples**

Copy two records in full:
```sh
subsetfx --input-fx-file in.fasta --output-fx-file out.fasta --file-type fasta \
    --subset-file names.txt
```

Extract positions 10–50 from every record:
```sh
subsetfx --input-fx-file in.fastq --output-fx-file out.fastq --file-type fastq \
    --start 10 --end 50
```

Select named records and trim to positions 10–50, preserving source order:
```sh
subsetfx --input-fx-file in.fasta --output-fx-file out.fasta --file-type fasta \
    --subset-file names.txt --start 10 --end 50 --keep-record-order
```

## Development

```sh
cargo build          # debug build
cargo test           # run all tests
cargo clippy         # lint
cargo fmt            # format
```
