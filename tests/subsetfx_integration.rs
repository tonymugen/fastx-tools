use assert_cmd::Command;
use std::io::Write;
use tempfile::NamedTempFile;

const SMALL_FASTA: &str = "tests/data/small.fasta";
const SMALL_FASTQ: &str = "tests/data/small.fastq";

fn subsetfx() -> Command {
    Command::cargo_bin("subsetfx").unwrap()
}

// ── FASTA happy-path tests ────────────────────────────────────────────────────

#[test]
fn test_fasta_subset_by_name_selects_correct_records() {
    let output = NamedTempFile::new().unwrap();
    let mut names = NamedTempFile::new().unwrap();
    writeln!(names, ">seq1").unwrap();
    writeln!(names, ">seq3").unwrap();

    subsetfx().args([
        "--input-fx-file",  SMALL_FASTA,
        "--output-fx-file", output.path().to_str().unwrap(),
        "--file-type",      "fasta",
        "--subset-file",    names.path().to_str().unwrap(),
    ]).assert().success();

    let content = std::fs::read_to_string(output.path()).unwrap();
    assert!(content.contains(">seq1"));
    assert!(content.contains(">seq3"));
    assert!(!content.contains(">seq2"));
    assert!(!content.contains(">seq4"));
}

#[test]
fn test_fasta_global_start_and_end_extracts_subsequence() {
    let output = NamedTempFile::new().unwrap();

    subsetfx().args([
        "--input-fx-file",  SMALL_FASTA,
        "--output-fx-file", output.path().to_str().unwrap(),
        "--file-type",      "fasta",
        "--start",          "4",
        "--end",            "6",
    ]).assert().success();

    let content = std::fs::read_to_string(output.path()).unwrap();
    // 1-based [4,6] → 0-based [3,6):
    // seq1 AAACCCGGGTTT → CCC
    // seq2 TTTGGGCCCAAA → GGG
    // seq3 CCCTTTAAAGGG → TTT
    // seq4 GGGAAATTTCCC → AAA
    assert!(content.contains("CCC"));
    assert!(content.contains("GGG"));
    assert!(content.contains("TTT"));
    assert!(content.contains("AAA"));
    assert!(!content.contains("AAACCCGGGTTT"));
}

#[test]
fn test_fasta_global_start_only_runs_to_end_of_sequence() {
    let output = NamedTempFile::new().unwrap();

    subsetfx().args([
        "--input-fx-file",  SMALL_FASTA,
        "--output-fx-file", output.path().to_str().unwrap(),
        "--file-type",      "fasta",
        "--start",          "7",
    ]).assert().success();

    let content = std::fs::read_to_string(output.path()).unwrap();
    // seq1 AAACCCGGGTTT → start 7 (1-based) → 0-based 6 → "GGGTTT"
    assert!(content.contains("GGGTTT"));
    assert!(!content.contains("AAACCC"));
}

#[test]
fn test_fasta_global_end_only_runs_from_start_of_sequence() {
    let output = NamedTempFile::new().unwrap();

    subsetfx().args([
        "--input-fx-file",  SMALL_FASTA,
        "--output-fx-file", output.path().to_str().unwrap(),
        "--file-type",      "fasta",
        "--end",            "3",
    ]).assert().success();

    let content = std::fs::read_to_string(output.path()).unwrap();
    // seq1 AAACCCGGGTTT → end 3 (1-based inclusive = 0-based exclusive 3) → "AAA"
    assert!(content.contains("AAA"));
    assert!(!content.contains("AAACCC"));
}

#[test]
fn test_fasta_per_record_ranges_in_subset_file() {
    let output = NamedTempFile::new().unwrap();
    let mut ranges = NamedTempFile::new().unwrap();
    writeln!(ranges, ">seq1 4 6").unwrap(); // AAACCCGGGTTT [4,6] → "CCC"
    writeln!(ranges, ">seq2 1 3").unwrap(); // TTTGGGCCCAAA [1,3] → "TTT"

    subsetfx().args([
        "--input-fx-file",  SMALL_FASTA,
        "--output-fx-file", output.path().to_str().unwrap(),
        "--file-type",      "fasta",
        "--subset-file",    ranges.path().to_str().unwrap(),
    ]).assert().success();

    let content = std::fs::read_to_string(output.path()).unwrap();
    assert!(content.contains("CCC"));
    assert!(content.contains("TTT"));
    assert!(!content.contains("AAACCCGGGTTT"));
    assert!(!content.contains("TTTGGGCCCAAA"));
}

#[test]
fn test_fasta_global_range_overrides_subset_file_ranges() {
    let output = NamedTempFile::new().unwrap();
    let mut subset = NamedTempFile::new().unwrap();
    writeln!(subset, ">seq1 1 3").unwrap(); // per-record range should be ignored
    writeln!(subset, ">seq2").unwrap();

    subsetfx().args([
        "--input-fx-file",  SMALL_FASTA,
        "--output-fx-file", output.path().to_str().unwrap(),
        "--file-type",      "fasta",
        "--subset-file",    subset.path().to_str().unwrap(),
        "--start",          "4",
        "--end",            "6",
    ]).assert().success();

    let content = std::fs::read_to_string(output.path()).unwrap();
    // global [4,6]: seq1→"CCC", seq2→"GGG"; per-record [1,3]→"AAA" was overridden
    assert!(content.contains("CCC"));
    assert!(content.contains("GGG"));
    assert!(!content.contains("AAA"));
}

#[test]
fn test_fasta_keep_record_order_preserves_source_order() {
    let output = NamedTempFile::new().unwrap();

    subsetfx().args([
        "--input-fx-file",  SMALL_FASTA,
        "--output-fx-file", output.path().to_str().unwrap(),
        "--file-type",      "fasta",
        "--keep-record-order",
    ]).assert().success();

    let content = std::fs::read_to_string(output.path()).unwrap();
    let headers: Vec<&str> = content.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers, vec![">seq1", ">seq2", ">seq3", ">seq4"]);
}

// ── FASTQ happy-path tests ────────────────────────────────────────────────────

#[test]
fn test_fastq_subset_by_name_selects_correct_records() {
    let output = NamedTempFile::new().unwrap();
    let mut names = NamedTempFile::new().unwrap();
    writeln!(names, "@record1").unwrap();
    writeln!(names, "@record3").unwrap();

    subsetfx().args([
        "--input-fx-file",  SMALL_FASTQ,
        "--output-fx-file", output.path().to_str().unwrap(),
        "--file-type",      "fastq",
        "--subset-file",    names.path().to_str().unwrap(),
    ]).assert().success();

    let content = std::fs::read_to_string(output.path()).unwrap();
    assert!(content.contains("@record1"));
    assert!(content.contains("@record3"));
    assert!(!content.contains("@record2"));
}

#[test]
fn test_fastq_global_start_and_end_extracts_subsequence() {
    let output = NamedTempFile::new().unwrap();

    subsetfx().args([
        "--input-fx-file",  SMALL_FASTQ,
        "--output-fx-file", output.path().to_str().unwrap(),
        "--file-type",      "fastq",
        "--start",          "1",
        "--end",            "4",
    ]).assert().success();

    let content = std::fs::read_to_string(output.path()).unwrap();
    // record1 ACGTACGTACGT → [1,4] → "ACGT"
    assert!(content.contains("ACGT"));
    assert!(!content.contains("ACGTACGTACGT"));
}

// ── Error / rejection tests ───────────────────────────────────────────────────

#[test]
fn test_missing_required_flag_exits_nonzero() {
    let output = NamedTempFile::new().unwrap();
    // --file-type is missing
    subsetfx().args([
        "--input-fx-file",  SMALL_FASTA,
        "--output-fx-file", output.path().to_str().unwrap(),
    ]).assert().failure();
}

#[test]
fn test_bad_input_file_exits_nonzero() {
    let output = NamedTempFile::new().unwrap();
    subsetfx().args([
        "--input-fx-file",  "/nonexistent/file.fasta",
        "--output-fx-file", output.path().to_str().unwrap(),
        "--file-type",      "fasta",
    ]).assert().failure();
}

#[test]
fn test_bad_subset_file_exits_nonzero() {
    let output  = NamedTempFile::new().unwrap();
    let gone    = { let t = NamedTempFile::new().unwrap(); t.path().to_str().unwrap().to_string() };
    subsetfx().args([
        "--input-fx-file",  SMALL_FASTA,
        "--output-fx-file", output.path().to_str().unwrap(),
        "--file-type",      "fasta",
        "--subset-file",    &gone,
    ]).assert().failure();
}
