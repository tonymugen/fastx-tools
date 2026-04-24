//! CLI tool to create subsets of FASTA and FASTQ files.

use std::env;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use fastx_tools::fastx;

/// Parse the input sequence list file.
///
/// # Arguments
///
/// * `name_file` path to the file that has sequence names and ranges to subset.
///
/// # Returns
///
/// A vector with `fastx::NameWithRange` objects. The sequence range reported is base-1 start and
/// end (inclusive). This will be changed before submitting for subsetting.
///
pub fn parse_name_file(name_file: &str) -> Result<Vec<fastx::NameWithRange>, String> {
    let mut names_with_ranges: Vec<fastx::NameWithRange> = Vec::new();
    let file = File::open(name_file).map_err(|e| e.to_string())?;
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|error| error.to_string())?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        match fields.len() {
            0 => continue,
            1 => names_with_ranges.push( fastx::NameWithRange {
                name: fields[0].to_string(),
                start: 0,
                end: 0,
            }),
            2 => return Err( format!("Ambiguous start/end on line: {line}") ),
            _ => {
                let start = fields[1].parse::<usize>()
                    .map_err( |error| format!("Invalid start '{}': {error}", fields[1]) )?;
                if start == 0 {
                    return Err( format!("Start value must be at least 1 on line: {line}") );
                }
                let end_field = fields[2].parse::<usize>()
                    .map_err( |error| format!("Invalid end '{}': {error}", fields[2]) )?;
                if end_field == 0 {
                    return Err( format!("End value must be at least 1 on line: {line}") );
                }
                names_with_ranges.push( fastx::NameWithRange {
                    name: fields[0].to_string(),
                    start,
                    end: end_field,
                });
            }
        }
    }
    Ok( names_with_ranges )
}

/// Split subsetting records.
///
/// Splits the input vector of `NameWithRange` object into a vector of sequence names to copy
/// completely and a vector of `NameWithRange` objects to copy subsequences. The latter will then
/// have base-0 start and end (exclusive).
///
/// # Returns
///
/// A tuple of (copy completely, copy subsequences).
///
pub fn split_subset_list(records: Vec<fastx::NameWithRange>) -> (Vec<String>, Vec<fastx::NameWithRange>) {
    let (with_range, name_only): (Vec<_>, Vec<_>) = records.into_iter()
        .partition(|r| r.start > 0 || r.end > 0);
    let names = name_only.into_iter().map(|r| r.name).collect();
    let ranges = with_range.into_iter()
        .map(|r| fastx::NameWithRange { start: r.start - 1, ..r })
        .collect();
    (names, ranges)
}

/// Holds parsed command line arguments.
#[derive(Debug)]
pub struct ParsedCLIFlags {
    string_values_: HashMap<String, String>,
    int_values_:    HashMap<String, usize>,
    bool_values_:   HashMap<String, bool>,
}

impl ParsedCLIFlags {
    fn parse_flag_pairs(args: &[String], bool_flags: &[&str]) -> Result<HashMap<String, String>, String> {
        let mut flags_and_fields: HashMap<String, String> = HashMap::new();
        let mut current_flag_name = String::new();
        let mut current_flag_value = String::new();
        for field in args {
            if field.starts_with("--") {
                if !current_flag_name.is_empty() {
                    flags_and_fields.insert( current_flag_name.clone(), current_flag_value.clone() );
                    current_flag_value.clear();
                }
                current_flag_name.clone_from(field);
                continue;
            }
            if current_flag_name.is_empty() {
                return Err( format!("Got value {field} before any flag name") );
            }
            if bool_flags.contains(&current_flag_name.as_str()) {
                return Err( format!("{} takes no value, got: {field}", current_flag_name) );
            }
            if !current_flag_value.is_empty() {
                return Err( format!("{} got multiple values", current_flag_name) );
            }
            current_flag_value.push_str(field);
        }
        if !current_flag_name.is_empty() {
            flags_and_fields.insert( current_flag_name.clone(), current_flag_value.clone() );
        }
        Ok( flags_and_fields )
    }

    /// Creates a new `ParsedCLIFlags` object.
    ///
    /// # Arguments
    ///
    /// * `command_line` - The command line arguments from `std::env::args()`.
    ///
    /// # Returns
    ///
    /// A `ParsedCLIFlags` object or parsing errors.
    pub fn new(command_line: Vec<String>) -> Result<ParsedCLIFlags, String> {

        let required_flags = [
            "--input-fx-file",
            "--output-fx-file",
            "--file-type"
        ];
        let bool_flags = ["--keep-record-order"];
        if command_line.len() < 2 {
            return Err( "No command line arguments.".to_string() )
        }
        let flags_and_fields = Self::parse_flag_pairs(&command_line[1..], &bool_flags)?;
        // an array good enough here because it is small
        // only larger data sets benefit from HashSet
        let known_flags = [
            "--input-fx-file",
            "--output-fx-file",
            "--file-type",
            "--subset-file",
            "--start",
            "--end",
            "--keep-record-order",
        ];
        for flag in flags_and_fields.keys() {
            if !known_flags.contains( &flag.as_str() ) {
                return Err( format!("Unknown flag: {flag}") );
            }
        }
        let mut missing_required_flags = String::new();
        for required_flag in required_flags {
            match flags_and_fields.get(required_flag) {
                None => {
                    missing_required_flags.push_str(required_flag);
                    missing_required_flags.push('\n');
                }
                Some(value) if value.is_empty() => {
                    return Err( format!("{required_flag} requires a value") );
                }
                _ => {}
            }
        }
        if !missing_required_flags.is_empty() {
            return Err( format!("Missing required flags:\n{missing_required_flags}") )
        }
        let file_type = flags_and_fields.get("--file-type").unwrap();
        if file_type != "fasta" && file_type != "fastq" {
            return Err( format!("--file-type must be 'fasta' or 'fastq', got: {file_type}") );
        }
        let string_flags = [
            "--input-fx-file",
            "--output-fx-file",
            "--file-type",
            "--subset-file"
        ];
        let mut string_flags_with_values: HashMap<String, String> = HashMap::new();
        for string_flag in string_flags {
            if let Some(value) = flags_and_fields.get(string_flag) {
                string_flags_with_values.insert( string_flag.to_string(), value.clone() );
            }
        }
        let int_flags = [
            "--start",
            "--end"
        ];
        let mut int_flags_with_values: HashMap<String, usize> = HashMap::new();
        for int_flag in int_flags {
            if let Some(value) = flags_and_fields.get(int_flag) {
                let int_value = value.parse::<usize>().map_err( |error| error.to_string() )?;
                int_flags_with_values.insert(int_flag.to_string(), int_value);
            }
        }
        let mut bool_flags_with_values: HashMap<String, bool> = HashMap::new();
        for bool_flag in bool_flags {
            if flags_and_fields.contains_key(bool_flag) {
                bool_flags_with_values.insert(bool_flag.to_string(), true);
            }
        }
        Ok(
            ParsedCLIFlags {
                string_values_: string_flags_with_values,
                int_values_:    int_flags_with_values,
                bool_values_:   bool_flags_with_values
            }
        )
     }

    /// Returns the value for a string-valued flag.
    ///
    /// # Arguments
    ///
    /// * `flag` - The flag name (e.g. `"--input-fx-file"`)
    ///
    /// # Returns
    ///
    /// The flag's string value, or an error if the flag was not provided.
    pub fn get_string_value(&self, flag: &str) -> Result<&str, String> {
        self.string_values_.get(flag)
            .map( |s| s.as_str() )
            .ok_or_else( || format!("{flag} not found") )
    }

    /// Returns the value for an integer-valued flag.
    ///
    /// # Arguments
    ///
    /// * `flag` - The flag name (e.g. `"--start"`)
    ///
    /// # Returns
    ///
    /// The flag's integer value, or an error if the flag was not provided.
    pub fn get_int_value(&self, flag: &str) -> Result<usize, String> {
        self.int_values_.get(flag)
            .copied()
            .ok_or_else( || format!("{flag} not found") )
    }

    /// Returns the value for a boolean flag.
    ///
    /// # Arguments
    ///
    /// * `flag` - The flag name (e.g. `"--keep-record-order"`)
    ///
    /// # Returns
    ///
    /// `true` if the flag was provided, or an error if it was absent.
    pub fn get_bool_value(&self, flag: &str) -> Result<bool, String> {
        self.bool_values_.get(flag)
            .copied()
            .ok_or_else( || format!("{flag} not found") )
    }
}

fn main() {
    let cli_flag_descritption = [
        "Command line flags (in any order):",
        "  --input-fx-file:      input FATSA or FASTQ file name (required)",
        "  --output-fx-file:     output FATSA or FASTQ file name (required)",
        "  --file-type:          file type (fasta or fastq) (required)",
        "  --subset-file:        file with names of records (optionally subsequence ranges) to save",
        "  --start:              subsequence start (base-1)",
        "  --end:                subsequence end (base-1 inclusive)",
        "  --keep-record-order:  if set (must have no value), will save the subset in the order of the original file",
        "",
        "NOTE: if either --start OR --end flags are present, ranges in the subset file will be overridden;",
        "      see the project README.md for details",
    ].join("\n");
    let cli_args = match ParsedCLIFlags::new( env::args().collect() ) {
        Ok(flags)  => flags,
        Err(error) => {
            eprintln!("ERROR: command line flag parsing failed: {error}");
            eprintln!("{cli_flag_descritption}");
            std::process::exit(1);
        }
    };
    let mut has_subset_file: bool    = true;
    let mut subset_file_name: String = String::new();
    match cli_args.get_string_value("--subset-file") {
        Ok(name) => { subset_file_name = name.to_string() },
        Err(_)   => {has_subset_file = false}
    }
    if has_subset_file {
        let records_to_save = match parse_name_file(&subset_file_name) {
            Ok(records) => records,
            Err(error)  => {
                eprintln!("Error parsing the subset file: {error}");
                std::process::exit(1);
            }
        };
        let mut has_start_flag: bool = true;
        let mut has_end_flag: bool = true;
        let subsequence_start = match cli_args.get_int_value("--start") {
            Ok(value) => value,
            Err(_)    => { has_start_flag = false; 0 }
        };
        let subsequence_end = match cli_args.get_int_value("--end") {
            Ok(value) => value,
            Err(_)    => { has_end_flag = false; 0 }
        };
    };

    todo!();
}

#[cfg(test)]
mod tests {
    use super::{ParsedCLIFlags, parse_name_file, split_subset_list};
    use fastx_tools::fastx;
    use tempfile::NamedTempFile;
    use std::io::Write;

    fn valid_args() -> Vec<String> {
        vec![
            "prog".to_string(),
            "--input-fx-file".to_string(), "in.fasta".to_string(),
            "--output-fx-file".to_string(), "out.fasta".to_string(),
            "--file-type".to_string(), "fasta".to_string(),
        ]
    }

    // ParsedCLIFlags::new tests

    #[test]
    fn test_new_no_arguments_returns_error() {
        let result = ParsedCLIFlags::new(vec!["prog".to_string()]);
        assert!(result.unwrap_err().contains("No command line arguments"));
    }

    #[test]
    fn test_new_value_before_flag_returns_error() {
        let result = ParsedCLIFlags::new(vec!["prog".to_string(), "orphan".to_string()]);
        assert!(result.unwrap_err().contains("before any flag name"));
    }

    #[test]
    fn test_new_multiple_values_for_flag_returns_error() {
        let result = ParsedCLIFlags::new(vec![
            "prog".to_string(),
            "--input-fx-file".to_string(), "a.fasta".to_string(), "b.fasta".to_string(),
        ]);
        assert!(result.unwrap_err().contains("got multiple values"));
    }

    #[test]
    fn test_new_unknown_flag_returns_error() {
        let mut args = valid_args();
        args.extend(["--unknown".to_string(), "val".to_string()]);
        assert!(ParsedCLIFlags::new(args).unwrap_err().contains("Unknown flag"));
    }

    #[test]
    fn test_new_missing_required_flag_returns_error() {
        let result = ParsedCLIFlags::new(vec![
            "prog".to_string(),
            "--output-fx-file".to_string(), "out.fasta".to_string(),
            "--file-type".to_string(), "fasta".to_string(),
        ]);
        assert!(result.unwrap_err().contains("--input-fx-file"));
    }

    #[test]
    fn test_new_multiple_missing_required_flags_all_reported() {
        let result = ParsedCLIFlags::new(vec![
            "prog".to_string(),
            "--file-type".to_string(), "fasta".to_string(),
        ]);
        let err = result.unwrap_err();
        assert!(err.contains("--input-fx-file"));
        assert!(err.contains("--output-fx-file"));
    }

    #[test]
    fn test_new_required_flag_with_empty_value_returns_error() {
        let result = ParsedCLIFlags::new(vec![
            "prog".to_string(),
            "--input-fx-file".to_string(),
            "--output-fx-file".to_string(), "out.fasta".to_string(),
            "--file-type".to_string(), "fasta".to_string(),
        ]);
        assert!(result.unwrap_err().contains("requires a value"));
    }

    #[test]
    fn test_new_flag_order_does_not_matter() {
        let args = vec![
            "prog".to_string(),
            "--file-type".to_string(), "fasta".to_string(),
            "--output-fx-file".to_string(), "out.fasta".to_string(),
            "--input-fx-file".to_string(), "in.fasta".to_string(),
        ];
        assert!(ParsedCLIFlags::new(args).is_ok());
    }

    #[test]
    fn test_new_invalid_file_type_returns_error() {
        let result = ParsedCLIFlags::new(vec![
            "prog".to_string(),
            "--input-fx-file".to_string(), "in.fasta".to_string(),
            "--output-fx-file".to_string(), "out.fasta".to_string(),
            "--file-type".to_string(), "bam".to_string(),
        ]);
        assert!(result.unwrap_err().contains("'fasta' or 'fastq'"));
    }

    #[test]
    fn test_new_invalid_start_value_returns_error() {
        let mut args = valid_args();
        args.extend(["--start".to_string(), "notanumber".to_string()]);
        assert!(ParsedCLIFlags::new(args).is_err());
    }

    #[test]
    fn test_new_valid_minimal_fasta_args_succeeds() {
        assert!(ParsedCLIFlags::new(valid_args()).is_ok());
    }

    #[test]
    fn test_new_valid_minimal_fastq_args_succeeds() {
        let args = vec![
            "prog".to_string(),
            "--input-fx-file".to_string(), "in.fastq".to_string(),
            "--output-fx-file".to_string(), "out.fastq".to_string(),
            "--file-type".to_string(), "fastq".to_string(),
        ];
        assert!(ParsedCLIFlags::new(args).is_ok());
    }

    #[test]
    fn test_new_all_optional_flags_accepted() {
        let mut args = valid_args();
        args.extend([
            "--subset-file".to_string(), "names.txt".to_string(),
            "--start".to_string(), "10".to_string(),
            "--end".to_string(), "100".to_string(),
            "--keep-record-order".to_string(),
        ]);
        assert!(ParsedCLIFlags::new(args).is_ok());
    }

    // get_string_value tests

    #[test]
    fn test_get_string_value_present_flag_returns_value() {
        let flags = ParsedCLIFlags::new(valid_args()).unwrap();
        assert_eq!(flags.get_string_value("--input-fx-file").unwrap(), "in.fasta");
    }

    #[test]
    fn test_get_string_value_absent_optional_flag_returns_error() {
        let flags = ParsedCLIFlags::new(valid_args()).unwrap();
        assert!(flags.get_string_value("--subset-file").is_err());
    }

    #[test]
    fn test_get_string_value_optional_flag_present_returns_value() {
        let mut args = valid_args();
        args.extend(["--subset-file".to_string(), "names.txt".to_string()]);
        let flags = ParsedCLIFlags::new(args).unwrap();
        assert_eq!(flags.get_string_value("--subset-file").unwrap(), "names.txt");
    }

    // get_int_value tests

    #[test]
    fn test_get_int_value_present_flag_returns_value() {
        let mut args = valid_args();
        args.extend(["--start".to_string(), "42".to_string()]);
        let flags = ParsedCLIFlags::new(args).unwrap();
        assert_eq!(flags.get_int_value("--start").unwrap(), 42);
    }

    #[test]
    fn test_get_int_value_absent_flag_returns_error() {
        let flags = ParsedCLIFlags::new(valid_args()).unwrap();
        assert!(flags.get_int_value("--start").is_err());
    }

    // get_bool_value tests

    #[test]
    fn test_get_bool_value_present_flag_returns_true() {
        let mut args = valid_args();
        args.push("--keep-record-order".to_string());
        let flags = ParsedCLIFlags::new(args).unwrap();
        assert_eq!(flags.get_bool_value("--keep-record-order").unwrap(), true);
    }

    #[test]
    fn test_new_bool_flag_with_value_returns_error() {
        let mut args = valid_args();
        args.extend(["--keep-record-order".to_string(), "false".to_string()]);
        assert!(ParsedCLIFlags::new(args).unwrap_err().contains("takes no value"));
    }

    #[test]
    fn test_get_bool_value_absent_flag_returns_error() {
        let flags = ParsedCLIFlags::new(valid_args()).unwrap();
        assert!(flags.get_bool_value("--keep-record-order").is_err());
    }

    // parse_name_file tests

    #[test]
    fn test_parse_name_file_nonexistent_file_returns_error() {
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path().to_str().unwrap().to_string();
        drop(tmp);
        assert!(parse_name_file(&path).is_err());
    }

    #[test]
    fn test_parse_name_file_empty_file_returns_empty_vec() {
        let tmp    = NamedTempFile::new().unwrap();
        let result = parse_name_file(tmp.path().to_str().unwrap()).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_parse_name_file_single_field_sets_zero_range() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "seq1").unwrap();
        let result = parse_name_file(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].name, "seq1");
        assert_eq!(result[0].start, 0);
        assert_eq!(result[0].end, 0);
    }

    #[test]
    fn test_parse_name_file_three_fields_parses_range() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "seq1 5 10").unwrap();
        let result = parse_name_file(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].name, "seq1");
        assert_eq!(result[0].start, 5);
        assert_eq!(result[0].end, 10);
    }

    #[test]
    fn test_parse_name_file_two_fields_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "seq1 5").unwrap();
        assert!(parse_name_file(tmp.path().to_str().unwrap())
            .unwrap_err().contains("Ambiguous start/end"));
    }

    #[test]
    fn test_parse_name_file_comment_lines_skipped() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "# this is a comment").unwrap();
        writeln!(tmp, "seq1").unwrap();
        let result = parse_name_file(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].name, "seq1");
    }

    #[test]
    fn test_parse_name_file_blank_lines_skipped() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "seq1").unwrap();
        writeln!(tmp, "").unwrap();
        writeln!(tmp, "seq2").unwrap();
        let result = parse_name_file(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_parse_name_file_whitespace_only_lines_skipped() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "seq1").unwrap();
        writeln!(tmp, "   ").unwrap();
        writeln!(tmp, "seq2").unwrap();
        let result = parse_name_file(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_parse_name_file_extra_fields_ignored() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "seq1 5 10 extra ignored").unwrap();
        let result = parse_name_file(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].start, 5);
        assert_eq!(result[0].end, 10);
    }

    #[test]
    fn test_parse_name_file_invalid_start_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "seq1 abc 10").unwrap();
        assert!(parse_name_file(tmp.path().to_str().unwrap())
            .unwrap_err().contains("Invalid start"));
    }

    #[test]
    fn test_parse_name_file_invalid_end_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "seq1 5 xyz").unwrap();
        assert!(parse_name_file(tmp.path().to_str().unwrap())
            .unwrap_err().contains("Invalid end"));
    }

    #[test]
    fn test_parse_name_file_zero_start_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "seq1 0 10").unwrap();
        assert!(parse_name_file(tmp.path().to_str().unwrap())
            .unwrap_err().contains("Start value must be at least 1"));
    }

    #[test]
    fn test_parse_name_file_zero_end_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "seq1 1 0").unwrap();
        assert!(parse_name_file(tmp.path().to_str().unwrap())
            .unwrap_err().contains("End value must be at least 1"));
    }

    #[test]
    fn test_parse_name_file_mixed_content() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "# header comment").unwrap();
        writeln!(tmp, "seq1").unwrap();
        writeln!(tmp, "").unwrap();
        writeln!(tmp, "seq2 3 8").unwrap();
        let result = parse_name_file(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].name, "seq1");
        assert_eq!(result[0].start, 0);
        assert_eq!(result[0].end, 0);
        assert_eq!(result[1].name, "seq2");
        assert_eq!(result[1].start, 3);
        assert_eq!(result[1].end, 8);
    }

    // split_subset_list tests

    fn make_name_only(name: &str) -> fastx::NameWithRange {
        fastx::NameWithRange { name: name.to_string(), start: 0, end: 0 }
    }

    fn make_with_range(name: &str, start: usize, end: usize) -> fastx::NameWithRange {
        fastx::NameWithRange { name: name.to_string(), start, end }
    }

    #[test]
    fn test_split_subset_list_empty_input() {
        let (names, ranges) = split_subset_list(vec![]);
        assert!(names.is_empty());
        assert!(ranges.is_empty());
    }

    #[test]
    fn test_split_subset_list_all_name_only() {
        let input = vec![make_name_only("seq1"), make_name_only("seq2")];
        let (names, ranges) = split_subset_list(input);
        assert_eq!(names, vec!["seq1", "seq2"]);
        assert!(ranges.is_empty());
    }

    #[test]
    fn test_split_subset_list_all_with_range() {
        let input = vec![make_with_range("seq1", 5, 10), make_with_range("seq2", 3, 8)];
        let (names, ranges) = split_subset_list(input);
        assert!(names.is_empty());
        assert_eq!(ranges.len(), 2);
    }

    #[test]
    fn test_split_subset_list_mixed() {
        let input = vec![
            make_name_only("seq1"),
            make_with_range("seq2", 5, 10),
            make_name_only("seq3"),
            make_with_range("seq4", 3, 8),
        ];
        let (names, ranges) = split_subset_list(input);
        assert_eq!(names, vec!["seq1", "seq3"]);
        assert_eq!(ranges.len(), 2);
    }

    #[test]
    fn test_split_subset_list_start_converted_to_base0() {
        let input = vec![make_with_range("seq1", 5, 10)];
        let (_, ranges) = split_subset_list(input);
        assert_eq!(ranges[0].start, 4);
    }

    #[test]
    fn test_split_subset_list_end_unchanged() {
        let input = vec![make_with_range("seq1", 5, 10)];
        let (_, ranges) = split_subset_list(input);
        assert_eq!(ranges[0].end, 10);
    }

    #[test]
    fn test_split_subset_list_name_preserved_in_ranges() {
        let input = vec![make_with_range("seq1", 5, 10)];
        let (_, ranges) = split_subset_list(input);
        assert_eq!(ranges[0].name, "seq1");
    }
}
