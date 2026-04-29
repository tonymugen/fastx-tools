//! FASTA and FASTQ record storage and manipulation.
//!
//! Types to store in memory and manipulate FASTA and FASTQ records.
//!


/// wrap it all in a module
pub mod fastx {
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufReader, Write};

    /// Holds FASTA/FASTQ record name and sequence range (for subsetting).
    #[derive(Debug)]
    pub struct NameWithRange {
        pub name:  String,
        pub start: usize,    // inclusive
        pub end:   usize,    // exclusive
    }

    fn clamp_slice(s: &str, start: usize, end: usize) -> &str {
        let real_start = start.min(end).min( s.len() );
        let real_end   = end.min( s.len() );
        &s[real_start..real_end]
    }

    /// Common interface for sequence record types stored in `FastxRecords`.
    pub trait FastxRecord: Clone {
        /// Returns the base-0 index of this record in the source file.
        fn get_index(&self) -> u32;
        /// Returns the sequence string.
        fn get_sequence(&self) -> &str;
        /// Returns a new record containing the subsequence between `start` (inclusive) and `end` (exclusive).
        fn subsequence(&self, start: usize, end: usize) -> Self;
        /// Returns the record formatted for writing to a file.
        fn format_output(&self, name: &str) -> String;
    }

    /// Holds the sequence and its base-0 index in the original FASTA file.
    #[derive(Clone, Debug)]
    pub struct IndexedSequence {
        original_index_: u32,
        sequence_:       String,
    }

    impl IndexedSequence {
        /// Creates a new `IndexedSequence` with the given index and sequence.
        ///
        /// # Arguments
        ///
        /// * `original_index` - The position of this sequence in the original FASTA file
        /// * `sequence` - The nucleotide or amino acid sequence string
        pub fn new(original_index: u32, sequence: &str) -> IndexedSequence {
            IndexedSequence {
                original_index_: original_index,
                sequence_: sequence.to_string(),
            }
        }

    }

    impl FastxRecord for IndexedSequence {
        fn get_index(&self) -> u32 { self.original_index_ }

        fn get_sequence(&self) -> &str { &self.sequence_ }

        fn subsequence(&self, start: usize, end: usize) -> IndexedSequence {
            IndexedSequence::new(self.original_index_, clamp_slice(&self.sequence_, start, end))
        }

        fn format_output(&self, name: &str) -> String {
            format!("{name}\n{}\n", self.sequence_)
        }
    }

    /// Holds the sequence, its quality scores, and its base-0 index in the original FASTQ file.
    #[derive(Clone, Debug)]
    pub struct IndexedSequenceWithQuality {
        original_index_: u32,
        quality_scores_: String,
        sequence_:       String,
    }

    impl IndexedSequenceWithQuality {
        /// Creates a new `IndexedSequenceWithQuality` with the given index, quality scores, and sequence.
        ///
        /// # Arguments
        ///
        /// * `original_index` - The position of this sequence in the original FASTQ file
        /// * `quality_scores` - The quality scores for this sequence
        /// * `sequence` - The nucleotide or amino acid sequence string
        ///
        /// # Returns
        ///
        /// An `IndexedSequenceWithQuality` object or an error if the quality and sequence strings
        /// are unequal in length.
        pub fn new(original_index: u32, quality_scores: &str, sequence: &str) -> Result<IndexedSequenceWithQuality, String> {
            if quality_scores.len() != sequence.len() {
                return Err( "Quality scores must be the same length as sequence".to_string() )
            }
            Ok(
                IndexedSequenceWithQuality {
                    original_index_: original_index,
                    quality_scores_: quality_scores.to_string(),
                    sequence_:       sequence.to_string(),
                }
            )
        }

        /// Returns the quality scores string.
        pub fn get_quality_scores(&self) -> &str { &self.quality_scores_ }
    }

    impl FastxRecord for IndexedSequenceWithQuality {
        fn get_index(&self) -> u32 { self.original_index_ }

        fn get_sequence(&self) -> &str { &self.sequence_ }

        fn subsequence(&self, start: usize, end: usize) -> IndexedSequenceWithQuality {
            IndexedSequenceWithQuality::new(
                self.original_index_,
                clamp_slice(&self.quality_scores_, start, end),
                clamp_slice(&self.sequence_, start, end),
            ).unwrap()
        }

        fn format_output(&self, name: &str) -> String {
            format!("{name}\n{}\n+\n{}\n", self.sequence_, self.quality_scores_)
        }
    }

    /// Generic collection of sequence records, indexed by name.
    #[derive(Debug)]
    pub struct FastxRecords<T> {
        records_:             HashMap<String, T>,
        max_sequence_length_: usize,
    }

    impl<T: FastxRecord> FastxRecords<T> {
        /// Return record count.
        pub fn num_records(&self) -> usize { self.records_.len() }
        /// Return length of the longest sequence
        pub fn get_max_length(&self) -> usize { self.max_sequence_length_ }

        /// Return a subset of records from a list of names.
        ///
        /// # Arguments
        ///
        /// * `names` - A list of record names, not necessarily all present in the current object.
        ///
        /// # Returns
        ///
        /// A tuple containing a `FastxRecords` object with any records present in the input list
        /// and a list of names not found. An empty list returns an empty object. Record indexes
        /// still refer to the original file. Duplicate record names are removed in the
        /// `FastxRecords` object but not in the list of absent records.
        ///
        pub fn records_by_name(&self, names: Vec<String>) -> (FastxRecords<T>, String) {
            let mut subset: HashMap<String, T> = HashMap::new();
            let mut absent_names               = String::new();
            let mut current_max_length: usize  = 0;
            for name in names {
                if let Some(record) = self.records_.get(&name) {
                    subset.insert( name, record.clone() );
                    current_max_length = current_max_length.max( record.get_sequence().len() );
                } else {
                    if !absent_names.is_empty() { absent_names.push('\n'); }
                    absent_names.push_str(&name);
                }
            }
            (
                FastxRecords {
                    records_: subset,
                    max_sequence_length_: current_max_length
                },
                absent_names
            )
        }

        /// Return subsequences of all records.
        ///
        /// # Arguments
        ///
        /// * `start` - Base-0 start position (inclusive)
        /// * `end` - Base-0 end position (exclusive)
        ///
        /// # Returns
        ///
        /// A `FastxRecords` object with all records, but with the subsequence between `start` and `end`.
        ///
        pub fn subsequences(&self, start: usize, end: usize) -> FastxRecords<T> {
            let mut subset: HashMap<String, T> = HashMap::new();
            // still must track because the actual lengths may not be (end - start) long
            let mut current_max_length: usize = 0;
            for (name, record) in &self.records_ {
                let local_subsequence = record.subsequence(start, end);
                current_max_length    = current_max_length.max( local_subsequence.get_sequence().len() );
                subset.insert(name.clone(), local_subsequence);
            }
            FastxRecords { records_: subset, max_sequence_length_: current_max_length }
        }

        /// Return subsequences of named records.
        ///
        /// # Arguments
        ///
        /// * `names_ranges` - Record names and ranges.
        ///
        /// # Returns
        ///
        /// A tuple containing a `FastxRecords` object with subsaets of any records present in the input list
        /// and a list of names not found. An empty list returns an empty object. Record indexes
        /// still refer to the original file. Duplicate record names are removed in the
        /// `FastxRecords` object but not in the list of absent records.
        ///
        pub fn subsequences_by_name(&self, names_ranges: Vec<NameWithRange>) -> (FastxRecords<T>, String) {
            let mut subset: HashMap<String, T> = HashMap::new();
            let mut absent_names               = String::new();
            let mut current_max_length: usize  = 0;
            for name_range in names_ranges {
                if let Some(record) = self.records_.get(&name_range.name) {
                    let local_subsequence = record.subsequence(name_range.start, name_range.end);
                    current_max_length    = current_max_length.max( local_subsequence.get_sequence().len() );
                    subset.insert(name_range.name, local_subsequence);
                } else {
                    if !absent_names.is_empty() { absent_names.push('\n'); }
                    absent_names.push_str(&name_range.name);
                }
            }
            (FastxRecords { records_: subset, max_sequence_length_: current_max_length }, absent_names)
        }

        /// Merge in a `FastxRecords` object.
        ///
        /// Any records also present in the current object are overwritten by new values.
        ///
        /// # Arguments
        ///
        /// * `from` - an object to merge
        ///
        /// # Returns
        ///
        /// Any errors.
        ///
        pub fn merge(&mut self, from: FastxRecords<T>) {
            self.max_sequence_length_ = self.max_sequence_length_.max(from.max_sequence_length_);
            self.records_.extend(from.records_);
        }

        /// Save the records to a file.
        ///
        /// # Arguments
        ///
        /// * `output_path` - Path to a file
        ///
        /// # Returns
        ///
        /// Any file write errors.
        ///
        pub fn save_records(&self, output_path: &str) -> Result<(), String> {
            let mut file = File::create(output_path)
                .map_err( |error| error.to_string() )?;
            for (name, record) in &self.records_ {
                file.write_all( record.format_output(name).as_bytes() )
                    .map_err( |error| error.to_string() )?;
            }
            Ok( () )
        }

        /// Save the records in the order they appeared in the original file.
        ///
        /// # Arguments
        ///
        /// * `output_path` - Path to a file
        ///
        /// # Returns
        ///
        /// Any file write errors.
        ///
        pub fn save_sorted_records(&self, output_path: &str) -> Result<(), String> {
            let mut file = File::create(output_path)
                .map_err( |error| error.to_string() )?;
            let mut sorted: Vec<(&String, &T)> = self.records_.iter().collect();
            sorted.sort_by_key( |(_, record)| record.get_index() );
            for (name, record) in sorted {
                file.write_all( record.format_output(name).as_bytes() )
                    .map_err( |error| error.to_string() )?;
            }
            Ok( () )
        }
    }

    /// Reads a FASTA file and returns a collection of records.
    ///
    /// # Arguments
    ///
    /// * `fasta_path` - Path to a FASTA file
    ///
    /// # Returns
    ///
    /// A `FastxRecords<IndexedSequence>` object or errors if the input file is inaccessible or invalid.
    ///
    pub fn read_fasta(fasta_path: &str) -> Result<FastxRecords<IndexedSequence>, String> {
        let mut local_records: HashMap<String, IndexedSequence> = HashMap::new();
        let mut current_header                                  = String::new();
        let mut current_sequence                                = String::new();
        let mut current_max_length: usize                       = 0;
        let file = File::open(fasta_path)
            .map_err( |error| error.to_string() )?;
        let mut record_idx: u32 = 0;
        for line in BufReader::new(file).lines() {
            let line = line.map_err( |error| error.to_string() )?;
            if line.starts_with('>') {
                if !current_header.is_empty() && !current_sequence.is_empty() {
                    local_records.insert( current_header.clone(), IndexedSequence::new(record_idx, &current_sequence) );
                    current_max_length = current_max_length.max( current_sequence.len() );
                    record_idx += 1;
                }
                current_header = line;
                current_sequence.clear();
                continue;
            }
            current_sequence.push_str(&line);
        }
        if !current_header.is_empty() && !current_sequence.is_empty() {
            local_records.insert(
                current_header.clone(),
                IndexedSequence::new(record_idx, &current_sequence)
            );
            current_max_length = current_max_length.max( current_sequence.len() );
        }
        if local_records.is_empty() {
            return Err( format!("No valid FASTA records in {fasta_path} file") )
        }
        Ok( FastxRecords { records_: local_records, max_sequence_length_: current_max_length } )
    }

    /// Reads a FASTQ file and returns a collection of records.
    ///
    /// # Arguments
    ///
    /// * `fastq_path` - Path to a FASTQ file
    ///
    /// # Returns
    ///
    /// A `FastxRecords<IndexedSequenceWithQuality>` object or errors if the input file is inaccessible or invalid.
    ///
    pub fn read_fastq(fastq_path: &str) -> Result<FastxRecords<IndexedSequenceWithQuality>, String> {
        let mut local_records: HashMap<String, IndexedSequenceWithQuality> = HashMap::new();
        let file = File::open(fastq_path)
            .map_err( |error| error.to_string() )?;
        let mut record_idx: u32           = 0;
        let mut current_max_length: usize = 0;
        let mut lines                     = BufReader::new(file).lines();
        loop {
            let header = match lines.next() {
                None    => break,
                Some(l) => l.map_err( |error| error.to_string() )?,
            };
            if header.is_empty() { continue; }
            if !header.starts_with('@') {
                return Err( format!("Expected '@' header line, got: {header}") );
            }
            let sequence = lines.next()
                .ok_or_else( || format!("Missing sequence line after header: {header}") )?
                .map_err( |error| error.to_string() )?;
            let plus = lines.next()
                .ok_or_else( || format!("Missing '+' line after sequence in record: {header}") )?
                .map_err( |error| error.to_string() )?;
            if !plus.starts_with('+') {
                return Err( format!("Expected '+' separator line, got: {plus} in record: {header}") );
            }
            let quality = lines.next()
                .ok_or_else( || format!("Missing quality line after '+' in record: {header}") )?
                .map_err( |error| error.to_string() )?;
            current_max_length = current_max_length.max( sequence.len() );
            local_records.insert(
                header.clone(),
                IndexedSequenceWithQuality::new(record_idx, &quality, &sequence)
                    .map_err( |error| format!("{error} in record: {header}") )?
            );
            record_idx += 1;
        }
        if local_records.is_empty() {
            return Err( format!("No valid FASTQ records in {fastq_path} file") );
        }
        Ok( FastxRecords { records_: local_records, max_sequence_length_: current_max_length } )
    }
}

#[cfg(test)]
mod tests {
    use super::fastx::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_standard_fasta() -> NamedTempFile {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">seq1").unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        writeln!(tmp, ">seq2").unwrap();
        writeln!(tmp, "TTGGCCAA").unwrap();
        writeln!(tmp, ">seq3").unwrap();
        writeln!(tmp, "GCGCGCGC").unwrap();
        tmp
    }

    fn make_standard_fastq() -> NamedTempFile {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "@record1").unwrap();
        writeln!(tmp, "ACGTACGTACGT").unwrap();
        writeln!(tmp, "+").unwrap();
        writeln!(tmp, "IIIIIIIIIIII").unwrap();
        writeln!(tmp, "@record2").unwrap();
        writeln!(tmp, "TTGGCCAATTGG").unwrap();
        writeln!(tmp, "+").unwrap();
        writeln!(tmp, "HHHHHHHHHHHH").unwrap();
        writeln!(tmp, "@record3").unwrap();
        writeln!(tmp, "GCGCGCGCGCGC").unwrap();
        writeln!(tmp, "+").unwrap();
        writeln!(tmp, "????????????").unwrap();
        tmp
    }

    // FastxRecords tests

    #[test]
    fn test_fasta_records_new_loads_all_records() {
        let fasta   = make_standard_fasta();
        let records = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        assert_eq!(records.num_records(), 3);
    }

    #[test]
    fn test_fasta_records_by_name_returns_present_records() {
        let fasta   = make_standard_fasta();
        let records = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names   = vec![">seq1".to_string(), ">seq2".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 2);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fasta_records_by_name_single_record() {
        let fasta            = make_standard_fasta();
        let records          = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names            = vec![">seq1".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fasta_records_by_name_absent_names_reported() {
        let fasta            = make_standard_fasta();
        let records          = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names            = vec!["not_a_real_record".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 0);
        assert!(absent.contains("not_a_real_record"));
    }

    #[test]
    fn test_fasta_records_by_name_empty_input_returns_empty() {
        let fasta            = make_standard_fasta();
        let records          = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let (subset, absent) = records.records_by_name(vec![]);
        assert_eq!(subset.num_records(), 0);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fasta_records_by_name_duplicate_present_name_deduplicated() {
        let fasta   = make_standard_fasta();
        let records = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names   = vec![">seq1".to_string(), ">seq1".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fasta_records_by_name_duplicate_absent_name_repeated_in_output() {
        let fasta            = make_standard_fasta();
        let records          = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names            = vec!["missing".to_string(), "missing".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 0);
        assert_eq!(absent, "missing\nmissing");
    }

    #[test]
    fn test_fasta_records_new_empty_file_returns_error() {
        let tmp    = NamedTempFile::new().unwrap();
        let result = read_fasta(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("No valid FASTA records in"));
    }

    #[test]
    fn test_fasta_records_new_no_headers_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        writeln!(tmp, "TTGGCCAA").unwrap();
        let result = read_fasta(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("No valid FASTA records in"));
    }

    #[test]
    fn test_fasta_records_new_header_without_sequence_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">solo_header").unwrap();
        let result = read_fasta(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("No valid FASTA records in"));
    }

    #[test]
    fn test_fasta_records_new_nonexistent_file_returns_error() {
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path().to_str().unwrap().to_string();
        drop(tmp);
        let result = read_fasta(&path);
        assert!(result.unwrap_err().contains("No such file or directory"));
    }

    // multi-line FASTA sequence tests
    #[test]
    fn test_fasta_multiline_sequence_is_concatenated() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">record1").unwrap();
        writeln!(tmp, "ACGT").unwrap();
        writeln!(tmp, "GCGC").unwrap();
        writeln!(tmp, "TTTT").unwrap();
        let records = read_fasta(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(records.num_records(), 1);
        let out = NamedTempFile::new().unwrap();
        records.save_sorted_records(out.path().to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(out.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines[1], "ACGTGCGCTTTT");
    }

    #[test]
    fn test_fasta_multiline_sequence_multiple_records() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">record1").unwrap();
        writeln!(tmp, "ACGT").unwrap();
        writeln!(tmp, "GCGC").unwrap();
        writeln!(tmp, ">record2").unwrap();
        writeln!(tmp, "TTTT").unwrap();
        writeln!(tmp, "AAAA").unwrap();
        let records = read_fasta(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(records.num_records(), 2);
    }

    // FastxRecords::subsequences tests
    #[test]
    fn test_fasta_records_subsequences_preserves_record_count() {
        let fasta   = make_standard_fasta();
        let records = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let sub     = records.subsequences(0, 10);
        assert_eq!(sub.num_records(), records.num_records());
    }

    #[test]
    fn test_fasta_records_subsequences_within_bounds() {
        let fasta       = make_standard_fasta();
        let records     = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let (first, _)  = records.records_by_name(vec![">seq1".to_string()]);
        let sub         = first.subsequences(2, 5);
        let (result, _) = sub.records_by_name(vec![">seq1".to_string()]);
        assert_eq!(result.num_records(), 1);
    }

    #[test]
    fn test_fasta_records_subsequences_full_sequence_preserves_length() {
        let fasta   = make_standard_fasta();
        let records = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let sub     = records.subsequences(0, usize::MAX);
        assert_eq!(sub.num_records(), 3);
    }

    #[test]
    fn test_fasta_records_subsequences_start_equals_end_returns_empty_sequences() {
        let fasta   = make_standard_fasta();
        let records = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let sub     = records.subsequences(5, 5);
        assert_eq!(sub.num_records(), 3);
    }

    #[test]
    fn test_fasta_records_subsequences_start_beyond_end_returns_empty_sequences() {
        let fasta   = make_standard_fasta();
        let records = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let sub     = records.subsequences(10, 2);
        assert_eq!(sub.num_records(), 3);
    }

    #[test]
    fn test_fasta_records_subsequences_start_beyond_sequence_returns_empty_sequences() {
        let fasta   = make_standard_fasta();
        let records = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let sub     = records.subsequences(usize::MAX - 1, usize::MAX);
        assert_eq!(sub.num_records(), 3);
    }

    // FastxRecords::subsequences_by_name tests
    #[test]
    fn test_subsequences_by_name_returns_present_records() {
        let fasta        = make_standard_fasta();
        let records      = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">seq1".to_string(), start: 0, end: 5 },
            NameWithRange { name: ">seq2".to_string(), start: 0, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 2);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_subsequences_by_name_single_record() {
        let fasta        = make_standard_fasta();
        let records      = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">seq1".to_string(), start: 2, end: 7 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_subsequences_by_name_absent_names_reported() {
        let fasta        = make_standard_fasta();
        let records      = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "not_a_real_record".to_string(), start: 0, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 0);
        assert_eq!(absent, "not_a_real_record");
    }

    #[test]
    fn test_subsequences_by_name_multiple_absent_names_separated_by_newline() {
        let fasta        = make_standard_fasta();
        let records      = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "missing_one".to_string(), start: 0, end: 5 },
            NameWithRange { name: "missing_two".to_string(), start: 0, end: 5 },
        ];
        let (_, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(absent, "missing_one\nmissing_two");
    }

    #[test]
    fn test_subsequences_by_name_mixed_present_and_absent() {
        let fasta        = make_standard_fasta();
        let records      = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">seq1".to_string(), start: 0, end: 5 },
            NameWithRange { name: "missing_record".to_string(), start: 0, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert_eq!(absent, "missing_record");
    }

    #[test]
    fn test_subsequences_by_name_empty_input_returns_empty() {
        let fasta            = make_standard_fasta();
        let records          = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let (subset, absent) = records.subsequences_by_name(vec![]);
        assert_eq!(subset.num_records(), 0);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_subsequences_by_name_start_equals_end() {
        let fasta        = make_standard_fasta();
        let records      = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">seq1".to_string(), start: 5, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_subsequences_by_name_start_beyond_end() {
        let fasta        = make_standard_fasta();
        let records      = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">seq1".to_string(), start: 10, end: 2 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_subsequences_by_name_start_beyond_sequence() {
        let fasta        = make_standard_fasta();
        let records      = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">seq1".to_string(), start: usize::MAX - 1, end: usize::MAX },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    // FastxRecords::save_records and save_sorted_records tests
    #[test]
    fn test_fasta_save_records_roundtrip() {
        let fasta   = make_standard_fasta();
        let records = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let tmp     = NamedTempFile::new().unwrap();
        records.save_records(tmp.path().to_str().unwrap()).unwrap();
        let saved = read_fasta(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(saved.num_records(), records.num_records());
    }

    #[test]
    fn test_fasta_save_sorted_records_roundtrip() {
        let fasta   = make_standard_fasta();
        let records = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let tmp     = NamedTempFile::new().unwrap();
        records.save_sorted_records(tmp.path().to_str().unwrap()).unwrap();
        let saved = read_fasta(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(saved.num_records(), records.num_records());
    }

    #[test]
    fn test_fasta_save_sorted_records_order() {
        let mut input_tmp = NamedTempFile::new().unwrap();
        writeln!(input_tmp, ">first").unwrap();
        writeln!(input_tmp, "AAAA").unwrap();
        writeln!(input_tmp, ">second").unwrap();
        writeln!(input_tmp, "CCCC").unwrap();
        writeln!(input_tmp, ">third").unwrap();
        writeln!(input_tmp, "GGGG").unwrap();
        let records = read_fasta(input_tmp.path().to_str().unwrap()).unwrap();
        let output_tmp = NamedTempFile::new().unwrap();
        records.save_sorted_records(output_tmp.path().to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(output_tmp.path()).unwrap();
        let headers: Vec<&str> = content.lines().filter( |l| l.starts_with('>') ).collect();
        assert_eq!(headers, vec![">first", ">second", ">third"]);
    }

    // FastqRecords tests

    #[test]
    fn test_fastq_records_new_loads_all_records() {
        let fastq   = make_standard_fastq();
        let records = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        assert_eq!(records.num_records(), 3);
    }

    #[test]
    fn test_fastq_records_by_name_returns_present_records() {
        let fastq            = make_standard_fastq();
        let records          = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names            = vec!["@record1".to_string(), "@record2".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 2);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_records_by_name_single_record() {
        let fastq            = make_standard_fastq();
        let records          = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names            = vec!["@record3".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_records_by_name_absent_names_reported() {
        let fastq            = make_standard_fastq();
        let records          = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names            = vec!["not_a_real_record".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 0);
        assert!(absent.contains("not_a_real_record"));
    }

    #[test]
    fn test_fastq_records_by_name_empty_input_returns_empty() {
        let fastq            = make_standard_fastq();
        let records          = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let (subset, absent) = records.records_by_name(vec![]);
        assert_eq!(subset.num_records(), 0);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_records_by_name_duplicate_present_name_deduplicated() {
        let fastq            = make_standard_fastq();
        let records          = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names            = vec!["@record1".to_string(), "@record1".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_records_by_name_duplicate_absent_name_repeated_in_output() {
        let fastq            = make_standard_fastq();
        let records          = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names            = vec!["missing".to_string(), "missing".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 0);
        assert_eq!(absent, "missing\nmissing");
    }

    #[test]
    fn test_fastq_records_new_empty_file_returns_error() {
        let tmp    = NamedTempFile::new().unwrap();
        let result = read_fastq(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("No valid FASTQ records in"));
    }

    #[test]
    fn test_fastq_records_new_no_at_header_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "record1").unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        writeln!(tmp, "+").unwrap();
        writeln!(tmp, "IIIIIIII").unwrap();
        let result = read_fastq(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("Expected '@' header line"));
    }

    #[test]
    fn test_fastq_records_new_missing_sequence_line_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "@record1").unwrap();
        let result = read_fastq(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("Missing sequence line after header"));
    }

    #[test]
    fn test_fastq_records_new_truncated_record_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "@record1").unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        let result = read_fastq(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("Missing '+' line"));
    }

    #[test]
    fn test_fastq_records_new_blank_lines_between_records_are_skipped() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "@record1").unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        writeln!(tmp, "+").unwrap();
        writeln!(tmp, "IIIIIIII").unwrap();
        writeln!(tmp, "").unwrap();
        writeln!(tmp, "@record2").unwrap();
        writeln!(tmp, "GCGCGCGC").unwrap();
        writeln!(tmp, "+").unwrap();
        writeln!(tmp, "HHHHHHHH").unwrap();
        let records = read_fastq(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(records.num_records(), 2);
    }

    #[test]
    fn test_fastq_records_new_invalid_separator_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "@record1").unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        writeln!(tmp, "GCGCGCGC").unwrap();
        writeln!(tmp, "IIIIIIII").unwrap();
        let result = read_fastq(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("Expected '+' separator line"));
    }

    #[test]
    fn test_fastq_records_new_missing_quality_line_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "@record1").unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        writeln!(tmp, "+").unwrap();
        let result = read_fastq(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("Missing quality line after '+'"));
    }

    #[test]
    fn test_fastq_records_new_mismatched_quality_length_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "@record1").unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        writeln!(tmp, "+").unwrap();
        writeln!(tmp, "IIII").unwrap();
        let result = read_fastq(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("Quality scores must be the same length"));
    }

    #[test]
    fn test_fastq_records_new_nonexistent_file_returns_error() {
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path().to_str().unwrap().to_string();
        drop(tmp);
        let result = read_fastq(&path);
        assert!(result.unwrap_err().contains("No such file or directory"));
    }

    // FastqRecords::subsequences tests
    #[test]
    fn test_fastq_records_subsequences_preserves_record_count() {
        let fastq   = make_standard_fastq();
        let records = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let sub     = records.subsequences(0, 5);
        assert_eq!(sub.num_records(), records.num_records());
    }

    #[test]
    fn test_fastq_records_subsequences_within_bounds() {
        let fastq       = make_standard_fastq();
        let records     = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let (first, _)  = records.records_by_name(vec!["@record1".to_string()]);
        let sub         = first.subsequences(2, 5);
        let (result, _) = sub.records_by_name(vec!["@record1".to_string()]);
        assert_eq!(result.num_records(), 1);
    }

    #[test]
    fn test_fastq_records_subsequences_full_sequence_preserves_count() {
        let fastq   = make_standard_fastq();
        let records = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let sub     = records.subsequences(0, usize::MAX);
        assert_eq!(sub.num_records(), 3);
    }

    #[test]
    fn test_fastq_records_subsequences_start_equals_end_returns_empty_sequences() {
        let fastq   = make_standard_fastq();
        let records = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let sub     = records.subsequences(5, 5);
        assert_eq!(sub.num_records(), 3);
    }

    #[test]
    fn test_fastq_records_subsequences_start_beyond_end_returns_empty_sequences() {
        let fastq   = make_standard_fastq();
        let records = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let sub     = records.subsequences(10, 2);
        assert_eq!(sub.num_records(), 3);
    }

    #[test]
    fn test_fastq_records_subsequences_start_beyond_sequence_returns_empty_sequences() {
        let fastq   = make_standard_fastq();
        let records = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let sub     = records.subsequences(usize::MAX - 1, usize::MAX);
        assert_eq!(sub.num_records(), 3);
    }

    // FastqRecords::subsequences_by_name tests
    #[test]
    fn test_fastq_subsequences_by_name_returns_present_records() {
        let fastq        = make_standard_fastq();
        let records      = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "@record1".to_string(), start: 0, end: 5 },
            NameWithRange { name: "@record2".to_string(), start: 0, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 2);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_subsequences_by_name_single_record() {
        let fastq        = make_standard_fastq();
        let records      = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "@record3".to_string(), start: 2, end: 7 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_subsequences_by_name_absent_names_reported() {
        let fastq        = make_standard_fastq();
        let records      = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "not_a_real_record".to_string(), start: 0, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 0);
        assert_eq!(absent, "not_a_real_record");
    }

    #[test]
    fn test_fastq_subsequences_by_name_multiple_absent_names_separated_by_newline() {
        let fastq        = make_standard_fastq();
        let records      = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "missing_one".to_string(), start: 0, end: 5 },
            NameWithRange { name: "missing_two".to_string(), start: 0, end: 5 },
        ];
        let (_, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(absent, "missing_one\nmissing_two");
    }

    #[test]
    fn test_fastq_subsequences_by_name_mixed_present_and_absent() {
        let fastq        = make_standard_fastq();
        let records      = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "@record1".to_string(), start: 0, end: 5 },
            NameWithRange { name: "missing_record".to_string(), start: 0, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert_eq!(absent, "missing_record");
    }

    #[test]
    fn test_fastq_subsequences_by_name_empty_input_returns_empty() {
        let fastq            = make_standard_fastq();
        let records          = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let (subset, absent) = records.subsequences_by_name(vec![]);
        assert_eq!(subset.num_records(), 0);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_subsequences_by_name_start_equals_end() {
        let fastq        = make_standard_fastq();
        let records      = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "@record1".to_string(), start: 5, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_subsequences_by_name_start_beyond_end() {
        let fastq        = make_standard_fastq();
        let records      = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "@record1".to_string(), start: 10, end: 2 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_subsequences_by_name_start_beyond_sequence() {
        let fastq        = make_standard_fastq();
        let records      = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "@record1".to_string(), start: usize::MAX - 1, end: usize::MAX },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    // FastqRecords::save_records and save_sorted_records tests
    #[test]
    fn test_fastq_save_records_roundtrip() {
        let fastq   = make_standard_fastq();
        let records = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let tmp     = NamedTempFile::new().unwrap();
        records.save_records(tmp.path().to_str().unwrap()).unwrap();
        let saved = read_fastq(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(saved.num_records(), records.num_records());
    }

    #[test]
    fn test_fastq_save_sorted_records_roundtrip() {
        let fastq   = make_standard_fastq();
        let records = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let tmp     = NamedTempFile::new().unwrap();
        records.save_sorted_records(tmp.path().to_str().unwrap()).unwrap();
        let saved = read_fastq(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(saved.num_records(), records.num_records());
    }

    #[test]
    fn test_fastq_save_sorted_records_order() {
        let fastq   = make_standard_fastq();
        let records = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        let tmp     = NamedTempFile::new().unwrap();
        records.save_sorted_records(tmp.path().to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        let headers: Vec<&str> = content.lines().filter( |l| l.starts_with('@') ).collect();
        assert_eq!(headers, vec!["@record1", "@record2", "@record3"]);
    }

    // FastxRecords::merge tests
    #[test]
    fn test_merge_disjoint_sets_combines_all_records() {
        let mut tmp1 = NamedTempFile::new().unwrap();
        writeln!(tmp1, ">r1").unwrap();
        writeln!(tmp1, "AAAA").unwrap();
        let mut tmp2 = NamedTempFile::new().unwrap();
        writeln!(tmp2, ">r2").unwrap();
        writeln!(tmp2, "CCCC").unwrap();
        let mut a = read_fasta(tmp1.path().to_str().unwrap()).unwrap();
        let b     = read_fasta(tmp2.path().to_str().unwrap()).unwrap();
        a.merge(b);
        assert_eq!(a.num_records(), 2);
    }

    #[test]
    fn test_merge_duplicate_key_overwrites_existing_record() {
        let mut tmp1 = NamedTempFile::new().unwrap();
        writeln!(tmp1, ">r1").unwrap();
        writeln!(tmp1, "AAAA").unwrap();
        let mut tmp2 = NamedTempFile::new().unwrap();
        writeln!(tmp2, ">r1").unwrap();
        writeln!(tmp2, "CCCC").unwrap();
        let mut a = read_fasta(tmp1.path().to_str().unwrap()).unwrap();
        let b     = read_fasta(tmp2.path().to_str().unwrap()).unwrap();
        a.merge(b);
        assert_eq!(a.num_records(), 1);
        let out = NamedTempFile::new().unwrap();
        a.save_records(out.path().to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(out.path()).unwrap();
        assert!(content.contains("CCCC"));
        assert!(!content.contains("AAAA"));
    }

    #[test]
    fn test_merge_into_empty_collection() {
        let mut tmp1 = NamedTempFile::new().unwrap();
        writeln!(tmp1, ">r1").unwrap();
        writeln!(tmp1, "AAAA").unwrap();
        let mut tmp2 = NamedTempFile::new().unwrap();
        writeln!(tmp2, ">r2").unwrap();
        writeln!(tmp2, "CCCC").unwrap();
        let a        = read_fasta(tmp1.path().to_str().unwrap()).unwrap();
        let b        = read_fasta(tmp2.path().to_str().unwrap()).unwrap();
        let (mut empty, _) = a.records_by_name(vec![]);
        empty.merge(b);
        assert_eq!(empty.num_records(), 1);
    }

    #[test]
    fn test_merge_empty_collection_into_existing() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">r1").unwrap();
        writeln!(tmp, "AAAA").unwrap();
        let mut a        = read_fasta(tmp.path().to_str().unwrap()).unwrap();
        let (empty, _)   = a.records_by_name(vec![]);
        a.merge(empty);
        assert_eq!(a.num_records(), 1);
    }

    #[test]
    fn test_merge_updates_max_length_when_incoming_is_longer() {
        let mut tmp1 = NamedTempFile::new().unwrap();
        writeln!(tmp1, ">r1").unwrap();
        writeln!(tmp1, "AAAA").unwrap();
        let mut tmp2 = NamedTempFile::new().unwrap();
        writeln!(tmp2, ">r2").unwrap();
        writeln!(tmp2, "CCCCCCCC").unwrap();
        let mut a = read_fasta(tmp1.path().to_str().unwrap()).unwrap();
        let b     = read_fasta(tmp2.path().to_str().unwrap()).unwrap();
        a.merge(b);
        assert_eq!(a.get_max_length(), 8);
    }

    #[test]
    fn test_merge_max_length_unchanged_when_incoming_is_shorter() {
        let mut tmp1 = NamedTempFile::new().unwrap();
        writeln!(tmp1, ">r1").unwrap();
        writeln!(tmp1, "CCCCCCCC").unwrap();
        let mut tmp2 = NamedTempFile::new().unwrap();
        writeln!(tmp2, ">r2").unwrap();
        writeln!(tmp2, "AAAA").unwrap();
        let mut a = read_fasta(tmp1.path().to_str().unwrap()).unwrap();
        let b     = read_fasta(tmp2.path().to_str().unwrap()).unwrap();
        a.merge(b);
        assert_eq!(a.get_max_length(), 8);
    }

    // get_max_length tests
    #[test]
    fn test_fasta_get_max_length_after_read() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">r1").unwrap();
        writeln!(tmp, "AAAA").unwrap();
        writeln!(tmp, ">r2").unwrap();
        writeln!(tmp, "CCCCCCCC").unwrap();
        writeln!(tmp, ">r3").unwrap();
        writeln!(tmp, "GGG").unwrap();
        let records = read_fasta(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(records.get_max_length(), 8);
    }

    #[test]
    fn test_fastq_get_max_length_after_read() {
        let fastq = make_standard_fastq();
        let records = read_fastq(fastq.path().to_str().unwrap()).unwrap();
        assert_eq!(records.get_max_length(), 12);
    }

    #[test]
    fn test_fasta_get_max_length_after_records_by_name() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">r1").unwrap();
        writeln!(tmp, "AAAA").unwrap();
        writeln!(tmp, ">r2").unwrap();
        writeln!(tmp, "CCCCCCCC").unwrap();
        writeln!(tmp, ">r3").unwrap();
        writeln!(tmp, "GGG").unwrap();
        let records = read_fasta(tmp.path().to_str().unwrap()).unwrap();
        let (subset, _) = records.records_by_name(vec![">r1".to_string(), ">r3".to_string()]);
        assert_eq!(subset.get_max_length(), 4);
    }

    #[test]
    fn test_fasta_get_max_length_after_subsequences() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">r1").unwrap();
        writeln!(tmp, "AAAA").unwrap();
        writeln!(tmp, ">r2").unwrap();
        writeln!(tmp, "CCCCCCCC").unwrap();
        writeln!(tmp, ">r3").unwrap();
        writeln!(tmp, "GGG").unwrap();
        let records = read_fasta(tmp.path().to_str().unwrap()).unwrap();
        let sub     = records.subsequences(0, 5);
        // r1=AAAA(4), r2=CCCCC(5), r3=GGG(3) — clamped at sequence length
        assert_eq!(sub.get_max_length(), 5);
    }

    #[test]
    fn test_fasta_get_max_length_after_subsequences_by_name() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">r1").unwrap();
        writeln!(tmp, "AAAA").unwrap();
        writeln!(tmp, ">r2").unwrap();
        writeln!(tmp, "CCCCCCCC").unwrap();
        writeln!(tmp, ">r3").unwrap();
        writeln!(tmp, "GGG").unwrap();
        let records = read_fasta(tmp.path().to_str().unwrap()).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">r2".to_string(), start: 0, end: 5 },
            NameWithRange { name: ">r1".to_string(), start: 0, end: 3 },
        ];
        let (subset, _) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.get_max_length(), 5);
    }

    #[test]
    fn test_fasta_get_max_length_empty_subset_is_zero() {
        let fasta       = make_standard_fasta();
        let records     = read_fasta(fasta.path().to_str().unwrap()).unwrap();
        let (subset, _) = records.records_by_name(vec![]);
        assert_eq!(subset.get_max_length(), 0);
    }

    // IndexedSequence tests
    #[test]
    fn test_new_and_get_index() {
        let seq = IndexedSequence::new( 3, "ACGT" );
        assert_eq!(seq.get_index(), 3);
    }

    #[test]
    fn test_slice_within_bounds() {
        let seq = IndexedSequence::new( 0, "ACGTACGT" );
        assert_eq!(seq.subsequence(2, 5).get_sequence(), "GTA");
    }

    #[test]
    fn test_slice_full_sequence() {
        let seq = IndexedSequence::new( 0, "ACGT" );
        assert_eq!(seq.subsequence(0, 4).get_sequence(), "ACGT");
    }

    #[test]
    fn test_slice_equal_start_and_end() {
        let seq = IndexedSequence::new( 0, "ACGT" );
        assert_eq!(seq.subsequence(2, 2).get_sequence(), "");
    }

    #[test]
    fn test_slice_end_beyond_sequence() {
        let seq = IndexedSequence::new( 0, "ACGT" );
        assert_eq!(seq.subsequence(1, 10).get_sequence(), "CGT");
    }

    #[test]
    fn test_slice_start_beyond_sequence() {
        let seq = IndexedSequence::new( 0, "ACGT" );
        assert_eq!(seq.subsequence(10, 20).get_sequence(), "");
    }

    #[test]
    fn test_slice_start_greater_than_end() {
        let seq = IndexedSequence::new( 0, "ACGT" );
        assert_eq!(seq.subsequence(3, 1).get_sequence(), "");
    }

    #[test]
    fn test_empty_sequence_zero_indices() {
        let seq = IndexedSequence::new( 0, "" );
        assert_eq!(seq.subsequence(0, 0).get_sequence(), "");
    }

    #[test]
    fn test_empty_sequence_nonzero_indices() {
        let seq = IndexedSequence::new( 0, "" );
        assert_eq!(seq.subsequence(2, 5).get_sequence(), "");
    }

    #[test]
    fn test_empty_sequence_start_greater_than_end() {
        let seq = IndexedSequence::new( 0, "" );
        assert_eq!(seq.subsequence(5, 2).get_sequence(), "");
    }

    #[test]
    fn test_clone_indexed_sequence() {
        let original = IndexedSequence::new( 7, "ACGT" );
        let cloned   = original.clone();
        assert_eq!(cloned.get_index(), original.get_index());
        assert_eq!(cloned.get_sequence(), original.get_sequence());
    }

    #[test]
    fn test_clone_indexed_sequence_is_independent() {
        let original = IndexedSequence::new( 2, "ACGT" );
        let cloned   = original.clone();
        assert_eq!(cloned.get_index(), 2);
        assert_eq!(cloned.get_sequence(), "ACGT");
    }

    // IndexedSequenceWithQuality tests
    #[test]
    fn test_quality_new_and_get_index() {
        let seq = IndexedSequenceWithQuality::new( 5, "IIII", "ACGT" ).unwrap();
        assert_eq!(seq.get_index(), 5);
    }

    #[test]
    fn test_quality_new_mismatched_lengths_returns_error() {
        let result = IndexedSequenceWithQuality::new( 0, "II", "ACGT" );
        assert!(result.is_err());
    }

    #[test]
    fn test_quality_subsequence_within_bounds() {
        let seq = IndexedSequenceWithQuality::new( 0, "IIHH????", "ACGTACGT" ).unwrap();
        let sub = seq.subsequence(2, 5);
        assert_eq!(sub.get_sequence(), "GTA");
        assert_eq!(sub.get_quality_scores(), "HH?");
    }

    #[test]
    fn test_quality_subsequence_full_sequence() {
        let seq = IndexedSequenceWithQuality::new( 0, "IIII", "ACGT" ).unwrap();
        let sub = seq.subsequence(0, 4);
        assert_eq!(sub.get_sequence(), "ACGT");
        assert_eq!(sub.get_quality_scores(), "IIII");
    }

    #[test]
    fn test_quality_subsequence_equal_start_and_end() {
        let seq = IndexedSequenceWithQuality::new( 0, "IIII", "ACGT" ).unwrap();
        let sub = seq.subsequence(2, 2);
        assert_eq!(sub.get_sequence(), "");
        assert_eq!(sub.get_quality_scores(), "");
    }

    #[test]
    fn test_quality_subsequence_end_beyond_sequence() {
        let seq = IndexedSequenceWithQuality::new( 0, "IIII", "ACGT" ).unwrap();
        let sub = seq.subsequence(1, 10);
        assert_eq!(sub.get_sequence(), "CGT");
        assert_eq!(sub.get_quality_scores(), "III");
    }

    #[test]
    fn test_quality_subsequence_start_beyond_sequence() {
        let seq = IndexedSequenceWithQuality::new( 0, "IIII", "ACGT" ).unwrap();
        let sub = seq.subsequence(10, 20);
        assert_eq!(sub.get_sequence(), "");
        assert_eq!(sub.get_quality_scores(), "");
    }

    #[test]
    fn test_quality_subsequence_start_greater_than_end() {
        let seq = IndexedSequenceWithQuality::new( 0, "IIII", "ACGT" ).unwrap();
        let sub = seq.subsequence(3, 1);
        assert_eq!(sub.get_sequence(), "");
        assert_eq!(sub.get_quality_scores(), "");
    }

    #[test]
    fn test_quality_empty_sequence_zero_indices() {
        let seq = IndexedSequenceWithQuality::new( 0, "", "" ).unwrap();
        let sub = seq.subsequence(0, 0);
        assert_eq!(sub.get_sequence(), "");
        assert_eq!(sub.get_quality_scores(), "");
    }

    #[test]
    fn test_quality_empty_sequence_nonzero_indices() {
        let seq = IndexedSequenceWithQuality::new( 0, "", "" ).unwrap();
        let sub = seq.subsequence(2, 5);
        assert_eq!(sub.get_sequence(), "");
        assert_eq!(sub.get_quality_scores(), "");
    }

    #[test]
    fn test_quality_empty_sequence_start_greater_than_end() {
        let seq = IndexedSequenceWithQuality::new( 0, "", "" ).unwrap();
        let sub = seq.subsequence(5, 2);
        assert_eq!(sub.get_sequence(), "");
        assert_eq!(sub.get_quality_scores(), "");
    }

    #[test]
    fn test_clone_indexed_sequence_with_quality() {
        let original = IndexedSequenceWithQuality::new( 3, "IIII", "ACGT" ).unwrap();
        let cloned   = original.clone();
        assert_eq!(cloned.get_index(), original.get_index());
        assert_eq!(cloned.get_sequence(), original.get_sequence());
        assert_eq!(cloned.get_quality_scores(), original.get_quality_scores());
    }

    #[test]
    fn test_clone_indexed_sequence_with_quality_is_independent() {
        let original = IndexedSequenceWithQuality::new( 1, "IIII", "ACGT" ).unwrap();
        let cloned   = original.clone();
        assert_eq!(cloned.get_index(), 1);
        assert_eq!(cloned.get_sequence(), "ACGT");
        assert_eq!(cloned.get_quality_scores(), "IIII");
    }

}
