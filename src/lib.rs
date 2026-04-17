//! FASTA and FATSQ record storage and manipulation.
//!
//! Classes to store in memeory and manupulate FASTA and FASTQ records.
//!


/// wrap it all in a module
pub mod fastx {
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    /// Holds FASTA/FASTQ record name and sequence range (for subsetting).
    pub struct NameWithRange {
        pub name: String,
        // inclusive
        pub start: usize,
        // exclusive
        pub end: usize,
    }
    /// Holds the sequence and its base-0 index in the original FASTA file.
    #[derive(Clone, Debug)]
    pub struct IndexedSequence {
        original_index_: u32,
        sequence_: String,
    }

    impl IndexedSequence {
        /// Creates a new `IndexedSequence` with the given index and sequence.
        ///
        /// # Arguments
        ///
        /// * `original_index` - The position of this sequence in the original FASTA file
        /// * `sequence` - The nucleotide or amino acid sequence string
        pub fn new(original_index: u32, sequence: &String) -> IndexedSequence {
            IndexedSequence {
                original_index_: original_index,
                sequence_: sequence.clone(),
            }
        }

        /// Returns the original (base-0) index of this sequence in the source FASTA file.
        pub fn get_index(&self) -> u32 { self.original_index_ }

        /// Returns the sequence string.
        pub fn get_sequence(&self) -> &str { &self.sequence_ }

        /// Returns a new `IndexedSequence` with the subsequence between `start` (inclusive) and `end` (exclusive).
        ///
        /// # Arguments
        ///
        /// * `start` - Base-0 start position (inclusive)
        /// * `end` - Base-0 end position (exclusive)
        ///
        /// # Returns
        ///
        /// If the `start` <= `end` and both are less than sequence length, returns new
        /// `IndexedSequence` object with the subsequence. If `start` is smaller than sequence length but end is larger,
        /// returns the sequence up to the end. If `start` is larger than `end` or sequence
        /// length, returns an empty slice.
        ///
        pub fn subsequence(&self, start: usize, end: usize) -> IndexedSequence {
            let real_start = start.min(end);
            let seq = &self.sequence_[real_start.min( self.sequence_.len() ) .. end.min( self.sequence_.len() )];
            IndexedSequence::new(self.original_index_, &seq.to_string())
        }
    }

    /// Holds the sequence, its quality scores, and its base-0 index in the original FASTQ file.
    #[derive(Clone, Debug)]
    pub struct IndexedSequenceWithQuality {
        original_index_: u32,
        quality_scores_: String,
        sequence_: String,
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
        /// An `IndexedSequenceWithQuality` object or an error if the length of the quality and
        /// sequence string are unequal in length.
        pub fn new(original_index: u32, quality_scores: &String, sequence: &String) -> Result<IndexedSequenceWithQuality, String> {
            if quality_scores.len() != sequence.len() {
                return Err( "Quality scores must be the same length as sequence".to_string() ) 
            }
            Ok(
                IndexedSequenceWithQuality {
                    original_index_: original_index,
                    quality_scores_: quality_scores.clone(),
                    sequence_: sequence.clone(),
                }
            )
        }

        /// Returns the original (base-0) index of this sequence in the source FASTQ file.
        pub fn get_index(&self) -> u32 { self.original_index_ }

        /// Returns the sequence string.
        pub fn get_sequence(&self) -> &str { &self.sequence_ }

        /// Returns the quality scores string.
        pub fn get_quality_scores(&self) -> &str { &self.quality_scores_ }

        /// Returns a new `IndexedSequenceWithQuality` with the subsequence between `start` (inclusive) and `end` (exclusive).
        ///
        /// # Arguments
        ///
        /// * `start` - Base-0 start position (inclusive)
        /// * `end` - Base-0 end position (exclusive)
        ///
        /// # Returns
        ///
        /// If the `start` <= `end` and both are less than sequence length, returns a new
        /// `IndexedSequenceWithQuality` object with the subsequence. If `start` is smaller than sequence length but end is larger,
        /// returns the sequence up to the end. If `start` is larger than `end` or sequence
        /// length, returns an empty object.
        ///
        pub fn subsequence(&self, start: usize, end: usize) -> IndexedSequenceWithQuality {
            let real_start = start.min(end);
            let seq = &self.sequence_[real_start.min( self.sequence_.len() ) .. end.min( self.sequence_.len() )];
            let qual = &self.quality_scores_[real_start.min( self.quality_scores_.len() ) .. end.min( self.quality_scores_.len() )];
            IndexedSequenceWithQuality::new(self.original_index_, &qual.to_string(), &seq.to_string()).unwrap()
        }
    }    

    /// Collection of FASTA records, indexed by name.
    #[derive(Debug)]
    pub struct FastaRecords {
        records_: HashMap<String, IndexedSequence>,
    }

    impl FastaRecords {
        /// Creates a new `FastaRecords` object from a FASTA file.
        ///
        /// # Arguments
        ///
        /// * `fasta_path` - Path to a FASTA file
        /// 
        /// # Returns
        ///
        /// A `FastaRecords` object or errors if the input file is inaccessible or invalid.
        ///
        pub fn new(fasta_path: &str) -> Result<FastaRecords, String> {
            let mut local_records: HashMap<String, IndexedSequence> = HashMap::new();
            let mut current_header = String::new();
            let mut current_sequence = String::new();
            let file = File::open(fasta_path)
                .map_err( |error| error.to_string() )?;
            let mut record_idx: u32 = 0;
            for line in BufReader::new(file).lines() {
                let line = line.map_err( |error| error.to_string() )?;
                if line.starts_with('>') {
                    if !current_header.is_empty() && !current_sequence.is_empty() {
                        local_records.insert( current_header.clone(), IndexedSequence::new(record_idx, &current_sequence) );
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
            }
            if local_records.is_empty() {
                return Err( format!("No valid FASTA records in {fasta_path} file") )
            }
            Ok(FastaRecords { records_: local_records })
        }

        /// Return record count.
        pub fn num_records(&self) -> usize { self.records_.len() }

        /// Return a subset of records from a list of names.
        ///
        /// # Arguments
        ///
        /// * `names` - A list of record names, not necessarily all present in the current object.
        ///
        /// # Returns 
        ///
        /// A tuple containing a `FastaRecords` object with any records present in the input list
        /// and a list of names not found. An empty list returns an empty object. Record indexes
        /// still refer to the original FASTA file. Duplicate record names are removed.
        ///
        pub fn records_by_name(&self, names: Vec<String>) -> (FastaRecords, String) {
            let mut subset: HashMap<String, IndexedSequence> = HashMap::new();
            let mut absent_names = String::new();
            for name in names {
                if self.records_.contains_key(&name) {
                    subset.insert( name.clone(), self.records_.get(&name).unwrap().clone() );
                    continue;
                }
                if !absent_names.is_empty() { absent_names.push('\n'); }
                absent_names.push_str(&name);
            }
            return (FastaRecords{ records_: subset }, absent_names)
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
        /// A `FastaRecords` object with all records, but with the subsequence
        /// between `start` and `end`
        ///
        pub fn subsequences(&self, start: usize, end: usize) -> FastaRecords {
            let mut subset: HashMap<String, IndexedSequence> = HashMap::new();
            for (name, record) in &self.records_ {
                subset.insert( name.clone(), record.subsequence(start, end) );
            }
            FastaRecords{ records_: subset }
        }

        /// Return subsequences of named records.
        ///
        /// # Arguments
        ///
        /// * `names_ranges` FASTA record names and ranges.
        ///
        /// # Returns
        ///
        /// A `FastaRecords` object with subsequences and a list of names that were not found.
        ///
        pub fn subsequences_by_name(&self, names_ranges: Vec<NameWithRange>) -> (FastaRecords, String) {
            let mut subset: HashMap<String, IndexedSequence> = HashMap::new();
            let mut absent_names = String::new();
            for name_range in names_ranges {
                if self.records_.contains_key(&name_range.name) {
                    subset.insert(
                        name_range.name.clone(),
                            self.records_.get(&name_range.name)
                                .unwrap()
                                    .subsequence(name_range.start, name_range.end)
                    );
                    continue;
                }
                if !absent_names.is_empty() { absent_names.push('\n'); }
                absent_names.push_str(&name_range.name);
            }
            return (FastaRecords{ records_: subset }, absent_names)
        }
    }
    /// Collection of FASTQ records, indexed by name.
    #[derive(Debug)]
    pub struct FastqRecords {
        records_: HashMap<String, IndexedSequenceWithQuality>,
    }

    impl FastqRecords {
        /// Creates a new `FastqRecords` object from a FASTQ file.
        ///
        /// # Arguments
        ///
        /// * `fastq_path` - Path to a FASTQ file
        ///
        /// # Returns
        ///
        /// A `FastqRecords` object or errors if the input file is inaccessible or invalid.
        ///
        pub fn new(fastq_path: &str) -> Result<FastqRecords, String> {
            let mut local_records: HashMap<String, IndexedSequenceWithQuality> = HashMap::new();
            let file = File::open(fastq_path)
                .map_err( |error| error.to_string() )?;
            let mut record_idx: u32 = 0;
            let mut lines = BufReader::new(file).lines();
            loop {
                let header = match lines.next() {
                    None => break,
                    Some(l) => l.map_err( |error| error.to_string() )?,
                };
                if !header.starts_with('@') {
                    return Err( format!("Expected '@' header line, got: {header}") );
                }
                let sequence = lines.next()
                    .ok_or_else( || format!("Missing sequence line after header: {header}") )?
                    .map_err( |error| error.to_string() )?;
                let _plus = lines.next()
                    .ok_or_else( || format!("Missing '+' line after sequence in record: {header}") )?
                    .map_err( |error| error.to_string() )?;
                let quality = lines.next()
                    .ok_or_else( || format!("Missing quality line after '+' in record: {header}") )?
                    .map_err( |error| error.to_string() )?;
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
            Ok(FastqRecords { records_: local_records })
        }

        /// Return record count.
        pub fn num_records(&self) -> usize { self.records_.len() }

        /// Return a subset of records from a list of names.
        ///
        /// # Arguments
        ///
        /// * `names` - A list of record names, not necessarily all present in the current object.
        ///
        /// # Returns
        ///
        /// A tuple containing a `FastqRecords` object with any records present in the input list
        /// and a list of names not found. An empty list returns an empty object. Record indexes
        /// still refer to the original FASTQ file. Duplicate record names are removed.
        ///
        pub fn records_by_name(&self, names: Vec<String>) -> (FastqRecords, String) {
            let mut subset: HashMap<String, IndexedSequenceWithQuality> = HashMap::new();
            let mut absent_names = String::new();
            for name in names {
                if self.records_.contains_key(&name) {
                    subset.insert( name.clone(), self.records_.get(&name).unwrap().clone() );
                    continue;
                }
                if !absent_names.is_empty() { absent_names.push('\n'); }
                absent_names.push_str(&name);
            }
            (FastqRecords{ records_: subset }, absent_names)
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
        /// A `FastqRecords` object with all records, but with the subsequence
        /// between `start` and `end`.
        ///
        pub fn subsequences(&self, start: usize, end: usize) -> FastqRecords {
            let mut subset: HashMap<String, IndexedSequenceWithQuality> = HashMap::new();
            for (name, record) in &self.records_ {
                subset.insert( name.clone(), record.subsequence(start, end) );
            }
            FastqRecords{ records_: subset }
        }

        /// Return subsequences of named records.
        ///
        /// # Arguments
        ///
        /// * `names_ranges` - FASTQ record names and ranges.
        ///
        /// # Returns
        ///
        /// A `FastqRecords` object with subsequences and a list of names that were not found.
        ///
        pub fn subsequences_by_name(&self, names_ranges: Vec<NameWithRange>) -> (FastqRecords, String) {
            let mut subset: HashMap<String, IndexedSequenceWithQuality> = HashMap::new();
            let mut absent_names = String::new();
            for name_range in names_ranges {
                if self.records_.contains_key(&name_range.name) {
                    subset.insert(
                        name_range.name.clone(),
                        self.records_.get(&name_range.name)
                            .unwrap()
                            .subsequence(name_range.start, name_range.end)
                    );
                    continue;
                }
                if !absent_names.is_empty() { absent_names.push('\n'); }
                absent_names.push_str(&name_range.name);
            }
            (FastqRecords{ records_: subset }, absent_names)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::fastx::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    // FastaRecords tests
    const CORRECT_FASTA: &str = "tests/data/correct.fasta";

    #[test]
    fn test_fasta_records_new_loads_all_records() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        assert_eq!(records.num_records(), 18);
    }

    #[test]
    fn test_fasta_records_by_name_returns_present_records() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let names = vec![
            ">B.FR.1983.IIIB_LAI.A04321".to_string(),
            ">C.IN.1993.93IN101.AB023804".to_string(),
        ];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 2);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fasta_records_by_name_single_record() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let names = vec![">B.US.1997.ARES2.AB078005".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fasta_records_by_name_absent_names_reported() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let names = vec!["not_a_real_record".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 0);
        assert!(absent.contains("not_a_real_record"));
    }

    #[test]
    fn test_fasta_records_by_name_empty_input_returns_empty() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let (subset, absent) = records.records_by_name(vec![]);
        assert_eq!(subset.num_records(), 0);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fasta_records_new_empty_file_returns_error() {
        let tmp = NamedTempFile::new().unwrap();
        let result = FastaRecords::new(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("No valid FASTA records in"));
    }

    #[test]
    fn test_fasta_records_new_no_headers_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        writeln!(tmp, "TTGGCCAA").unwrap();
        let result = FastaRecords::new(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("No valid FASTA records in"));
    }

    #[test]
    fn test_fasta_records_new_header_without_sequence_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, ">solo_header").unwrap();
        let result = FastaRecords::new(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("No valid FASTA records in"));
    }

    #[test]
    fn test_fasta_records_new_nonexistent_file_returns_error() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path().to_str().unwrap().to_string();
        drop(tmp);
        let result = FastaRecords::new(&path);
        assert!(result.unwrap_err().contains("No such file or directory"));
    }

    // FastaRecords::subsequences tests
    #[test]
    fn test_fasta_records_subsequences_preserves_record_count() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let sub = records.subsequences(0, 10);
        assert_eq!(sub.num_records(), records.num_records());
    }

    #[test]
    fn test_fasta_records_subsequences_within_bounds() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let (first, _) = records.records_by_name(vec![">B.FR.1983.IIIB_LAI.A04321".to_string()]);
        let sub = first.subsequences(2, 5);
        let (result, _) = sub.records_by_name(vec![">B.FR.1983.IIIB_LAI.A04321".to_string()]);
        assert_eq!(result.num_records(), 1);
    }

    #[test]
    fn test_fasta_records_subsequences_full_sequence_preserves_length() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let sub = records.subsequences(0, usize::MAX);
        assert_eq!(sub.num_records(), 18);
    }

    #[test]
    fn test_fasta_records_subsequences_start_equals_end_returns_empty_sequences() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let sub = records.subsequences(5, 5);
        assert_eq!(sub.num_records(), 18);
    }

    #[test]
    fn test_fasta_records_subsequences_start_beyond_end_returns_empty_sequences() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let sub = records.subsequences(10, 2);
        assert_eq!(sub.num_records(), 18);
    }

    #[test]
    fn test_fasta_records_subsequences_start_beyond_sequence_returns_empty_sequences() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let sub = records.subsequences(usize::MAX - 1, usize::MAX);
        assert_eq!(sub.num_records(), 18);
    }

    // FastaRecords::subsequences_by_name tests
    #[test]
    fn test_subsequences_by_name_returns_present_records() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">B.FR.1983.IIIB_LAI.A04321".to_string(), start: 0, end: 10 },
            NameWithRange { name: ">C.IN.1993.93IN101.AB023804".to_string(), start: 0, end: 10 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 2);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_subsequences_by_name_single_record() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">B.US.1997.ARES2.AB078005".to_string(), start: 2, end: 7 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_subsequences_by_name_absent_names_reported() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "not_a_real_record".to_string(), start: 0, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 0);
        assert_eq!(absent, "not_a_real_record");
    }

    #[test]
    fn test_subsequences_by_name_multiple_absent_names_separated_by_newline() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "missing_one".to_string(), start: 0, end: 5 },
            NameWithRange { name: "missing_two".to_string(), start: 0, end: 5 },
        ];
        let (_, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(absent, "missing_one\nmissing_two");
    }

    #[test]
    fn test_subsequences_by_name_mixed_present_and_absent() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">B.FR.1983.IIIB_LAI.A04321".to_string(), start: 0, end: 5 },
            NameWithRange { name: "missing_record".to_string(), start: 0, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert_eq!(absent, "missing_record");
    }

    #[test]
    fn test_subsequences_by_name_empty_input_returns_empty() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let (subset, absent) = records.subsequences_by_name(vec![]);
        assert_eq!(subset.num_records(), 0);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_subsequences_by_name_start_equals_end() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">B.FR.1983.IIIB_LAI.A04321".to_string(), start: 5, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_subsequences_by_name_start_beyond_end() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">B.FR.1983.IIIB_LAI.A04321".to_string(), start: 10, end: 2 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_subsequences_by_name_start_beyond_sequence() {
        let records = FastaRecords::new(CORRECT_FASTA).unwrap();
        let names_ranges = vec![
            NameWithRange { name: ">B.FR.1983.IIIB_LAI.A04321".to_string(), start: usize::MAX - 1, end: usize::MAX },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    // FastqRecords tests
    const CORRECT_FASTQ: &str = "tests/data/correct.fastq";

    #[test]
    fn test_fastq_records_new_loads_all_records() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        assert_eq!(records.num_records(), 3);
    }

    #[test]
    fn test_fastq_records_by_name_returns_present_records() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let names = vec!["@record1".to_string(), "@record2".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 2);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_records_by_name_single_record() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let names = vec!["@record3".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_records_by_name_absent_names_reported() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let names = vec!["not_a_real_record".to_string()];
        let (subset, absent) = records.records_by_name(names);
        assert_eq!(subset.num_records(), 0);
        assert!(absent.contains("not_a_real_record"));
    }

    #[test]
    fn test_fastq_records_by_name_empty_input_returns_empty() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let (subset, absent) = records.records_by_name(vec![]);
        assert_eq!(subset.num_records(), 0);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_records_new_empty_file_returns_error() {
        let tmp = NamedTempFile::new().unwrap();
        let result = FastqRecords::new(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("No valid FASTQ records in"));
    }

    #[test]
    fn test_fastq_records_new_no_at_header_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "record1").unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        writeln!(tmp, "+").unwrap();
        writeln!(tmp, "IIIIIIII").unwrap();
        let result = FastqRecords::new(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("Expected '@' header line"));
    }

    #[test]
    fn test_fastq_records_new_truncated_record_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "@record1").unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        let result = FastqRecords::new(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("Missing '+' line"));
    }

    #[test]
    fn test_fastq_records_new_mismatched_quality_length_returns_error() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "@record1").unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        writeln!(tmp, "+").unwrap();
        writeln!(tmp, "IIII").unwrap();
        let result = FastqRecords::new(tmp.path().to_str().unwrap());
        assert!(result.unwrap_err().contains("Quality scores must be the same length"));
    }

    #[test]
    fn test_fastq_records_new_nonexistent_file_returns_error() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path().to_str().unwrap().to_string();
        drop(tmp);
        let result = FastqRecords::new(&path);
        assert!(result.unwrap_err().contains("No such file or directory"));
    }

    // FastqRecords::subsequences tests
    #[test]
    fn test_fastq_records_subsequences_preserves_record_count() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let sub = records.subsequences(0, 5);
        assert_eq!(sub.num_records(), records.num_records());
    }

    #[test]
    fn test_fastq_records_subsequences_within_bounds() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let (first, _) = records.records_by_name(vec!["@record1".to_string()]);
        let sub = first.subsequences(2, 5);
        let (result, _) = sub.records_by_name(vec!["@record1".to_string()]);
        assert_eq!(result.num_records(), 1);
    }

    #[test]
    fn test_fastq_records_subsequences_full_sequence_preserves_count() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let sub = records.subsequences(0, usize::MAX);
        assert_eq!(sub.num_records(), 3);
    }

    #[test]
    fn test_fastq_records_subsequences_start_equals_end_returns_empty_sequences() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let sub = records.subsequences(5, 5);
        assert_eq!(sub.num_records(), 3);
    }

    #[test]
    fn test_fastq_records_subsequences_start_beyond_end_returns_empty_sequences() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let sub = records.subsequences(10, 2);
        assert_eq!(sub.num_records(), 3);
    }

    #[test]
    fn test_fastq_records_subsequences_start_beyond_sequence_returns_empty_sequences() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let sub = records.subsequences(usize::MAX - 1, usize::MAX);
        assert_eq!(sub.num_records(), 3);
    }

    // FastqRecords::subsequences_by_name tests
    #[test]
    fn test_fastq_subsequences_by_name_returns_present_records() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
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
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "@record3".to_string(), start: 2, end: 7 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_subsequences_by_name_absent_names_reported() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "not_a_real_record".to_string(), start: 0, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 0);
        assert_eq!(absent, "not_a_real_record");
    }

    #[test]
    fn test_fastq_subsequences_by_name_multiple_absent_names_separated_by_newline() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "missing_one".to_string(), start: 0, end: 5 },
            NameWithRange { name: "missing_two".to_string(), start: 0, end: 5 },
        ];
        let (_, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(absent, "missing_one\nmissing_two");
    }

    #[test]
    fn test_fastq_subsequences_by_name_mixed_present_and_absent() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
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
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let (subset, absent) = records.subsequences_by_name(vec![]);
        assert_eq!(subset.num_records(), 0);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_subsequences_by_name_start_equals_end() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "@record1".to_string(), start: 5, end: 5 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_subsequences_by_name_start_beyond_end() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "@record1".to_string(), start: 10, end: 2 },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    #[test]
    fn test_fastq_subsequences_by_name_start_beyond_sequence() {
        let records = FastqRecords::new(CORRECT_FASTQ).unwrap();
        let names_ranges = vec![
            NameWithRange { name: "@record1".to_string(), start: usize::MAX - 1, end: usize::MAX },
        ];
        let (subset, absent) = records.subsequences_by_name(names_ranges);
        assert_eq!(subset.num_records(), 1);
        assert!(absent.is_empty());
    }

    // IndexedSequence tests
    #[test]
    fn test_new_and_get_index() {
        let seq = IndexedSequence::new( 3, &"ACGT".to_string() );
        assert_eq!(seq.get_index(), 3);
    }

    #[test]
    fn test_slice_within_bounds() {
        let seq = IndexedSequence::new( 0, &"ACGTACGT".to_string() );
        assert_eq!(seq.subsequence(2, 5).get_sequence(), "GTA");
    }

    #[test]
    fn test_slice_full_sequence() {
        let seq = IndexedSequence::new( 0, &"ACGT".to_string() );
        assert_eq!(seq.subsequence(0, 4).get_sequence(), "ACGT");
    }

    #[test]
    fn test_slice_equal_start_and_end() {
        let seq = IndexedSequence::new( 0, &"ACGT".to_string() );
        assert_eq!(seq.subsequence(2, 2).get_sequence(), "");
    }

    #[test]
    fn test_slice_end_beyond_sequence() {
        let seq = IndexedSequence::new( 0, &"ACGT".to_string() );
        assert_eq!(seq.subsequence(1, 10).get_sequence(), "CGT");
    }

    #[test]
    fn test_slice_start_beyond_sequence() {
        let seq = IndexedSequence::new( 0, &"ACGT".to_string() );
        assert_eq!(seq.subsequence(10, 20).get_sequence(), "");
    }

    #[test]
    fn test_slice_start_greater_than_end() {
        let seq = IndexedSequence::new( 0, &"ACGT".to_string() );
        assert_eq!(seq.subsequence(3, 1).get_sequence(), "");
    }

    #[test]
    fn test_empty_sequence_zero_indices() {
        let seq = IndexedSequence::new( 0, &"".to_string() );
        assert_eq!(seq.subsequence(0, 0).get_sequence(), "");
    }

    #[test]
    fn test_empty_sequence_nonzero_indices() {
        let seq = IndexedSequence::new( 0, &"".to_string() );
        assert_eq!(seq.subsequence(2, 5).get_sequence(), "");
    }

    #[test]
    fn test_empty_sequence_start_greater_than_end() {
        let seq = IndexedSequence::new( 0, &"".to_string() );
        assert_eq!(seq.subsequence(5, 2).get_sequence(), "");
    }

    #[test]
    fn test_clone_indexed_sequence() {
        let original = IndexedSequence::new( 7, &"ACGT".to_string() );
        let cloned = original.clone();
        assert_eq!(cloned.get_index(), original.get_index());
        assert_eq!(cloned.get_sequence(), original.get_sequence());
    }

    #[test]
    fn test_clone_indexed_sequence_is_independent() {
        let original = IndexedSequence::new( 2, &"ACGT".to_string() );
        let cloned = original.clone();
        assert_eq!(cloned.get_index(), 2);
        assert_eq!(cloned.get_sequence(), "ACGT");
    }

    // IndexedSequenceWithQuality tests
    #[test]
    fn test_quality_new_and_get_index() {
        let seq = IndexedSequenceWithQuality::new( 5, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        assert_eq!(seq.get_index(), 5);
    }

    #[test]
    fn test_quality_new_mismatched_lengths_returns_error() {
        let result = IndexedSequenceWithQuality::new( 0, &"II".to_string(), &"ACGT".to_string() );
        assert!(result.is_err());
    }

    #[test]
    fn test_quality_subsequence_within_bounds() {
        let seq = IndexedSequenceWithQuality::new( 0, &"IIHH????".to_string(), &"ACGTACGT".to_string() ).unwrap();
        let sub = seq.subsequence(2, 5);
        assert_eq!(sub.get_sequence(), "GTA");
        assert_eq!(sub.get_quality_scores(), "HH?");
    }

    #[test]
    fn test_quality_subsequence_full_sequence() {
        let seq = IndexedSequenceWithQuality::new( 0, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        let sub = seq.subsequence(0, 4);
        assert_eq!(sub.get_sequence(), "ACGT");
        assert_eq!(sub.get_quality_scores(), "IIII");
    }

    #[test]
    fn test_quality_subsequence_equal_start_and_end() {
        let seq = IndexedSequenceWithQuality::new( 0, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        let sub = seq.subsequence(2, 2);
        assert_eq!(sub.get_sequence(), "");
        assert_eq!(sub.get_quality_scores(), "");
    }

    #[test]
    fn test_quality_subsequence_end_beyond_sequence() {
        let seq = IndexedSequenceWithQuality::new( 0, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        let sub = seq.subsequence(1, 10);
        assert_eq!(sub.get_sequence(), "CGT");
        assert_eq!(sub.get_quality_scores(), "III");
    }

    #[test]
    fn test_quality_subsequence_start_beyond_sequence() {
        let seq = IndexedSequenceWithQuality::new( 0, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        let sub = seq.subsequence(10, 20);
        assert_eq!(sub.get_sequence(), "");
        assert_eq!(sub.get_quality_scores(), "");
    }

    #[test]
    fn test_quality_subsequence_start_greater_than_end() {
        let seq = IndexedSequenceWithQuality::new( 0, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        let sub = seq.subsequence(3, 1);
        assert_eq!(sub.get_sequence(), "");
        assert_eq!(sub.get_quality_scores(), "");
    }

    #[test]
    fn test_quality_empty_sequence_zero_indices() {
        let seq = IndexedSequenceWithQuality::new( 0, &"".to_string(), &"".to_string() ).unwrap();
        let sub = seq.subsequence(0, 0);
        assert_eq!(sub.get_sequence(), "");
        assert_eq!(sub.get_quality_scores(), "");
    }

    #[test]
    fn test_quality_empty_sequence_nonzero_indices() {
        let seq = IndexedSequenceWithQuality::new( 0, &"".to_string(), &"".to_string() ).unwrap();
        let sub = seq.subsequence(2, 5);
        assert_eq!(sub.get_sequence(), "");
        assert_eq!(sub.get_quality_scores(), "");
    }

    #[test]
    fn test_quality_empty_sequence_start_greater_than_end() {
        let seq = IndexedSequenceWithQuality::new( 0, &"".to_string(), &"".to_string() ).unwrap();
        let sub = seq.subsequence(5, 2);
        assert_eq!(sub.get_sequence(), "");
        assert_eq!(sub.get_quality_scores(), "");
    }

    #[test]
    fn test_clone_indexed_sequence_with_quality() {
        let original = IndexedSequenceWithQuality::new( 3, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        let cloned = original.clone();
        assert_eq!(cloned.get_index(), original.get_index());
        assert_eq!(cloned.get_sequence(), original.get_sequence());
        assert_eq!(cloned.get_quality_scores(), original.get_quality_scores());
    }

    #[test]
    fn test_clone_indexed_sequence_with_quality_is_independent() {
        let original = IndexedSequenceWithQuality::new( 1, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        let cloned = original.clone();
        assert_eq!(cloned.get_index(), 1);
        assert_eq!(cloned.get_sequence(), "ACGT");
        assert_eq!(cloned.get_quality_scores(), "IIII");
    }

}
