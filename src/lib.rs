//! FASTA and FATSQ record storage and manipulation.
//!
//! Classes to store in memeory and manupulate FASTA and FASTQ records.
//!


/// wrap it all in a module
pub mod fastx {
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    /// Holds the sequence and its base-0 index in the original FASTA file.
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

        /// Returns a slice of the sequence between `start` (inclusive) and `end` (exclusive).
        ///
        /// # Arguments
        ///
        /// * `start` - Base-0 start position (inclusive)
        /// * `end` - Base-0 end position (exclusive)
        /// 
        /// # Returns
        ///
        /// If the `start` <= `end` and both are less than sequence length, returns the
        /// subsequence. If `start` is smaller than sequence length but end is larger,
        /// returns the sequence up to the end. If `start` is larger than `end` or sequence
        /// length, returns an empty slice.
        ///
        pub fn subsequence(&self, start: usize, end: usize) -> &str {
            let real_start = start.min(end);
            &self.sequence_[real_start.min( self.sequence_.len() ) .. end.min( self.sequence_.len() )]
        }
    }

    /// Holds the sequence, its quality scores, and its base-0 index in the original FASTQ file.
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

        /// Returns the original (base-0) index of this sequence in the source FASTA file.
        pub fn get_index(&self) -> u32 { self.original_index_ }

        /// Returns a slice of the sequence and quality between `start` (inclusive) and `end` (exclusive).
        ///
        /// # Arguments
        ///
        /// * `start` - Base-0 start position (inclusive)
        /// * `end` - Base-0 end position (exclusive)
        /// 
        /// # Returns
        ///
        /// If the `start` <= `end` and both are less than sequence length, returns the
        /// subsequence and quality score substring, in that order. If `start` is smaller than sequence length but end is larger,
        /// returns the sequence up to the end. If `start` is larger than `end` or sequence
        /// length, returns an empty slice.
        ///
        pub fn subsequence(&self, start: usize, end: usize) -> (&str, &str) {
            let real_start = start.min(end);
            (
                &self.sequence_[real_start.min( self.sequence_.len() ) .. end.min( self.sequence_.len() )],
                &self.quality_scores_[real_start.min( self.quality_scores_.len() ) .. end.min( self.quality_scores_.len() )]
            )
        }
    }    

    /// Collection of FASTA records, indexed by name.
    pub struct FastaRecords {
        records_: HashMap<String, IndexedSequence>,
    }

    impl FastaRecords {
        /// Creates a new `FastaRecords` object from a FASTA file.
        fn new(fasta_path: &str) -> Result<FastaRecords, String> {
            let mut current_header = String::new();
            let mut current_sequence = String::new();
            let file = File::open(fasta_path)
                .map_err( |error| error.to_string() )?;
            for line in BufReader::new(file).lines() {
                let line = line.map_err( |error| error.to_string() )?;
                if line.starts_with('>') {
                    current_header = line;
                }
            }
            todo!()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::fastx::*;

    #[test]
    fn test_new_and_get_index() {
        let seq = IndexedSequence::new( 3, &"ACGT".to_string() );
        assert_eq!(seq.get_index(), 3);
    }

    #[test]
    fn test_slice_within_bounds() {
        let seq = IndexedSequence::new( 0, &"ACGTACGT".to_string() );
        assert_eq!(seq.subsequence(2, 5), "GTA");
    }

    #[test]
    fn test_slice_full_sequence() {
        let seq = IndexedSequence::new( 0, &"ACGT".to_string() );
        assert_eq!(seq.subsequence(0, 4), "ACGT");
    }

    #[test]
    fn test_slice_equal_start_and_end() {
        let seq = IndexedSequence::new( 0, &"ACGT".to_string() );
        assert_eq!(seq.subsequence(2, 2), "");
    }

    #[test]
    fn test_slice_end_beyond_sequence() {
        let seq = IndexedSequence::new( 0, &"ACGT".to_string() );
        assert_eq!(seq.subsequence(1, 10), "CGT");
    }

    #[test]
    fn test_slice_start_beyond_sequence() {
        let seq = IndexedSequence::new( 0, &"ACGT".to_string() );
        assert_eq!(seq.subsequence(10, 20), "");
    }

    #[test]
    fn test_slice_start_greater_than_end() {
        let seq = IndexedSequence::new( 0, &"ACGT".to_string() );
        assert_eq!(seq.subsequence(3, 1), "");
    }

    #[test]
    fn test_empty_sequence_zero_indices() {
        let seq = IndexedSequence::new( 0, &"".to_string() );
        assert_eq!(seq.subsequence(0, 0), "");
    }

    #[test]
    fn test_empty_sequence_nonzero_indices() {
        let seq = IndexedSequence::new( 0, &"".to_string() );
        assert_eq!(seq.subsequence(2, 5), "");
    }

    #[test]
    fn test_empty_sequence_start_greater_than_end() {
        let seq = IndexedSequence::new( 0, &"".to_string() );
        assert_eq!(seq.subsequence(5, 2), "");
    }

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
        assert_eq!(seq.subsequence(2, 5), ("GTA", "HH?"));
    }

    #[test]
    fn test_quality_subsequence_full_sequence() {
        let seq = IndexedSequenceWithQuality::new( 0, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        assert_eq!(seq.subsequence(0, 4), ("ACGT", "IIII"));
    }

    #[test]
    fn test_quality_subsequence_equal_start_and_end() {
        let seq = IndexedSequenceWithQuality::new( 0, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        assert_eq!(seq.subsequence(2, 2), ("", ""));
    }

    #[test]
    fn test_quality_subsequence_end_beyond_sequence() {
        let seq = IndexedSequenceWithQuality::new( 0, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        assert_eq!(seq.subsequence(1, 10), ("CGT", "III"));
    }

    #[test]
    fn test_quality_subsequence_start_beyond_sequence() {
        let seq = IndexedSequenceWithQuality::new( 0, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        assert_eq!(seq.subsequence(10, 20), ("", ""));
    }

    #[test]
    fn test_quality_subsequence_start_greater_than_end() {
        let seq = IndexedSequenceWithQuality::new( 0, &"IIII".to_string(), &"ACGT".to_string() ).unwrap();
        assert_eq!(seq.subsequence(3, 1), ("", ""));
    }

    #[test]
    fn test_quality_empty_sequence_zero_indices() {
        let seq = IndexedSequenceWithQuality::new( 0, &"".to_string(), &"".to_string() ).unwrap();
        assert_eq!(seq.subsequence(0, 0), ("", ""));
    }

    #[test]
    fn test_quality_empty_sequence_nonzero_indices() {
        let seq = IndexedSequenceWithQuality::new( 0, &"".to_string(), &"".to_string() ).unwrap();
        assert_eq!(seq.subsequence(2, 5), ("", ""));
    }

    #[test]
    fn test_quality_empty_sequence_start_greater_than_end() {
        let seq = IndexedSequenceWithQuality::new( 0, &"".to_string(), &"".to_string() ).unwrap();
        assert_eq!(seq.subsequence(5, 2), ("", ""));
    }
}
