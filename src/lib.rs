//! FASTA and FATSQ record storage and manipulation.
//!
//! Classes to store in memeory and manupulate FASTA and FASTQ records.
//!

/// wrap it all in a module
pub mod fastx {

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

}

#[cfg(test)]
mod tests {
    use super::fastx::IndexedSequence;

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
}
