use bitvec::bitvec;
use bitvec::order::Lsb0;
use bitvec::vec::BitVec;

/// A growable k-mer represented as a bit vector.
#[derive(Debug)]
struct Kmer {
    inner: BitVec<usize, Lsb0>,
}

impl Kmer {
    pub fn new(k: usize) -> Self {
        Self {
            inner: bitvec!(usize, Lsb0; 0; k * 2),
        }
    }

    pub fn size(&self) -> usize {
        self.inner.len()
    }
}
