#![cfg(feature = "bitvec")]

use std::iter::FusedIterator;

use bitvec::bitbox;
use bitvec::boxed::BitBox;
use bitvec::field::BitField;
use bitvec::order::Msb0;
use bitvec::slice::BitSlice;
use bitvec::view::BitView;

use crate::Base;

/// A fixed-size k-mer represented as a bit vector.
///
/// Stores 1 to `usize::MAX / 2` bases.
#[derive(Debug)]
pub struct Kmer {
    inner: BitBox<usize, Msb0>,
    /// Index of the first base, in bits.
    start: usize,
}

impl Kmer {
    /// Creates a new k-mer with the given size.
    pub fn new(k: usize) -> Self {
        Self {
            inner: bitbox![usize, Msb0; 0; k * 2],
            start: 0,
        }
    }

    pub fn from_bytes(bytes: &[u8]) -> Self {
        let mut inner = bitbox!(usize, Msb0; 0; bytes.len() * 8);
        for (byte, chunk) in bytes.iter().zip(inner.chunks_exact_mut(8)) {
            chunk.store(*byte);
        }
        Self { inner, start: 0 }
    }

    pub fn push(&mut self, base: Base) {
        self.inner
            .get_mut(self.start..self.start + 2)
            .expect("slice is too small")
            .copy_from_bitslice(
                &(base as usize).view_bits::<Msb0>()[(<usize>::BITS - 2) as usize..],
            );

        self.start = (self.start + 2) % self.inner.len();
    }

    pub fn size(&self) -> usize {
        self.inner.len() / 2
    }

    pub fn bases(&self) -> Bases<'_> {
        Bases {
            inner: self,
            num_read: 0,
        }
    }
}

/// An iterator over the bases in a k-mer.
#[derive(Debug)]
pub struct Bases<'a> {
    inner: &'a Kmer,
    /// Number of bases that have already been read.
    num_read: usize,
}

impl<'a> Iterator for Bases<'a> {
    type Item = Base;

    fn next(&mut self) -> Option<Self::Item> {
        if self.num_read >= self.inner.size() {
            return None;
        }

        let bit_pos = (self.num_read * 2 + self.inner.start) % self.inner.inner.len();

        let base = unsafe {
            Base::from_u8_unchecked(self.inner.inner.get(bit_pos..bit_pos + 2)?.load::<u8>())
        };
        self.num_read += 1;
        Some(base)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.inner.size() - self.num_read;
        (remaining, Some(remaining))
    }
}

impl<'a> FusedIterator for Bases<'a> {}

impl<'a> ExactSizeIterator for Bases<'a> {}

#[cfg(test)]
mod tests {
    use bitvec::field::BitField;

    use super::*;

    #[test]
    fn new() {
        let kmer = Kmer::new(47);
        dbg!(kmer.inner.as_raw_slice());
        assert!(kmer.inner.not_any());
        assert_eq!(kmer.inner.len(), 47 * 2);
    }

    #[test]
    fn push() {
        let mut kmer = Kmer::new(1);
        kmer.push(Base::A);
        kmer.push(Base::T);

        assert_eq!(kmer.inner.len(), 2);
        assert_eq!(
            unsafe {
                Base::from_u8_unchecked(kmer.inner.chunks_exact(2).next().unwrap().load::<u8>())
            },
            Base::T
        );
    }

    #[test]
    fn bases() {
        let mut kmer = Kmer::new(7);
        for _ in 0..3 {
            kmer.push(Base::A);
            kmer.push(Base::T);
        }
        kmer.push(Base::C);
        kmer.push(Base::G);

        let bases: Vec<Base> = kmer.bases().collect();
        assert_eq!(
            bases,
            vec![
                Base::T,
                Base::A,
                Base::T,
                Base::A,
                Base::T,
                Base::C,
                Base::G
            ]
        );
    }

    #[test]
    fn from_bytes() {
        let bytes = [0b00011011u8];
        let kmer = Kmer::from_bytes(&bytes);
        assert_eq!(kmer.inner.len(), 8);
        assert_eq!(
            kmer.bases().collect::<Vec<_>>(),
            vec![Base::C, Base::A, Base::T, Base::G]
        );
    }
}
