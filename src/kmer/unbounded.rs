use std::iter::FusedIterator;

use bitvec::bitbox;
use bitvec::boxed::BitBox;
use bitvec::field::BitField;
use bitvec::order::Lsb0;
use bitvec::view::BitView;

use crate::Base;

/// A fixed-size k-mer represented as a bit vector.
///
/// Stores 1 to `usize::MAX / 2` bases.
#[derive(Debug)]
pub struct Kmer {
    store: BitBox<usize, Lsb0>,
    /// Index of the first base, in bits.
    start: usize,
}

impl Kmer {
    /// Creates a new k-mer with the given size.
    pub fn new(k: usize) -> Self {
        Self {
            store: bitbox![usize, Lsb0; 0; k * 2],
            start: k * 2,
        }
    }

    pub fn from_bytes(bytes: &[u8]) -> Self {
        let mut inner = bitbox!(usize, Lsb0; 0; bytes.len() * 8);
        for (byte, chunk) in bytes.iter().rev().zip(inner.chunks_exact_mut(8)) {
            chunk.store(*byte);
        }
        Self {
            store: inner,
            start: bytes.len() * 8,
        }
    }

    pub fn push(&mut self, base: Base) {
        self.store
            .get_mut(self.start - 2..self.start)
            .expect("slice is too small")
            .copy_from_bitslice(&(base as usize).view_bits()[..2]);

        self.start = (self.start - 2) % self.store.len();
        self.start += (self.start == 0) as usize * self.store.len(); // wrap around at 0
    }

    pub fn size(&self) -> usize {
        self.store.len() / 2
    }

    pub fn bases(&self) -> Bases<'_> {
        Bases {
            kmer: self,
            num_read: 0,
        }
    }
}

/// An iterator over the bases in a k-mer.
#[derive(Debug)]
pub struct Bases<'a> {
    kmer: &'a Kmer,
    /// Number of bases that have already been read.
    num_read: usize,
}

impl<'a> Iterator for Bases<'a> {
    type Item = Base;

    fn next(&mut self) -> Option<Self::Item> {
        if self.num_read >= self.kmer.size() {
            return None;
        }

        let bit_pos = match (self.kmer.start as isize - self.num_read as isize * 2)
            .rem_euclid(self.kmer.store.len() as isize)
        {
            0 => self.kmer.store.len(),
            p => p as usize,
        };

        let base = unsafe {
            Base::from_u8_unchecked(self.kmer.store.get(bit_pos - 2..bit_pos)?.load::<u8>())
        };
        self.num_read += 1;
        Some(base)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.kmer.size() - self.num_read;
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
        dbg!(kmer.store.as_raw_slice());
        assert!(kmer.store.not_any());
        assert_eq!(kmer.store.len(), 47 * 2);
    }

    #[test]
    fn push() {
        let mut kmer = Kmer::new(1);
        kmer.push(Base::A);
        kmer.push(Base::T);

        assert_eq!(kmer.store.len(), 2);
        assert_eq!(
            unsafe {
                Base::from_u8_unchecked(kmer.store.chunks_exact(2).next().unwrap().load::<u8>())
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
        let bytes = [0x1B, 0xAA, 0xF0, 0x0F, 0xCC, 0xFF, 0x00, 0x3C, 0xCF];
        let kmer = Kmer::from_bytes(&bytes);
        assert_eq!(kmer.store.len(), 9 * <u8>::BITS as usize);
        assert_eq!(kmer.store.as_raw_slice(), &[0xAAF00FCCFF003CCF, 0x1B]);
        assert_eq!(
            kmer.bases().collect::<Vec<_>>()[..8],
            vec![
                Base::C,
                Base::A,
                Base::T,
                Base::G,
                Base::T,
                Base::T,
                Base::T,
                Base::T
            ][..]
        );
        assert_eq!(
            kmer.bases().collect::<Vec<_>>()[28..],
            vec![
                Base::C,
                Base::G,
                Base::G,
                Base::C,
                Base::G,
                Base::C,
                Base::G,
                Base::G
            ][..]
        );
    }
}
