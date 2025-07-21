use std::iter::FusedIterator;

use bitvec::field::BitField as _;
use bitvec::order::Lsb0;
use bitvec::slice::ChunksExact;
use bitvec::store::BitStore;
use bitvec::vec::BitVec;

use crate::Base;

pub struct Sequence<B: BitStore> {
    store: BitVec<B, Lsb0>,
}

impl<B: BitStore> Sequence<B> {
    pub fn new() -> Self {
        Self {
            store: BitVec::new(),
        }
    }

    pub fn push(&mut self, base: Base) {
        let bits = match base {
            Base::C => (false, false),
            Base::A => (false, true),
            Base::T => (true, false),
            Base::G => (true, true),
        };
        self.store.push(bits.1);
        self.store.push(bits.0);
    }

    pub fn kmers<const K: usize>(&self) -> SmallKmerIter<'_, K, B> {
        let mut kmer = crate::small::Kmer::<K>::new();
        let mut bases = self.store.chunks_exact(2);
        for _ in 0..K - 1 {
            if let Some(chunk) = bases.next() {
                // SAFETY: 2 bit bases are always valid.
                kmer.push(unsafe { Base::from_u8_unchecked(chunk.load::<u8>()) });
            }
        }

        SmallKmerIter { bases, kmer }
    }
}

pub struct SmallKmerIter<'a, const K: usize, B: BitStore> {
    bases: ChunksExact<'a, B, Lsb0>,
    kmer: crate::small::Kmer<K>,
}

impl<'a, const K: usize, B: BitStore> Iterator for SmallKmerIter<'a, K, B> {
    type Item = crate::small::Kmer<K>;

    fn next(&mut self) -> Option<Self::Item> {
        self.kmer
            .push(unsafe { Base::from_u8_unchecked(self.bases.next()?.load::<u8>()) });
        Some(self.kmer)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.bases.size_hint()
    }
}

impl<'a, const K: usize, B: BitStore> FusedIterator for SmallKmerIter<'a, K, B> {}

impl<'a, const K: usize, B: BitStore> ExactSizeIterator for SmallKmerIter<'a, K, B> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_kmer() {
        let mut bases = Vec::new();
        let mut seq = Sequence::<usize>::new();
        for _ in 0..11 {
            seq.push(Base::A);
            bases.push(Base::A);
            seq.push(Base::T);
            bases.push(Base::T);
        }
        seq.push(Base::C);
        bases.push(Base::C);

        let kmers = seq.kmers();
        assert_eq!(kmers.len(), 1);
        let kmers: Vec<crate::small::Kmer<23>> = kmers.collect();
        assert_eq!(kmers.len(), 1);
        eprintln!("{}", kmers[0]);

        let iter_bases: Vec<Base> = kmers[0].bases().collect();
        assert_eq!(iter_bases, bases);
    }
}
