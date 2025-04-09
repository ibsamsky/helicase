#![cfg(feature = "bitvec")]

use bitvec::field::BitField;
use bitvec::order::Msb0;
use bitvec::vec::BitVec;
use bitvec::view::BitView;

use crate::Base;

pub struct Sequence {
    bases: BitVec<usize, Msb0>,
}

impl Sequence {
    pub fn new() -> Self {
        Self {
            bases: BitVec::new(),
        }
    }

    pub fn push(&mut self, base: Base) {
        self.bases
            .extend_from_bitslice(&(base as usize).view_bits::<Msb0>()[usize::BITS as usize - 2..]);
    }

    pub fn kmers<const K: usize>(&self) -> SmallKmerIter<'_, K> {
        let kmer = crate::small::Kmer::<K>::new();
        self.bases.chunks_exact(2).take(K - 1).for_each(|base| {
            // SAFETY: 2-bit chunks are always in the range `0..4`.
            kmer.push(unsafe { Base::from_u8_unchecked(base.load::<u8>()) });
        });

        SmallKmerIter {
            seq: self,
            kmer,
            pos: K - 1,
        }
    }
}

pub struct SmallKmerIter<'a, const K: usize> {
    seq: &'a Sequence,
    kmer: crate::small::Kmer<K>,
    pos: usize,
}

impl<'a, const K: usize> Iterator for SmallKmerIter<'a, K> {
    type Item = crate::small::Kmer<K>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos * 2 >= self.seq.bases.len() {
            return None;
        }

        assert!(self.pos * 2 + 1 < self.seq.bases.len());

        let bit_pos = self.pos * 2;

        self.kmer.push(unsafe {
            Base::from_u8_unchecked(self.seq.bases[bit_pos..bit_pos + 2].load::<u8>())
        });

        self.pos += 1;
        Some(self.kmer.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_kmer() {
        let mut seq = Sequence::new();
        for _ in 0..23 {
            seq.push(Base::A);
        }

        let kmers: Vec<crate::small::Kmer<23>> = seq.kmers().collect();
        assert_eq!(kmers.len(), 1);
        eprintln!("{}", kmers[0]);

        let bases: Vec<Base> = kmers[0].bases().collect();
        assert_eq!(bases, vec![Base::A; 23]);
    }
}
