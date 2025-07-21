use std::iter::FusedIterator;

use crate::Base;

pub struct Sequence {
    bases: Vec<Base>,
}

impl Sequence {
    pub fn new() -> Self {
        Self { bases: vec![] }
    }

    pub fn push(&mut self, base: Base) {
        self.bases.push(base);
    }

    pub fn kmers<const K: usize>(&self) -> SmallKmerIter<'_, K> {
        let mut kmer = crate::small::Kmer::<K>::new();
        self.bases.iter().take(K - 1).copied().for_each(|base| {
            kmer.push(base);
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
        self.kmer.push(*self.seq.bases.get(self.pos)?);
        self.pos += 1;
        Some(self.kmer)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.seq.bases.len() - self.pos;
        (remaining, Some(remaining))
    }
}

impl<'a, const K: usize> FusedIterator for SmallKmerIter<'a, K> {}

impl<'a, const K: usize> ExactSizeIterator for SmallKmerIter<'a, K> {}

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
