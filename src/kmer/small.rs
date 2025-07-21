//! Fixed-size k-mers up to 32 bases.
//!
//! # Example
//!
//! ```
//! use helicase::small::Kmer;
//! use helicase::Base;
//!
//! let mut kmer = Kmer::<5>::new();
//! kmer.push(Base::A).push(Base::A).push(Base::G).push(Base::T);
//! assert_eq!(kmer.as_masked(), 0b00_01_01_11_10);
//! ```
//!
//! # Limitations
//!
//! This implementation is not suitable for k-mers with more than 32 bases, as
//! it uses a `u64` to store the k-mer.

use std::fmt::Display;
use std::iter::FusedIterator;

use crate::base::Base;
use crate::utils;

/// A fixed-size k-mer represented as a 64-bit integer.
///
/// Stores between 1 and 32 bases. For larger k-mers, see
#[cfg_attr(feature = "unstable_nightly", doc = "[`unbounded`] and")]
/// [`growable`].
///
#[cfg_attr(
    feature = "unstable_nightly",
    doc = "[`unbounded`]: crate::kmer::unbounded"
)]
/// [`growable`]: crate::kmer::growable
#[derive(Debug, Clone, Copy)]
pub struct Kmer<const K: usize> {
    inner: u64,
}

impl<const K: usize> Display for Kmer<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.bases().map(|b| b.to_string()).collect::<String>()
        )
    }
}

impl<const K: usize> Default for Kmer<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const K: usize> From<u64> for Kmer<K> {
    fn from(value: u64) -> Self {
        Self { inner: value }
    }
}

impl<const K: usize> Kmer<K> {
    /// Creates a new k-mer.
    ///
    /// All bases are initialized to `Base::C`.
    pub const fn new() -> Self {
        utils::const_eval::assert_less::<0, K>();
        utils::const_eval::assert_leq::<K, 32>();
        Self { inner: 0 }
    }

    pub fn from_bases(bases: [Base; K]) -> Self {
        let inner = bases.iter().fold(0, |acc, base| acc << 2 | *base as u64);
        Self { inner }
    }

    /// Pushes a base onto the k-mer.
    ///
    /// Bases are pushed to the end of the k-mer, and the bases are shifted to the left, removing the first base.
    pub fn push(&mut self, base: Base) -> &mut Self {
        self.inner = (self.inner << 2) | (base as u64);
        self
    }

    /// Returns an iterator over the bases in the k-mer.
    pub const fn bases(&self) -> Bases<'_, K> {
        Bases {
            inner: self,
            pos: 0,
        }
    }

    /// Convert the k-mer into its inner value, masked to `K` bases.
    pub const fn as_masked(&self) -> u64 {
        bitfrob::u64_get_region(0, if K > 31 { 63 } else { K as u32 * 2 - 1 }, self.inner)
    }

    /// Shrinks the k-mer to a new size.
    ///
    /// # Panics
    ///
    /// Panics if `L` is greater than `K`.
    ///
    /// # Examples
    ///
    /// ```
    /// use helicase::small::Kmer;
    /// use helicase::Base;
    ///
    /// let mut kmer = Kmer::<5>::new();
    /// kmer.push(Base::A).push(Base::A).push(Base::G).push(Base::T);
    /// let resized = kmer.shrink_to::<3>();
    /// assert_eq!(resized.as_masked(), 0b01_11_10);
    /// ```
    pub const fn shrink_to<const L: usize>(self) -> Kmer<L> {
        utils::const_eval::assert_leq::<L, K>();
        Kmer { inner: self.inner }
    }

    #[cfg(feature = "unstable_nightly")]
    /// Joins two k-mers into a new k-mer.
    ///
    /// # Panics
    ///
    /// Panics if the combined size of the k-mers is greater than 32.
    pub fn join<const L: usize>(self, other: Kmer<L>) -> Kmer<{ K + L }> {
        utils::const_eval::assert_sum_leq::<K, L, 32>();
        let mut kmer = Kmer::<{ K + L }>::new();
        kmer.inner = (self.inner << (L * 2)) | other.inner;
        kmer
    }
}

/// An iterator over the bases in a k-mer.
#[derive(Debug)]
pub struct Bases<'a, const K: usize> {
    inner: &'a Kmer<K>,
    pos: usize,
}

impl<'a, const K: usize> Iterator for Bases<'a, K> {
    type Item = Base;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= K {
            return None;
        }

        let base = (self.inner.inner >> ((K - self.pos - 1) * 2)) as u8 & 3;
        self.pos += 1;

        // SAFETY: `base` is always in the range `0..4`.
        Some(unsafe { Base::from_u8_unchecked(base) })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = K - self.pos;
        (remaining, Some(remaining))
    }
}

impl<'a, const K: usize> FusedIterator for Bases<'a, K> {}

impl<'a, const K: usize> ExactSizeIterator for Bases<'a, K> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new() {
        macro_rules! test_new {
            ($k:expr $(,)?) => {
                assert_eq!(Kmer::<{ $k }>::new().inner, 0);
            };
            ($($k:expr),*) => {
                $(
                    assert_eq!(Kmer::<{ $k }>::new().inner, 0);
                )*
            };
        }

        test_new!(
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
            25, 26, 27, 28, 29, 30, 31, 32
        );
    }

    #[test]
    fn bases_simple() {
        let mut kmer = Kmer::<5>::new();
        kmer.push(Base::A)
            .push(Base::C)
            .push(Base::G)
            .push(Base::T)
            .push(Base::A);

        let bases: Vec<Base> = kmer.bases().collect();
        assert_eq!(bases, vec![Base::A, Base::C, Base::G, Base::T, Base::A]);
    }

    #[test]
    fn bases_large() {
        let mut kmer = Kmer::<23>::new();
        for _ in 0..10 {
            kmer.push(Base::A).push(Base::C);
        }

        for _ in 0..4 {
            kmer.push(Base::T);
        }

        // (A)CACACACACACACACACACTTTT

        let bases: Vec<Base> = kmer.bases().collect();

        for i in 0..9 {
            assert_eq!(bases[i * 2], Base::C);
            assert_eq!(bases[i * 2 + 1], Base::A);
        }

        assert_eq!(bases[18], Base::C);

        assert_eq!(bases[19..23], vec![Base::T, Base::T, Base::T, Base::T]);
    }

    #[test]
    fn convert() {
        let mut kmer = Kmer::<5>::new();
        kmer.push(Base::A)
            .push(Base::C)
            .push(Base::G)
            .push(Base::T)
            .push(Base::A);

        let small = kmer.shrink_to::<2>();
        let bases: Vec<Base> = small.bases().collect();
        assert_eq!(bases, vec![Base::T, Base::A]);
    }

    #[test]
    fn into_masked() {
        let mut kmer: Kmer<2> = Kmer::new();
        kmer.push(Base::A).push(Base::A);
        assert_eq!(kmer.as_masked(), 0x05);
    }

    #[cfg(feature = "unstable_nightly")]
    #[test]
    fn join() {
        let mut kmer1 = Kmer::<2>::new();
        kmer1.push(Base::A).push(Base::C);
        let mut kmer2 = Kmer::<2>::new();
        kmer2.push(Base::G).push(Base::T);
        let kmer3 = kmer1.join(kmer2);
        let bases: Vec<Base> = kmer3.bases().collect();
        assert_eq!(bases, vec![Base::A, Base::C, Base::G, Base::T]);
    }
}
