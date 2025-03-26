//! kmer types

#![cfg_attr(feature = "unstable_nightly", feature(generic_const_exprs))]
#![cfg_attr(docsrs, feature(doc_auto_cfg))]
#![warn(clippy::all, missing_docs, rust_2018_idioms, unreachable_pub)]

mod base;
mod kmer;

pub use base::Base;
pub use kmer::small;
#[cfg(feature = "bitvec")]
pub use kmer::{growable, unbounded};

pub(crate) mod utils {
    pub(crate) const fn saturating_bitmask(bits: u32) -> u64 {
        u64::MAX.unbounded_shr(64u32.saturating_sub(bits))
    }

    pub(crate) const fn bitmask_range(m: u32, n: u32) -> u64 {
        saturating_bitmask(n.saturating_sub(m)).unbounded_shl(m)
    }

    pub(crate) mod const_eval {
        pub(crate) const fn assert_less<const L: usize, const K: usize>() {
            assert!(L < K);
        }

        pub(crate) const fn assert_leq<const L: usize, const K: usize>() {
            assert!(L <= K);
        }

        pub(crate) const fn assert_sum_less<const L: usize, const K: usize, const M: usize>() {
            assert!(L + K < M);
        }

        pub(crate) const fn assert_sum_leq<const L: usize, const K: usize, const M: usize>() {
            assert!(L + K <= M);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bitmask() {
        assert_eq!(utils::saturating_bitmask(0), 0);
        assert_eq!(utils::saturating_bitmask(1), 1);
        assert_eq!(utils::saturating_bitmask(2), 3);
        assert_eq!(utils::saturating_bitmask(3), 7);

        for i in 0..64 {
            assert_eq!(utils::saturating_bitmask(i), (1 << i) - 1);
        }

        for i in 64..128 {
            assert_eq!(utils::saturating_bitmask(i), u64::MAX);
        }
    }

    #[test]
    fn bitmask_range() {
        assert_eq!(utils::bitmask_range(0, 63), 0x7fffffffffffffff);
        assert_eq!(utils::bitmask_range(1, 63), 0x7ffffffffffffffe);
        assert_eq!(utils::bitmask_range(0, 64), 0xffffffffffffffff);
        assert_eq!(utils::bitmask_range(1, 64), 0xfffffffffffffffe);
    }
}
