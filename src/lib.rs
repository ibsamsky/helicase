//! kmer types

#![cfg_attr(feature = "unstable_nightly", feature(generic_const_exprs))]
#![cfg_attr(docsrs, feature(doc_auto_cfg))]
#![warn(clippy::all, missing_docs, rust_2018_idioms, unreachable_pub)]

mod base;
mod kmer;
mod sequence;

pub use base::Base;
pub use kmer::small;
#[cfg(feature = "bitvec")]
pub use kmer::{growable, unbounded};
pub use sequence::Sequence;

pub(crate) mod utils {
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
