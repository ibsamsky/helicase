#![cfg(feature = "bitvec")]

use bitvec::order::Msb0;
use bitvec::vec::BitVec;

pub struct Sequence<const K: usize> {
    bases: BitVec<usize, Msb0>,
}

impl<const K: usize> Sequence<K> {
    pub fn new() -> Self {
        Self {
            bases: BitVec::new(),
        }
    }
}
