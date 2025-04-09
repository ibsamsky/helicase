use std::fmt::Display;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]

/// A nucleotide base.
pub enum Base {
    /// Cytosine.
    C,
    /// Adenine.
    A,
    /// Thymine.
    T,
    /// Guanine.
    G,
}

impl Display for Base {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Base::C => write!(f, "C"),
            Base::A => write!(f, "A"),
            Base::T => write!(f, "T"),
            Base::G => write!(f, "G"),
        }
    }
}

impl TryFrom<u8> for Base {
    type Error = ();

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Base::C),
            1 => Ok(Base::A),
            2 => Ok(Base::T),
            3 => Ok(Base::G),
            _ => Err(()),
        }
    }
}

impl Base {
    /// Converts a u8 to a Base without checking bounds.
    ///
    /// # Safety
    /// The value must be in the range `0..4`.
    pub const unsafe fn from_u8_unchecked(value: u8) -> Self {
        match value {
            0 => Base::C,
            1 => Base::A,
            2 => Base::T,
            3 => Base::G,
            _ => unsafe { std::hint::unreachable_unchecked() },
        }
    }

    /// Converts an ASCII character to a Base.
    ///
    /// Returns `None` if the character is not a valid base.
    pub const fn from_ascii(value: u8) -> Option<Self> {
        match value {
            b'C' | b'c' => Some(Base::C),
            b'A' | b'a' => Some(Base::A),
            b'T' | b't' => Some(Base::T),
            b'G' | b'g' => Some(Base::G),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_u8() {
        assert_eq!(unsafe { Base::from_u8_unchecked(0) }, Base::C);
        assert_eq!(unsafe { Base::from_u8_unchecked(1) }, Base::A);
        assert_eq!(unsafe { Base::from_u8_unchecked(2) }, Base::T);
        assert_eq!(unsafe { Base::from_u8_unchecked(3) }, Base::G);
    }

    #[test]
    fn from_ascii() {
        assert_eq!(Base::from_ascii(b'C'), Some(Base::C));
        assert_eq!(Base::from_ascii(b'c'), Some(Base::C));
        assert_eq!(Base::from_ascii(b'A'), Some(Base::A));
        assert_eq!(Base::from_ascii(b'a'), Some(Base::A));
        assert_eq!(Base::from_ascii(b'T'), Some(Base::T));
        assert_eq!(Base::from_ascii(b't'), Some(Base::T));
        assert_eq!(Base::from_ascii(b'G'), Some(Base::G));
        assert_eq!(Base::from_ascii(b'g'), Some(Base::G));
        assert_eq!(Base::from_ascii(b'X'), None);
    }
}
