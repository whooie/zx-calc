use std::fmt;
use num_complex::Complex64 as C64;

/// Like [`Ord`], but without the expectation that the ordering is based on the
/// whole value; this trait should produce orderings based *only* on quantum
/// state information.
pub trait StateOrd {
    fn state_cmp(&self, other: &Self) -> std::cmp::Ordering;
}

/// The state on a single wire.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum State {
    /// ∣0⟩ = ∣+⟩ + ∣–⟩
    Zero,
    /// ∣1⟩ = ∣+⟩ – ∣–⟩
    One,
    /// ∣+⟩ = ∣0⟩ + ∣1⟩
    Plus,
    /// ∣–⟩ = ∣0⟩ - ∣1⟩
    Minus,
}

impl State {
    /// Calculate the dot product `⟨self|rhs⟩`.
    pub fn dot(&self, other: &Self) -> C64 {
        use std::f64::consts::FRAC_1_SQRT_2;
        match (self, other) {
            (Self::Zero,  Self::Zero ) => 1.0.into(),
            (Self::Zero,  Self::One  ) => 0.0.into(),
            (Self::Zero,  Self::Plus ) => FRAC_1_SQRT_2.into(),
            (Self::Zero,  Self::Minus) => FRAC_1_SQRT_2.into(),
            (Self::One,   Self::Zero ) => 0.0.into(),
            (Self::One,   Self::One  ) => 1.0.into(),
            (Self::One,   Self::Plus ) => FRAC_1_SQRT_2.into(),
            (Self::One,   Self::Minus) => (-FRAC_1_SQRT_2).into(),
            (Self::Plus,  Self::Zero ) => FRAC_1_SQRT_2.into(),
            (Self::Plus,  Self::One  ) => FRAC_1_SQRT_2.into(),
            (Self::Plus,  Self::Plus ) => 1.0.into(),
            (Self::Plus,  Self::Minus) => 0.0.into(),
            (Self::Minus, Self::Zero ) => FRAC_1_SQRT_2.into(),
            (Self::Minus, Self::One  ) => (-FRAC_1_SQRT_2).into(),
            (Self::Minus, Self::Plus ) => 0.0.into(),
            (Self::Minus, Self::Minus) => 1.0.into(),
        }
    }

    /// Return `true` if `self` is a member of `basis`.
    pub fn is_basis(&self, basis: Basis) -> bool {
        matches!(
            (self, basis),
            (Self::Zero,  Basis::Z)
            | (Self::One,   Basis::Z)
            | (Self::Plus,  Basis::X)
            | (Self::Minus, Basis::X)
        )
    }

    /// Return `self` decomposed in `basis` if `self` is not already a member of
    /// `basis`.
    pub fn decomp(&self, basis: Basis) -> Option<((C64, Self), (C64, Self))> {
        if self.is_basis(basis) {
            None
        } else {
            let up = basis.up();
            let down = basis.down();
           Some( ( (self.dot(&up), up), (self.dot(&down), down) ) )
        }
    }

    /// Apply a Hadamard transformation to `self` in place, swapping `0 ↔ +` and
    /// `1 ↔ -`.
    pub fn hadamard_mut(&mut self) {
        match *self {
            Self::Zero  => { *self = Self::Plus  },
            Self::One   => { *self = Self::Minus },
            Self::Plus  => { *self = Self::Zero  },
            Self::Minus => { *self = Self::One   },
        }
    }

    /// Apply a Hadamard transformation `self`, swapping `0 ↔ +` and `1 ↔ -`.
    pub fn hadamard(&self) -> Self {
        match self {
            Self::Zero  => Self::Plus,
            Self::One   => Self::Minus,
            Self::Plus  => Self::Zero,
            Self::Minus => Self::One,
        }
    }
}

impl StateOrd for State {
    fn state_cmp(&self, other: &Self) -> std::cmp::Ordering { self.cmp(other) }
}

impl fmt::Display for State {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::Zero  => "0".fmt(f),
            Self::One   => "1".fmt(f),
            Self::Plus  => "+".fmt(f),
            Self::Minus => "-".fmt(f),
        }
    }
}

/// A specific basis for representing states.
///
/// The *Z* basis refers to [`State`]s `Zero` and `One`, while the *X* basis
/// refers to `Plus` and `Minus`. The *Z* basis is chosen by default for actual
/// computations.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Basis {
    Z,
    X,
}

impl Basis {
    // get the "up" (0/+) state for the basis
    pub(crate) fn up(&self) -> State {
        match self {
            Self::Z => State::Zero,
            Self::X => State::Plus,
        }
    }

    // get the "down" (1/-) state for the basis
    pub(crate) fn down(&self) -> State {
        match self {
            Self::Z => State::One,
            Self::X => State::Minus,
        }
    }
}

/// Describes a non-continguously indexed collection of [`State`]s, to be used
/// as one "half" of a ket-bra.
pub type States = crate::indexmap::IndexMap<State>;

/// Iterator over references to the states on one "half" of a ket-bra.
///
/// The iterator item type is `(usize, &`[`State`]`)`.
pub type StatesIter<'a> = crate::indexmap::Iter<'a, State>;

/// Iterator over mutable references to the states on one "half" of a ket-bra.
///
/// The iterator item type is `(usize, &'a mut `[`State`]`)`.
pub type StatesIterMut<'a> = crate::indexmap::IterMut<'a, State>;

/// Iterator over the states on one "half" of a ket-bra.
///
/// The iterator item type is `(usize, `[`State`]`)`.
pub type StatesIntoIter<'a> = crate::indexmap::IntoIter<State>;

