//! Like [`phase`][crate::phase] but limited to integer multiples of π/4,
//! corresponding to only phases accessible by combinations of *T* gates.
//!
//! All phases and arithmetic operations thereof are automatically performed
//! modulo 2*π*.

use num_complex::Complex64 as C64;
use crate::{
    graph2::clifft::Complex,
    phase::Phase,
};

/// A real phase, constrained to positive multiples of *π*/4, modulo 2*π*.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum TPhase {
    /// 0
    T0,
    /// *π*/4
    T1,
    /// *π*/2
    T2,
    /// 3*π*/4
    T3,
    /// *π*
    T4,
    /// 5*π*/4
    T5,
    /// 3*π*/2
    T6,
    /// 7*π*/4
    T7,
}

impl From<TPhase> for Phase {
    fn from(tph: TPhase) -> Self {
        match tph {
            TPhase::T0 => Self::zero(),
            TPhase::T1 => Self::pi4(),
            TPhase::T2 => Self::pi2(),
            TPhase::T4 => Self::pi(),
            _ => Self::new(tph as i64, 8),
        }
    }
}

impl From<TPhase> for Complex {
    fn from(tph: TPhase) -> Self {
        match tph {
            TPhase::T0 => Self { div2: 0, re:  1, im:  0, ph_pos:  0, ph_neg:  0 },
            TPhase::T1 => Self { div2: 0, re:  0, im:  0, ph_pos:  1, ph_neg:  0 },
            TPhase::T2 => Self { div2: 0, re:  0, im:  1, ph_pos:  0, ph_neg:  0 },
            TPhase::T3 => Self { div2: 0, re:  0, im:  0, ph_pos:  0, ph_neg: -1 },
            TPhase::T4 => Self { div2: 0, re: -1, im:  0, ph_pos:  0, ph_neg:  0 },
            TPhase::T5 => Self { div2: 0, re:  0, im:  0, ph_pos: -1, ph_neg:  0 },
            TPhase::T6 => Self { div2: 0, re:  0, im: -1, ph_pos:  0, ph_neg:  0 },
            TPhase::T7 => Self { div2: 0, re:  0, im:  0, ph_pos:  0, ph_neg:  1 },
        }
    }
}

macro_rules! impl_tphase_from_int {
    ( $int:ty ) => {
        impl From<$int> for TPhase {
            fn from(n: $int) -> Self {
                match n.rem_euclid(8) {
                    0 => Self::T0,
                    1 => Self::T1,
                    2 => Self::T2,
                    3 => Self::T3,
                    4 => Self::T4,
                    5 => Self::T5,
                    6 => Self::T6,
                    7 => Self::T7,
                    _ => unreachable!(),
                }
            }
        }
    }
}
impl_tphase_from_int!(u8);
impl_tphase_from_int!(u16);
impl_tphase_from_int!(u32);
impl_tphase_from_int!(u64);
impl_tphase_from_int!(u128);
impl_tphase_from_int!(usize);
impl_tphase_from_int!(i8);
impl_tphase_from_int!(i16);
impl_tphase_from_int!(i32);
impl_tphase_from_int!(i64);
impl_tphase_from_int!(i128);
impl_tphase_from_int!(isize);

impl TPhase {
    /// Convert to a complex number with modulus 1 and argument equal to `self`.
    pub fn cis(self) -> C64 { Phase::from(self).cis() }

    /// Convert to a complex number with modulus `r` and argument equal to
    /// `self`.
    pub fn as_polar(self, r: f64) -> C64 { Phase::from(self).as_polar(r) }

    pub(crate) fn label(&self) -> String {
        match self {
            Self::T0 => "".to_string(),
            Self::T1 => "π/4".to_string(),
            Self::T2 => "π/2".to_string(),
            Self::T3 => "3π/4".to_string(),
            Self::T4 => "π".to_string(),
            Self::T5 => "5π/4".to_string(),
            Self::T6 => "3π/2".to_string(),
            Self::T7 => "7π/4".to_string(),
        }
    }
}

impl std::ops::Neg for TPhase {
    type Output = Self;

    fn neg(self) -> Self { (-(self as i8)).into() }
}

macro_rules! impl_arith_tphase {
    (
        $trait:ident,
        $fun:ident,
        $op:tt,
        $trait_assign:ident,
        $fun_assign:ident,
        $op_assign:tt
    ) => {
        impl std::ops::$trait<TPhase> for TPhase {
            type Output = TPhase;

            fn $fun(self, rhs: TPhase) -> Self::Output {
                Self::from(self as i8 $op rhs as i8)
            }
        }

        impl std::ops::$trait_assign<TPhase> for TPhase {
            fn $fun_assign(&mut self, rhs: TPhase) {
                *self = *self $op rhs;
            }
        }
    }
}
impl_arith_tphase!(Add, add, +, AddAssign, add_assign, +=);
impl_arith_tphase!(Sub, sub, -, SubAssign, sub_assign, -=);
impl_arith_tphase!(Rem, rem, %, RemAssign, rem_assign, %=);

impl std::iter::Sum for TPhase {
    fn sum<I>(iter: I) -> Self
    where I: IntoIterator<Item = Self>
    {
        let mut acc = Self::T0;
        for tph in iter.into_iter() { acc += tph; }
        acc
    }
}

macro_rules! impl_muldiv_int {
    (
        $int:ty
    ) => {
        impl std::ops::Mul<$int> for TPhase {
            type Output = TPhase;

            fn mul(self, rhs: $int) -> Self::Output {
                Self::from(self as $int * rhs)
            }
        }

        impl std::ops::MulAssign<$int> for TPhase {
            fn mul_assign(&mut self, rhs: $int) {
                *self = *self * rhs;
            }
        }

        impl std::ops::Mul<TPhase> for $int {
            type Output = TPhase;

            fn mul(self, rhs: TPhase) -> Self::Output {
                TPhase::from(self * rhs as $int)
            }
        }

        impl std::ops::Div<$int> for TPhase {
            type Output = TPhase;

            fn div(self, rhs: $int) -> Self::Output {
                Self::from(self as $int / rhs)
            }
        }

        impl std::ops::DivAssign<$int> for TPhase {
            fn div_assign(&mut self, rhs: $int) {
                *self = *self / rhs;
            }
        }
    }
}
impl_muldiv_int!(u8);
impl_muldiv_int!(u16);
impl_muldiv_int!(u32);
impl_muldiv_int!(u64);
impl_muldiv_int!(u128);
impl_muldiv_int!(usize);
impl_muldiv_int!(i8);
impl_muldiv_int!(i16);
impl_muldiv_int!(i32);
impl_muldiv_int!(i64);
impl_muldiv_int!(i128);
impl_muldiv_int!(isize);



