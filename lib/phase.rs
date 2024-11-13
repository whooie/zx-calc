//! Numerically exact, real phases backed by rational numbers.
//!
//! All phases and arithmetic operations thereof are automatically performed
//! modulo 2*π*.

use std::f64::consts::TAU;
use num_complex::Complex64 as C64;
use num_rational::Rational64 as R64;
use num_traits::{ One, Zero };

// via Euclid's algorithm
fn gcd(mut a: i64, mut b: i64) -> i64 {
    let mut t: i64;
    while b != 0 {
        t = b;
        b = a % b;
        a = t;
    }
    a.abs()
}

fn lcm(a: i64, b: i64) -> i64 { (a / gcd(a, b)) * b }

// return the reduction of `a` modulo `m`, constrained to positive values
pub(crate) fn rempos(a: R64, m: R64) -> R64 {
    let d = lcm(*a.denom(), *m.denom());
    let b = (*(a * d).numer()).rem_euclid(*(m * d).numer());
    R64::new(b, d)
}

// return `true` if `a` and `b` are equal, modulo 1.
pub(crate) fn phase_eq(a: R64, b: R64) -> bool {
    rempos(a - b, R64::new(1, 1)) == R64::new(0, 1)
}

// convert a rational number to a floating-point number.
pub(crate) fn r2f(a: R64) -> f64 { *a.numer() as f64 / *a.denom() as f64 }

/// A description of a phase.
///
/// This type relies on rational approximation, holding an inner [`R64`]
/// representing the number *φ* such that the phase represented by a `Phase` as
/// a whole is 2*π* × *φ*.
///
/// Note that *φ* is constrained to positive values modulo 2*π* in all
/// operations, which has implications for multiplication and division.
///
/// ```
/// # use zx_calc::phase::Phase;
/// assert_eq!(  Phase::new(3, 4),      -Phase::new(1, 4) );
/// assert_eq!(  Phase::new(3, 4) / 2,   Phase::new(3, 8) );
/// assert_eq!( -Phase::new(1, 4) / 2,  -Phase::new(5, 8) );
/// assert_eq!( -Phase::new(5, 8),       Phase::new(3, 8) );
/// ```
#[derive(Copy, Clone, Debug)]
pub struct Phase(pub R64);

impl From<f64> for Phase {
    /// *Panics if the original floating-point number is non-normal.*
    fn from(f: f64) -> Self {
        let ph =
            R64::approximate_float(f / TAU)
            .expect("error converting to phase: unrepresentable float");
        Self(ph)
    }
}

impl From<Phase> for f64 {
    fn from(ph: Phase) -> Self {
        TAU * (*ph.0.numer() as f64 / *ph.0.denom() as f64)
    }
}

impl PartialEq for Phase {
    fn eq(&self, other: &Self) -> bool {
        rempos(self.0 - other.0, R64::one()) == R64::zero()
    }
}

impl Eq for Phase { }

impl PartialOrd for Phase {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Phase {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        rempos(self.0, R64::one()).cmp(&rempos(other.0, R64::one()))
    }
}

impl Phase {
    /// Construct a new `Phase` as `(numer / denom) × 2π`.
    pub fn new(numer: i64, denom: i64) -> Self {
        Self(rempos(R64::new(numer, denom), R64::one()))
    }

    /// Convert from a floating-point number.
    ///
    /// *Panics if the original number is non-normal.*
    pub fn from_float(f: f64) -> Self { f.into() }

    /// Convert to a floating-point number.
    pub fn into_float(self) -> f64 { self.into() }

    /// Return the `Phase` representation of 0 ≡ 2π mod 2π.
    pub fn zero() -> Self { Self(R64::zero()) }

    // /// Return the `Phase` representation of τ = 2π.
    // pub fn tau() -> Self { Self(R64::new(1, 1)) }

    /// Return the `Phase` representation of π.
    pub fn pi() -> Self { Self(R64::new(1, 2)) }

    /// Return the `Phase` representation of π/2.
    pub fn pi2() -> Self { Self(R64::new(1, 4)) }

    /// Return the `Phase` representation of π/4.
    pub fn pi4() -> Self { Self(R64::new(1, 8)) }

    /// Return the `Phase` representation of π/8.
    pub fn pi8() -> Self { Self(R64::new(1, 16)) }

    /// Return the `Phase` representation of 2π/`n`.
    pub fn frac(n: i64) -> Self { Self(rempos(R64::new(1, n), R64::one())) }

    /// Return a copy of `self` reduced modulo 2π.
    pub fn reduced(self) -> Self { Self(rempos(self.0, R64::one())) }

    /// Return `true` if `self` is an integer multiple of 2π/`n`.
    pub fn is_mult(self, n: i64) -> bool {
        self % Self::frac(n) == Self::zero()
    }

    /// Convert to a complex number with modulus 1 and argument equal to `self`.
    pub fn cis(self) -> C64 { C64::cis(self.into()) }

    /// Convert to a complex number with modulus `r` and argument equal to
    /// `self`.
    pub fn as_polar(self, r: f64) -> C64 { C64::from_polar(r, self.into()) }

    pub(crate) fn label(&self) -> String {
        if *self == Self::zero() {
            return "".to_string();
        } else if *self == Self::pi() {
            return "π".to_string();
        }
        let modpi = 2 * *self;
        if *modpi.0.numer() == 1 {
            format!("π/{}", modpi.0.denom())
        } else if *modpi.0.denom() <= 1000 {
            format!("({})π", modpi.0)
        } else {
            format!("{}π", r2f(modpi.0))
        }
    }
}

impl std::ops::Neg for Phase {
    type Output = Phase;

    fn neg(self) -> Self::Output {
        Self(rempos(-self.0, R64::one()))
    }
}

macro_rules! impl_addsubrem_phase {
    (
        $trait:ident,
        $fun:ident,
        $op:tt,
        $trait_assign:ident,
        $fun_assign:ident,
        $op_assign:tt
    ) => {
        impl std::ops::$trait<Phase> for Phase {
            type Output = Phase;

            fn $fun(self, rhs: Phase) -> Self::Output {
                Self(rempos(self.0 $op rhs.0, R64::one()))
            }
        }

        impl std::ops::$trait_assign<Phase> for Phase {
            fn $fun_assign(&mut self, rhs: Phase) {
                *self = *self $op rhs;
            }
        }
    }
}
impl_addsubrem_phase!(Add, add, +, AddAssign, add_assign, +=);
impl_addsubrem_phase!(Sub, sub, -, SubAssign, sub_assign, -=);
impl_addsubrem_phase!(Rem, rem, %, RemAssign, rem_assign, %=);

impl std::iter::Sum for Phase {
    fn sum<I>(iter: I) -> Self
    where I: IntoIterator<Item = Self>
    {
        let mut acc = Self::zero();
        for ph in iter.into_iter() { acc += ph; }
        acc
    }
}

impl std::ops::Mul<f64> for Phase {
    type Output = Phase;

    fn mul(self, rhs: f64) -> Self::Output {
        let rhs =
            R64::approximate_float(rhs)
            .expect("error in phase mul: unrepresentable float");
        Self(rempos(self.0 * rhs, R64::one()))
    }
}

impl std::ops::MulAssign<f64> for Phase {
    fn mul_assign(&mut self, rhs: f64) {
        let rhs =
            R64::approximate_float(rhs)
            .expect("error in phase mul_assign: unrepresentable float");
        *self = Phase(rempos(self.0 * rhs, R64::one()));
    }
}

impl std::ops::Mul<Phase> for f64 {
    type Output = Phase;

    fn mul(self, rhs: Phase) -> Self::Output {
        let lhs =
            R64::approximate_float(self)
            .expect("error in phase mul: unrepresentable float");
        Phase(rempos(rhs.0 * lhs, R64::one()))
    }
}

impl std::ops::Mul<i64> for Phase {
    type Output = Phase;

    fn mul(self, rhs: i64) -> Self::Output {
        Self(rempos(self.0 * rhs, R64::one()))
    }
}

impl std::ops::MulAssign<i64> for Phase {
    fn mul_assign(&mut self, rhs: i64) {
        *self = Phase(rempos(self.0 * rhs, R64::one()));
    }
}

impl std::ops::Mul<Phase> for i64 {
    type Output = Phase;

    fn mul(self, rhs: Phase) -> Self::Output {
        Phase(rempos(rhs.0 * self, R64::one()))
    }
}

impl std::ops::Div<f64> for Phase {
    type Output = Phase;

    fn div(self, rhs: f64) -> Self::Output {
        let rhs =
            R64::approximate_float(rhs)
            .expect("error in phase div: unrepresentable float");
        Self(rempos(self.0 / rhs, R64::one()))
    }
}

impl std::ops::DivAssign<f64> for Phase {
    fn div_assign(&mut self, rhs: f64) {
        let rhs =
            R64::approximate_float(rhs)
            .expect("error in phase div_assign: unrepresentable float");
        *self = Self(rempos(self.0 / rhs, R64::one()));
    }
}

impl std::ops::Div<i64> for Phase {
    type Output = Phase;

    fn div(self, rhs: i64) -> Self::Output {
        Self(rempos(self.0 / rhs, R64::one()))
    }
}

impl std::ops::DivAssign<i64> for Phase {
    fn div_assign(&mut self, rhs: i64) {
        *self = Self(rempos(self.0 / rhs, R64::one()));
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn init_mod() {
        assert_eq!(Phase::new(5, 3), Phase(R64::new(2, 3)));
        assert_eq!(Phase::new(4, 3), Phase::new(1, 3));
        assert_eq!(Phase::new(-1, 3), Phase::new(2, 3));
        assert_eq!(Phase::new(1, -3), Phase::new(2, 3));
        assert_eq!(Phase::zero(), Phase::zero());
        assert_eq!(Phase::pi(), Phase::new(1, 2));
        assert_eq!(Phase::pi2(), Phase::new(1, 4));
        assert_eq!(Phase::pi4(), Phase::new(1, 8));
        assert_eq!(Phase::pi8(), Phase::new(1, 16));
        assert_eq!(Phase::frac(1), Phase::zero());
        assert_eq!(Phase::frac(2), Phase::pi());
        assert_eq!(Phase::frac(3), Phase::new(1, 3));
    }

    #[test]
    fn float_conv() {
        fn approx_eq(f1: f64, f2: f64) -> bool { (f1 - f2).abs() < 1e-15 }

        assert_eq!(Phase::from_float(TAU), Phase::zero());
        assert_eq!(Phase::from_float(TAU / 2.0), Phase::pi());
        assert_eq!(Phase::from_float(TAU / 3.0), Phase::new(1, 3));
        assert!(approx_eq(Phase::zero().into_float(), 0.0));
        assert!(approx_eq(Phase::pi().into_float(), TAU / 2.0));
        assert!(approx_eq(Phase::frac(3).into_float(), TAU / 3.0));
    }

    #[test]
    fn add() {
        assert_eq!(Phase::zero() + Phase::zero(), Phase::zero());
        assert_eq!(Phase::pi() + Phase::pi(), Phase::zero());
        assert_eq!(Phase::new(1, 3) + Phase::new(2, 3), Phase::zero());
        assert_eq!(Phase::new(1, 3) + Phase::zero(), Phase::new(1, 3));
        assert_eq!(Phase::new(2, 3) + Phase::new(2, 3), Phase::new(1, 3));
    }

    #[test]
    fn sub() {
        assert_eq!(Phase::zero() - Phase::zero(), Phase::zero());
        assert_eq!(Phase::pi() - Phase::pi(), Phase::zero());
        assert_eq!(Phase::new(1, 3) - Phase::new(2, 3), Phase::new(2, 3));
        assert_eq!(Phase::new(1, 3) - Phase::zero(), Phase::new(1, 3));
        assert_eq!(Phase::new(2, 3) - Phase::new(2, 3), Phase::zero());
    }

    #[test]
    fn rem() {
        assert_eq!(Phase::pi() % Phase::pi(), Phase::zero());
        assert_eq!(Phase::new(1, 3) % Phase::new(2, 3), Phase::new(1, 3));
        assert_eq!(Phase::new(2, 3) % Phase::pi(), Phase::new(1, 6));
        assert_eq!(Phase::new(-2, 3) % Phase::pi(), Phase::new(1, 3));
        assert_eq!(Phase::new(2, 3) % Phase::new(2, 3), Phase::zero());
    }

    #[test]
    fn mul() {
        assert_eq!(Phase::new(1, 3) * 2, Phase::new(2, 3));
        assert_eq!(2 * Phase::new(1, 3), Phase::new(2, 3));
        assert_eq!(Phase::new(1, 3) * 3, Phase::zero());
        assert_eq!(Phase::new(1, 3) * 5, Phase::new(2, 3));
        assert_eq!(Phase::new(1, 3) * 0.5, Phase::new(1, 6));
        assert_eq!(Phase::new(1, 3) * (1.0 / 3.0), Phase::new(1, 9));
    }

    #[test]
    fn div() {
        assert_eq!(Phase::new(1, 3) / 2, Phase::new(1, 6));
        assert_eq!(Phase::new(1, 3) / 3, Phase::new(1, 9));
        assert_eq!(Phase::new(1, 3) / 5, Phase::new(1, 15));
        assert_eq!(Phase::new(1, 3) / 0.5, Phase::new(2, 3));
        assert_eq!(Phase::new(1, 3) / (1.0 / 3.0), Phase::zero());
    }
}


