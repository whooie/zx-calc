//! Complex numbers represented as members of the ring
//! **D**[exp(*iπ*/4)], where **D** is the ring of dyadic rational numbers.
//!
//! All such elements can be represented using five integers
//! (*k*, *a*, *b*, *c*, *d*), giving a complex number as
//! (*a* + *b* *i* + *c* exp(*iπ*/4) + *d* exp(–*iπ*/4)). These are
//! exact representations of all possible scalars that can arise from
//! Clifford+*T* quantum circuits.

use std::ops::{ Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Neg };
use num_complex::Complex64 as C64;
use num_traits::{ Zero, One };
use crate::graph::ComplexRing;

/// An element of the ring **D**[exp(*i* *π* / 4)], where **D** is the ring of
/// dyadic rational numbers.
///
/// All elements of this ring are complex numbers that can be written in the
/// form
///
/// <blockquote>
///   <p style="font-size:20px">
///     (
///         <i>a</i>
///         + <i>b i</i>
///         + <i>c</i> exp(<i>i π</i> / 4)
///         + <i>d</i> exp(–<i>i π</i> / 4)
///     ) / 2<sup><i>k</i></sup>
///   </p>
/// </blockquote>
///
/// Further, every possible scalar that can arise from Clifford+T quantum
/// circuits can be represented exactly in this form.
#[derive(Copy, Clone, Debug)]
pub struct Complex {
    /// Exponent on the outer factor of 1/2.
    pub div2: i32,
    /// Purely real part.
    pub re: i32,
    /// Purely imaginary part.
    pub im: i32,
    /// exp(<i>i π</i> / 4) component.
    pub ph_pos: i32,
    /// exp(–<i>i π</i> / 4) component.
    pub ph_neg: i32,
}

#[derive(Copy, Clone, Debug)]
enum IsPow2 {
    Pow(u8),
    Zero,
    Nil,
}

impl From<Option<u8>> for IsPow2 {
    fn from(mb_pow2: Option<u8>) -> Self {
        match mb_pow2 {
            Some(pow) => IsPow2::Pow(pow),
            None => IsPow2::Nil,
        }
    }
}

fn is_pow_2(n: i32) -> IsPow2 {
    if n == 0 { return IsPow2::Zero; }
    if n == i32::MIN { return IsPow2::Pow(31); }
    let nabs = n.abs();
    if nabs == 1 { return IsPow2::Pow(0); }
    let mut check: i32 = 1;
    (1..=30_u8).find(|_| {
        check <<= 1;
        check == nabs
    }).into()
}

impl Complex {
    /// The constant value 0.
    pub const ZERO: Self =
        Self { div2: 0, re: 0, im: 0, ph_pos: 0, ph_neg: 0 }; 

    /// The real unit 1.
    pub const ONE: Self =
        Self { div2: 0, re: 1, im: 0, ph_pos: 0, ph_neg: 0 }; 

    /// The imaginary unit *i*.
    pub const I: Self =
        Self { div2: 0, re: 0, im: 1, ph_pos: 0, ph_neg: 0 }; 

    /// The positive-phase unit exp(*i* *π* / 4).
    pub const PH_POS: Self =
        Self { div2: 0, re: 0, im: 0, ph_pos: 1, ph_neg: 0 }; 

    /// The negative-phase unit exp(–*i* *π* / 4).
    pub const PH_NEG: Self =
        Self { div2: 0, re: 0, im: 0, ph_pos: 0, ph_neg: 1 }; 

    /// The real number √2.
    pub const SQRT_2: Self =
        Self { div2: 0, re: 0, im: 0, ph_pos: 1, ph_neg: 1 }; 

    /// The imaginary number *i*√2.
    pub const ISQRT_2: Self =
        Self { div2: 0, re: 0, im: 0, ph_pos: 1, ph_neg: -1 }; 

    /// The real number 1/√2.
    pub const FRAC_1_SQRT_2: Self =
        Self { div2: 1, re: 0, im: 0, ph_pos: 1, ph_neg: 1 }; 

    /// The imaginary unmber *i*/√2.
    pub const FRAC_I_SQRT_2: Self =
        Self { div2: 1, re: 0, im: 0, ph_pos: 1, ph_neg: -1 }; 

    /// the constant value 0.
    pub fn zero() -> Self { Self::ZERO }

    /// the real unit 1.
    pub fn one() -> Self { Self::ONE }

    /// the imaginary unit *i*.
    pub fn i() -> Self { Self::I }

    /// the positive-phase unit exp(*i* *π* / 4).
    pub fn ph_pos() -> Self { Self::PH_POS }

    /// the negative-phase unit exp(–*i* *π* / 4).
    pub fn ph_neg() -> Self { Self::PH_NEG }

    /// the real number √2.
    pub fn sqrt_2() -> Self { Self::SQRT_2 }

    /// the imaginary number *i*√2.
    pub fn isqrt_2() -> Self { Self::ISQRT_2 }

    /// the real number 1/√2.
    pub fn frac_1_sqrt_2() -> Self { Self::FRAC_1_SQRT_2 }

    /// the imaginary unmber *i*/√2.
    pub fn frac_i_sqrt_2() -> Self { Self::FRAC_I_SQRT_2 }

    /// Create a new `Complex`.
    pub fn new(div2: i32, re: i32, im: i32, ph_pos: i32, ph_neg: i32) -> Self {
        let mut z = Self { div2, re, im, ph_pos, ph_neg };
        z.reduce();
        z
    }

    /// Create a new, purely real `Complex`.
    pub fn from_re(re: i32, div2: i32) -> Self {
        let mut z = Self { div2, re, im: 0, ph_pos: 0, ph_neg: 0 };
        z.reduce();
        z
    }

    /// Create a new, purely imaginary `Complex`.
    pub fn from_im(im: i32, div2: i32) -> Self {
        let mut z = Self { div2, re: 0, im, ph_pos: 0, ph_neg: 0 };
        z.reduce();
        z
    }

    /// Create a new `Complex` proportional to exp(*i* *π* / 4).
    pub fn from_ph_pos(ph_pos: i32, div2: i32) -> Self {
        let mut z = Self { div2, re: 0, im: 0, ph_pos, ph_neg: 0 };
        z.reduce();
        z
    }

    /// Create a new `Complex` proportional to exp(–*i* *π* / 4).
    pub fn from_ph_neg(ph_neg: i32, div2: i32) -> Self {
        let mut z = Self { div2, re: 0, im: 0, ph_pos: 0, ph_neg };
        z.reduce();
        z
    }

    pub(crate) fn reduce(&mut self) {
        use IsPow2::*;
        let mb_pow_re = is_pow_2(self.re);
        let mb_pow_im = is_pow_2(self.im);
        let mb_pow_ph_pos = is_pow_2(self.ph_pos);
        let mb_pow_ph_neg = is_pow_2(self.ph_neg);
        match (mb_pow_re, mb_pow_im, mb_pow_ph_pos, mb_pow_ph_neg) {
            (Pow(pow_re), Pow(pow_im), Pow(pow_pos), Pow(pow_neg)) => {
                let pow = pow_re.min(pow_im).min(pow_pos).min(pow_neg);
                self.div2 -= i32::from(pow);
                self.re >>= pow;
                self.im >>= pow;
                self.ph_pos >>= pow;
                self.ph_neg >>= pow;
            },
            (Pow(pow_re), Pow(pow_im), Pow(pow_pos), Zero) => {
                let pow = pow_re.min(pow_im).min(pow_pos);
                self.div2 -= i32::from(pow);
                self.re >>= pow;
                self.im >>= pow;
                self.ph_pos >>= pow;
            },
            (Pow(pow_re), Pow(pow_im), Zero, Pow(pow_neg)) => {
                let pow = pow_re.min(pow_im).min(pow_neg);
                self.div2 -= i32::from(pow);
                self.re >>= pow;
                self.im >>= pow;
                self.ph_neg >>= pow;
            },
            (Pow(pow_re), Zero, Pow(pow_pos), Pow(pow_neg)) => {
                let pow = pow_re.min(pow_pos).min(pow_neg);
                self.div2 -= i32::from(pow);
                self.re >>= pow;
                self.ph_pos >>= pow;
                self.ph_neg >>= pow;
            },
            (Zero, Pow(pow_im), Pow(pow_pos), Pow(pow_neg)) => {
                let pow = pow_im.min(pow_pos).min(pow_neg);
                self.div2 -= i32::from(pow);
                self.im >>= pow;
                self.ph_pos >>= pow;
                self.ph_neg >>= pow;
            },
            (Pow(pow_re), Pow(pow_im), Zero, Zero) => {
                let pow = pow_re.min(pow_im);
                self.div2 -= i32::from(pow);
                self.re >>= pow;
                self.im >>= pow;
            },
            (Pow(pow_re), Zero, Pow(pow_pos), Zero) => {
                let pow = pow_re.min(pow_pos);
                self.div2 -= i32::from(pow);
                self.re >>= pow;
                self.ph_pos >>= pow;
            },
            (Pow(pow_re), Zero, Zero, Pow(pow_neg)) => {
                let pow = pow_re.min(pow_neg);
                self.div2 -= i32::from(pow);
                self.re >>= pow;
                self.ph_neg >>= pow;
            },
            (Zero, Pow(pow_im), Pow(pow_pos), Zero) => {
                let pow = pow_im.min(pow_pos);
                self.div2 -= i32::from(pow);
                self.im >>= pow;
                self.ph_pos >>= pow;
            },
            (Zero, Pow(pow_im), Zero, Pow(pow_neg)) => {
                let pow = pow_im.min(pow_neg);
                self.div2 -= i32::from(pow);
                self.im >>= pow;
                self.ph_neg >>= pow;
            },
            (Zero, Zero, Pow(pow_pos), Pow(pow_neg)) => {
                let pow = pow_pos.min(pow_neg);
                self.div2 -= i32::from(pow);
                self.ph_pos >>= pow;
                self.ph_neg >>= pow;
            },
            (Pow(pow_re), Zero, Zero, Zero) => {
                self.div2 -= i32::from(pow_re);
                self.re >>= pow_re;
            },
            (Zero, Pow(pow_im), Zero, Zero) => {
                self.div2 -= i32::from(pow_im);
                self.im >>= pow_im;
            },
            (Zero, Zero, Pow(pow_pos), Zero) => {
                self.div2 -= i32::from(pow_pos);
                self.ph_pos >>= pow_pos;
            },
            (Zero, Zero, Zero, Pow(pow_neg)) => {
                self.div2 -= i32::from(pow_neg);
                self.ph_neg >>= pow_neg;
            },
            (Zero, Zero, Zero, Zero) => {
                self.div2 = 0;
            },
            _ => { },
        }
    }

    pub(crate) fn reduced(&self) -> Self {
        let mut new = *self;
        new.reduce();
        new
    }

    /// Raise `self` to a non-negative integer power.
    pub fn pow(self, pow: u32) -> Self {
        (0..pow).fold(Self::ONE, |mut acc, _| { acc *= self; acc })
    }

    /// Multiply by a power of two.
    pub fn mul2(mut self, pow: i32) -> Self {
        self.div2 -= pow;
        self
    }

    /// Alias for `self.mul2(1)`.
    pub fn double(self) -> Self { self.mul2(1) }

    /// Alias for `self.mul2(-1)`.
    pub fn half(self) -> Self { self.mul2(-1) }

    /// Multiply by an integer power of exp(*i* *π* / 4).
    pub fn rot(mut self, pow: i32) -> Self {
        use std::mem::swap;
        match pow.rem_euclid(8) {
            0 => { },
            1 => {
                let mut temp: i32 = 0;
                swap(&mut temp, &mut self.re);
                swap(&mut temp, &mut self.ph_pos);
                swap(&mut temp, &mut self.im);
                swap(&mut temp, &mut self.ph_neg);
                self.ph_neg *= -1;
                swap(&mut temp, &mut self.re);
            },
            2 => {
                swap(&mut self.re, &mut self.im);
                self.re *= -1;
                swap(&mut self.ph_pos, &mut self.ph_neg);
                self.ph_neg *= -1;
            },
            3 => {
                let mut temp: i32 = 0;
                swap(&mut temp, &mut self.re);
                swap(&mut temp, &mut self.ph_neg);
                self.ph_neg *= -1;
                swap(&mut temp, &mut self.im);
                swap(&mut temp, &mut self.ph_pos);
                self.ph_pos *= -1;
                swap(&mut temp, &mut self.re);
                self.re *= -1;
            },
            4 => {
                self.re *= -1;
                self.im *= -1;
                self.ph_pos *= -1;
                self.ph_neg *= -1;
            },
            5 => {
                let mut temp: i32 = 0;
                swap(&mut temp, &mut self.re);
                swap(&mut temp, &mut self.ph_pos);
                self.ph_pos *= -1;
                swap(&mut temp, &mut self.im);
                self.im *= -1;
                swap(&mut temp, &mut self.ph_neg);
                swap(&mut temp, &mut self.re);
                self.re *= -1;
            },
            6 => {
                swap(&mut self.re, &mut self.im);
                self.im *= -1;
                swap(&mut self.ph_pos, &mut self.ph_neg);
                self.ph_pos *= -1;
            },
            7 => {
                let mut temp: i32 = 0;
                swap(&mut temp, &mut self.re);
                swap(&mut temp, &mut self.ph_neg);
                swap(&mut temp, &mut self.im);
                self.im *= -1;
                swap(&mut temp, &mut self.ph_pos);
                swap(&mut temp, &mut self.re);
            },
            _ => unreachable!(),
        }
        self
    }
}

impl PartialEq for Complex {
    fn eq(&self, other: &Self) -> bool {
        let l = self.reduced();
        let r = other.reduced();
        l.div2 == r.div2
            && l.re == r.re
            && l.im == r.im
            && l.ph_pos == r.ph_pos
            && l.ph_neg == r.ph_neg
    }
}

impl Eq for Complex { }

impl From<i32> for Complex {
    fn from(re: i32) -> Self {
        Self { div2: 0, re, im: 0, ph_pos: 0, ph_neg: 0 }
    }
}

impl From<Complex> for C64 {
    fn from(z: Complex) -> Self {
        use std::f64::consts::FRAC_1_SQRT_2 as ONRT2;
        const PH_POS_PI4: C64 = C64 { re: ONRT2, im:  ONRT2 };
        const PH_NEG_PI4: C64 = C64 { re: ONRT2, im: -ONRT2 };
        const PH_POS_PI2: C64 = C64 { re: 0.0,   im: 1.0    };
        let a = f64::from(z.re);
        let b = f64::from(z.im);
        let c = f64::from(z.ph_pos);
        let d = f64::from(z.ph_neg);
        (a + b * PH_POS_PI2 + c * PH_POS_PI4 + d * PH_NEG_PI4)
            / 2.0_f64.powi(z.div2)
    }
}

impl Neg for Complex {
    type Output = Self;

    fn neg(mut self) -> Self {
        self.re = -self.re;
        self.im = -self.im;
        self.ph_pos = -self.ph_pos;
        self.ph_neg = -self.ph_neg;
        self
    }
}

impl AddAssign<Complex> for Complex {
    fn add_assign(&mut self, mut rhs: Complex) {
        match self.div2.cmp(&rhs.div2) {
            std::cmp::Ordering::Less => {
                let powdiff = 2_i32.pow(rhs.div2.abs_diff(self.div2));
                self.div2 += (rhs.div2 - self.div2).abs();
                self.re *= powdiff;
                self.im *= powdiff;
                self.ph_pos *= powdiff;
                self.ph_neg *= powdiff;
                self.re += rhs.re;
                self.im += rhs.im;
                self.ph_pos += rhs.ph_pos;
                self.ph_neg += rhs.ph_neg;
            },
            std::cmp::Ordering::Greater => {
                let powdiff = 2_i32.pow(self.div2.abs_diff(rhs.div2));
                self.div2 += (self.div2 - rhs.div2).abs();
                rhs.re *= powdiff;
                rhs.im *= powdiff;
                rhs.ph_pos *= powdiff;
                rhs.ph_neg *= powdiff;
                self.re += rhs.re;
                self.im += rhs.im;
                self.ph_pos += rhs.ph_pos;
                self.ph_neg += rhs.ph_neg;
            },
            std::cmp::Ordering::Equal => {
                self.re += rhs.re;
                self.im += rhs.im;
                self.ph_pos += rhs.ph_pos;
                self.ph_neg += rhs.ph_neg;
            },
        }
        self.reduce();
    }
}

impl Add<Complex> for Complex {
    type Output = Self;

    fn add(mut self, rhs: Complex) -> Self {
        self += rhs;
        self
    }
}

impl AddAssign<i32> for Complex {
    fn add_assign(&mut self, rhs: i32) {
        match self.div2.cmp(&0) {
            std::cmp::Ordering::Less => {
                let pow = 2_i32.pow(self.div2.abs_diff(0));
                self.div2 = 0;
                self.re *= pow;
                self.im *= pow;
                self.ph_pos *= pow;
                self.ph_neg *= pow;
                self.re += rhs;
            },
            std::cmp::Ordering::Greater => {
                let pow = 2_i32.pow(self.div2.abs_diff(0));
                self.re += pow * rhs;
            },
            std::cmp::Ordering::Equal => {
                self.re += rhs;
            },
        }
        self.reduce();
    }
}

impl Add<i32> for Complex {
    type Output = Self;

    fn add(mut self, rhs: i32) -> Self {
        self += rhs;
        self
    }
}

impl Add<Complex> for i32 {
    type Output = Complex;

    fn add(self, mut rhs: Complex) -> Complex {
        rhs += self;
        rhs
    }
}

impl SubAssign<Complex> for Complex {
    fn sub_assign(&mut self, mut rhs: Complex) {
        match self.div2.cmp(&rhs.div2) {
            std::cmp::Ordering::Less => {
                let powdiff = 2_i32.pow(rhs.div2.abs_diff(self.div2));
                self.div2 -= (rhs.div2 - self.div2).abs();
                self.re *= powdiff;
                self.im *= powdiff;
                self.ph_pos *= powdiff;
                self.ph_neg *= powdiff;
                self.re -= rhs.re;
                self.im -= rhs.im;
                self.ph_pos -= rhs.ph_pos;
                self.ph_neg -= rhs.ph_neg;
            },
            std::cmp::Ordering::Greater => {
                let powdiff = 2_i32.pow(self.div2.abs_diff(rhs.div2));
                self.div2 -= (self.div2 - rhs.div2).abs();
                rhs.re *= powdiff;
                rhs.im *= powdiff;
                rhs.ph_pos *= powdiff;
                rhs.ph_neg *= powdiff;
                self.re -= rhs.re;
                self.im -= rhs.im;
                self.ph_pos -= rhs.ph_pos;
                self.ph_neg -= rhs.ph_neg;
            },
            std::cmp::Ordering::Equal => {
                self.re -= rhs.re;
                self.im -= rhs.im;
                self.ph_pos -= rhs.ph_pos;
                self.ph_neg -= rhs.ph_neg;
            },
        }
        self.reduce();
    }
}

impl Sub<Complex> for Complex {
    type Output = Self;

    fn sub(mut self, rhs: Complex) -> Self {
        self -= rhs;
        self
    }
}

impl SubAssign<i32> for Complex {
    fn sub_assign(&mut self, rhs: i32) {
        match self.div2.cmp(&0) {
            std::cmp::Ordering::Less => {
                let pow = 2_i32.pow(self.div2.abs_diff(0));
                self.div2 = 0;
                self.re *= pow;
                self.im *= pow;
                self.ph_pos *= pow;
                self.ph_neg *= pow;
                self.re -= rhs;
            },
            std::cmp::Ordering::Greater => {
                let pow = 2_i32.pow(self.div2.abs_diff(0));
                self.re -= pow * rhs;
            },
            std::cmp::Ordering::Equal => {
                self.re -= rhs;
            },
        }
        self.reduce();
    }
}

impl Sub<i32> for Complex {
    type Output = Self;

    fn sub(mut self, rhs: i32) -> Self {
        self -= rhs;
        self
    }
}

impl Sub<Complex> for i32 {
    type Output = Complex;

    fn sub(self, mut rhs: Complex) -> Complex {
        rhs -= -self;
        rhs
    }
}

impl MulAssign<Complex> for Complex {
    fn mul_assign(&mut self, rhs: Complex) {
        let Self { div2: k, re: a, im: b, ph_pos: c, ph_neg: d } = *self;
        let Self { div2: j, re: e, im: f, ph_pos: g, ph_neg: h } = rhs;
        self.div2 = k + j;
        self.re = a * e - b * f + c * h + d * g;
        self.im = a * f + b * e + c * g - d * h;
        self.ph_pos = a * g + b * h + c * e + d * f;
        self.ph_neg = a * h - b * g - c * f + d * e;
        self.reduce();
    }
}

impl Mul<Complex> for Complex {
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self {
        self *= rhs;
        self
    }
}

impl MulAssign<i32> for Complex {
    fn mul_assign(&mut self, rhs: i32) {
        match is_pow_2(rhs) {
            IsPow2::Pow(k) => {
                self.div2 -= i32::from(k);
                self.re *= rhs.signum();
                self.im *= rhs.signum();
                self.ph_pos *= rhs.signum();
                self.ph_neg *= rhs.signum();
            },
            IsPow2::Zero => {
                self.div2 = 0;
                self.re = 0;
                self.im = 0;
                self.ph_pos = 0;
                self.ph_neg = 0;
            },
            IsPow2::Nil => {
                self.re *= rhs;
                self.im *= rhs;
                self.ph_pos *= rhs;
                self.ph_neg *= rhs;
            },
        }
    }
}

impl Mul<i32> for Complex {
    type Output = Self;

    fn mul(mut self, rhs: i32) -> Self {
        self *= rhs;
        self
    }
}

impl Mul<Complex> for i32 {
    type Output = Complex;

    fn mul(self, mut rhs: Complex) -> Complex {
        rhs *= self;
        rhs
    }
}

impl Zero for Complex {
    fn zero() -> Self { Self::ZERO }

    fn is_zero(&self) -> bool { *self == Self::ZERO }
}

impl One for Complex {
    fn one() -> Self { Self::ONE }
}

impl ComplexRing for Complex {
    fn conj(mut self) -> Self {
        self.im = -self.im;
        std::mem::swap(&mut self.ph_pos, &mut self.ph_neg);
        self
    }
}

impl std::fmt::Display for Complex {
    #[allow(clippy::comparison_chain)]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.re == 0 && self.im == 0
            && self.ph_pos == 0 && self.ph_neg == 0
        {
            write!(f, "0")?;
            return Ok(());
        }
        if self.div2 < 0 {
            write!(f, "{:.0}*(", 2.0_f64.powi(-self.div2))?;
        } else if self.div2 > 0 {
            write!(f, "(")?;
        }
        let mut prev = false;
        if self.re != 0 {
            write!(f, "{}", self.re)?;
            prev = true;
        }
        if prev {
            if self.im > 0 {
                write!(f, " + {}i", self.im)?;
            } else if self.im < 0 {
                write!(f, " - {}i", self.im.abs())?;
            }
        } else if self.im != 0 {
            write!(f, "{}i", self.im)?;
            prev = true;
        }
        if prev {
            if self.ph_pos > 0 {
                write!(f, " + {}*exp(+iπ/4)", self.ph_pos)?;
            } else if self.ph_pos < 0 {
                write!(f, " - {}*exp(+iπ/4)", self.ph_pos.abs())?;
            }
        } else if self.ph_pos != 0 {
            write!(f, "{}*exp(+iπ/4)", self.ph_pos)?;
            prev = true;
        }
        if prev {
            if self.ph_neg > 0 {
                write!(f, " + {}*exp(-iπ/4)", self.ph_neg)?;
            } else if self.ph_neg < 0 {
                write!(f, " - {}*exp(-iπ/4)", self.ph_neg.abs())?;
            }
        } else if self.ph_neg != 0 {
            write!(f, "{}*exp(-iπ/4)", self.ph_neg)?;
        }
        if self.div2 > 0 {
            write!(f, ")/{:.0}", 2.0_f64.powi(self.div2))?;
        } else if self.div2 < 0 {
            write!(f, ")")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn direct_eq(lhs: Complex, rhs: Complex) -> bool {
        lhs.div2 == rhs.div2
            && lhs.re == rhs.re
            && lhs.im == rhs.im
            && lhs.ph_pos == rhs.ph_pos
            && lhs.ph_neg == rhs.ph_neg
    }

    #[test]
    fn reduce() {
        let z = Complex { div2: 0, re: 2, im: 2, ph_pos: -2, ph_neg: 2 };
        let w = Complex { div2: -1, re: 1, im: 1, ph_pos: -1, ph_neg: 1 };
        assert!(direct_eq(z.reduced(), w));

        let z = Complex { div2: 0, re: 2, im: 2, ph_pos: 2, ph_neg: 1 };
        let w = Complex { div2: 0, re: 2, im: 2, ph_pos: 2, ph_neg: 1 };
        assert!(direct_eq(z.reduced(), w));

        let z = Complex { div2: 0, re: 4, im: 0, ph_pos: 0, ph_neg: 4 };
        let w = Complex { div2: -2, re: 1, im: 0, ph_pos: 0, ph_neg: 1 };
        assert!(direct_eq(z.reduced(), w));

        let z = Complex { div2: 0, re: 0, im: 0, ph_pos: 0, ph_neg: 0 };
        let w = Complex { div2: 0, re: 0, im: 0, ph_pos: 0, ph_neg: 0 };
        assert!(direct_eq(z.reduced(), w));

        let z = Complex::new(2, -8, 0, 8, 0);
        let w = Complex { div2: -1, re: -1, im: 0, ph_pos: 1, ph_neg: 0 };
        assert!(direct_eq(z.reduced(), w));

        let z = Complex { div2: 0, re: 4, im: 0, ph_pos: 2, ph_neg: 0 };
        let w = Complex { div2: -1, re: 2, im: 0, ph_pos: 1, ph_neg: 0 };
        assert!(direct_eq(z.reduced(), w));

        let z = Complex { div2: 0, re: 4, im: 0, ph_pos: -3, ph_neg: 0 };
        let w = Complex { div2: 0, re: 4, im: 0, ph_pos: -3, ph_neg: 0 };
        assert!(direct_eq(z.reduced(), w));
    }

    #[test]
    fn add() {
        let z0 = Complex { div2: 0, re: 4, im: 0, ph_pos: 3, ph_neg: 0 };
        let z1 = Complex { div2: 0, re: 0, im: 0, ph_pos: 1, ph_neg: 0 };
        let w = Complex { div2: -2, re: 1, im: 0, ph_pos: 1, ph_neg: 0 };
        assert!(direct_eq(z0 + z1, w));

        let z0 = Complex { div2: 1, re: 1, im: 0, ph_pos: 3, ph_neg: 0 };
        let z1 = 1;
        let w = Complex { div2: 1, re: 3, im: 0, ph_pos: 3, ph_neg: 0 };
        assert!(direct_eq(z0 + z1, w));

        let z0 = Complex { div2: -1, re: 1, im: 0, ph_pos: 3, ph_neg: 0 };
        let z1 = 1;
        let w = Complex { div2: 0, re: 3, im: 0, ph_pos: 6, ph_neg: 0 };
        assert!(direct_eq(z0 + z1, w));
    }

    #[test]
    fn sub() {
        let z0 = Complex { div2: 0, re: 4, im: 0, ph_pos: 3, ph_neg: 0 };
        let z1 = Complex { div2: 0, re: 0, im: 0, ph_pos: 1, ph_neg: 0 };
        let w = Complex { div2: -1, re: 2, im: 0, ph_pos: 1, ph_neg: 0 };
        assert!(direct_eq(z0 - z1, w));

        let z0 = Complex { div2: 1, re: 1, im: 0, ph_pos: 3, ph_neg: 0 };
        let z1 = 1;
        let w = Complex { div2: 1, re: -1, im: 0, ph_pos: 3, ph_neg: 0 };
        assert!(direct_eq(z0 - z1, w));

        let z0 = Complex { div2: -1, re: 1, im: 0, ph_pos: 3, ph_neg: 0 };
        let z1 = 1;
        let w = Complex { div2: 0, re: 1, im: 0, ph_pos: 6, ph_neg: 0 };
        assert!(direct_eq(z0 - z1, w));
    }

    #[test]
    fn mul() {
        let z0 = Complex { div2: 0, re: 1, im: 0, ph_pos: 0, ph_neg: 0 };
        let z1 = Complex { div2: 0, re: 4, im: 0, ph_pos: 0, ph_neg: 0 };
        let w = Complex { div2: -2, re: 1, im: 0, ph_pos: 0, ph_neg: 0 };
        assert!(direct_eq(z0 * z1, w));

        let z0 = Complex { div2: 0, re: 3, im: 0, ph_pos: 0, ph_neg: 0 };
        let z1 = Complex { div2: 0, re: 1, im: 2, ph_pos: 0, ph_neg: 0 };
        let w = Complex { div2: 0, re: 3, im: 6, ph_pos: 0, ph_neg: 0 };
        assert!(direct_eq(z0 * z1, w));

        let z0 = Complex { div2: -1, re: 0, im: 3, ph_pos: 0, ph_neg: 0 };
        let z1 = Complex { div2: 1, re: 1, im: 2, ph_pos: 0, ph_neg: 0 };
        let w = Complex { div2: 0, re: -6, im: 3, ph_pos: 0, ph_neg: 0 };
        assert!(direct_eq(z0 * z1, w));

        let z0 = Complex { div2: 0, re: 0, im: 0, ph_pos: 3, ph_neg: 0 };
        let z1 = Complex { div2: 0, re: 1, im: 0, ph_pos: 0, ph_neg: 2 };
        let w = Complex { div2: 0, re: 6, im: 0, ph_pos: 3, ph_neg: 0 };
        assert!(direct_eq(z0 * z1, w));

        let z0 = Complex { div2: 0, re: 0, im: 0, ph_pos: 0, ph_neg: 3 };
        let z1 = Complex { div2: 0, re: 0, im: 2, ph_pos: 0, ph_neg: 3 };
        let w = Complex { div2: 0, re: 0, im: -9, ph_pos: 6, ph_neg: 0 };
        assert!(direct_eq(z0 * z1, w));

        let z0 = Complex { div2: -1, re: 7, im: -3, ph_pos: 1, ph_neg: 5 };
        let z1 = -6;
        let w = Complex { div2: -1, re: -42, im: 18, ph_pos: -6, ph_neg: -30 };
        assert!(direct_eq(z0 * z1, w));

        let z0 = Complex { div2: 1, re: 7, im: 3, ph_pos: 1, ph_neg: 5 };
        let z1 = 8;
        let w = Complex { div2: -2, re: 7, im: 3, ph_pos: 1, ph_neg: 5 };
        assert!(direct_eq(z0 * z1, w));

        let z0 = Complex { div2: 1, re: 7, im: 3, ph_pos: 1, ph_neg: 5 };
        let z1 = 0;
        let w = Complex { div2: 0, re: 0, im: 0, ph_pos: 0, ph_neg: 0 };
        assert!(direct_eq(z0 * z1, w));
    }

    #[test]
    fn rot() {
        let z = Complex { div2: -5, re: 7, im: 3, ph_pos: -8, ph_neg: 1 };
        assert_eq!(z.rot(0), z * Complex::ONE);
        assert_eq!(z.rot(1), z * Complex::PH_POS);
        assert_eq!(z.rot(2), z * Complex::I);
        assert_eq!(z.rot(3), z * (-Complex::PH_NEG));
        assert_eq!(z.rot(4), -z);
        assert_eq!(z.rot(5), z * (-Complex::PH_POS));
        assert_eq!(z.rot(6), z * (-Complex::I));
        assert_eq!(z.rot(7), z * Complex::PH_NEG);
        assert_eq!(z.rot(8), z);
        assert_eq!(z.rot(-1), z * Complex::PH_NEG);
        assert_eq!(z.rot(-14), z.rot(2));
    }
}

