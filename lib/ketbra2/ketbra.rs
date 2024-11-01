use std::fmt;
use itertools::Itertools;
use num_complex::Complex64 as C64;
use crate::ketbra2::{
    Basis,
    Element,
    KBError::*,
    KBResult,
    State,
    States,
    StatesIter,
    StateOrd
};

/// Represents a single matrix element as a ketbra with an associated
/// amplitude.
#[derive(Clone, Debug)]
pub struct KetBra {
    pub(crate) ampl: C64,
    pub(crate) ket: States,
    pub(crate) bra: States,
}

impl PartialEq for KetBra {
    fn eq(&self, other: &Self) -> bool {
        const EPSILON: f64 = 1e-12;
        (self.ampl - other.ampl).norm() < EPSILON
            && self.ket == other.ket
            && self.bra == other.bra
    }
}

impl StateOrd for KetBra {
    fn state_cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.bra.cmp(&other.bra).then_with(|| self.ket.cmp(&other.ket))
    }
}

impl KetBra {
    /// Create a new `KetBra`.
    pub fn new<I, J>(ampl: C64, ket: I, bra: J) -> Self
    where
        I: IntoIterator<Item = (usize, State)>,
        J: IntoIterator<Item = (usize, State)>,
    {
        let ket: States = ket.into_iter().collect();
        let bra: States = bra.into_iter().collect();
        Self { ampl, ket, bra }
    }

    /// Return the amplitude of the ketbra.
    pub fn ampl(&self) -> C64 { self.ampl }

    /// Return `true` if `self` has all empty ket and bra states.
    pub fn is_scalar(&self) -> bool {
        self.ket.is_empty() && self.bra.is_empty()
    }

    /// If `self` has all empty ket and bra states, return the amplitude as a
    /// single scalar.
    pub fn as_scalar(&self) -> Option<C64> {
        (self.ket.is_empty() && self.bra.is_empty()).then_some(self.ampl)
    }

    /// Return `true` if `self` has all empty bra states and non-empty ket
    /// states.
    pub fn is_ket(&self) -> bool {
        !self.ket.is_empty() && self.bra.is_empty()
    }

    /// Return a reference to the collection of ket states.
    pub fn ket(&self) -> &States { &self.ket }

    pub(crate) fn ket_mut(&mut self) -> &mut States { &mut self.ket }

    /// Return an iterator over all ket states.
    pub fn ket_iter(&self) -> StatesIter<'_> { self.ket.iter() }

    /// Return `true` if `self` has all empty ket states and non-empty bra
    /// states.
    pub fn is_bra(&self) -> bool {
        self.ket.is_empty() && !self.bra.is_empty()
    }

    /// Return a reference to the collection of bra states.
    pub fn bra(&self) -> &States { &self.bra }

    pub(crate) fn bra_mut(&mut self) -> &mut States { &mut self.bra }

    /// Return an iterator over all bra states.
    pub fn bra_iter(&self) -> StatesIter<'_> { self.bra.iter() }

    /// Return `true` if `self` and `other` are identical, up to a scalar.
    pub fn eq_states(&self, other: &Self) -> bool {
        self.ket == other.ket && self.bra == other.bra
    }

    /// Compute the dot-product `self · rhs`.
    ///
    /// Unmatched kets and bras are left untouched and "passed through" the
    /// product. Fails if an unmatched ket or bra leaves multiple states on a
    /// single index.
    pub fn dot(&self, rhs: &Self) -> KBResult<Self> {
        // check for duplicate kets in result
        if let Some(rep) =
            rhs.ket.iter()
            .find_map(|(id, _)| {
                (self.ket.contains_key(id) && !self.bra.contains_key(id))
                    .then_some(id)
            })
        {
            return Err(DotDuplicateKetKey(rep));
        }
        // check for duplicate bras in result
        if let Some(rep) =
            self.bra.iter()
            .find_map(|(id, _)| {
                (rhs.bra.contains_key(id) && !rhs.ket.contains_key(id))
                    .then_some(id)
            })
        {
            return Err(DotDuplicateBraKey(rep));
        }

        let mut ampl = self.ampl * rhs.ampl;
        let mut ket = self.ket.clone();
        let mut bra = rhs.bra.clone();
        let matched: Vec<usize> =
            self.bra.iter()
            .filter_map(|(id, s_l)| {
                if let Some(s_r) = rhs.ket.get(id) {
                    ampl *= s_l.dot(s_r);
                    Some(id)
                } else {
                    bra.insert(id, *s_l);
                    None
                }
            })
            .collect();
        rhs.ket.iter()
            .for_each(|(id, s_r)| {
                if !matched.contains(&id) { ket.insert(id, *s_r); }
            });
        Ok(Self { ampl, ket, bra })
    }

    /// Compute the dot-product `self · rhs`, consuming both.
    ///
    /// Unmatched kets and bras are left untouched and "passed through" the
    /// product. Fails if an unmatched ket or bra leaves multiple states on a
    /// single index.
    pub fn into_dot(self, rhs: Self) -> KBResult<Self> {
        // check for duplicate kets in result
        for (id, _) in rhs.ket.iter() {
            if self.ket.contains_key(id) && !self.bra.contains_key(id) {
                return Err(DotDuplicateKetKey(id));
            }
        }
        // check for duplicate bras in result
        for (id, _) in self.bra.iter() {
            if rhs.bra.contains_key(id) && !rhs.ket.contains_key(id) {
                return Err(DotDuplicateBraKey(id));
            }
        }
        let mut ampl = self.ampl * rhs.ampl;
        let mut ket = self.ket;
        let mut bra = rhs.bra;
        let matched: Vec<usize> =
            self.bra.into_iter()
            .filter_map(|(id, s_l)| {
                if let Some(s_r) = rhs.ket.get(id) {
                    ampl *= s_l.dot(s_r);
                    Some(id)
                } else {
                    bra.insert(id, s_l);
                    None
                }
            })
            .collect();
        rhs.ket.into_iter()
            .for_each(|(id, s_r)| {
                if !matched.contains(&id) { ket.insert(id, s_r); }
            });
        Ok(Self { ampl, ket, bra })
    }

    /// Change the basis of operation in place by applying a Hadamard transform
    /// to each qubit index, swapping `0 ↔ +` and `1 ↔ -`.
    pub fn hadamard_mut(&mut self) {
        self.ket.iter_mut().for_each(|(_, s)| { s.hadamard_mut(); });
        self.bra.iter_mut().for_each(|(_, s)| { s.hadamard_mut(); });
    }

    /// Change the basis of operation by applying a Hadamard transform to each
    /// qubit index, swapping `0 ↔ +` and `1 ↔ -`.
    pub fn into_hadamard(mut self) -> Self {
        self.hadamard_mut();
        self
    }

    /// Change the basis of operation by applying a Hadamard transform to each
    /// qubit index, swapping `0 ↔ +` and `1 ↔ -`. Returns a copy of `self` with
    /// the transformation applied.
    pub fn hadamard(&self) -> Self {
        let mut new = self.clone();
        new.hadamard_mut();
        new
    }

    /// Conjugate `self` in place, swapping all kets and bras, and conjugating
    /// the amplitude.
    pub fn adjoint_mut(&mut self) {
        std::mem::swap(&mut self.ket, &mut self.bra);
        self.ampl = self.ampl.conj();
    }

    /// Conjugate `self`, swapping all kets and bras, and conjugating the
    /// amplitude.
    pub fn into_adjoint(mut self) -> Self {
        self.adjoint_mut();
        self
    }

    /// Return the complex conjugate of `self`, with all kets swapped with bras,
    /// and the amplitude conjugated.
    pub fn adjoint(&self) -> Self {
        let mut new = self.clone();
        new.adjoint_mut();
        new
    }

    /// Express `self` in `basis`, consuming `self`.
    ///
    /// This returns an [`Element`] because conversion to a different basis
    /// typically introduces factors expressed as sums, resulting in the need
    /// for multiple ket-bras in the final expression. In particular, the
    /// returned `Element` will always be [non-atomic][Element::is_atomic].
    pub fn into_basis(self, basis: Basis) -> Element {
        if self.ket.iter().all(|(_, s)| s.is_basis(basis))
            && self.bra.iter().all(|(_, s)| s.is_basis(basis))
        {
            self.into()
        } else {
            let Self { ampl, ket, bra } = self;
            let ket_len = ket.len();
            let ket_bounds = ket.idx_bounds();
            let bra_bounds = bra.idx_bounds();
            let ket_iter =
                ket.into_iter()
                .map(|(id, s)| {
                    if let Some(((a_up, s_up), (a_dn, s_dn))) =
                        s.decomp(basis)
                    {
                        DecompIter::Two {
                            yielded: None,
                            item0: (id, a_up, s_up),
                            item1: (id, a_dn, s_dn),
                        }
                    } else {
                        DecompIter::One {
                            yielded: false,
                            item: (id, C64::from(1.0), s),
                        }
                    }
                });
            let bra_iter =
                bra.into_iter()
                .map(|(id, s)| {
                    if let Some(((a_up, s_up), (a_dn, s_dn))) =
                        s.decomp(basis)
                    {
                        DecompIter::Two {
                            yielded: None,
                            item0: (id, a_up, s_up),
                            item1: (id, a_dn, s_dn),
                        }
                    } else {
                        DecompIter::One {
                            yielded: false,
                            item: (id, C64::from(1.0), s),
                        }
                    }
                });
            let terms =
                ket_iter.chain(bra_iter)
                .multi_cartesian_product()
                .map(|mut states| {
                    let mut term_ampl = ampl;
                    let mut ket =
                        ket_bounds
                        .map(|(min, max)| States::with_capacity(min, max))
                        .unwrap_or_default();
                    let mut bra =
                        bra_bounds
                        .map(|(min, max)| States::with_capacity(min, max))
                        .unwrap_or_default();
                    states.drain(..ket_len)
                        .for_each(|(id, a, s)| {
                            println!("ket: {} {} {}", id, s, a);
                            term_ampl *= a;
                            ket.insert(id, s);
                        });
                    states.drain(..)
                        .for_each(|(id, a, s)| {
                            println!("bra: {} {} {}", id, s, a);
                            term_ampl *= a;
                            bra.insert(id, s);
                        });
                    KetBra { ampl: term_ampl, ket, bra }
                });
            Element::from_terms(terms)
        }
    }

}

fn subscript_str(n: usize) -> String {
    format!("{}", n).chars()
        .map(|c| {
            match c {
                '0' => '₀',
                '1' => '₁',
                '2' => '₂',
                '3' => '₃',
                '4' => '₄',
                '5' => '₅',
                '6' => '₆',
                '7' => '₇',
                '8' => '₈',
                '9' => '₉',
                _ => unreachable!(),
            }
        })
        .collect()
}

impl fmt::Display for KetBra {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "(")?;
        self.ampl.fmt(f)?;
        write!(f, ")")?;
        if !self.ket.is_empty() {
            write!(f, "∣")?;
            for (idx, s) in self.ket.iter() {
                write!(f, "{}{}", s, subscript_str(idx))?;
            }
            write!(f, "⟩")?;
        }
        if !self.bra.is_empty() {
            write!(f, "⟨")?;
            for (idx, s) in self.bra.iter() {
                write!(f, "{}{}", s, subscript_str(idx))?;
            }
            write!(f, "∣")?;
        }
        Ok(())
    }
}

// I just really don't want to allocate a Vec
#[derive(Copy, Clone, Debug)]
enum DecompIter {
    One {
        yielded: bool,
        item: (usize, C64, State),
    },
    Two {
        yielded: Option<bool>, // ≅ 1 + 2 = 3
        item0: (usize, C64, State),
        item1: (usize, C64, State),
    },
}

impl Iterator for DecompIter {
    type Item = (usize, C64, State);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::One { yielded, item } => {
                if *yielded {
                    None
                } else {
                    *yielded = true;
                    Some(*item)
                }
            },
            Self::Two { yielded, item0, item1 } => {
                if yielded.is_some_and(|b| b) {
                    None
                } else if let Some(b) = yielded.as_mut() {
                    *b = true;
                    Some(*item1)
                } else {
                    *yielded = Some(false);
                    Some(*item0)
                }
            },
        }
    }
}

impl std::ops::MulAssign<C64> for KetBra {
    fn mul_assign(&mut self, z: C64) { self.ampl *= z; }
}

impl std::ops::Mul<C64> for KetBra {
    type Output = KetBra;

    fn mul(mut self, z: C64) -> Self::Output {
        self *= z;
        self
    }
}

impl std::ops::Mul<KetBra> for C64 {
    type Output = KetBra;

    fn mul(self, mut rhs: KetBra) -> Self::Output {
        rhs *= self;
        rhs
    }
}

impl std::ops::MulAssign<f64> for KetBra {
    fn mul_assign(&mut self, z: f64) { self.ampl *= z; }
}

impl std::ops::Mul<f64> for KetBra {
    type Output = KetBra;

    fn mul(mut self, z: f64) -> Self::Output {
        self *= z;
        self
    }
}

impl std::ops::Mul<KetBra> for f64 {
    type Output = KetBra;

    fn mul(self, mut rhs: KetBra) -> Self::Output {
        rhs *= self;
        rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ketbra2::Basis;
    use crate::ketbra2::State;

    #[test]
    fn scalars() {
        let scalar = KetBra::new(C64::from(1.0), [], []);
        let nonscalar =
            KetBra::new(C64::from(1.0), [(0, State::Zero)], [(0, State::Plus)]);
        assert!(scalar.is_scalar());
        assert!(!nonscalar.is_scalar());
        assert_eq!(scalar.as_scalar(), Some(C64::from(1.0)));
        assert_eq!(nonscalar.as_scalar(), None);
    }

    #[test]
    fn ket() {
        let scalar = KetBra::new(C64::from(1.0), [], []);
        let ket = KetBra::new(C64::from(1.0), [(0, State::Zero)], []);
        let nonket =
            KetBra::new(C64::from(1.0), [(0, State::Zero)], [(0, State::Plus)]);
        assert!(!scalar.is_ket());
        assert!(ket.is_ket());
        assert!(!nonket.is_ket());
    }

    #[test]
    fn bra() {
        let scalar = KetBra::new(C64::from(1.0), [], []);
        let bra = KetBra::new(C64::from(1.0), [], [(0, State::Plus)]);
        let nonbra =
            KetBra::new(C64::from(1.0), [(0, State::Zero)], [(0, State::Plus)]);
        assert!(!scalar.is_bra());
        assert!(bra.is_bra());
        assert!(!nonbra.is_bra());
    }

    #[test]
    fn dot() {
        let onrt2 = C64::from(std::f64::consts::FRAC_1_SQRT_2);

        let s0 = KetBra::new(C64::from(1.0), [], [(0, State::Zero)]);
        let s1 = KetBra::new(C64::i(), [(0, State::Minus)], []);
        let dot01 = s0.dot(&s1);
        let dot01_expected =
            KetBra::new(C64::i() * onrt2, [], []);
        let dot10 = s1.dot(&s0);
        let dot10_expected =
            KetBra::new(C64::i(), [(0, State::Minus)], [(0, State::Zero)]);
        assert!(dot01.is_ok());
        assert_eq!(dot01.unwrap(), dot01_expected);
        assert!(dot10.is_ok());
        assert_eq!(dot10.unwrap(), dot10_expected);

        let s2 = KetBra::new(C64::from(1.0), [], [(0, State::Zero)]);
        let s3 = KetBra::new(C64::from(1.0), [], [(1, State::One)]);
        let s4 = KetBra::new(C64::from(1.0), [], [(0, State::Plus)]);
        let dot23 = s2.clone().into_dot(s3);
        let dot23_expected =
            KetBra::new(C64::from(1.0), [], [(0, State::Zero), (1, State::One)]);
        assert!(dot23.is_ok());
        assert_eq!(dot23.unwrap(), dot23_expected);
        assert!(s2.dot(&s4).is_err());
    }

    #[test]
    fn hadamard() {
        let mut s = KetBra::new(
            C64::i(),
            [(0, State::Zero), (2, State::Minus)],
            [(1, State::Plus), (3, State::One)],
        );
        s.hadamard_mut();
        let s_expected = KetBra::new(
            C64::i(),
            [(0, State::Plus), (2, State::One)],
            [(1, State::Zero), (3, State::Minus)],
        );
        assert_eq!(s, s_expected);
    }

    #[test]
    fn adjoint() {
        let mut s = KetBra::new(
            C64::i(),
            [(0, State::Zero), (2, State::Minus)],
            [(1, State::Plus), (3, State::One)],
        );
        s.adjoint_mut();
        let s_expected = KetBra::new(
            -C64::i(),
            [(1, State::Plus), (3, State::One)],
            [(0, State::Zero), (2, State::Minus)],
        );
        assert_eq!(s, s_expected);

        let s0 = KetBra::new(
            C64::i(),
            [],
            [(1, State::Plus), (3, State::One)],
        );
        let mut s1 = s0.adjoint();
        let s1_expected = KetBra::new(
            -C64::i(),
            [(1, State::Plus), (3, State::One)],
            [],
        );
        assert_eq!(s1, s1_expected);
        s1.adjoint_mut();
        assert_eq!(s0, s1);
    }

    #[test]
    fn into_basis() {
        use State::*;
        let onrt2 = C64::from(std::f64::consts::FRAC_1_SQRT_2);

        let s = KetBra::new(C64::i(), [(0, Plus), (1, Minus)], []);
        let elem = s.into_basis(Basis::Z);
        let mb_terms = elem.terms();
        assert!(mb_terms.is_some());
        let terms = mb_terms.unwrap();
        let terms_expected: Vec<KetBra> =
            vec![
                KetBra::new( C64::i() / 2.0, [(0, Zero), (1, Zero)], []),
                KetBra::new(-C64::i() / 2.0, [(0, Zero), (1, One )], []),
                KetBra::new( C64::i() / 2.0, [(0, One ), (1, Zero)], []),
                KetBra::new(-C64::i() / 2.0, [(0, One ), (1, One )], []),
            ];
        assert_eq!(terms.len(), terms_expected.len());
        assert!(terms.iter().all(|kb| terms_expected.contains(kb)));

        let s = KetBra::new(C64::from(1.0), [(0, One)], [(0, Plus)]);
        let elem = s.into_basis(Basis::X);
        let terms = elem.terms().unwrap();
        let terms_expected: Vec<KetBra> =
            vec![
                KetBra::new( onrt2, [(0, Plus )], [(0, Plus)]),
                KetBra::new(-onrt2, [(0, Minus)], [(0, Plus)]),
            ];
        assert_eq!(terms.len(), terms_expected.len());
        assert!(terms.iter().all(|kb| terms_expected.contains(kb)));
    }

}

