use std::fmt;
use itertools::Itertools;
use num_complex::Complex64 as C64;
use crate::{
    ketbra2::{ Basis, KBError, KBResult, KetBra, State, StateOrd },
    phase::Phase,
};
use KBError::*;

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum ElementData {
    Z { ket: Vec<usize>, bra: Vec<usize>, ph: Phase },
    X { ket: Vec<usize>, bra: Vec<usize>, ph: Phase },
    H { ket: Vec<usize>, bra: Vec<usize>, arg: C64 },
    Swap { a: usize, b: usize },
    Cup { a: usize, b: usize },
    Cap { a: usize, b: usize },
    Terms { terms: Vec<KetBra> },
}

impl ElementData {
    pub(crate) fn is_z(&self) -> bool { matches!(self, Self::Z { .. }) }

    pub(crate) fn is_x(&self) -> bool { matches!(self, Self::X { .. }) }

    pub(crate) fn is_h(&self) -> bool { matches!(self, Self::H { .. }) }

    pub(crate) fn is_cup(&self) -> bool { matches!(self, Self::Cup { .. }) }

    pub(crate) fn is_cap(&self) -> bool { matches!(self, Self::Cap { .. }) }

    pub(crate) fn is_spider(&self) -> bool {
        matches!(self, Self::Z { .. } | Self::X { .. })
    }

    pub(crate) fn is_atomic(&self) -> bool {
        !matches!(self, Self::Terms { .. })
    }

    pub(crate) fn is_scalar(&self) -> bool {
        match self {
            Self::Z { ket, bra, ph: _ } => ket.is_empty() && bra.is_empty(),
            Self::X { ket, bra, ph: _ } => ket.is_empty() && bra.is_empty(),
            Self::H { ket, bra, arg: _ } => ket.is_empty() && bra.is_empty(),
            Self::Swap { .. } => false,
            Self::Cup { .. } => false,
            Self::Cap { .. } => false,
            Self::Terms { terms } => terms.iter().all(|kb| kb.is_scalar()),
        }
    }

    pub(crate) fn as_scalar(&self) -> Option<C64> {
        match self {
            Self::Z { ket, bra, ph } if ket.is_empty() && bra.is_empty() => {
                Some(ph.cis() + 1.0)
            },
            Self::X { ket, bra, ph } if ket.is_empty() && bra.is_empty() => {
                Some(ph.cis() + 1.0)
            },
            Self::H { ket, bra, arg } if ket.is_empty() && bra.is_empty() => {
                Some(*arg)
            },
            Self::Terms { terms } => {
                terms.iter()
                    .try_fold(C64::from(0.0), |mut acc, kb| {
                        kb.as_scalar().map(|a| { acc += a; acc })
                    })
            },
            _ => None,
        }
    }

    pub(crate) fn z<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let mut ket: Vec<usize> = Vec::new();
        outs.into_iter()
            .for_each(|id| { if !ket.contains(&id) { ket.push(id); } });
        ket.sort_unstable();
        let mut bra: Vec<usize> = Vec::new();
        ins.into_iter()
            .for_each(|id| { if !bra.contains(&id) { bra.push(id); } });
        bra.sort_unstable();
        Self::Z { ket, bra, ph: phase.unwrap_or(Phase::zero()) }
    }

    pub(crate) fn x<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let mut ket: Vec<usize> = Vec::new();
        outs.into_iter()
            .for_each(|id| { if !ket.contains(&id) { ket.push(id); } });
        ket.sort_unstable();
        let mut bra: Vec<usize> = Vec::new();
        ins.into_iter()
            .for_each(|id| { if !bra.contains(&id) { bra.push(id); } });
        bra.sort_unstable();
        Self::X { ket, bra, ph: phase.unwrap_or(Phase::zero()) }
    }

    pub(crate) fn h<I, J>(ins: I, outs: J, arg: Option<C64>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let mut ket: Vec<usize> = Vec::new();
        outs.into_iter()
            .for_each(|id| { if !ket.contains(&id) { ket.push(id); } });
        ket.sort_unstable();
        let mut bra: Vec<usize> = Vec::new();
        ins.into_iter()
            .for_each(|id| { if !bra.contains(&id) { bra.push(id); } });
        bra.sort_unstable();
        Self::H { ket, bra, arg: arg.unwrap_or((-1.0).into()) }
    }

    pub(crate) fn cup(a: usize, b: usize) -> Self {
        (a != b).then_some(())
            .expect("Element::cup: output wires cannot be equal");
        Self::Cup { a, b }
    }

    pub(crate) fn cap(a: usize, b: usize) -> Self {
        (a != b).then_some(())
            .expect("Element::cap: input wires cannot be equal");
        Self::Cap { a, b }
    }

    pub(crate) fn swap(a: usize, b: usize) -> Self {
        (a != b).then_some(())
            .expect("Element::swap: wires cannot be equal");
        Self::Swap { a, b }
    }

    pub(crate) fn from_terms<I>(terms: I) -> Self
    where I: IntoIterator<Item = KetBra>
    {
        Self::Terms { terms: terms.into_iter().collect() }
    }

    pub(crate) fn terms(&self) -> Option<&Vec<KetBra>> {
        match &self {
            Self::Terms { terms } => Some(terms),
            _ => None,
        }
    }

    pub(crate) fn terms_mut(&mut self) -> Option<&mut Vec<KetBra>> {
        match self {
            Self::Terms { terms } => Some(terms),
            _ => None,
        }
    }

    pub(crate) fn contains_ket(&self, k: usize) -> bool {
        match self {
            Self::Z { ket, bra: _, ph: _ } => ket.contains(&k),
            Self::X { ket, bra: _, ph: _ } => ket.contains(&k),
            Self::H { ket, bra: _, arg: _ } => ket.contains(&k),
            Self::Swap { a, b } => *a == k || *b == k,
            Self::Cup { a, b } => *a == k || *b == k,
            Self::Cap { .. } => false,
            Self::Terms { terms } => {
                terms.iter().any(|kb| kb.ket().contains_key(k))
            },
        }
    }

    pub(crate) fn contains_bra(&self, k: usize) -> bool {
        match self {
            Self::Z { ket: _, bra, ph: _ } => bra.contains(&k),
            Self::X { ket: _, bra, ph: _ } => bra.contains(&k),
            Self::H { ket: _, bra, arg: _ } => bra.contains(&k),
            Self::Swap { a, b } => *a == k || *b == k,
            Self::Cup { .. } => false,
            Self::Cap { a, b } => *a == k || *b == k,
            Self::Terms { terms } => {
                terms.iter().any(|kb| kb.bra().contains_key(k))
            },
        }
    }

    pub(crate) fn ins(&self) -> Vec<usize> {
        match self {
            Self::Z { ket: _, bra, ph: _ } => bra.clone(),
            Self::X { ket: _, bra, ph: _ } => bra.clone(),
            Self::H { ket: _, bra, arg: _ } => bra.clone(),
            Self::Swap { a, b } => vec![*a, *b],
            Self::Cup { .. } => Vec::new(),
            Self::Cap { a, b } => vec![*a, *b],
            Self::Terms { terms } => {
                let mut acc: Vec<usize> = Vec::new();
                terms.iter()
                    .for_each(|kb| {
                        kb.bra_iter().for_each(|(k, _)| {
                            if !acc.contains(&k) { acc.push(k); }
                        });
                    });
                acc
            },
        }
    }

    pub(crate) fn outs(&self) -> Vec<usize> {
        match self {
            Self::Z { ket, bra: _, ph: _ } => ket.clone(),
            Self::X { ket, bra: _, ph: _ } => ket.clone(),
            Self::H { ket, bra: _, arg: _ } => ket.clone(),
            Self::Swap { a, b } => vec![*a, *b],
            Self::Cup { a, b } => vec![*a, *b],
            Self::Cap { .. } => Vec::new(),
            Self::Terms { terms } => {
                let mut acc: Vec<usize> = Vec::new();
                terms.iter()
                    .for_each(|kb| {
                        kb.ket_iter().for_each(|(k, _)| {
                            if !acc.contains(&k) { acc.push(k); }
                        });
                    });
                acc
            },
        }
    }

    #[allow(clippy::uninit_assumed_init, invalid_value)]
    pub(crate) fn make_terms(&mut self) {
        let temp =
            unsafe { std::mem::MaybeUninit::uninit().assume_init() };
        let temp = std::mem::replace(self, temp);
        let terms: Vec<KetBra> =
            match temp {
                Self::Z { ket, bra, ph } => vec![
                    KetBra::new(
                        C64::from(1.0),
                        ket.iter().map(|id| (*id, State::Zero)),
                        bra.iter().map(|id| (*id, State::Zero)),
                    ),
                    KetBra::new(
                        ph.cis(),
                        ket.into_iter().map(|id| (id, State::One)),
                        bra.into_iter().map(|id| (id, State::One)),
                    ),
                ],
                Self::X { ket, bra, ph } => vec![
                    KetBra::new(
                        C64::from(1.0),
                        ket.iter().map(|id| (*id, State::Plus)),
                        bra.iter().map(|id| (*id, State::Plus)),
                    ),
                    KetBra::new(
                        ph.cis(),
                        ket.into_iter().map(|id| (id, State::Minus)),
                        bra.into_iter().map(|id| (id, State::Minus)),
                    ),
                ],
                Self::H { ket, bra, arg } => {
                    let ket_len = ket.len();
                    ket.into_iter().chain(bra)
                        .map(|id| [(id, State::Zero), (id, State::One)])
                        .multi_cartesian_product()
                        .map(|states| {
                            let ampl =
                                if states.iter().all(|ks| ks.1 == State::One) {
                                    arg
                                } else {
                                    C64::from(1.0)
                                };
                            KetBra::new(
                                ampl,
                                states.iter().take(ket_len).copied(),
                                states.iter().skip(ket_len).copied(),
                            )
                        })
                        .collect()
                },
                Self::Swap { a, b } => vec![
                    KetBra::new(
                        C64::from(1.0),
                        [(a, State::Zero), (b, State::Zero)],
                        [(a, State::Zero), (b, State::Zero)],
                    ),
                    KetBra::new(
                        C64::from(1.0),
                        [(a, State::Zero), (b, State::One )],
                        [(a, State::One ), (b, State::Zero)],
                    ),
                    KetBra::new(
                        C64::from(1.0),
                        [(a, State::One ), (b, State::Zero)],
                        [(a, State::Zero), (b, State::One )],
                    ),
                    KetBra::new(
                        C64::from(1.0),
                        [(a, State::One ), (b, State::One )],
                        [(a, State::One ), (b, State::One )],
                    ),
                ],
                Self::Cup { a, b } => vec![
                    KetBra::new(
                        C64::from(1.0),
                        [(a, State::Zero), (b, State::Zero)],
                        [],
                    ),
                    KetBra::new(
                        C64::from(1.0),
                        [(a, State::One ), (b, State::One )],
                        [],
                    ),
                ],
                Self::Cap { a, b } => vec![
                    KetBra::new(
                        C64::from(1.0),
                        [],
                        [(a, State::Zero), (b, State::Zero)],
                    ),
                    KetBra::new(
                        C64::from(1.0),
                        [],
                        [(a, State::One ), (b, State::One )],
                    ),
                ],
                Self::Terms { terms } => terms,
            };
        let new = Self::Terms { terms };
        let temp = std::mem::replace(self, new);
        std::mem::forget(temp);
    }

    pub(crate) fn as_terms(&self) -> Vec<KetBra> {
        let mut new = self.clone();
        new.make_terms();
        let Self::Terms { terms } = new else { unreachable!() };
        terms
    }

    pub(crate) fn into_terms(mut self) -> Vec<KetBra> {
        self.make_terms();
        let Self::Terms { terms } = self else { unreachable!() };
        terms
    }

    /// Combine like terms in place, if `self` is non-atomic.
    pub(crate) fn simplify(&mut self) {
        const EPSILON: f64 = 1e-12;
        if let Some(terms) = self.terms_mut() {
            let n = terms.len();
            let mut remove: Vec<bool> = vec![false; n];
            let mut term0: &KetBra;
            let mut da: C64;
            for k in 0..n {
                if remove[k] { continue; }
                term0 = &terms[k];
                da =
                    terms.iter().zip(remove.iter_mut())
                    .skip(k + 1)
                    .filter_map(|(term1, v)| {
                        term0.eq_states(term1)
                            .then(|| {
                                *v = true;
                                term1.ampl
                            })
                    })
                    .sum();
                terms[k].ampl += da;
                remove[k] = terms[k].ampl.norm() < EPSILON;
            }
            remove.into_iter().enumerate().rev()
                .for_each(|(k, r)| { if r { terms.swap_remove(k); } });
        }
    }

    pub(crate) fn iso(&self, other: &Self) -> bool {
        const EPSILON: f64 = 1e-12;
        match (&self, &other) {
            (
                Self::Z { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::Z { ket: ket_r, bra: bra_r, ph: ph_r },
            ) => {
                ph_l == ph_r
                    && ket_l.len() == ket_r.len()
                    && ket_l.iter().all(|id| ket_r.contains(id))
                    && bra_l.len() == bra_r.len()
                    && bra_l.iter().all(|id| bra_r.contains(id))
            },
            (
                Self::X { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::X { ket: ket_r, bra: bra_r, ph: ph_r },
            ) => {
                ph_l == ph_r
                    && ket_l.len() == ket_r.len()
                    && ket_l.iter().all(|id| ket_r.contains(id))
                    && bra_l.len() == bra_r.len()
                    && bra_l.iter().all(|id| bra_r.contains(id))
            },
            (
                Self::Z { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::X { ket: ket_r, bra: bra_r, ph: ph_r },
            ) => {
                (
                    ket_l.len() == 1 && ket_l == ket_r
                    && bra_l.len() == 1 && bra_l == bra_r
                    && *ph_l == Phase::zero()
                    && *ph_r == Phase::zero()
                ) || (
                    ket_l.len() == 1 && bra_l.is_empty()
                    && ket_r.len() == 1 && bra_r.is_empty()
                    && ket_l == ket_r
                    && (
                        (*ph_l == Phase::pi2() && *ph_r == -Phase::pi2())
                        || (*ph_l == -Phase::pi2() && *ph_r == Phase::pi2())
                    )
                ) || (
                    ket_l.is_empty() && bra_l.len() == 1
                    && ket_r.is_empty() && bra_r.len() == 1
                    && bra_l == bra_r
                    && (
                        (*ph_l == Phase::pi2() && *ph_r == -Phase::pi2())
                        || (*ph_l == -Phase::pi2() && *ph_r == Phase::pi2())
                    )
                ) || (
                    ket_l.is_empty() && bra_l.is_empty()
                    && ket_r.is_empty() && bra_r.is_empty()
                    && *ph_l == *ph_r
                )
            },
            (
                Self::X { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::Z { ket: ket_r, bra: bra_r, ph: ph_r },
            ) => {
                (
                    ket_l.len() == 1 && ket_l == ket_r
                    && bra_l.len() == 1 && bra_l == bra_r
                    && *ph_l == Phase::zero()
                    && *ph_r == Phase::zero()
                ) || (
                    ket_l.len() == 1 && bra_l.is_empty()
                    && ket_r.len() == 1 && bra_r.is_empty()
                    && ket_l == ket_r
                    && (
                        (*ph_l == Phase::pi2() && *ph_r == -Phase::pi2())
                        || (*ph_l == -Phase::pi2() && *ph_r == Phase::pi2())
                    )
                ) || (
                    ket_l.is_empty() && bra_l.len() == 1
                    && ket_r.is_empty() && bra_r.len() == 1
                    && bra_l == bra_r
                    && (
                        (*ph_l == Phase::pi2() && *ph_r == -Phase::pi2())
                        || (*ph_l == -Phase::pi2() && *ph_r == Phase::pi2())
                    )
                ) || (
                    ket_l.is_empty() && bra_l.is_empty()
                    && ket_r.is_empty() && bra_r.is_empty()
                    && *ph_l == *ph_r
                )
            },
            (
                Self::H { ket: ket_l, bra: bra_l, arg: arg_l },
                Self::H { ket: ket_r, bra: bra_r, arg: arg_r },
            ) => {
                ket_l == ket_r && bra_l == bra_r
                    && (arg_l - arg_r).norm() < EPSILON
            },
            (
                Self::H { ket: ket_l, bra: bra_l, arg: arg_l },
                Self::Z { ket: ket_r, bra: bra_r, ph: ph_r },
            ) => {
                (
                    (
                        ket_l.len() == 1 && bra_l.is_empty()
                        && ket_r.len() == 1 && bra_r.is_empty()
                        && ket_l == ket_r
                    ) || (
                        ket_l.is_empty() && bra_l.len() == 1
                        && ket_r.is_empty() && bra_r.len() == 1
                        && bra_l == bra_r
                    ) || (
                        ket_l.is_empty() && bra_l.is_empty()
                        && ket_r.is_empty() && bra_r.is_empty()
                    )
                ) && (arg_l - ph_r.cis()).norm() < EPSILON
            },
            (
                Self::Z { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::H { ket: ket_r, bra: bra_r, arg: arg_r },
            ) => {
                (
                    (
                        ket_l.len() == 1 && bra_l.is_empty()
                        && ket_r.len() == 1 && bra_r.is_empty()
                        && ket_l == ket_r
                    ) || (
                        ket_l.is_empty() && bra_l.len() == 1
                        && ket_r.is_empty() && bra_r.len() == 1
                        && bra_l == bra_r
                    ) || (
                        ket_l.is_empty() && bra_l.is_empty()
                        && ket_r.is_empty() && bra_r.is_empty()
                    )
                ) && (ph_l.cis() - arg_r).norm() < EPSILON
            },
            (
                Self::H { ket: ket_l, bra: bra_l, arg: arg_l },
                Self::X { ket: ket_r, bra: bra_r, ph: ph_r },
            ) => {
                ket_l.is_empty() && bra_l.is_empty()
                    && ket_r.is_empty() && bra_r.is_empty()
                    && (arg_l - ph_r.cis()).norm() < EPSILON
            },
            (
                Self::X { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::H { ket: ket_r, bra: bra_r, arg: arg_r },
            ) => {
                ket_l.is_empty() && bra_l.is_empty()
                    && ket_r.is_empty() && bra_r.is_empty()
                    && (ph_l.cis() - arg_r).norm() < EPSILON
            },
            (
                Self::Swap { a: a_l, b: b_l },
                Self::Swap { a: a_r, b: b_r },
            ) => (a_l == a_r && b_l == b_r) || (a_l == b_r && b_l == a_r),
            (
                Self::Cup { a: a_l, b: b_l },
                Self::Cup { a: a_r, b: b_r },
            ) => (a_l == a_r && b_l == b_r) || (a_l == b_r && b_l == a_r),
            (
                Self::Cup { a: a_l, b: b_l },
                Self::Z { ket: ket_r, bra: bra_r, ph: ph_r },
            ) => {
                *ph_r == Phase::zero()
                    && bra_r.is_empty()
                    && ket_r.len() == 2
                    && ket_r.contains(a_l) && ket_r.contains(b_l)
            },
            (
                Self::Z { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::Cup { a: a_r, b: b_r },
            ) => {
                *ph_l == Phase::zero()
                    && bra_l.is_empty()
                    && ket_l.len() == 2
                    && ket_l.contains(a_r) && ket_l.contains(b_r)
            },
            (
                Self::Cup { a: a_l, b: b_l },
                Self::X { ket: ket_r, bra: bra_r, ph: ph_r },
            ) => {
                *ph_r == Phase::zero()
                    && bra_r.is_empty()
                    && ket_r.len() == 2
                    && ket_r.contains(a_l) && ket_r.contains(b_l)
            },
            (
                Self::X { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::Cup { a: a_r, b: b_r },
            ) => {
                *ph_l == Phase::zero()
                    && bra_l.is_empty()
                    && ket_l.len() == 2
                    && ket_l.contains(a_r) && ket_l.contains(b_r)
            },
            (
                Self::Cap { a: a_l, b: b_l },
                Self::Cap { a: a_r, b: b_r },
            ) => (a_l == a_r && b_l == b_r) || (a_l == b_r && b_l == a_r),
            (
                Self::Cap { a: a_l, b: b_l },
                Self::Z { ket: ket_r, bra: bra_r, ph: ph_r },
            ) => {
                *ph_r == Phase::zero()
                    && ket_r.is_empty()
                    && bra_r.len() == 2
                    && bra_r.contains(a_l) && bra_r.contains(b_l)
            },
            (
                Self::Z { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::Cap { a: a_r, b: b_r },
            ) => {
                *ph_l == Phase::zero()
                    && ket_l.is_empty()
                    && bra_l.len() == 2
                    && bra_l.contains(a_r) && bra_l.contains(b_r)
            },
            (
                Self::Cap { a: a_l, b: b_l },
                Self::X { ket: ket_r, bra: bra_r, ph: ph_r },
            ) => {
                *ph_r == Phase::zero()
                    && ket_r.is_empty()
                    && bra_r.len() == 2
                    && bra_r.contains(a_l) && bra_r.contains(b_l)
            },
            (
                Self::X { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::Cap { a: a_r, b: b_r },
            ) => {
                *ph_l == Phase::zero()
                    && ket_l.is_empty()
                    && bra_l.len() == 2
                    && bra_l.contains(a_r) && bra_l.contains(b_r)
            },
            (
                Self::Terms { terms: terms_l },
                Self::Terms { terms: terms_r },
            ) => {
                terms_l.len() == terms_r.len()
                    && terms_l.iter().all(|term_l| terms_r.contains(term_l))
            },
            _ => false,
        }
    }

    pub(crate) fn iso_mat(&self, other: &Self) -> bool {
        const EPSILON: f64 = 1e-12;
        if self.iso(other) { return true; }
        if let (Some(l), Some(r)) = (self.as_scalar(), other.as_scalar()) {
            return (l.norm() < EPSILON && r.norm() < EPSILON)
                || (l.norm() >= EPSILON && r.norm() >= EPSILON)
        }

        todo!()
    }

    pub(crate) fn dot(&self, rhs: &Self) -> KBResult<Self> {
        const EPSILON: f64 = 1e-12;
        match (&self, &rhs) {
            (
                Self::Z { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::Z { ket: ket_r, bra: bra_r, ph: ph_r },
            ) => {
                let mut common: Vec<usize> =
                    Vec::with_capacity(bra_l.len().max(ket_r.len()));
                // check for duplicate kets in result
                for id in ket_r.iter() {
                    if bra_l.contains(id) {
                        common.push(*id);
                    } else if ket_l.contains(id) {
                        return Err(DotDuplicateKetKey(*id));
                    }
                }
                // check for duplicate bras in result
                for id in bra_l.iter() {
                    if !ket_r.contains(id) && bra_r.contains(id) {
                        return Err(DotDuplicateBraKey(*id));
                    }
                }
                // do the dot product
                if common.is_empty() {
                    let mut new =
                        Self::from_terms(
                            self.as_terms().into_iter()
                                .cartesian_product(rhs.as_terms())
                                .map(|(l, r)| l.into_dot(r).unwrap())
                        );
                    new.simplify();
                    Ok(new)
                } else {
                    let mut ket = ket_l.clone();
                    let mut bra = bra_r.clone();
                    ket_r.iter()
                        .filter(|id| !common.contains(id))
                        .for_each(|id| { ket.push(*id); });
                    bra_l.iter()
                        .filter(|id| !common.contains(id))
                        .for_each(|id| { bra.push(*id); });
                    ket.sort_unstable();
                    bra.sort_unstable();
                    let ph = *ph_l + *ph_r;
                    Ok(Self::Z { ket, bra, ph })
                }
            },
            (
                Self::X { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::X { ket: ket_r, bra: bra_r, ph: ph_r },
            ) => {
                let mut common: Vec<usize> =
                    Vec::with_capacity(bra_l.len().max(ket_r.len()));
                // check for duplicate kets in result
                for id in ket_r.iter() {
                    if bra_l.contains(id) {
                        common.push(*id);
                    } else if ket_l.contains(id) {
                        return Err(DotDuplicateKetKey(*id));
                    }
                }
                // check for duplicate bras in result
                for id in bra_l.iter() {
                    if !ket_r.contains(id) && bra_r.contains(id) {
                        return Err(DotDuplicateBraKey(*id));
                    }
                }
                // do the dot product
                if common.is_empty() {
                    let mut new =
                        Self::from_terms(
                            self.as_terms().into_iter()
                                .cartesian_product(rhs.as_terms())
                                .map(|(l, r)| l.into_dot(r).unwrap())
                        );
                    new.simplify();
                    Ok(new)
                } else {
                    let mut ket = ket_l.clone();
                    let mut bra = bra_r.clone();
                    ket_r.iter()
                        .filter(|id| !common.contains(id))
                        .for_each(|id| { ket.push(*id); });
                    bra_l.iter()
                        .filter(|id| !common.contains(id))
                        .for_each(|id| { bra.push(*id); });
                    ket.sort_unstable();
                    bra.sort_unstable();
                    let ph = *ph_l + *ph_r;
                    Ok(Self::X { ket, bra, ph })
                }
            },
            (
                Self::H { ket: ket_l, bra: bra_l, arg: arg_l },
                Self::H { ket: ket_r, bra: bra_r, arg: arg_r },
            ) if ket_l.len() == 1 && bra_l.len() == 1
                && (*arg_l + 1.0).norm() < EPSILON
                && ket_r.len() == 1 && bra_r.len() == 1
                && (*arg_r + 1.0).norm() < EPSILON
                && bra_l == ket_r
            => {
                let new = Self::Z {
                    ket: ket_l.clone(),
                    bra: bra_r.clone(),
                    ph: Phase::zero(),
                };
                Ok(new)
            },
            (
                Self::H { ket: ket_l, bra: bra_l, arg: arg_l },
                Self::Z { ket: ket_r, bra: bra_r, ph: _ },
            ) if ket_r.len() == 1 && bra_r.is_empty()
                && ket_l.len() == 1 && bra_l.len() == 1
                && (*arg_l + 1.0).norm() < EPSILON
                && bra_l.contains(&ket_r[0])
            => {
                Ok(rhs.hadamard())
            },
            (
                Self::Z { ket: ket_l, bra: bra_l, ph: _ },
                Self::H { ket: ket_r, bra: bra_r, arg: arg_r },
            ) if ket_l.is_empty() && bra_l.len() == 1
                && ket_r.len() == 1 && bra_r.len() == 1
                && (*arg_r + 1.0).norm() < EPSILON
                && ket_r.contains(&bra_l[0])
            => {
                Ok(self.hadamard())
            },
            (
                Self::H { ket: ket_l, bra: bra_l, arg: arg_l },
                Self::X { ket: ket_r, bra: bra_r, ph: _ },
            ) if ket_r.len() == 1 && bra_r.is_empty()
                && ket_l.len() == 1 && bra_l.len() == 1
                && (*arg_l + 1.0).norm() < EPSILON
                && bra_l.contains(&ket_r[0])
            => {
                Ok(rhs.hadamard())
            },
            (
                Self::X { ket: ket_l, bra: bra_l, ph: _ },
                Self::H { ket: ket_r, bra: bra_r, arg: arg_r },
            ) if ket_l.is_empty() && bra_l.len() == 1
                && ket_r.len() == 1 && bra_r.len() == 1
                && (*arg_r + 1.0).norm() < EPSILON
                && ket_r.contains(&bra_l[0])
            => {
                Ok(self.hadamard())
            },
            (
                l,
                Self::H { ket: ket_r, bra: bra_r, arg: arg_r },
            ) if ket_r.len() == 1 && bra_r.len() == 1
                && (*arg_r + 1.0).norm() < EPSILON
                && l.contains_bra(ket_r[0])
            => {
                let id = ket_r[0];
                let mut terms: Vec<KetBra> = l.as_terms();
                terms.iter_mut()
                    .for_each(|kb| {
                        kb.bra_mut().get_mut(id).unwrap().hadamard_mut();
                    });
                let mut new = Self::Terms { terms };
                new.simplify();
                Ok(new)
            },
            (
                Self::H { ket: ket_l, bra: bra_l, arg: arg_l },
                r,
            ) if ket_l.len() == 1 && bra_l.len() == 1
                && (*arg_l + 1.0).norm() < EPSILON
                && r.contains_ket(bra_l[0])
            => {
                let id = bra_l[0];
                let mut terms: Vec<KetBra> = r.as_terms();
                terms.iter_mut()
                    .for_each(|kb| {
                        kb.ket_mut().get_mut(id).unwrap().hadamard_mut();
                    });
                let mut new = Self::Terms { terms };
                new.simplify();
                Ok(new)
            },
            (
                Self::Terms { terms: terms_l },
                Self::Terms { terms: terms_r },
            ) => {
                let terms: Vec<KetBra> =
                    terms_l.iter()
                    .cartesian_product(terms_r.iter())
                    .map(|(l, r)| l.dot(r))
                    .collect::<KBResult<_>>()?;
                let mut new = Self::Terms { terms };
                new.simplify();
                Ok(new)
            },
            _ => {
                let terms: Vec<KetBra> =
                    self.as_terms().into_iter()
                    .cartesian_product(rhs.as_terms())
                    .map(|(l, r)| l.into_dot(r))
                    .collect::<KBResult<_>>()?;
                let mut new = Self::Terms { terms };
                new.simplify();
                Ok(new)
            },
        }
    }

    pub(crate) fn into_dot(self, rhs: Self) -> KBResult<Self> {
        const EPSILON: f64 = 1e-12;
        match (self, rhs) {
            (
                Self::Z { ket: mut ket_l, bra: bra_l, ph: mut ph_l },
                Self::Z { ket: ket_r, bra: mut bra_r, ph: ph_r },
            ) => {
                let mut common: Vec<usize> =
                    Vec::with_capacity(bra_l.len().max(ket_r.len()));
                // check for duplicate kets in result
                for id in ket_r.iter() {
                    if bra_l.contains(id) {
                        common.push(*id);
                    } else if ket_l.contains(id) {
                        return Err(DotDuplicateKetKey(*id));
                    }
                }
                // check for duplicate bras in result
                for id in bra_l.iter() {
                    if !ket_r.contains(id) && bra_r.contains(id) {
                        return Err(DotDuplicateBraKey(*id));
                    }
                }
                // do the dot product
                if common.is_empty() {
                    // hoping zero-cost abstractions make this okay
                    let l = Self::Z { ket: ket_l, bra: bra_l, ph: ph_l };
                    let r = Self::Z { ket: ket_r, bra: bra_r, ph: ph_r };
                    let mut new =
                        Self::from_terms(
                            l.into_terms().into_iter()
                                .cartesian_product(r.into_terms())
                                .map(|(l, r)| l.into_dot(r).unwrap())
                        );
                    new.simplify();
                    Ok(new)
                } else {
                    ket_r.into_iter()
                        .filter(|id| !common.contains(id))
                        .for_each(|id| { ket_l.push(id); });
                    bra_l.into_iter()
                        .filter(|id| !common.contains(id))
                        .for_each(|id| { bra_r.push(id); });
                    ket_l.sort_unstable();
                    bra_r.sort_unstable();
                    ph_l += ph_r;
                    Ok(Self::Z { ket: ket_l, bra: bra_r, ph: ph_l })
                }
            },
            (
                Self::X { ket: mut ket_l, bra: bra_l, ph: mut ph_l },
                Self::X { ket: ket_r, bra: mut bra_r, ph: ph_r },
            ) => {
                let mut common: Vec<usize> =
                    Vec::with_capacity(bra_l.len().max(ket_r.len()));
                // check for duplicate kets in result
                for id in ket_r.iter() {
                    if bra_l.contains(id) {
                        common.push(*id);
                    } else if ket_l.contains(id) {
                        return Err(DotDuplicateKetKey(*id));
                    }
                }
                // check for duplicate bras in result
                for id in bra_l.iter() {
                    if !ket_r.contains(id) && bra_r.contains(id) {
                        return Err(DotDuplicateBraKey(*id));
                    }
                }
                // do the dot product
                if common.is_empty() {
                    // hoping zero-cost abstractions make this okay
                    let l = Self::X { ket: ket_l, bra: bra_l, ph: ph_l };
                    let r = Self::X { ket: ket_r, bra: bra_r, ph: ph_r };
                    let mut new =
                        Self::from_terms(
                            l.into_terms().into_iter()
                                .cartesian_product(r.into_terms())
                                .map(|(l, r)| l.into_dot(r).unwrap())
                        );
                    new.simplify();
                    Ok(new)
                } else {
                    ket_r.into_iter()
                        .filter(|id| !common.contains(id))
                        .for_each(|id| { ket_l.push(id); });
                    bra_l.into_iter()
                        .filter(|id| !common.contains(id))
                        .for_each(|id| { bra_r.push(id); });
                    ket_l.sort_unstable();
                    bra_r.sort_unstable();
                    ph_l += ph_r;
                    Ok(Self::X { ket: ket_l, bra: bra_r, ph: ph_l })
                }
            },
            (
                Self::H { ket: ket_l, bra: bra_l, arg: arg_l },
                Self::H { ket: ket_r, bra: bra_r, arg: arg_r },
            ) if ket_l.len() == 1 && bra_l.len() == 1
                && (arg_l + 1.0).norm() < EPSILON
                && ket_r.len() == 1 && bra_r.len() == 1
                && (arg_r + 1.0).norm() < EPSILON
                && bra_l == ket_r
            => {
                let new = Self::Z { ket: ket_l, bra: bra_r, ph: Phase::zero() };
                Ok(new)
            },
            (
                Self::H { ket: ket_l, bra: bra_l, arg: arg_l },
                Self::Z { ket: ket_r, bra: bra_r, ph: ph_r },
            ) if ket_r.len() == 1 && bra_r.is_empty()
                && ket_l.len() == 1 && bra_l.len() == 1
                && (arg_l + 1.0).norm() < EPSILON
                && bra_l.contains(&ket_r[0])
            => {
                Ok(Self::Z { ket: ket_r, bra: bra_r, ph: ph_r }.into_hadamard())
            },
            (
                Self::Z { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::H { ket: ket_r, bra: bra_r, arg: arg_r },
            ) if ket_l.is_empty() && bra_l.len() == 1
                && ket_r.len() == 1 && bra_r.len() == 1
                && (arg_r + 1.0).norm() < EPSILON
                && ket_r.contains(&bra_l[0])
            => {
                Ok(Self::Z { ket: ket_l, bra: bra_l, ph: ph_l }.into_hadamard())
            },
            (
                Self::H { ket: ket_l, bra: bra_l, arg: arg_l },
                Self::X { ket: ket_r, bra: bra_r, ph: ph_r },
            ) if ket_r.len() == 1 && bra_r.is_empty()
                && ket_l.len() == 1 && bra_l.len() == 1
                && (arg_l + 1.0).norm() < EPSILON
                && bra_l.contains(&ket_r[0])
            => {
                Ok(Self::X { ket: ket_r, bra: bra_r, ph: ph_r }.into_hadamard())
            },
            (
                Self::X { ket: ket_l, bra: bra_l, ph: ph_l },
                Self::H { ket: ket_r, bra: bra_r, arg: arg_r },
            ) if ket_l.is_empty() && bra_l.len() == 1
                && ket_r.len() == 1 && bra_r.len() == 1
                && (arg_r + 1.0).norm() < EPSILON
                && ket_r.contains(&bra_l[0])
            => {
                Ok(Self::X { ket: ket_l, bra: bra_l, ph: ph_l }.into_hadamard())
            },
            (
                l,
                Self::H { ket: ket_r, bra: bra_r, arg: arg_r },
            ) if ket_r.len() == 1 && bra_r.len() == 1
                && (arg_r + 1.0).norm() < EPSILON
                && l.contains_bra(ket_r[0])
            => {
                let id = ket_r[0];
                let mut terms: Vec<KetBra> = l.into_terms();
                terms.iter_mut()
                    .for_each(|kb| {
                        kb.bra_mut().get_mut(id).unwrap().hadamard_mut();
                    });
                let mut new = Self::Terms { terms };
                new.simplify();
                Ok(new)
            },
            (
                Self::H { ket: ket_l, bra: bra_l, arg: arg_l },
                r,
            ) if ket_l.len() == 1 && bra_l.len() == 1
                && (arg_l + 1.0).norm() < EPSILON
                && r.contains_ket(bra_l[0])
            => {
                let id = bra_l[0];
                let mut terms: Vec<KetBra> = r.into_terms();
                terms.iter_mut()
                    .for_each(|kb| {
                        kb.ket_mut().get_mut(id).unwrap().hadamard_mut();
                    });
                let mut new = Self::Terms { terms };
                new.simplify();
                Ok(new)
            },
            (
                Self::Terms { terms: terms_l },
                Self::Terms { terms: terms_r },
            ) => {
                let terms: Vec<KetBra> =
                    terms_l.into_iter()
                    .cartesian_product(terms_r)
                    .map(|(l, r)| l.into_dot(r))
                    .collect::<KBResult<_>>()?;
                let mut new = Self::Terms { terms };
                new.simplify();
                Ok(new)
            },
            (l, r) => {
                let terms: Vec<KetBra> =
                    l.into_terms().into_iter()
                    .cartesian_product(r.into_terms())
                    .map(|(l, r)| l.into_dot(r))
                    .collect::<KBResult<_>>()?;
                let mut new = Self::Terms { terms };
                new.simplify();
                Ok(new)
            },
        }
    }

    #[allow(clippy::uninit_assumed_init, invalid_value)]
    pub(crate) fn hadamard_mut(&mut self) {
        const EPSILON: f64 = 1e-12;
        match self {
            Self::Z { .. } => {
                let temp =
                    unsafe { std::mem::MaybeUninit::uninit().assume_init() };
                let Self::Z { ket, bra, ph } = std::mem::replace(self, temp)
                    else { unreachable!() };
                let new = Self::X { ket, bra, ph };
                let temp = std::mem::replace(self, new);
                std::mem::forget(temp);
            },
            Self::X { .. } => {
                let temp =
                    unsafe { std::mem::MaybeUninit::uninit().assume_init() };
                let Self::X { ket, bra, ph } = std::mem::replace(self, temp)
                    else { unreachable!() };
                let new = Self::Z { ket, bra, ph };
                let temp = std::mem::replace(self, new);
                std::mem::forget(temp);
            },
            Self::H { ket, bra, arg } => {
                if ket.len() == 1 && bra.len() == 1
                    && (*arg + 1.0).norm() < EPSILON
                { return; }
                let temp =
                    unsafe { std::mem::MaybeUninit::uninit().assume_init() };
                let mut data = std::mem::replace(self, temp);
                data.make_terms();
                let Self::Terms { mut terms } = data else { unreachable!() };
                terms.iter_mut().for_each(|term| { term.hadamard_mut(); });
                let new = Self::Terms { terms };
                let temp = std::mem::replace(self, new);
                std::mem::forget(temp);
            },
            Self::Terms { terms } => {
                terms.iter_mut().for_each(|term| { term.hadamard_mut(); });
            },
            _ => { },
        }
    }

    pub(crate) fn into_hadamard(mut self) -> Self {
        self.hadamard_mut();
        self
    }

    pub(crate) fn hadamard(&self) -> Self {
        let mut new = self.clone();
        new.hadamard_mut();
        new
    }

    pub(crate) fn adjoint_mut(&mut self) {
        match self {
            Self::Z { ket, bra, ph } => {
                std::mem::swap(ket, bra);
                *ph = -*ph;
            },
            Self::X { ket, bra, ph } => {
                std::mem::swap(ket, bra);
                *ph = -*ph;
            },
            Self::H { ket, bra, arg } => {
                std::mem::swap(ket, bra);
                *arg = arg.conj();
            },
            Self::Swap { .. } => { },
            Self::Cup { a, b } => {
                let a = *a;
                let b = *b;
                let _ = std::mem::replace(self, Self::Cap { a, b });
            },
            Self::Cap { a, b } => {
                let a = *a;
                let b = *b;
                let _ = std::mem::replace(self, Self::Cup { a, b });
            },
            Self::Terms { terms } => {
                terms.iter_mut().for_each(|term| { term.adjoint_mut(); });
            },
        }
    }

    pub(crate) fn into_adjoint(mut self) -> Self {
        self.adjoint_mut();
        self
    }

    pub(crate) fn adjoint(&self) -> Self {
        let mut new = self.clone();
        new.adjoint_mut();
        new
    }

    pub(crate) fn into_basis(self, basis: Basis) -> Self {
        match self {
            Self::Z { ket, bra, ph } if ket.is_empty() && bra.is_empty() =>
                Self::Z { ket, bra, ph },
            Self::X { ket, bra, ph } if ket.is_empty() && bra.is_empty() =>
                Self::X { ket, bra, ph },
            Self::H { ket, bra, arg }
                if (ket.is_empty() && bra.is_empty())
                    || (ket.len() == 1 && bra.len() == 1) =>
                Self::H { ket, bra, arg },
            Self::Swap { a, b } => Self::Swap { a, b },
            Self::Cup { a, b } => Self::Cup { a, b },
            Self::Cap { a, b } => Self::Cap { a, b },
            data => {
                let mut new = Self::from_terms(
                    data.into_terms()
                    .into_iter()
                    .flat_map(|term| term.into_basis(basis).into_terms())
                );
                new.simplify();
                new
            },
        }
    }

}

impl fmt::Display for ElementData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut terms: Vec<KetBra> = self.as_terms();
        terms.sort_by(|l, r| l.state_cmp(r));
        let n: usize = terms.len();
        for (k, term) in terms.iter().enumerate() {
            term.fmt(f)?;
            if k < n - 1 { write!(f, " + ")?; }
        }
        Ok(())
    }
}

/// Represents a single element of a diagram as a sum over [`KetBra`]s.
///
/// Note that while every named object in the ZX(H)-calculus (i.e., spiders,
/// cups, caps, and H-boxes) is an `Element` and can be created as such, not
/// every `Element` will be a named object.
#[derive(Clone, Debug)]
pub struct Element {
    pub(crate) data: ElementData,
}

impl From<ElementData> for Element {
    fn from(data: ElementData) -> Self { Self { data } }
}

impl From<KetBra> for Element {
    fn from(ketbra: KetBra) -> Self {
        ElementData::Terms { terms: vec![ketbra] }.into()
    }
}

impl Element {
    /// Return `true` if `self` is a Z-spider.
    pub fn is_z(&self) -> bool { self.data.is_z() }

    /// Return `true` if `self` is an X-spider.
    pub fn is_x(&self) -> bool { self.data.is_x() }

    /// Return `true` if `self` is an H-box.
    pub fn is_h(&self) -> bool { self.data.is_h() }

    /// Return `true` if `self` is a cup (Bell state).
    pub fn is_cup(&self) -> bool { self.data.is_cup() }

    /// Return `true` if `self` is a cap (Bell effect).
    pub fn is_cap(&self) -> bool { self.data.is_cap() }

    /// Return `true` if `self` is a spider.
    pub fn is_spider(&self) -> bool { self.data.is_spider() }

    /// Return `true` if `self` is an atomic form (i.e. a spider, H-box, or Bell
    /// state/effect).
    pub fn is_atomic(&self) -> bool { self.data.is_atomic() }

    /// Return `true` if `self` has all empty ket and bra states.
    pub fn is_scalar(&self) -> bool { self.data.is_scalar() }

    /// If `self` has all empty ket and bra states, return the amplitude as a
    /// single scalar.
    pub fn as_scalar(&self) -> Option<C64> { self.data.as_scalar() }

    /// Create a new Z-spider.
    pub fn z<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        ElementData::z(ins, outs, phase).into()
    }

    /// Create a new X-spider.
    pub fn x<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        ElementData::x(ins, outs, phase).into()
    }

    /// Create a new H-box.
    pub fn h<I, J>(ins: I, outs: J, arg: Option<C64>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        ElementData::h(ins, outs, arg).into()
    }

    /// Create a new SWAP on wires `a` and `b`.
    ///
    /// *Panics if `a == b`*.
    pub fn swap(a: usize, b: usize) -> Self { ElementData::swap(a, b).into() }

    /// Create a new cup (Bell state) on wires `a` and `b`.
    ///
    /// *Panics if `a == b`*.
    pub fn cup(a: usize, b: usize) -> Self { ElementData::cup(a, b).into() }

    /// Create a new cap (Bell effect) on wires `a` and `b`.
    ///
    /// *Panics if `a == b`*.
    pub fn cap(a: usize, b: usize) -> Self { ElementData::cap(a, b).into() }

    /// Create a new `Element` from a series of [`KetBra`]s.
    pub fn from_terms<I>(terms: I) -> Self
    where I: IntoIterator<Item = KetBra>
    {
        ElementData::from_terms(terms).into()
    }

    pub(crate) fn terms(&self) -> Option<&Vec<KetBra>> {
        self.data.terms()
    }

    pub(crate) fn terms_mut(&mut self) -> Option<&mut Vec<KetBra>> {
        self.data.terms_mut()
    }

    /// Return `true` if `self` contains a ket on wire `k`.
    pub fn contains_ket(&self, k: usize) -> bool { self.data.contains_ket(k) }

    /// Return `true` if `self` contains a bra on wire `k`.
    pub fn contains_bra(&self, k: usize) -> bool { self.data.contains_bra(k) }

    /// Return a list of all input wire indices.
    pub fn ins(&self) -> Vec<usize> { self.data.ins() }

    /// Return a list of all output wire indices.
    pub fn outs(&self) -> Vec<usize> { self.data.outs() }

    pub(crate) fn make_terms(&mut self) { self.data.make_terms(); }

    pub(crate) fn as_terms(&self) -> Vec<KetBra> { self.data.as_terms() }

    pub(crate) fn into_terms(self) -> Vec<KetBra> { self.data.into_terms() }

    /// Combine like terms in place, if `self` is non-atomic (no-op otherwise).
    pub fn simplify(&mut self) { self.data.simplify() }

    /// Non-exhaustive equality checking: Returns `true` if two `Element`s are
    /// literally the same, plus a few computationally easy edge cases where the
    /// two correspond to "effectively" the same mapping.
    ///
    /// This method will likely not be able to accurately check equality for two
    /// non-atomic elements unless they comprise exactly the same `KetBra`s
    /// expressed in the same basis; this is to avoid making equality checks
    /// computationally expensive. See also [`iso_mat`][Self::iso_mat].
    pub fn iso(&self, other: &Self) -> bool { self.data.iso(&other.data) }

    /// Exhaustive equality checking: Returns `true` if two `Element`s have all
    /// equal matrix elements.
    ///
    /// This method will compute *all* relevant matrix elements for *all* states
    /// on *all* wires, and thus has exponential overhead. See also
    /// [`iso`][Self::iso], which is called before checking matrix elements to
    /// catch easy cases.
    pub fn iso_mat(&self, other: &Self) -> bool {
        self.data.iso_mat(&other.data)
    }

    /// Compute the dot product `self · rhs`, combining like terms in the
    /// result.
    pub fn dot(&self, rhs: &Self) -> KBResult<Self> {
        self.data.dot(&rhs.data).map(Self::from)
    }

    /// Compute the dot product `self · rhs`, consuming both and combining like
    /// terms in the result.
    pub fn into_dot(self, rhs: Self) -> KBResult<Self> {
        self.data.into_dot(rhs.data).map(Self::from)
    }

    /// Thin wrapper around [`dot`][Self::dot] but with operands reversed (i.e.
    /// computing `rhs · self`) in order to conform to a left-to-right style of
    /// composition.
    pub fn then(&self, rhs: &Self) -> KBResult<Self> { rhs.dot(self) }

    /// Thin wrapper around [`into_dot`][Self::into_dot] but with operands
    /// reversed (i.e. computing `rhs · self`) in order to conform to a
    /// left-to-right style of composition.
    pub fn into_then(self, rhs: Self) -> KBResult<Self> { rhs.into_dot(self) }

    /// Change the basis of operation in place by applying a Hadamard transform
    /// to each qubit index, swapping `0 ↔ +` and `1 ↔ -`. Note that this may
    /// transform `self` into a non-atomic form.
    pub fn hadamard_mut(&mut self) { self.data.hadamard_mut(); }

    /// Change the basis of operation by applying a Hadamard transform to each
    /// qubit index, swapping `0 ↔ +` and `1 ↔ -`. Note that this may transform
    /// `self` into a non-atomic form.
    pub fn into_hadamard(self) -> Self { self.data.into_hadamard().into() }

    /// Change the basis of operation by applying a Hadamard transform to each
    /// qubit index, swapping `0 ↔ +` and `1 ↔ -`. Returns a copy of `self` with
    /// the transformation applied. Note that this may transform `self` into a
    /// non-atomic form.
    pub fn hadamard(&self) -> Self { self.data.hadamard().into() }

    /// Conjugate `self` in place, swapping all kets and bras, and conjugating
    /// the any arguments.
    pub fn adjoint_mut(&mut self) { self.data.adjoint(); }

    /// Conjugate `self`, swapping all kets and bras, and conjugating any
    /// arguments.
    pub fn into_adjoint(self) -> Self { self.data.into_adjoint().into() }

    /// Return the complex conjugate of `self`, with all kets swapped with bras,
    /// and any arguments conjugated.
    pub fn adjoint(&self) -> Self { self.data.adjoint().into() }

    /// Express `self` in `basis`, consuming `self`.
    ///
    /// This will result in a non-atomic form unless `self` is a scalar,
    /// unary and empty H-box, swap, cup, or cap.
    pub fn into_basis(self, basis: Basis) -> Element {
        self.data.into_basis(basis).into()
    }

}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.data.fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalars() {
        const EPSILON: f64 = 1e-12;
        let z_scalar = Element::z([], [], None);
        let z_nonscalar = Element::z([0, 1], [2], None);
        assert_eq!(z_scalar.as_scalar(), Some(C64::from(2.0)));
        assert_eq!(z_nonscalar.as_scalar(), None);

        let x_scalar = Element::x([], [], Some(Phase::pi()));
        let x_nonscalar = Element::x([0], [], None);
        assert!(x_scalar.is_scalar());
        assert!(x_scalar.as_scalar().unwrap().norm() < EPSILON);
        assert_eq!(x_nonscalar.as_scalar(), None);

        let h_scalar = Element::h([], [], Some(C64::i()));
        let h_nonscalar = Element::h([], [0], None);
        assert_eq!(h_scalar.as_scalar(), Some(C64::i()));
        assert_eq!(h_nonscalar.as_scalar(), None);

        let terms_scalar = Element::from_terms([
            KetBra::new(C64::i(),       [], []),
            KetBra::new(C64::from(1.0), [], []),
        ]);
        let terms_nonscalar = Element::from_terms([
            KetBra::new(C64::from(1.0), [(0, State::Zero)], []),
            KetBra::new(C64::from(1.0), [(0, State::One )], []),
        ]);
        let terms_empty = Element::from_terms([]);
        assert_eq!(terms_scalar.as_scalar(), Some(C64::from(1.0) + C64::i()));
        assert_eq!(terms_nonscalar.as_scalar(), None);
        assert_eq!(terms_empty.as_scalar(), Some(C64::from(0.0)));
    }

    #[test]
    #[should_panic]
    fn bad_swap() {
        let _ = Element::swap(0, 0);
    }

    #[test]
    #[should_panic]
    fn bad_cup() {
        let _ = Element::cup(1, 1);
    }

    #[test]
    #[should_panic]
    fn bad_cap() {
        let _ = Element::cap(2, 2);
    }

    #[test]
    fn contains() {
        let z_ket = Element::z([], [0, 1, 2], None);
        let z_bra = Element::z([0, 1, 2], [], None);
        assert!(z_ket.contains_ket(0));
        assert!(!z_ket.contains_ket(3));
        assert!(!z_ket.contains_bra(0));
        assert!(!z_bra.contains_ket(0));
        assert!(!z_bra.contains_bra(3));
        assert!(z_bra.contains_bra(0));

        let x_ket = Element::z([], [0, 1, 2], None);
        let x_bra = Element::z([0, 1, 2], [], None);
        assert!(x_ket.contains_ket(0));
        assert!(!x_ket.contains_ket(3));
        assert!(!x_ket.contains_bra(0));
        assert!(!x_bra.contains_ket(0));
        assert!(!x_bra.contains_bra(3));
        assert!(x_bra.contains_bra(0));

        let h_ket = Element::z([], [0, 1, 2], None);
        let h_bra = Element::z([0, 1, 2], [], None);
        assert!(h_ket.contains_ket(0));
        assert!(!h_ket.contains_ket(3));
        assert!(!h_ket.contains_bra(0));
        assert!(!h_bra.contains_ket(0));
        assert!(!h_bra.contains_bra(3));
        assert!(h_bra.contains_bra(0));

        let swap = Element::swap(0, 1);
        assert!(swap.contains_ket(0));
        assert!(swap.contains_ket(1));
        assert!(!swap.contains_ket(2));
        assert!(swap.contains_bra(0));
        assert!(swap.contains_bra(1));
        assert!(!swap.contains_bra(2));

        let cup = Element::cup(0, 1);
        assert!(cup.contains_ket(0));
        assert!(cup.contains_ket(1));
        assert!(!cup.contains_ket(2));
        assert!(!cup.contains_bra(0));
        assert!(!cup.contains_bra(1));
        assert!(!cup.contains_bra(2));

        let cap = Element::cap(0, 1);
        assert!(!cap.contains_ket(0));
        assert!(!cap.contains_ket(1));
        assert!(!cap.contains_ket(2));
        assert!(cap.contains_bra(0));
        assert!(cap.contains_bra(1));
        assert!(!cap.contains_bra(2));
    }

    #[test]
    fn terms() {
        use State::*;
        let z = Element::z([0], [0, 1], Some(Phase::pi2()));
        assert!(z.terms().is_none());
        let terms = z.as_terms();
        let terms_expected: Vec<KetBra> =
            vec![
                KetBra::new(C64::from(1.0), [(0, Zero), (1, Zero)], [(0, Zero)]),
                KetBra::new(C64::i(),       [(0, One ), (1, One )], [(0, One )]),
            ];
        assert_eq!(terms.len(), terms_expected.len());
        assert!(terms.iter().all(|term| terms_expected.contains(&term)));
        
        let x = Element::x([0, 1], [0], Some(-Phase::pi2()));
        assert!(x.terms().is_none());
        let terms = x.as_terms();
        let terms_expected: Vec<KetBra> =
            vec![
                KetBra::new(C64::from(1.0), [(0, Plus )], [(0, Plus ), (1, Plus )]),
                KetBra::new(-C64::i(),      [(0, Minus)], [(0, Minus), (1, Minus)]),
            ];
        assert_eq!(terms.len(), terms_expected.len());
        assert!(terms.iter().all(|term| terms_expected.contains(&term)));

        let h = Element::h([0, 1], [0], Some(C64::i()));
        assert!(h.terms().is_none());
        let terms = h.as_terms();
        let terms_expected: Vec<KetBra> =
            vec![
                KetBra::new(C64::from(1.0), [(0, Zero)], [(0, Zero), (1, Zero)]),
                KetBra::new(C64::from(1.0), [(0, Zero)], [(0, Zero), (1, One )]),
                KetBra::new(C64::from(1.0), [(0, Zero)], [(0, One ), (1, Zero)]),
                KetBra::new(C64::from(1.0), [(0, Zero)], [(0, One ), (1, One )]),
                KetBra::new(C64::from(1.0), [(0, One )], [(0, Zero), (1, Zero)]),
                KetBra::new(C64::from(1.0), [(0, One )], [(0, Zero), (1, One )]),
                KetBra::new(C64::from(1.0), [(0, One )], [(0, One ), (1, Zero)]),
                KetBra::new(C64::i(),       [(0, One )], [(0, One ), (1, One )]),
            ];
        assert_eq!(terms.len(), terms_expected.len());
        assert!(terms.iter().all(|term| terms_expected.contains(&term)));

        let swap = Element::swap(0, 1);
        assert!(swap.terms().is_none());
        let terms = swap.as_terms();
        let terms_expected: Vec<KetBra> =
            vec![
                KetBra::new(C64::from(1.0), [(0, Zero), (1, Zero)], [(0, Zero), (1, Zero)]),
                KetBra::new(C64::from(1.0), [(0, Zero), (1, One )], [(0, One ), (1, Zero)]),
                KetBra::new(C64::from(1.0), [(0, One ), (1, Zero)], [(0, Zero), (1, One )]),
                KetBra::new(C64::from(1.0), [(0, One ), (1, One )], [(0, One ), (1, One )]),
            ];
        assert_eq!(terms.len(), terms_expected.len());
        assert!(terms.iter().all(|term| terms_expected.contains(&term)));

        let cup = Element::cup(0, 1);
        assert!(cup.terms().is_none());
        let terms = cup.as_terms();
        let terms_expected: Vec<KetBra> =
            vec![
                KetBra::new(C64::from(1.0), [(0, Zero), (1, Zero)], []),
                KetBra::new(C64::from(1.0), [(0, One ), (1, One )], []),
            ];
        assert_eq!(terms.len(), terms_expected.len());
        assert!(terms.iter().all(|term| terms_expected.contains(&term)));

        let cap = Element::cap(0, 1);
        assert!(cap.terms().is_none());
        let terms = cap.as_terms();
        let terms_expected: Vec<KetBra> =
            vec![
                KetBra::new(C64::from(1.0), [], [(0, Zero), (1, Zero)]),
                KetBra::new(C64::from(1.0), [], [(0, One ), (1, One )]),
            ];
        assert_eq!(terms.len(), terms_expected.len());
        assert!(terms.iter().all(|term| terms_expected.contains(&term)));

        let nonatomic = Element::from_terms([
            KetBra::new( C64::i(), [(0, Zero)], [(2, Plus)]),
            KetBra::new(-C64::i(), [(5, Minus)], []),
        ]);
        assert!(nonatomic.terms().is_some());
        let terms = nonatomic.as_terms();
        let terms_expected: Vec<KetBra> =
            vec![
                KetBra::new( C64::i(), [(0, Zero)], [(2, Plus)]),
                KetBra::new(-C64::i(), [(5, Minus)], []),
            ];
        assert_eq!(terms.len(), terms_expected.len());
        assert!(terms.iter().all(|term| terms_expected.contains(&term)));
    }

    #[test]
    fn simplify() {
        use State::*;
        let mut nonatomic = Element::from_terms([
            KetBra::new(Phase::new(0, 4).cis(), [(0, Zero)], [(0, Zero)]),
            KetBra::new(Phase::new(1, 4).cis(), [(0, Zero)], [(0, Zero)]),
            KetBra::new(C64::from(1.0),         [(0, Zero)], [(1, Zero)]),
            KetBra::new(Phase::new(2, 4).cis(), [(0, Zero)], [(0, Zero)]),
            KetBra::new(Phase::new(3, 4).cis(), [(0, Zero)], [(0, Zero)]),
            KetBra::new(C64::from(1.0),         [(0, Zero)], [(0, One )]),
            KetBra::new(C64::from(-0.5),        [(0, Zero)], [(0, One )]),
        ]);
        nonatomic.simplify();
        let terms = nonatomic.as_terms();
        let terms_expected: Vec<KetBra> =
            vec![
                KetBra::new(C64::from(1.0),         [(0, Zero)], [(1, Zero)]),
                KetBra::new(C64::from(0.5),         [(0, Zero)], [(0, One )]),
            ];
        assert_eq!(terms.len(), terms_expected.len());
        assert!(terms.iter().all(|term| terms_expected.contains(&term)));
    }

    #[test]
    fn iso() {
        let z0 = Element::z([1, 2], [3, 4], Some(Phase::pi2()));
        let z1 = Element::z([1, 2], [4, 3], Some(Phase::pi2()));
        let z2 = Element::z([1, 2], [3, 4], Some(Phase::pi()));
        let z3 = Element::z([1, 4], [3, 4], Some(Phase::pi()));
        assert!(z0.iso(&z1));
        assert!(z1.iso(&z0));
        assert!(!z0.iso(&z2));
        assert!(!z2.iso(&z0));
        assert!(!z0.iso(&z3));
        assert!(!z3.iso(&z0));

        let x0 = Element::x([1, 2], [3, 4], Some(Phase::pi2()));
        let x1 = Element::x([1, 2], [3, 4], Some(Phase::pi2()));
        let x2 = Element::x([1, 2], [3, 4], Some(Phase::pi()));
        let x3 = Element::x([1, 4], [3, 4], Some(Phase::pi()));
        assert!(x0.iso(&x1));
        assert!(x1.iso(&x0));
        assert!(!x0.iso(&x2));
        assert!(!x2.iso(&x0));
        assert!(!x0.iso(&x3));
        assert!(!x3.iso(&x0));

        let z_id = Element::z([0], [0], None);
        let x_id = Element::x([0], [0], None);
        let x_id_not = Element::x([1], [0], None);
        assert!(z_id.iso(&x_id));
        assert!(x_id.iso(&z_id));
        assert!(!z_id.iso(&x_id_not));
        assert!(!x_id_not.iso(&z_id));

        let z_ystate_p = Element::z([], [0], Some( Phase::pi2()));
        let x_ystate_p = Element::x([], [0], Some(-Phase::pi2()));
        let z_yeffect_m = Element::z([0], [], Some(-Phase::pi2()));
        let x_yeffect_m = Element::x([0], [], Some( Phase::pi2()));
        assert!(z_ystate_p.iso(&x_ystate_p));
        assert!(x_ystate_p.iso(&z_ystate_p));
        assert!(z_yeffect_m.iso(&x_yeffect_m));
        assert!(x_yeffect_m.iso(&z_yeffect_m));
        assert!(!z_ystate_p.iso(&x_yeffect_m));
        assert!(!x_yeffect_m.iso(&z_ystate_p));

        let z_scalar = Element::z([], [], Some(Phase::pi4()));
        let x_scalar = Element::x([], [], Some(Phase::pi4()));
        let x_nonscalar = Element::x([0], [], Some(Phase::pi4()));
        let z_diffscalar = Element::z([], [], Some(Phase::pi2()));
        let h_scalar = Element::h([], [], Some(Phase::pi4().cis()));
        assert!(z_scalar.iso(&x_scalar));
        assert!(x_scalar.iso(&z_scalar));
        assert!(!z_scalar.iso(&x_nonscalar));
        assert!(!x_nonscalar.iso(&z_scalar));
        assert!(!z_scalar.iso(&z_diffscalar));
        assert!(!z_diffscalar.iso(&z_scalar));
        assert!(z_scalar.iso(&h_scalar));
        assert!(h_scalar.iso(&z_scalar));
        assert!(x_scalar.iso(&h_scalar));
        assert!(h_scalar.iso(&x_scalar));

        let h_state = Element::h([], [0], None);
        let h_effect = Element::h([0], [], Some(Phase::zero().cis()));
        let z_state = Element::z([], [0], Some(Phase::pi()));
        let z_effect = Element::z([0], [], Some(Phase::zero()));
        assert!(h_state.iso(&z_state));
        assert!(z_state.iso(&h_state));
        assert!(h_effect.iso(&z_effect));
        assert!(z_effect.iso(&h_effect));

        let swap1 = Element::swap(1, 2);
        let swap2 = Element::swap(2, 1);
        assert!(swap1.iso(&swap2));
        assert!(swap2.iso(&swap1));

        let cup1 = Element::cup(1, 2);
        let cup2 = Element::cup(2, 1);
        let z_cup = Element::z([], [1, 2], None);
        let x_cup = Element::x([], [1, 2], None);
        assert!(cup1.iso(&cup2));
        assert!(cup2.iso(&cup1));
        assert!(z_cup.iso(&cup1));
        assert!(cup1.iso(&z_cup));
        assert!(cup2.iso(&x_cup));
        assert!(x_cup.iso(&cup2));

        let cap1 = Element::cap(1, 2);
        let cap2 = Element::cap(2, 1);
        let z_cap = Element::z([1, 2], [], None);
        let x_cap = Element::x([2, 1], [], None);
        assert!(cap1.iso(&cap2));
        assert!(cap2.iso(&cap1));
        assert!(cap2.iso(&z_cap));
        assert!(z_cap.iso(&cap2));
        assert!(x_cap.iso(&cap1));
        assert!(cap1.iso(&x_cap));

        let terms1 = Element::from_terms([
            KetBra::new(C64::i(), [(0, State::Zero), (1, State::One)], []),
            KetBra::new(-C64::i(), [(2, State::Plus)], [(2, State::Minus)]),
        ]);
        let terms2 = Element::from_terms([
            KetBra::new(-C64::i(), [(2, State::Plus)], [(2, State::Minus)]),
            KetBra::new(C64::i(), [(0, State::Zero), (1, State::One)], []),
        ]);
        assert!(terms1.iso(&terms2));
        assert!(terms2.iso(&terms1));
    }

    // TODO: iso_mat

    #[test]
    fn dot() {
        let z0 = Element::z([1, 2], [3, 4], Some(Phase::pi2()));
        let z1 = Element::z([5, 6], [1, 2, 7], Some(Phase::pi2()));
        let dot = z0.dot(&z1);
        assert!(dot.is_ok());
        let dot = dot.unwrap();
        assert!(dot.is_z());
        assert!(dot.iso(&Element::z([5, 6], [3, 4, 7], Some(Phase::pi()))));

        let x0 = Element::x([1, 2], [3, 4], Some(Phase::pi2()));
        let x1 = Element::x([5, 6], [1, 2, 7], Some(Phase::pi2()));
        let dot = x0.dot(&x1);
        assert!(dot.is_ok());
        let dot = dot.unwrap();
        assert!(dot.is_x());
        assert!(dot.iso(&Element::x([5, 6], [3, 4, 7], Some(Phase::pi()))));

        let z0 = Element::z([1, 2], [3, 4], Some(Phase::pi2()));
        let z1 = Element::z([5, 6], [8, 9, 7], Some(Phase::pi2()));
        let dot = z0.dot(&z1);
        assert!(dot.is_ok());
        let dot = dot.unwrap();
        assert!(!dot.is_atomic());

        let h0 = Element::h([0], [0], None);
        let h1 = Element::h([0], [0], None);
        let h2 = Element::h([1], [0], None);
        let h3 = Element::h([0], [1], None);
        let h4 = Element::h([0], [0], Some(Phase::zero().cis()));
        let h5 = Element::h([0, 1], [0], None);
        let dot01 = h0.dot(&h1);
        assert!(dot01.is_ok());
        assert!(dot01.unwrap().is_z());
        let dot02 = h0.dot(&h2);
        assert!(dot02.is_ok());
        assert!(dot02.unwrap().is_z());
        let dot03 = h0.dot(&h3);
        assert!(dot03.is_err());
        let dot04 = h0.dot(&h4);
        assert!(dot04.is_ok());
        assert!(!dot04.unwrap().is_atomic());
        let dot05 = h0.dot(&h5);
        assert!(dot05.is_ok());
        assert!(!dot05.unwrap().is_atomic());

        let h = Element::h([0], [0], None);
        let z0 = Element::z([], [0], Some(Phase::pi()));
        let x0 = Element::x([], [0], Some(Phase::pi()));
        let dot0 = h.dot(&z0);
        assert!(dot0.is_ok());
        let dot0 = dot0.unwrap();
        assert!(dot0.iso(&x0));
        let z1 = Element::z([0], [], Some(Phase::pi2()));
        let x1 = Element::x([0], [], Some(Phase::pi2()));
        let dot1 = x1.dot(&h);
        assert!(dot1.is_ok());
        let dot1 = dot1.unwrap();
        assert!(dot1.iso(&z1));
    }

    #[test]
    fn hadamard() {
        let z = Element::z([0], [1, 2], Some(Phase::pi2() + Phase::pi()));
        let hzh = z.hadamard();
        let expected = Element::x([0], [1, 2], Some(Phase::pi2() + Phase::pi()));
        assert!(hzh.iso(&expected));

        let x = Element::x([0], [1, 2], Some(Phase::pi2() + Phase::pi()));
        let hxh = x.hadamard();
        let expected = Element::z([0], [1, 2], Some(Phase::pi2() + Phase::pi()));
        assert!(hxh.iso(&expected));

        let h = Element::h([0], [1], None);
        assert!(h.hadamard().iso(&h));
        let h = Element::h([0], [1], Some(C64::i()));
        assert!(!h.hadamard().iso(&h));

        let terms = Element::from_terms([
            KetBra::new(C64::i(), [(0, State::Zero)], [(1, State::Minus)]),
            KetBra::new(C64::i(), [(0, State::One), (2, State::Plus)], []),
        ]);
        let htermsh = terms.hadamard();
        let expected = Element::from_terms([
            KetBra::new(C64::i(), [(0, State::Plus)], [(1, State::One)]),
            KetBra::new(C64::i(), [(0, State::Minus), (2, State::Zero)], []),
        ]);
        assert!(htermsh.iso(&expected));
    }

    #[test]
    fn adjoint() {
        let z = Element::z([0, 3, 2], [1, 5, 7], Some(Phase::pi2()));
        let z_dag = z.adjoint();
        let expected = Element::z([1, 7, 5], [0, 2, 3], Some(-Phase::pi2()));
        assert!(z_dag.iso(&expected));

        let x = Element::x([0, 3, 2], [1, 5, 7], Some(Phase::pi2()));
        let x_dag = x.adjoint();
        let expected = Element::x([1, 7, 5], [0, 2, 3], Some(-Phase::pi2()));
        assert!(x_dag.iso(&expected));

        let h = Element::h([0, 1], [0, 1], None);
        assert!(h.adjoint().iso(&h));
        let hi = Element::h([2, 1], [0], Some(C64::i()));
        let expected = Element::h([0], [2, 1], Some(-C64::i()));
        assert!(hi.adjoint().iso(&expected));

        let swap = Element::swap(0, 1);
        assert!(swap.adjoint().iso(&swap));

        let cup = Element::cup(0, 1);
        let cap = Element::cap(1, 0);
        assert!(cup.adjoint().iso(&cap));
        assert!(cap.adjoint().iso(&cup));

        let terms = Element::from_terms([
            KetBra::new(C64::i(), [(0, State::Zero)], [(1, State::Minus)]),
            KetBra::new(C64::i(), [(0, State::One), (2, State::Plus)], []),
        ]);
        let terms_dag = terms.adjoint();
        let expected = Element::from_terms([
            KetBra::new(-C64::i(), [(1, State::Minus)], [(0, State::Zero)]),
            KetBra::new(-C64::i(), [], [(0, State::One), (2, State::Plus)]),
        ]);
        assert!(terms_dag.iso(&expected));
    }

    #[test]
    fn into_basis() {
        let z_scalar = Element::z([], [], Some(Phase::pi2()));
        let z_scalar_x = z_scalar.clone().into_basis(Basis::X);
        assert!(z_scalar.iso(&z_scalar_x));

        let x_scalar = Element::x([], [], Some(Phase::pi4()));
        let x_scalar_z = x_scalar.clone().into_basis(Basis::Z);
        assert!(x_scalar.iso(&x_scalar_z));

        let h = Element::h([0], [0], None);
        let h_z = h.clone().into_basis(Basis::Z);
        let h_x = h.clone().into_basis(Basis::X);
        assert!(h.iso(&h_z));
        assert!(h.iso(&h_x));

        let swap = Element::swap(0, 1);
        let swap_z = swap.clone().into_basis(Basis::Z);
        let swap_x = swap.clone().into_basis(Basis::X);
        assert!(swap.iso(&swap_z));
        assert!(swap.iso(&swap_x));

        let cup = Element::cup(0, 1);
        let cup_z = cup.clone().into_basis(Basis::Z);
        let cup_x = cup.clone().into_basis(Basis::X);
        assert!(cup.iso(&cup_z));
        assert!(cup.iso(&cup_x));

        let cap = Element::cap(0, 1);
        let cap_z = cap.clone().into_basis(Basis::Z);
        let cap_x = cap.clone().into_basis(Basis::X);
        assert!(cap.iso(&cap_z));
        assert!(cap.iso(&cap_x));

    }
}

