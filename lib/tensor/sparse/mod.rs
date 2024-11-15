//! Sparse tensor support for [`Element`][super::Element].

use itertools::Itertools;
use ndarray as nd;
use num_complex::Complex64 as C64;
use thiserror::Error;
use crate::{ c64_eq, phase::Phase };
use super::ElementData;

#[derive(Debug, Error)]
pub enum SpError {
    #[error("duplicate ket key in dot {0}")]
    DuplicateKetKey(usize),

    #[error("duplicate bra key in dot {0}")]
    DuplicateBraKey(usize),
}
pub type SpResult<T> = Result<T, SpError>;

pub(crate) mod state;
pub use state::*;

pub(crate) mod ketbra;
pub use ketbra::*;

/// Sparse storage type for use with [`Element`][super::Element].
#[derive(Clone, Debug)]
pub struct Sp(pub(crate) Vec<KetBra>);

impl From<KetBra> for Sp {
    fn from(data: KetBra) -> Self { Self(vec![data]) }
}

impl From<C64> for Sp {
    fn from(a: C64) -> Self {
        if c64_eq(a, 0.0) { Self(vec![]) } else { Self(vec![a.into()]) }
    }
}

impl ElementData for Sp {
    type Error = SpError;
    type InputIter<'a> = BraIndices<'a>;
    type OutputIter<'a> = KetIndices<'a>;

    fn input_iter(&self) -> Self::InputIter<'_> {
        if let Some(first) = self.0.first() {
            if first.is_scalar() {
                BraIndicesData::Scalar.into()
            } else {
                BraIndicesData::Elem(first.bra_iter()).into()
            }
        } else {
            BraIndicesData::Scalar.into()
        }
    }

    fn shift_indices(&mut self, sh: usize) {
        self.0.iter_mut()
            .for_each(|term| {
                term.ket_mut().shift_indices_r(sh);
                term.bra_mut().shift_indices_r(sh);
            });
    }

    fn output_iter(&self) -> Self::OutputIter<'_> {
        if let Some(first) = self.0.first() {
            if first.is_scalar() {
                KetIndicesData::Scalar.into()
            } else {
                KetIndicesData::Elem(first.ket_iter()).into()
            }
        } else {
            KetIndicesData::Scalar.into()
        }
    }

    fn has_input(&self, k: usize) -> bool {
        self.0.first().is_some_and(|kb| kb.bra().contains_key(k))
    }

    fn has_output(&self, k: usize) -> bool {
        self.0.first().is_some_and(|kb| kb.ket().contains_key(k))
    }

    fn scalar<T>(val: T) -> Self
    where T: Into<C64>
    {
        Self::from(val.into())
    }

    fn as_scalar(&self) -> Option<C64> {
        self.0.iter()
            .try_fold(
                C64::from(0.0),
                |acc, term| term.as_scalar().map(|a| acc + a),
            )
    }

    fn is_scalar(&self) -> bool {
        if let Some(first) = self.0.first() {
            first.is_scalar()
        } else {
            true
        }
    }

    fn id(i: usize) -> Self {
        let terms =
            vec![
                KetBra::new(1.0, [(i, State::Zero)], [(i, State::Zero)]),
                KetBra::new(1.0, [(i, State::One )], [(i, State::One )]),
            ];
        Self(terms)
    }

    fn z<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let ph = phase.unwrap_or_else(Phase::zero);
        let ins: Vec<usize> = collect_unique(ins);
        let outs: Vec<usize> = collect_unique(outs);
        if ins.is_empty() && outs.is_empty() {
            return Self::scalar(1.0 + ph.cis());
        }
        let ket0 = outs.iter().copied().zip(std::iter::repeat(State::Zero));
        let ket1 = outs.iter().copied().zip(std::iter::repeat(State::One ));
        let bra0 = ins.iter().copied().zip(std::iter::repeat(State::Zero));
        let bra1 = ins.iter().copied().zip(std::iter::repeat(State::One ));
        let terms =
            vec![
                KetBra::new(1.0,      ket0, bra0),
                KetBra::new(ph.cis(), ket1, bra1),
            ];
        Self(terms)
    }

    fn x<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let ph = phase.unwrap_or_else(Phase::zero);
        let ins: Vec<usize> = collect_unique(ins);
        let outs: Vec<usize> = collect_unique(outs);
        if ins.is_empty() && outs.is_empty() {
            return Self::scalar(1.0 + ph.cis());
        }
        let ketp = outs.iter().copied().zip(std::iter::repeat(State::Plus ));
        let ketm = outs.iter().copied().zip(std::iter::repeat(State::Minus));
        let brap = ins.iter().copied().zip(std::iter::repeat(State::Plus ));
        let bram = ins.iter().copied().zip(std::iter::repeat(State::Minus));
        let terms =
            vec![
                KetBra::new(1.0,      ketp, brap),
                KetBra::new(ph.cis(), ketm, bram),
            ];
        Self(terms)
    }


    fn h<I, J>(ins: I, outs: J, a: Option<C64>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        const STATES: [State; 2] = [State::Zero, State::One];
        let a = a.unwrap_or_else(|| -C64::from(1.0));
        let ket_idx: Vec<usize> = collect_unique(outs);
        let bra_idx: Vec<usize> = collect_unique(ins);

        let empty_ket = ket_idx.is_empty();
        let empty_bra = bra_idx.is_empty();
        if empty_ket && empty_bra {
            Self::from(a)
        } else if empty_ket {
            let terms: Vec<KetBra> =
                bra_idx.iter()
                .map(|_| STATES)
                .multi_cartesian_product()
                .map(|bra_states| {
                    let ampl =
                        if bra_states.iter().all(|s| *s == State::One)
                        { a } else { 1.0.into() };
                    let bra = bra_idx.iter().copied().zip(bra_states);
                    KetBra::new(ampl, [], bra)
                })
                .collect();
            Self(terms)
        } else if empty_bra {
            let terms: Vec<KetBra> =
                ket_idx.iter()
                .map(|_| STATES)
                .multi_cartesian_product()
                .map(|ket_states| {
                    let ampl =
                        if ket_states.iter().all(|s| *s == State::One)
                        { a } else { 1.0.into() };
                    let ket = ket_idx.iter().copied().zip(ket_states);
                    KetBra::new(ampl, ket, [])
                })
                .collect();
            Self(terms)
        } else if ket_idx.len() + bra_idx.len() == 2 && c64_eq(a, -1.0) {
            use State::*;
            let ampl = C64::from(std::f64::consts::FRAC_1_SQRT_2);
            let terms: Vec<KetBra> =
                match (&ket_idx[..], &bra_idx[..]) {
                    (&[i, j], &[]) => vec![
                        KetBra::new( ampl, [(i, Zero), (j, Zero)], []),
                        KetBra::new( ampl, [(i, Zero), (j, One )], []),
                        KetBra::new( ampl, [(i, One ), (j, Zero)], []),
                        KetBra::new(-ampl, [(i, One ), (j, One )], []),
                    ],
                    (&[i], &[j]) => vec![
                        KetBra::new( ampl, [(i, Zero)], [(j, Zero)]),
                        KetBra::new( ampl, [(i, Zero)], [(j, One )]),
                        KetBra::new( ampl, [(i, One )], [(j, Zero)]),
                        KetBra::new(-ampl, [(i, One )], [(j, One )]),
                    ],
                    (&[], &[i, j]) => vec![
                        KetBra::new( ampl, [], [(i, Zero), (j, Zero)]),
                        KetBra::new( ampl, [], [(i, Zero), (j, One )]),
                        KetBra::new( ampl, [], [(i, One ), (j, Zero)]),
                        KetBra::new(-ampl, [], [(i, One ), (j, One )]),
                    ],
                    _ => unreachable!(),
                };
            Self(terms)
        } else {
            let ket_iter =
                ket_idx.iter().map(|_| STATES).multi_cartesian_product();
            let bra_iter =
                bra_idx.iter().map(|_| STATES).multi_cartesian_product();
            let terms: Vec<KetBra> =
                ket_iter.cartesian_product(bra_iter)
                .map(|(ket_states, bra_states)| {
                    let ampl =
                        if ket_states.iter().all(|s| *s == State::One)
                            && bra_states.iter().all(|s| *s == State::One)
                        { a } else { 1.0.into() };
                    let ket = ket_idx.iter().copied().zip(ket_states);
                    let bra = bra_idx.iter().copied().zip(bra_states);
                    KetBra::new(ampl, ket, bra)
                })
                .collect();
            Self(terms)
        }
    }

    fn swap(i0: usize, i1: usize) -> Self {
        use State::*;
        if i0 == i1 { return Self::id(i0); }
        let terms: Vec<KetBra> =
            vec![
                KetBra::new(1.0, [(i0, Zero), (i1, Zero)], [(i0, Zero), (i1, Zero)]),
                KetBra::new(1.0, [(i0, Zero), (i1, One )], [(i0, One ), (i1, Zero)]),
                KetBra::new(1.0, [(i0, One ), (i1, Zero)], [(i0, Zero), (i1, One )]),
                KetBra::new(1.0, [(i0, One ), (i1, One )], [(i0, One ), (i1, One )]),
            ];
        Self(terms)
    }

    fn cup(i0: usize, i1: usize) -> Self {
        use State::*;
        if i0 == i1 { return Self::z([], [i0], None); }
        let terms: Vec<KetBra> =
            vec![
                KetBra::new(1.0, [(i0, Zero), (i1, Zero)], []),
                KetBra::new(1.0, [(i0, One ), (i1, One )], []),
            ];
        Self(terms)
    }

    fn cap(i0: usize, i1: usize) -> Self {
        use State::*;
        if i0 == i1 { return Self::z([i0], [], None); }
        let terms: Vec<KetBra> =
            vec![
                KetBra::new(1.0, [], [(i0, Zero), (i1, Zero)]),
                KetBra::new(1.0, [], [(i0, One ), (i1, One )]),
            ];
        Self(terms)
    }

    fn dot(self, rhs: Self) -> Result<Self, Self::Error> {
        let mut terms: Vec<KetBra> =
            self.0.into_iter()
            .cartesian_product(rhs.0)
            .map(|(l, r)| l.into_dot(r))
            .collect::<SpResult<_>>()?;
        if terms.first().is_some_and(|kb| kb.is_scalar()) {
            let a: C64 = terms.drain(..).map(|kb| kb.ampl).sum();
            terms.push(KetBra::new(a, [], []));
            Ok(Self(terms))
        } else {
            let mut new = Self(terms);
            new.simplify();
            Ok(new)
        }
    }

    fn scalar_mul<C>(&mut self, scalar: C)
    where C: Into<C64>
    {
        let scalar = scalar.into();
        if c64_eq(scalar, 1.0) { return; }
        self.0.iter_mut().for_each(|term| { term.ampl *= scalar; });
    }

    fn conjugate(&mut self) {
        self.0.iter_mut().for_each(|term| { term.adjoint_mut(); });
    }
}

impl Sp {
    /// Create a new `Sp` from a function over states on input/output indices.
    ///
    /// This function is meant to generate a decomposition from an iterator over
    /// all possible [`State`]s over all (unique) indices passed to this
    /// function. Ket (output) states are passed to the left argument of the
    /// function and bra (input) states are passed to the right; the positions
    /// of the states in the arguments correspond to the *unique* indices passed
    /// to this function *in the order in which they're passed*. Specification
    /// of the basis limits the iterator to the states of that basis.
    ///
    /// # Example
    /// TODO
    pub fn new_terms<I, J, F>(
        ins: I,
        outs: J,
        mut ampl: F,
        basis: Basis,
    ) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
        F: FnMut(&[State], &[State]) -> Option<C64>,
    {
        const ZSTATES: [State; 2] = [State::Zero, State::One  ];
        const XSTATES: [State; 2] = [State::Plus, State::Minus];

        let states =
            match basis {
                Basis::Z => ZSTATES,
                Basis::X => XSTATES,
            };

        let ket_idx: Vec<usize> = collect_unique(outs);
        let bra_idx: Vec<usize> = collect_unique(ins);

        let empty_ket = ket_idx.is_empty();
        let empty_bra = bra_idx.is_empty();
        if empty_ket && empty_bra {
            ampl(&[], &[])
                .unwrap_or(0.0.into())
                .into()
        } else if empty_ket {
            let terms: Vec<KetBra> =
                bra_idx.iter()
                .map(|_| states)
                .multi_cartesian_product()
                .filter_map(|bra_states| {
                    ampl(&[], &bra_states)
                        .map(|a| {
                            let bra = bra_idx.iter().copied().zip(bra_states);
                            KetBra::new(a, [], bra)
                        })
                })
                .collect();
            Self(terms)
        } else if empty_bra {
            let terms: Vec<KetBra> =
                ket_idx.iter()
                .map(|_| states)
                .multi_cartesian_product()
                .filter_map(|ket_states| {
                    ampl(&ket_states, &[])
                        .map(|a| {
                            let ket = ket_idx.iter().copied().zip(ket_states);
                            KetBra::new(a, ket, [])
                        })
                })
                .collect();
            Self(terms)
        } else {
            let ket_iter =
                ket_idx.iter().map(|_| states).multi_cartesian_product();
            let bra_iter =
                bra_idx.iter().map(|_| states).multi_cartesian_product();
            let terms: Vec<KetBra> =
                ket_iter.cartesian_product(bra_iter)
                .filter_map(|(ket_states, bra_states)| {
                    ampl(&ket_states, &bra_states)
                        .map(|a| {
                            let ket = ket_idx.iter().copied().zip(ket_states);
                            let bra = bra_idx.iter().copied().zip(bra_states);
                            KetBra::new(a, ket, bra)
                        })
                })
                .collect();
            Self(terms)
        }
    }

    /// Return an iterator over all ket-bra terms.
    ///
    /// All terms are guaranteed to have states on the same set of input and
    /// output wire indices.
    ///
    /// The iterator item type is `&`[`KetBra`].
    pub fn terms(&self) -> Terms<'_> { Terms { iter: self.0.iter() } }

    /// Sort terms in place by their indices and states lexicographically.
    pub fn sort_terms(&mut self) {
        self.0.sort_unstable_by(|l, r| l.state_cmp(r));
    }

    /// Combine like terms in place, ensuring minimal storage without changing
    /// basis.
    pub fn simplify(&mut self) {
        let n = self.0.len();
        let mut remove: Vec<bool> = vec![false; n];
        let mut term0: &KetBra;
        let mut da: C64;
        for k in 0..n {
            if remove[k] { continue; }
            term0 = &self.0[k];
            da =
                self.0.iter().zip(remove.iter_mut())
                .skip(k + 1)
                .filter_map(|(term1, v)| {
                    term0.eq_states(term1)
                        .then(|| { *v = true; term1.ampl })
                })
                .sum();
            self.0[k].ampl += da;
            remove[k] = c64_eq(self.0[k].ampl, 0.0);
        }
        remove.into_iter().enumerate().rev()
            .for_each(|(k, r)| { if r { self.0.swap_remove(k); } });
    }

    /// Re-express all terms in a single basis, in place.
    ///
    /// Terms are [simplified][Self::simplify] afterward.
    pub fn make_basis(&mut self, basis: Basis) {
        let mut new_terms: Vec<KetBra> =
            self.0.drain(..)
            .flat_map(|term| term.into_basis_terms(basis))
            .collect();
        self.0.append(&mut new_terms);
        self.simplify();
    }

    /// Return `true` if `self` and `other` denote the same tensor to within an
    /// optional threshold, which defaults to `1e-12`.
    ///
    /// Specifically, `self` and `other` must have identical states on identical
    /// indices, and the modulus of the difference between any two corresponding
    /// terms must be less than `thresh`. `self` and `other` will both be
    /// [re-expressed in the *Z*-basis][Self::make_basis], which requires
    /// mutable references to both.
    pub fn approx_eq(&mut self, other: &mut Self, thresh: Option<f64>)
        -> bool
    {
        let eps = thresh.unwrap_or(1e-12);
        self.make_basis(Basis::Z);
        self.sort_terms();
        other.make_basis(Basis::Z);
        other.sort_terms();
        self.0.iter().zip(&other.0)
            .all(|(l, r)| {
                l.ket() == r.ket() && l.bra() == r.bra()
                    && (l.ampl - r.ampl).norm() < eps
            })
    }

    /// Convert to a corresponding dense representation.
    pub fn to_dense(&self) -> super::dense::De {
        use super::dense as de;
        if let Some(first) = self.0.first() {
            let n_kets = first.ket().len();
            let n_bras = first.bra().len();
            if n_kets + n_bras == 0 {
                de::De::scalar(first.ampl)
            } else {
                let mut indices: Vec<de::Q> = Vec::new();
                self.output_iter()
                    .for_each(|idx| { indices.push(de::Q::Ket(idx)); });
                self.input_iter()
                    .for_each(|idx| { indices.push(de::Q::Bra(idx)); });
                let shape: Vec<usize> = vec![2; n_kets + n_bras];
                let mut array: nd::ArrayD<C64> = nd::ArrayD::zeros(shape);
                for term in self.0.iter() {
                    let decomp_iter =
                        term.ket_iter().chain(term.bra_iter())
                        .map(|(id, s)| {
                            DecompIter::new(id, *s, Basis::Z)
                                .map(|(_, a, s)| (a, s))
                        })
                        .multi_cartesian_product();
                    for decomp in decomp_iter {
                        let decomp_ampl: C64 =
                            decomp.iter()
                            .map(|(a, _)| *a)
                            .product();
                        let pos: Vec<usize> =
                            decomp.iter()
                            .map(|(_, s)| *s as usize)
                            .collect();
                        array[pos.as_slice()] += term.ampl * decomp_ampl;
                    }
                }
                let data =
                    de::Tensor::from_array_unchecked(indices, array);
                de::De(data)
            }
        } else {
            de::De::scalar(0.0)
        }
    }
}

impl std::fmt::Display for Sp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let n: usize = self.0.len();
        for (k, term) in self.0.iter().enumerate() {
            term.fmt(f)?;
            if k < n - 1 { write!(f, " + ")?; }
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
enum KetIndicesData<'a> {
    Scalar,
    Elem(StatesIter<'a>),
}

/// Iterator over all ket-side wire indices of a tensor.
///
/// The iterator item type is `usize`.
#[derive(Clone, Debug)]
pub struct KetIndices<'a>(KetIndicesData<'a>);

impl<'a> From<KetIndicesData<'a>> for KetIndices<'a> {
    fn from(data: KetIndicesData<'a>) -> Self { Self(data) }
}

impl<'a> Iterator for KetIndices<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        match &mut self.0 {
            KetIndicesData::Scalar => None,
            KetIndicesData::Elem(ref mut iter) => {
                iter.next().map(|(idx, _)| idx)
            },
        }
    }
}

impl<'a> DoubleEndedIterator for KetIndices<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        match &mut self.0 {
            KetIndicesData::Scalar => None,
            KetIndicesData::Elem(ref mut iter) => {
                iter.next_back().map(|(idx, _)| idx)
            },
        }
    }
}

impl<'a> ExactSizeIterator for KetIndices<'a> {
    fn len(&self) -> usize {
        match &self.0 {
            KetIndicesData::Scalar => 0,
            KetIndicesData::Elem(ref iter) => iter.len(),
        }
    }
}

impl<'a> std::iter::FusedIterator for KetIndices<'a> { }

#[derive(Clone, Debug)]
enum BraIndicesData<'a> {
    Scalar,
    Elem(StatesIter<'a>),
}

/// Iterator over all bra-side wire indices of a tensor.
///
/// The iterator item type is `usize`.
#[derive(Clone, Debug)]
pub struct BraIndices<'a>(BraIndicesData<'a>);

impl<'a> From<BraIndicesData<'a>> for BraIndices<'a> {
    fn from(data: BraIndicesData<'a>) -> Self { Self(data) }
}

impl<'a> Iterator for BraIndices<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        match &mut self.0 {
            BraIndicesData::Scalar => None,
            BraIndicesData::Elem(ref mut iter) => {
                iter.next().map(|(idx, _)| idx)
            },
        }
    }
}

impl<'a> DoubleEndedIterator for BraIndices<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        match &mut self.0 {
            BraIndicesData::Scalar => None,
            BraIndicesData::Elem(ref mut iter) => {
                iter.next_back().map(|(idx, _)| idx)
            },
        }
    }
}

impl<'a> ExactSizeIterator for BraIndices<'a> {
    fn len(&self) -> usize {
        match &self.0 {
            BraIndicesData::Scalar => 0,
            BraIndicesData::Elem(ref iter) => iter.len(),
        }
    }
}

impl<'a> std::iter::FusedIterator for BraIndices<'a> { }

/// Iterator over the ket-bra terms of an [`Element`][super::Element].
///
/// The iterator item type is `&`[`KetBra`].
#[derive(Clone, Debug)]
pub struct Terms<'a> {
    iter: std::slice::Iter<'a, KetBra>,
}

impl<'a> Iterator for Terms<'a> {
    type Item = &'a KetBra;

    fn next(&mut self) -> Option<Self::Item> { self.iter.next() }
}

impl<'a> DoubleEndedIterator for Terms<'a> {
    fn next_back(&mut self) -> Option<Self::Item> { self.iter.next_back() }
}

impl<'a> ExactSizeIterator for Terms<'a> {
    fn len(&self) -> usize { self.iter.len() }
}

impl<'a> std::iter::FusedIterator for Terms<'a> { }

fn collect_unique<I, T>(iter: I) -> Vec<T>
where
    I: IntoIterator<Item = T>,
    T: PartialEq,
{
    let mut items: Vec<T> = Vec::new();
    iter.into_iter()
        .for_each(|it| { if !items.contains(&it) { items.push(it); } });
    items
}

