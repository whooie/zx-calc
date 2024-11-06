use std::fmt;
use itertools::Itertools;
use ndarray as nd;
use num_complex::Complex64 as C64;
use crate::{
    ketbra2::{
        Basis,
        DecompIter,
        KBResult,
        KetBra,
        State,
        StatesIter,
        StateOrd,
    },
    phase::Phase,
};

#[derive(Copy, Clone, Debug)]
pub(crate) enum Kind {
    Id(usize),
    Scalar(C64),
    Z(Phase),
    X(Phase),
    H(C64),
    Swap(usize, usize),
    Cup(usize, usize),
    Cap(usize, usize),
    Unknown,
}

/// Represents a single element of a diagram as a sum over ket-bra terms.
///
/// Note that while every named object in the ZX(H)-calculus (i.e. spiders,
/// cups, caps, and H-boxes) is an `Element` and can be created as such, not
/// every `Element` will be a named object.
#[derive(Clone, Debug)]
pub struct Element {
    pub(crate) kind: Kind, // used only for rendering purposes
    pub(crate) terms: Vec<KetBra>,
}

impl From<KetBra> for Element {
    fn from(data: KetBra) -> Self {
        if let Some(a) = data.as_scalar() {
            Self { kind: Kind::Scalar(a), terms: vec![data] }
        } else {
            Self { kind: Kind::Unknown, terms: vec![data] }
        }
    }
}

impl From<C64> for Element {
    fn from(val: C64) -> Self {
        if val.norm() < 1e-12 {
            Self { kind: Kind::Scalar(0.0.into()), terms: vec![] }
        } else {
            Self {
                kind: Kind::Scalar(val),
                terms: vec![KetBra::new(val, [], [])]
            }
        }
    }
}

impl Element {
    /// Create a new `Element` from a function over states on input/output
    /// indices.
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
            Self { kind: Kind::Unknown, terms }
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
            Self { kind: Kind::Unknown, terms }
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
            Self { kind: Kind::Unknown, terms }
        }
    }

    /// Create a new `Element` from an arbitrary scalar.
    pub fn new_scalar<T>(val: T) -> Self
    where T: Into<C64>
    {
        Self::from(val.into())
    }

    /// Return `true` if `self` comprises only [scalar
    /// terms][KetBra::is_scalar].
    pub fn is_scalar(&self) -> bool {
        self.terms.iter().all(|term| term.is_scalar())
    }

    /// If `self` comprises only [scalar terms][KetBra::is_scalar], extract
    /// their total value as a single scalar.
    pub fn as_scalar(&self) -> Option<C64> {
        self.terms.iter()
            .try_fold(
                C64::from(0.0),
                |acc, term| term.as_scalar().map(|a| acc + a),
            )
    }

    /// Return an iterator over all input wire indices.
    pub fn ins(&self) -> BraIndices<'_> {
        if let Some(first) = self.terms.first() {
            if first.is_scalar() {
                BraIndicesData::Scalar.into()
            } else {
                BraIndicesData::Elem(first.bra_iter()).into()
            }
        } else {
            BraIndicesData::Scalar.into()
        }
    }

    /// Return an iterator over all output wire indices.
    pub fn outs(&self) -> KetIndices<'_> {
        if let Some(first) = self.terms.first() {
            if first.is_scalar() {
                KetIndicesData::Scalar.into()
            } else {
                KetIndicesData::Elem(first.ket_iter()).into()
            }
        } else {
            KetIndicesData::Scalar.into()
        }
    }

    /// Return an identity element, corresponding to an empty wire.
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     id(<i>i</i>) =
    ///       ∣0<sub><i>i</i></sub>⟩⟨0<sub><i>i</i></sub>∣
    ///       + ∣1<sub><i>i</i></sub>⟩⟨1<sub><i>i</i></sub>∣
    ///   </p>
    /// </blockquote>
    pub fn id(i: usize) -> Self {
        let kind = Kind::Id(i);
        let terms =
            vec![
                KetBra::new(1.0, [(i, State::Zero)], [(i, State::Zero)]),
                KetBra::new(1.0, [(i, State::One )], [(i, State::One )]),
            ];
        Self { kind, terms }
    }

    /// Create a Z-spider.
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     <i>Z</i>([<i>i</i><sub>0</sub>, ..., <i>i</i><sub><i>n</i></sub>], [<i>j</i><sub>0</sub>, ..., <i>i</i><sub><i>m</i></sub>], <i>α</i>) =
    ///       ∣0<sub><i>j</i><sub>0</sub></sub>, ..., 0<sub><i>j</i><sub><i>m</i></sub></sub>⟩⟨0<sub><i>i</i><sub>0</sub></sub>, ..., 0<sub><i>i</i><sub><i>n</i></sub></sub>∣
    ///       + <i>e</i><sup><i>iα</i></sup>
    ///         ∣1<sub><i>j</i><sub>0</sub></sub>, ..., 1<sub><i>j</i><sub><i>m</i></sub></sub>⟩⟨1<sub><i>i</i><sub>0</sub></sub>, ..., 1<sub><i>i</i><sub><i>n</i></sub></sub>∣
    ///   </p>
    /// </blockquote>
    ///
    /// The phase defaults to zero, and duplicate wire indices are removed.
    pub fn z<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let ph = phase.unwrap_or_else(Phase::zero);
        let kind = Kind::Z(ph);
        let ins: Vec<usize> = collect_unique(ins);
        let outs: Vec<usize> = collect_unique(outs);
        let ket0 = outs.iter().copied().zip(std::iter::repeat(State::Zero));
        let ket1 = outs.iter().copied().zip(std::iter::repeat(State::One ));
        let bra0 = ins.iter().copied().zip(std::iter::repeat(State::Zero));
        let bra1 = ins.iter().copied().zip(std::iter::repeat(State::One ));
        let terms =
            vec![
                KetBra::new(1.0,      ket0, bra0),
                KetBra::new(ph.cis(), ket1, bra1),
            ];
        Self { kind, terms }
    }

    /// Create an X-spider.
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     <i>X</i>([<i>i</i><sub>0</sub>, ..., <i>i</i><sub><i>n</i></sub>], [<i>j</i><sub>0</sub>, ..., <i>i</i><sub><i>m</i></sub>], <i>α</i>) =
    ///       ∣+<sub><i>j</i><sub>0</sub></sub>, ..., +<sub><i>j</i><sub><i>m</i></sub></sub>⟩⟨+<sub><i>i</i><sub>0</sub></sub>, ..., +<sub><i>i</i><sub><i>n</i></sub></sub>∣
    ///       + <i>e</i><sup><i>iα</i></sup>
    ///         ∣–<sub><i>j</i><sub>0</sub></sub>, ..., –<sub><i>j</i><sub><i>m</i></sub></sub>⟩⟨–<sub><i>i</i><sub>0</sub></sub>, ..., –<sub><i>i</i><sub><i>n</i></sub></sub>∣
    ///   </p>
    /// </blockquote>
    ///
    /// The phase defaults to zero, and duplicate wire indices are removed.
    pub fn x<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let ph = phase.unwrap_or_else(Phase::zero);
        let kind = Kind::X(ph);
        let ins: Vec<usize> = collect_unique(ins);
        let outs: Vec<usize> = collect_unique(outs);
        let ketp = outs.iter().copied().zip(std::iter::repeat(State::Plus ));
        let ketm = outs.iter().copied().zip(std::iter::repeat(State::Minus));
        let brap = ins.iter().copied().zip(std::iter::repeat(State::Plus ));
        let bram = ins.iter().copied().zip(std::iter::repeat(State::Minus));
        let terms =
            vec![
                KetBra::new(1.0,      ketp, brap),
                KetBra::new(ph.cis(), ketm, bram),
            ];
        Self { kind, terms }
    }

    /// Create an H-box.
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     <i>H</i>([<i>i</i><sub>0</sub>, ..., <i>i</i><sub><i>n</i></sub>], [<i>j</i><sub>0</sub>, ..., <i>i</i><sub><i>m</i></sub>], <i>a</i>) =
    ///       Σ<sub><i>s</i><sub><i>k</i></sub> ∊ {0, 1}</sub>
    ///         <i>a</i><sup>Π<sub><i>k</i></sub> <i>s</i><sub><i>k</i></sub></sup>
    ///           ∣<i>s</i><sub><i>j</i><sub>0</sub></sub>, ..., <i>s</i><sub><i>j</i><sub><i>m</i></sub></sub>⟩⟨<i>s</i><sub><i>i</i><sub>0</sub></sub>, ..., <i>s</i><sub><i>i</i><sub><i>n</i></sub></sub>∣
    ///   </p>
    /// </blockquote>
    ///
    /// The parameter defaults to –1. Duplicate wire indices are removed. An
    /// extra, global factor of 1/√2 is added if the total number of wires is 2
    /// and `a` is –1 so that the tensor will coincide with the usual definition
    /// of the Hadamard gate.
    pub fn h<I, J>(ins: I, outs: J, a: Option<C64>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        const EPSILON: f64 = 1e-12;
        const STATES: [State; 2] = [State::Zero, State::One];
        let a = a.unwrap_or_else(|| -C64::from(1.0));
        let kind = Kind::H(a);
        let ket_idx: Vec<usize> = outs.into_iter().collect();
        let bra_idx: Vec<usize> = ins.into_iter().collect();

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
            Self { kind, terms }
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
            Self { kind, terms }
        } else if ket_idx.len() + bra_idx.len() == 2
            && (a + 1.0).norm() < EPSILON
        {
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
            Self { kind, terms }
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
            Self { kind, terms }
        }
    }

    /// Create a swap.
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     SWAP(<i>i</i><sub>0</sub>, <i>i</i><sub>1</sub>) =
    ///       ∣0<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>⟩⟨0<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>∣
    ///       + ∣0<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>⟩⟨1<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>∣
    ///       + ∣1<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>⟩⟨0<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>∣
    ///       + ∣1<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>⟩⟨1<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>∣
    ///   </p>
    /// </blockquote>
    ///
    /// If `i0 == i1`, the returned element is an [identity][Self::id].
    pub fn swap(i0: usize, i1: usize) -> Self {
        use State::*;
        if i0 == i1 { return  Self::id(i0); }
        let kind = Kind::Swap(i0, i1);
        let terms: Vec<KetBra> =
            vec![
                KetBra::new(1.0, [(i0, Zero), (i1, Zero)], [(i0, Zero), (i1, Zero)]),
                KetBra::new(1.0, [(i0, Zero), (i1, One )], [(i0, One ), (i1, Zero)]),
                KetBra::new(1.0, [(i0, One ), (i1, Zero)], [(i0, Zero), (i1, One )]),
                KetBra::new(1.0, [(i0, One ), (i1, One )], [(i0, One ), (i1, One )]),
            ];
        Self { kind, terms }
    }

    /// Create a cup (Bell state).
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     cup(<i>i</i><sub>0</sub>, <i>i</i><sub>1</sub>) =
    ///       ∣0<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>⟩
    ///       + ∣1<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>⟩
    ///   </p>
    /// </blockquote>
    ///
    /// *Panics if `i0 == i1`.*
    pub fn cup(i0: usize, i1: usize) -> Self {
        use State::*;
        if i0 == i1 { panic!("Element::cup: wires cannot be equal"); }
        let kind = Kind::Cup(i0, i1);
        let terms: Vec<KetBra> =
            vec![
                KetBra::new(1.0, [(i0, Zero), (i1, Zero)], []),
                KetBra::new(1.0, [(i0, One ), (i1, One )], []),
            ];
        Self { kind, terms }
    }

    /// Create a cap (Bell effect).
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     cap(<i>i</i><sub>0</sub>, <i>i</i><sub>1</sub>) =
    ///       ⟨0<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>∣
    ///       + ⟨1<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>∣
    ///   </p>
    /// </blockquote>
    ///
    /// *Panics if `i0 == i1`.*
    pub fn cap(i0: usize, i1: usize) -> Self {
        use State::*;
        if i0 == i1 { panic!("Element::cap: wires cannot be equal"); }
        let kind = Kind::Cap(i0, i1);
        let terms: Vec<KetBra> =
            vec![
                KetBra::new(1.0, [], [(i0, Zero), (i1, Zero)]),
                KetBra::new(1.0, [], [(i0, One ), (i1, One )]),
            ];
        Self { kind, terms }
    }

    /// Return an iterator over all ket-bra terms.
    ///
    /// All terms are guaranteed to have states on the same set of input and
    /// output wire indices.
    ///
    /// The iterator item type is `&`[`KetBra`].
    pub fn terms(&self) -> Terms<'_> { Terms { iter: self.terms.iter() } }

    /// Return `true` if `self` has ket terms involving the `k`-th wire index.
    pub fn contains_ket(&self, k: usize) -> bool {
        self.terms.first().is_some_and(|kb| kb.ket().contains_key(k))
    }

    /// Return `true` if `self` has bra terms involving the `k`-th wire index.
    pub fn contains_bra(&self, k: usize) -> bool {
        self.terms.first().is_some_and(|kb| kb.bra().contains_key(k))
    }

    /// Sort terms in place by their indices and states lexicographically.
    pub fn sort_terms(&mut self) {
        self.terms.sort_unstable_by(|l, r| l.state_cmp(r));
    }

    /// Sort terms by their indices and states lexicographically, consuming
    /// `self`.
    pub fn sorted_terms(mut self) -> Self {
        self.sort_terms();
        self
    }

    /// Combine like terms in place, ensuring minimal storage without changing
    /// basis.
    pub fn simplify(&mut self) {
        const EPSILON: f64 = 1e-12;
        let n = self.terms.len();
        let mut remove: Vec<bool> = vec![false; n];
        let mut term0: &KetBra;
        let mut da: C64;
        for k in 0..n {
            if remove[k] { continue; }
            term0 = &self.terms[k];
            da =
                self.terms.iter().zip(remove.iter_mut())
                .skip(k + 1)
                .filter_map(|(term1, v)| {
                    term0.eq_states(term1)
                        .then(|| { *v = true; term1.ampl })
                })
                .sum();
            self.terms[k].ampl += da;
            remove[k] = self.terms[k].ampl.norm() < EPSILON;
        }
        remove.into_iter().enumerate().rev()
            .for_each(|(k, r)| { if r { self.terms.swap_remove(k); } });
    }

    /// Combine like terms, ensuring minimal storage without changing basis,
    /// consuming `self`.
    pub fn simplified(mut self) -> Self {
        self.simplify();
        self
    }

    /// Re-express all terms in a single basis, in place.
    ///
    /// Terms are [simplified][Self::simplify] afterward.
    pub fn make_basis(&mut self, basis: Basis) {
        let mut new_terms: Vec<KetBra> =
            self.terms.drain(..)
            .flat_map(|term| term.into_basis_terms(basis))
            .collect();
        self.terms.append(&mut new_terms);
        self.simplify();
    }

    /// Re-express all terms in a single basis, consuming `self`.
    ///
    /// Terms are [simplified][Self::simplify] afterward.
    pub fn into_basis(mut self, basis: Basis) -> Self {
        self.make_basis(basis);
        self
    }

    /// Compute the dot product `self · rhs`, consuming both.
    ///
    /// Terms are [simplified][Self::simplify] afterward.
    pub fn dot(self, rhs: Self) -> KBResult<Self> {
        let mut terms: Vec<KetBra> =
            self.terms.into_iter()
            .cartesian_product(rhs.terms)
            .map(|(l, r)| l.into_dot(r))
            .collect::<KBResult<_>>()?;
        if terms.first().is_some_and(|kb| kb.is_scalar()) {
            let a: C64 = terms.drain(..).map(|kb| kb.ampl).sum();
            terms.push(KetBra::new(a, [], []));
            let kind = Kind::Scalar(a);
            Ok(Self { kind, terms })
        } else {
            let kind = Kind::Unknown;
            let mut new = Self { kind, terms };
            new.simplify();
            Ok(new)
        }
    }

    /// Thin wrapper around [`dot`][Self::dot] but with the operands reversed
    /// (i.e. computing `rhs · self`) in order to conform to a left-to-right
    /// style of composition.
    pub fn then(self, rhs: Self) -> KBResult<Self> { rhs.dot(self) }

    /// Multiply by a scalar, modifying `self` in place.
    pub fn scalar_mul_mut(&mut self, scalar: C64) {
        const EPSILON: f64 = 1e-12;
        if (scalar + 1.0).norm() < EPSILON { return; }
        self.terms.iter_mut().for_each(|term| { term.ampl *= scalar; });
        self.kind = Kind::Unknown;
    }

    /// Multiply by a scalar, consuming `self`.
    pub fn scalar_mul(mut self, scalar: C64) -> Self {
        self.scalar_mul_mut(scalar);
        self
    }

    /// Conjugate `self` in place, converting kets to bras and bras to kets, and
    /// conjugating the amplitudes of all terms.
    pub fn adjoint_mut(&mut self) {
        match &mut self.kind {
            Kind::Scalar(a) => { *a = a.conj(); },
            Kind::Z(ph) => { *ph = -*ph; },
            Kind::X(ph) => { *ph = -*ph; },
            Kind::H(a) => { *a = a.conj(); },
            Kind::Cup(a, b) => { self.kind = Kind::Cap(*a, *b); },
            Kind::Cap(a, b) => { self.kind = Kind::Cup(*a, *b); },
            _ => { },
        }
        self.terms.iter_mut().for_each(|term| { term.adjoint_mut(); });
    }

    /// Conjugate `self`, consuming `self` to convert kets to bras and bras to
    /// kets, and replace each element of the inner tensor with its complex
    /// conjugate.
    pub fn adjoint(mut self) -> Self {
        self.adjoint_mut();
        self
    }

    /// Return `true` if `self` and `other` denote the same tensor to within an
    /// optional threshold, which defaults to `1e-12`.
    ///
    /// Specifically, `self` and `other` must have identical states on identical
    /// indices, and the modulus of the difference between any two corresponding
    /// terms must be less than `thresh`. `self` and `other` will both be
    /// [re-expressed in the *Z*-basis][Self::make_basis], which requires
    /// mutable references to both.
    pub fn approx_eq(&mut self, other: &mut Self, thresh: Option<f64>) -> bool {
        let eps = thresh.unwrap_or(1e-12);
        self.make_basis(Basis::Z);
        self.sort_terms();
        other.make_basis(Basis::Z);
        other.sort_terms();
        self.terms.iter().zip(&other.terms)
            .all(|(l, r)| {
                l.ket() == r.ket() && l.bra() == r.bra()
                    && (l.ampl - r.ampl).norm() < eps
            })
    }

    /// Convert to a [`tensor::Element`][crate::tensor::Element].
    pub fn as_tensor(&self) -> crate::tensor::Element {
        use crate::tensor as tn;
        match self.kind {
            Kind::Scalar(a) => tn::Element::new_scalar(a),
            Kind::Id(i) => tn::Element::id(i),
            Kind::Z(ph) => tn::Element::z(self.ins(), self.outs(), Some(ph)),
            Kind::X(ph) => tn::Element::x(self.ins(), self.outs(), Some(ph)),
            Kind::H(a) => tn::Element::h(self.ins(), self.outs(), Some(a)),
            Kind::Swap(a, b) => tn::Element::swap(a, b),
            Kind::Cup(a, b) => tn::Element::cup(a, b),
            Kind::Cap(a, b) => tn::Element::cap(a, b),
            Kind::Unknown => {
                if let Some(first) = self.terms.first() {
                    let n_kets = first.ket().len();
                    let n_bras = first.bra().len();
                    let mut indices: Vec<tn::Q> = Vec::new();
                    self.outs()
                        .for_each(|idx| { indices.push(tn::Q::Ket(idx)); });
                    self.ins()
                        .for_each(|idx| { indices.push(tn::Q::Bra(idx)); });
                    if indices.is_empty() {
                        tn::Element::new_scalar(first.ampl)
                    } else {
                        let shape: Vec<usize> = vec![2; n_kets + n_bras];
                        let mut array: nd::ArrayD<C64> =
                            nd::ArrayD::zeros(shape);
                        for term in self.terms.iter() {
                            let iter =
                                term.ket_iter().chain(term.bra_iter())
                                .map(|(id, s)| {
                                    DecompIter::new(id, *s, Basis::Z)
                                        .map(|(_, a, s)| (a, s))
                                })
                                .multi_cartesian_product();
                            for decomp in iter {
                                let decomp_ampl: C64 =
                                    decomp.iter()
                                    .map(|(a, _)| *a)
                                    .product();
                                let pos: Vec<usize> =
                                    decomp.iter()
                                    .map(|(_, s)| *s as usize)
                                    .collect();
                                array[pos.as_slice()] +=
                                    term.ampl * decomp_ampl;
                            }
                        }
                        let data =
                            tn::Tensor::from_array_unchecked(indices, array);
                        tn::Element::new(data)
                    }
                } else {
                    tn::Element::new_scalar(0.0)
                }
            },
        }
    }

}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let n: usize = self.terms.len();
        for (k, term) in self.terms.iter().enumerate() {
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

/// Iterator over the ket-bra terms of an [`Element`].
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

