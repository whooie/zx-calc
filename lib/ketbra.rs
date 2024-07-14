//! Naive tools to create and compute end products of diagrams in the
//! ZX(H)-calculus.
//!
//! You will most likely only want to deal directly with [`Element`] and
//! [`Diagram`].

use std::{
    cmp::Ordering,
    // collections::{ HashMap, HashSet },
    fmt,
    fs,
    io::Write,
    ops::{ Deref, DerefMut },
    path::Path,
};
pub use std::f64::consts::{
    SQRT_2 as RT2,
    FRAC_1_SQRT_2 as ONRT2,
    FRAC_PI_2 as PI2,
    FRAC_PI_3 as PI3,
    FRAC_PI_4 as PI4,
    FRAC_PI_6 as PI6,
    FRAC_PI_8 as PI8,
    PI,
    TAU,
};
use itertools::Itertools;
use num_complex::Complex64 as C64;
use rustc_hash::{ FxHashMap as HashMap, FxHashSet as HashSet };
use thiserror::Error;
use crate::{ c, graph::{ self, format_phase } };

#[derive(Debug, Error)]
pub enum KBError {
    /// `element: all kets must have the same keys`
    #[error("element: all kets must have the same keys")]
    ElementKetSameKeys,

    /// `element: all bras must have the same keys`
    #[error("element: all bras must have the same keys")]
    ElementBraSameKeys,

    /// `ketbra: ket keys must be disjoint in kron`
    #[error("ketbra: ket keys must be disjoint in kron;\nlhs: {0:?}\nrhs: {1:?}")]
    KronDisjointKetKeys(Box<[usize]>, Box<[usize]>),

    /// `ketbra: bra keys must be disjoint in kron`
    #[error("ketbra: bra keys must be disjoint in kron;\nlhs: {0:?}\nrhs: {1:?}")]
    KronDisjointBraKeys(Box<[usize]>, Box<[usize]>),

    /// `ketbra: duplicate ket key in dot resultant`
    #[error("ketbra: duplicate ket key in dot resultant: {0}")]
    DotDuplicateKetKey(usize),

    /// `ketbra: duplicate bra key in dot resultant`
    #[error("ketbra: duplicate bra key in dot resultant: {0}")]
    DotDuplicateBraKey(usize),

    /// `diagram: collected elements must have pair-wise matching input/output wires`
    #[error("diagram: collected elements must have pair-wise matching input/output wires")]
    DiagramSliceMatchingWires,

    /// `diagram: encountered a non-generator element`
    #[error("diagram: encountered a non-generator element")]
    DiagramConversionNonGenerator,

    #[error("error constructing GraphViz representation: {0}")]
    GraphVizError(String),

    #[error("I/O error: {0}")]
    IOError(#[from] std::io::Error),
}
pub type KBResult<T> = Result<T, KBError>;

/// A state of a single wire.
///
/// States of different bases are defined relative to each other according to
/// ```math
/// \begin{aligned}
///     \ket{0} &= \frac{1}{\sqrt{2}} \Big( \ket{+} + \ket{-} \Big)
///     \\
///     \ket{1} &= \frac{1}{\sqrt{2}} \Big( \ket{+} - \ket{-} \Big)
///     \\
///     \ket{+} &= \frac{1}{\sqrt{2}} \Big( \ket{0} + \ket{1} \Big)
///     \\
///     \ket{-} &= \frac{1}{\sqrt{2}} \Big( \ket{0} - \ket{1} \Big)
/// \end{aligned}
/// ```
/// with
/// ```math
/// \braket{0|1} = \braket{+|-} = 0
/// ```
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum State {
    Zero,
    One,
    Plus,
    Minus,
}
use State::*;

impl fmt::Display for State {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Zero  => write!(f, "0"),
            One   => write!(f, "1"),
            Plus  => write!(f, "+"),
            Minus => write!(f, "-"),
        }
    }
}

impl State {
    /// Convert `self` to the appropriate basis as an [`Element`] representing a
    /// state (ket) on a single wire.
    pub fn to_ket(self, idx: usize, basis: Basis) -> Element {
        let ketbras: Vec<KetBra>
            = match (basis, self) {
                (Basis::Z, Zero ) | (Basis::Z, One  ) => vec![
                    KetBra::new(c!(1.0), [(idx, self)], []),
                ],
                (Basis::Z, Plus ) => vec![
                    KetBra::new(c!( ONRT2), [(idx, Zero )], []),
                    KetBra::new(c!( ONRT2), [(idx, One  )], []),
                ],
                (Basis::Z, Minus) => vec![
                    KetBra::new(c!( ONRT2), [(idx, Zero )], []),
                    KetBra::new(c!(-ONRT2), [(idx, One  )], []),
                ],
                (Basis::X, Plus ) | (Basis::X, Minus) => vec![
                    KetBra::new(c!(1.0), [(idx, self)], []),
                ],
                (Basis::X, Zero ) => vec![
                    KetBra::new(c!( ONRT2), [(idx, Plus )], []),
                    KetBra::new(c!( ONRT2), [(idx, Minus)], []),
                ],
                (Basis::X, One  ) => vec![
                    KetBra::new(c!( ONRT2), [(idx, Plus )], []),
                    KetBra::new(c!(-ONRT2), [(idx, Minus)], []),
                ],
            };
        Element { kind: ElementKind::Unknown, terms: ketbras }
    }

    /// Convert `self` to the appropriate basis as in [`Element`] representing
    /// an effect (bra) on a single wire.
    pub fn to_bra(self, idx: usize, basis: Basis) -> Element {
        let ketbras: Vec<KetBra>
            = match (basis, self) {
                (Basis::Z, Zero ) | (Basis::Z, One  ) => vec![
                    KetBra::new(c!(1.0), [], [(idx, self)])
                ],
                (Basis::Z, Plus ) => vec![
                    KetBra::new(c!( ONRT2), [], [(idx, Zero )]),
                    KetBra::new(c!( ONRT2), [], [(idx, One  )]),
                ],
                (Basis::Z, Minus) => vec![
                    KetBra::new(c!( ONRT2), [], [(idx, Zero )]),
                    KetBra::new(c!(-ONRT2), [], [(idx, One  )]),
                ],
                (Basis::X, Plus ) | (Basis::X, Minus) => vec![
                    KetBra::new(c!(1.0), [], [(idx, self)])
                ],
                (Basis::X, Zero ) => vec![
                    KetBra::new(c!( ONRT2), [], [(idx, Plus )]),
                    KetBra::new(c!( ONRT2), [], [(idx, Minus)]),
                ],
                (Basis::X, One  ) => vec![
                    KetBra::new(c!( ONRT2), [], [(idx, Plus )]),
                    KetBra::new(c!(-ONRT2), [], [(idx, Minus)]),
                ],
            };
        Element { kind: ElementKind::Unknown, terms: ketbras }
    }
}

/// A specific basis for representing states.
///
/// The `$Z$` basis refers to [`State`]s `Zero` and `One`, while the `$X$` basis
/// refers to `Plus` and `Minus`. The `$Z$` basis is chosen by default for
/// actual computations.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Basis {
    Z,
    X,
}

/// Represents a single matrix element as a ket-bra with an associated
/// amplitude.
///
/// Mathematically this is an expression fo the form
/// ```math
/// a \ket{s_{j_0}, \dots, s_{j_n}} \bra{s_{i_0}, \dots, s_{i_m}}
/// ```
/// where `$a$` is the amplitude and `$s_k \in \{0, 1, +, -\}$`. `$i_0, \dots,
/// i_n$` and `$j_0, \dots, j_m$` index incoming and outgoing wires.
#[derive(Clone, Debug, PartialEq)]
pub struct KetBra {
    ampl: C64,
    ket: HashMap<usize, State>,
    bra: HashMap<usize, State>,
}

impl KetBra {
    /// Create a new `KetBra`.
    pub fn new<I, J>(ampl: C64, ket: I, bra: J) -> Self
    where
        I: IntoIterator<Item = (usize, State)>,
        J: IntoIterator<Item = (usize, State)>,
    {
        Self {
            ampl,
            ket: ket.into_iter().collect(),
            bra: bra.into_iter().collect(),
        }
    }

    /// If `self` is empty then return the amplitude as a single scalar.
    pub fn as_scalar(&self) -> Option<C64> {
        (self.ket.is_empty() && self.bra.is_empty()).then_some(self.ampl)
    }

    /// Express `self` in `basis`.
    ///
    /// This returns an [`Element`] because conversion to a different basis
    /// typically introduces factors expressed as sums, resulting in the need
    /// for multiple ket-bras in the final expression.
    pub fn to_basis(&self, basis: Basis) -> Element {
        if self.ket.is_empty() && self.bra.is_empty() {
            self.clone().into()
        } else {
            let KetBra { ampl, ket, bra } = self;
            let terms: Vec<KetBra>
                = Iterator::chain(
                    ket.iter().map(|(idx, s)| (*s).to_ket(*idx, basis)),
                    bra.iter().map(|(idx, s)| (*s).to_bra(*idx, basis)),
                )
                .multi_cartesian_product()
                .map(|ketbras| {
                    let mut term = Self::into_multikron_unchecked(ketbras);
                    term.ampl *= ampl;
                    term
                })
                .collect();
            Element { kind: ElementKind::Unknown, terms }.simplified()
        }
    }

    /// Express `self` in `basis`, consuming `self`.
    ///
    /// This returns an [`Element`] because conversion to a different basis
    /// typically introduces factors expressed as sums, resulting in the need
    /// for multiple ket-bras in the final expression.
    pub fn into_basis(self, basis: Basis) -> Element {
        if self.ket.is_empty() && self.bra.is_empty() {
            self.into()
        } else {
            let KetBra { ampl, ket, bra } = self;
            let terms: Vec<KetBra>
                = Iterator::chain(
                    ket.into_iter().map(|(idx, s)| s.to_ket(idx, basis)),
                    bra.into_iter().map(|(idx, s)| s.to_bra(idx, basis)),
                )
                .multi_cartesian_product()
                .map(|ketbras| {
                    let mut term = Self::into_multikron_unchecked(ketbras);
                    term.ampl *= ampl;
                    term
                })
                .collect();
            Element { kind: ElementKind::Unknown, terms }.simplified()
        }
    }

    /// Compute the Kronecker product of `self` and `rhs`.
    ///
    /// No inner products are evaluated; the set of wire indices contained
    /// within `self` and `rhs`'s kets and bras must be disjoint, otherwise this
    /// operation fails.
    pub fn kron(&self, rhs: &KetBra) -> KBResult<KetBra> {
        let self_ket_keys: HashSet<&usize> = self.ket.keys().collect();
        let rhs_ket_keys: HashSet<&usize> = rhs.ket.keys().collect();
        self_ket_keys.is_disjoint(&rhs_ket_keys).then_some(())
            .ok_or(KBError::KronDisjointKetKeys(
                self_ket_keys.into_iter().copied().collect(),
                rhs_ket_keys.into_iter().copied().collect(),
            ))?;
        let self_bra_keys: HashSet<&usize> = self.bra.keys().collect();
        let rhs_bra_keys: HashSet<&usize> = rhs.bra.keys().collect();
        self_bra_keys.is_disjoint(&rhs_bra_keys).then_some(())
            .ok_or(KBError::KronDisjointBraKeys(
                self_bra_keys.into_iter().copied().collect(),
                rhs_bra_keys.into_iter().copied().collect(),
            ))?;
        Ok(
            Self::new(
                self.ampl * rhs.ampl,
                self.ket.iter().map(|(k, sk)| (*k, *sk))
                    .chain(rhs.ket.iter().map(|(k, sk)| (*k, *sk))),
                self.bra.iter().map(|(k, sk)| (*k, *sk))
                    .chain(rhs.bra.iter().map(|(k, sk)| (*k, *sk))),
            )
        )
    }

    /// Compute the Kronecker product of `self` and `rhs`, consuming both.
    ///
    /// No inner products are evaluated; the set of wire indices contained
    /// within `self` and `rhs`'s kets and bras must be disjoint, otherwise this
    /// operation fails.
    pub fn into_kron(self, rhs: KetBra) -> KBResult<KetBra> {
        let self_ket_keys: HashSet<&usize> = self.ket.keys().collect();
        let rhs_ket_keys: HashSet<&usize> = rhs.ket.keys().collect();
        self_ket_keys.is_disjoint(&rhs_ket_keys).then_some(())
            .ok_or(KBError::KronDisjointKetKeys(
                self_ket_keys.into_iter().copied().collect(),
                rhs_ket_keys.into_iter().copied().collect(),
            ))?;
        let self_bra_keys: HashSet<&usize> = self.bra.keys().collect();
        let rhs_bra_keys: HashSet<&usize> = rhs.bra.keys().collect();
        self_bra_keys.is_disjoint(&rhs_bra_keys).then_some(())
            .ok_or(KBError::KronDisjointBraKeys(
                self_bra_keys.into_iter().copied().collect(),
                rhs_bra_keys.into_iter().copied().collect(),
            ))?;
        Ok(
            Self::new(
                self.ampl * rhs.ampl,
                self.ket.into_iter().chain(rhs.ket),
                self.bra.into_iter().chain(rhs.bra),
            )
        )
    }

    /// Compute the Kronecker product of multiple ket-bras.
    ///
    /// No inner products are evaluated; the set of wire indices contained
    /// within `self` and `rhs`'s kets and bras must be disjoint, otherwise this
    /// operation fails.
    pub fn multikron<'a, I>(ketbras: I) -> KBResult<KetBra>
    where I: IntoIterator<Item = &'a KetBra>
    {

        let mut acc = Self::new(c!(1.0), [], []);
        for kb in ketbras.into_iter() {
            acc = acc.kron(kb)?;
        }
        Ok(acc)
    }

    /// Compute the Kronecker product of multiple ket-bras, consuming all
    /// inputs.
    ///
    /// No inner products are evaluated; the set of wire indices contained
    /// within `self` and `rhs`'s kets and bras must be disjoint, otherwise this
    /// operation fails.
    pub fn into_multikron<I>(ketbras: I) -> KBResult<KetBra>
    where I: IntoIterator<Item = KetBra>
    {
        let mut acc = Self::new(c!(1.0), [], []);
        for kb in ketbras.into_iter() {
            acc = acc.into_kron(kb)?;
        }
        Ok(acc)
    }

    fn multikron_unckecked<'a, I>(ketbras: I) -> KetBra
    where I: IntoIterator<Item = &'a KetBra>
    {
        let mut ampl: C64 = c!(1.0);
        let mut ket: HashMap<usize, State> = HashMap::default();
        let mut bra: HashMap<usize, State> = HashMap::default();
        for kb in ketbras.into_iter() {
            let KetBra { ampl: ampl_k, ket: ket_k, bra: bra_k } = kb;
            ampl *= ampl_k;
            ket_k.iter().for_each(|(idx, s)| { ket.insert(*idx, *s); });
            bra_k.iter().for_each(|(idx, s)| { bra.insert(*idx, *s); });
        }
        Self { ampl, ket, bra }
    }

    fn into_multikron_unchecked<I>(ketbras: I) -> KetBra
    where I: IntoIterator<Item = KetBra>
    {
        let mut ampl: C64 = c!(1.0);
        let mut ket: HashMap<usize, State> = HashMap::default();
        let mut bra: HashMap<usize, State> = HashMap::default();
        for kb in ketbras.into_iter() {
            let KetBra { ampl: ampl_k, ket: ket_k, bra: bra_k } = kb;
            ampl *= ampl_k;
            ket_k.into_iter().for_each(|(idx, s)| { ket.insert(idx, s); });
            bra_k.into_iter().for_each(|(idx, s)| { bra.insert(idx, s); });
        }
        Self { ampl, ket, bra }
    }

    /// Compute the dot-product of `self` and `rhs`.
    ///
    /// Unmatched kets and bras are left untouched and "passed through" the
    /// product; e.g.
    /// ```math
    /// \ket{0_0 +_1 1_2} \bra{0_0 -_2}
    ///     \cdot \ket{-_2 1_3} \bra{+_1 1_3}
    ///     = \ket{0_0 +_1 1_2 1_3} \bra{0_0 +_1 1_3}
    /// ```
    /// The product must not leave multiple states on the same wire.
    pub fn dot(&self, rhs: &Self) -> KBResult<Element> {
        // check for duplicate kets in result
        let mut wire_counts: HashMap<usize, isize>
            = self.ket.keys().map(|idx| (*idx, 1)).collect();
        self.bra.keys()
            .for_each(|idx| {
                wire_counts.entry(*idx)
                    .and_modify(|count| *count -= 1)
                    .or_insert(-1);
            });
        rhs.ket.keys()
            .for_each(|idx| {
                wire_counts.entry(*idx)
                    .and_modify(|count| *count += 1)
                    .or_insert(1);
            });
        if let Some((idx, _))
            = wire_counts.iter().find(|(_idx, count)| **count > 1)
        {
            return Err(KBError::DotDuplicateKetKey(*idx));
        }

        // check for duplicate bras in result
        wire_counts = self.bra.keys().map(|idx| (*idx, 1)).collect();
        rhs.ket.keys()
            .for_each(|idx| {
                wire_counts.entry(*idx)
                    .and_modify(|count| *count -= 1)
                    .or_insert(-1);
            });
        rhs.bra.keys()
            .for_each(|idx| {
                wire_counts.entry(*idx)
                    .and_modify(|count| *count += 1)
                    .or_insert(1);
            });
        if let Some((idx, _))
            = wire_counts.iter().find(|(_idx, count)| count > &&1)
        {
            return Err(KBError::DotDuplicateBraKey(*idx));
        }

        let terms: Vec<KetBra>
            = Itertools::cartesian_product(
                self.to_basis(Basis::Z).into_iter(),
                rhs.to_basis(Basis::Z).into_iter(),
            )
            .filter_map(|(mut l, mut r)| {
                let common: HashSet<usize>
                    = l.bra.keys()
                    .filter_map(|l| r.ket.contains_key(l).then_some(*l))
                    .chain(
                        r.ket.keys()
                        .filter_map(|r| l.bra.contains_key(r).then_some(*r))
                    )
                    .collect();
                let dot: bool
                    = common.iter()
                    .all(|idx| {
                        let ls: &State = l.bra.get(idx).unwrap();
                        let rs: &State = r.bra.get(idx).unwrap();
                        matches!(
                            (ls, rs),
                            (&Zero, &Zero) | (&One , &One )
                        )
                    });
                if dot && !common.is_empty() {
                    None
                } else {
                    common.iter().for_each(|idx| { l.bra.remove(idx); });
                    common.iter().for_each(|idx| { r.ket.remove(idx); });
                    Some(l.into_kron(r))
                }
            })
            .collect::<KBResult<Vec<KetBra>>>()?;
        Ok(Element::new(terms)?.simplified())
    }

    /// Compute the dot-product of `self` and `rhs`, consuming both.
    ///
    /// Unmatched kets and bras are left untouched and "passed through" the
    /// product; e.g.
    /// ```math
    /// \ket{0_0 +_1 1_2} \bra{0_0 -_2}
    ///     \cdot \ket{-_2 1_3} \bra{+_1 1_3}
    ///     = \ket{0_0 +_1 1_2 1_3} \bra{0_0 +_1 1_3}
    /// ```
    /// The product must not leave multiple states on the same wire.
    pub fn into_dot(self, rhs: Self) -> KBResult<Element> { self.dot(&rhs) }
}

fn subscript_str(n: usize) -> String {
    format!("{}", n)
        .replace('0', "₀")
        .replace('1', "₁")
        .replace('2', "₂")
        .replace('3', "₃")
        .replace('4', "₄")
        .replace('5', "₅")
        .replace('6', "₆")
        .replace('7', "₇")
        .replace('8', "₈")
        .replace('9', "₉")
}

impl fmt::Display for KetBra {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "(")?;
        self.ampl.fmt(f)?;
        write!(f, ")")?;
        if !self.ket.is_empty() {
            write!(f, "∣")?;
            let mut ket: Vec<(&usize, &State)> = self.ket.iter().collect();
            ket.sort_by(|(idx_l, _), (idx_r, _)| idx_l.cmp(idx_r));
            for (idx, s) in ket.into_iter() {
                write!(f, "{}{}", s, subscript_str(*idx))?;
            }
            write!(f, "⟩")?;
        }
        if !self.bra.is_empty() {
            write!(f, "⟨")?;
            let mut bra: Vec<(&usize, &State)> = self.bra.iter().collect();
            bra.sort_by(|(idx_l, _), (idx_r, _)| idx_l.cmp(idx_r));
            for (idx, s) in bra.into_iter() {
                write!(f, "{}{}", s, subscript_str(*idx))?;
            }
            write!(f, "∣")?;
        }
        Ok(())
    }
}

#[derive(Copy, Clone, Debug)]
pub(crate) enum ElementKind {
    Z(f64),
    X(f64),
    H(C64),
    Swap(usize, usize),
    Cup(usize, usize),
    Cap(usize, usize),
    Unknown,
}

/// Represents a single element of a diagram as a sum over [`KetBra`].
///
/// Note that while every named object in the ZX(H)-calculus (i.e., spiders,
/// cups, caps, and H-boxes) is an `Element` and can be created as such, not
/// every `Element` will be a named object.
#[derive(Clone, Debug)]
pub struct Element {
    kind: ElementKind,
    terms: Vec<KetBra>
}

impl Deref for Element {
    type Target = Vec<KetBra>;

    fn deref(&self) -> &Self::Target { &self.terms }
}

impl DerefMut for Element {
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.terms }
}

/// Iterator type for [`Element`]s over each [`KetBra`] in the sum.
#[derive(Clone)]
pub struct ElementIntoIter {
    iter: std::vec::IntoIter<KetBra>
}

impl Iterator for ElementIntoIter {
    type Item = KetBra;

    fn next(&mut self) -> Option<Self::Item> { self.iter.next() }
}

impl IntoIterator for Element {
    type Item = KetBra;
    type IntoIter = ElementIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        ElementIntoIter { iter: self.terms.into_iter() }
    }
}

impl From<KetBra> for Element {
    fn from(kb: KetBra) -> Self {
        Self { kind: ElementKind::Unknown, terms: vec![kb] }
    }
}

trait LogicalStateOrd {
    fn state_ord(&self, rhs: &Self) -> Ordering;
}

impl LogicalStateOrd for HashMap<usize, State> {
    fn state_ord(&self, rhs: &Self) -> Ordering {
        let mut common: Vec<&usize>
            = self.keys().collect::<HashSet<&usize>>()
            .union(&rhs.keys().collect::<HashSet<&usize>>())
            .copied()
            .collect();
        common.sort();
        for idx in common.into_iter() {
            match (self.get(idx), rhs.get(idx)) {
                (Some(ls), Some(rs)) => match ls.cmp(rs) {
                    Ordering::Equal => { continue; },
                    x => { return x; },
                },
                (Some(ls), None) => match *ls {
                    Zero => { continue; }
                    _ => { return Ordering::Greater; },
                },
                (None, Some(rs)) => match *rs {
                    Zero => { continue; }
                    _ => { return Ordering::Less; },
                },
                (None, None) => { continue; }
            }
        }
        Ordering::Equal
    }
}

impl LogicalStateOrd for KetBra {
    fn state_ord(&self, rhs: &Self) -> Ordering {
        self.bra.state_ord(&rhs.bra)
            .then(self.ket.state_ord(&rhs.ket))
    }
}

/// A single spider.
#[derive(Copy, Clone, Debug)]
pub enum Spider {
    Z(f64),
    X(f64),
}

impl Element {
    /// Create a new `Element` from a list of [`KetBra`]s.
    ///
    /// Each term in the list must have its ket and bra defined over the same
    /// wires as every other, otherwise this function fails.
    pub fn new<I>(ketbras: I) -> KBResult<Self>
    where I: IntoIterator<Item = KetBra>
    {
        let ketbras: Vec<KetBra> = ketbras.into_iter().collect();
        let ket_idx: HashSet<&usize>
            = ketbras.iter().flat_map(|kb| kb.ket.keys()).collect();
        ketbras.iter()
            .all(|kb| {
                let ket_keys: HashSet<&usize> = kb.ket.keys().collect();
                ket_keys == ket_idx
            })
            .then_some(())
            .ok_or(KBError::ElementKetSameKeys)?;
        let bra_idx: HashSet<&usize>
            = ketbras.iter().flat_map(|kb| kb.bra.keys()).collect();
        ketbras.iter()
            .all(|kb| {
                let bra_keys: HashSet<&usize> = kb.bra.keys().collect();
                bra_keys == bra_idx
            })
            .then_some(())
            .ok_or(KBError::ElementBraSameKeys)?;
        Ok(Self { kind: ElementKind::Unknown, terms: ketbras }.simplified())
    }

    /// If `self` contains only one term with an empty ket and bra, extract its
    /// amplitude as a single scalar.
    pub fn as_scalar(&self) -> Option<C64> {
        (self.terms.len() == 1)
            .then(|| self.terms[0].as_scalar())
            .flatten()
    }

    fn collect_unchecked<I>(ketbras: I) -> Self
    where I: IntoIterator<Item = KetBra>
    {
        Self {
            kind: ElementKind::Unknown,
            terms: ketbras.into_iter().collect(),
        }
    }

    /// Simplify `self` by folding terms with identical kets and bras into each
    /// other, summing their amplitudes.
    ///
    /// Terms with zero amplitude are removed at the end.
    pub fn simplify(&mut self) {
        fn has_term<'a>(acc: &'a mut [KetBra], term: &KetBra)
            -> Option<&'a mut KetBra>
        {
            acc.iter_mut()
                .find(|term_acc| {
                    term.ket == term_acc.ket && term.bra == term_acc.bra
                })
        }

        let mut acc: Vec<KetBra> = Vec::new();
        for term in self.terms.drain(..) {
            if let Some(term_acc) = has_term(&mut acc, &term) {
                term_acc.ampl += term.ampl;
            } else {
                acc.push(term);
            }
        }
        self.terms
            = acc.into_iter().filter(|kb| kb.ampl.norm() > 1e-9).collect();
    }

    /// Simplify `self` by folding terms with identical kets and braw into each
    /// other, summing their amplitudes, and return `self` at the end.
    ///
    /// Terms with zero amplitude are removed.
    pub fn simplified(mut self) -> Self { self.simplify(); self }

    /// Create a scalar (i.e. empty ket and bra) with amplitude `0`.
    pub fn zero() -> Self {
        Self {
            kind: ElementKind::Unknown,
            terms: vec![KetBra::new(c!(0.0), [], [])],
        }
    }

    /// Create a scalar (i.e. empty ket and bra) with amplitude `1`.
    pub fn one() -> Self {
        Self {
            kind: ElementKind::Unknown,
            terms: vec![KetBra::new(c!(1.0), [], [])],
        }
    }

    /// Create a scalar (i.e. empty ket and bra) with arbitrary amplitude.
    pub fn scalar(ampl: C64) -> Self {
        Self {
            kind: ElementKind::Unknown,
            terms: vec![KetBra::new(ampl, [], [])],
        }
    }

    /// Return a set of all input wire indices.
    pub fn ins(&self) -> HashSet<usize> {
        match self.kind {
            ElementKind::Swap(w1, w2) => [w1, w2].into_iter().collect(),
            ElementKind::Cup(_, _) => HashSet::default(),
            ElementKind::Cap(w1, w2) => [w1, w2].into_iter().collect(),
            _ => {
                self.terms.iter()
                    .flat_map(|term| term.bra.keys())
                    .copied()
                    .collect()
            },
        }
    }

    /// Return a set of all output wire indices.
    pub fn outs(&self) -> HashSet<usize> {
        match self.kind {
            ElementKind::Swap(w1, w2) => [w1, w2].into_iter().collect(),
            ElementKind::Cup(w1, w2) => [w1, w2].into_iter().collect(),
            ElementKind::Cap(_, _) => HashSet::default(),
            _ => {
                self.terms.iter()
                    .flat_map(|term| term.ket.keys())
                    .copied()
                    .collect()
            },
        }
    }

    /// Express all of `self`'s terms in `basis`, simplifying the results.
    pub fn to_basis(&self, basis: Basis) -> Self {
        let terms: Vec<KetBra>
            = self.terms.iter()
            .flat_map(|kb| kb.to_basis(basis))
            .collect();
        Self {
            kind: ElementKind::Unknown,
            terms,
        }.simplified()
    }

    /// Create the identity element on wire `idx`.
    pub fn id(idx: usize) -> Self { Self::Z([idx], [idx], Some(0.0)) }

    /// Create a Z-spider.
    ///
    /// ```math
    /// Z([i_0, \dots, i_n], [j_0, \dots, j_m], \alpha)
    ///     = \ket{0_{j_0}, \dots, 0_{j_n}} \bra{0_{i_0}, \dots, 0_{i_m}}
    ///     + e^{i \alpha}
    ///         \ket{1_{j_0}, \dots, 1_{j_n}} \bra{1_{i_0}, \dots, 1_{i_m}}
    /// ```
    ///
    /// The phase defaults to `0.0`. Duplicate wire indices are removed.
    pub fn Z<I, J>(ins: I, outs: J, phase: Option<f64>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let ins: HashSet<usize> = ins.into_iter().collect();
        let outs: HashSet<usize> = outs.into_iter().collect();
        Self {
            kind: ElementKind::Z(phase.unwrap_or(0.0)),
            terms: vec![
                KetBra::new(
                    c!(1.0),
                    outs.iter().map(|idx| (*idx, Zero)),
                    ins.iter().map(|idx| (*idx, Zero)),
                ),
                KetBra::new(
                    c!(e phase.unwrap_or(0.0)),
                    outs.iter().map(|idx| (*idx, One)),
                    ins.iter().map(|idx| (*idx, One)),
                ),
            ],
        }
    }

    /// Create an X-spider.
    ///
    /// ```math
    /// X([i_0, \dots, i_n], [j_0, \dots, j_m], \alpha)
    ///     = \ket{+_{j_0}, \dots, +_{j_n}} \bra{+_{i_0}, \dots, +_{i_m}}
    ///     + e^{i \alpha}
    ///         \ket{-_{j_0}, \dots, -_{j_n}} \bra{-_{i_0}, \dots, -_{i_m}}
    /// ```
    ///
    /// The phase defaults to `0.0`. Duplicate wire indices are removed.
    pub fn X<I, J>(ins: I, outs: J, phase: Option<f64>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let ins: HashSet<usize> = ins.into_iter().collect();
        let outs: HashSet<usize> = outs.into_iter().collect();
        Self {
            kind: ElementKind::X(phase.unwrap_or(0.0)),
            terms: vec![
                KetBra::new(
                    c!(1.0),
                    outs.iter().map(|idx| (*idx, Plus)),
                    ins.iter().map(|idx| (*idx, Plus)),
                ),
                KetBra::new(
                    c!(e phase.unwrap_or(0.0)),
                    outs.iter().map(|idx| (*idx, Minus)),
                    ins.iter().map(|idx| (*idx, Minus)),
                ),
            ],
        }
    }

    /// Create an H-box.
    ///
    /// ```math
    /// H([i_0, \dots, i_n], [j_0, \dots, j_m], a)
    ///     = \sum_{s_k \in \{0, 1\}}
    ///         a^{\prod s_k}
    ///             \ket{s_{i_0}, \dots, s_{i_n}} \bra{s_{j_0}, \dots, s_{j_m}}
    /// ```
    ///
    /// The parameter defaults to `-1.0`. Duplicate wire indices are removed.
    pub fn H<I, J>(ins: I, outs: J, a: Option<C64>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let a: C64 = a.unwrap_or(c!(-1.0));
        let ins: HashSet<usize> = ins.into_iter().collect();
        let outs: HashSet<usize> = outs.into_iter().collect();
        let terms: Vec<KetBra>
            = if ins.is_empty() {
                outs.into_iter()
                    .map(|idx| [(idx, Zero), (idx, One)])
                    .multi_cartesian_product()
                    .map(|outs| {
                        let ampl: C64
                            = outs.iter()
                            .any(|(_, s)| *s == Zero)
                            .then_some(c!(1.0))
                            .unwrap_or(a);
                        KetBra::new(ampl, outs, [])
                    })
                .collect()
            } else if outs.is_empty() {
                ins.into_iter()
                    .map(|idx| [(idx, Zero), (idx, One)])
                    .multi_cartesian_product()
                    .map(|ins| {
                        let ampl: C64
                            = ins.iter()
                            .any(|(_, s)| *s == Zero)
                            .then_some(c!(1.0))
                            .unwrap_or(a);
                        KetBra::new(ampl, [], ins)
                    })
                .collect()
            } else {
                Itertools::cartesian_product(
                    ins.into_iter()
                        .map(|idx| [(idx, Zero), (idx, One)])
                        .multi_cartesian_product(),
                    outs.into_iter()
                        .map(|idx| [(idx, Zero), (idx, One)])
                        .multi_cartesian_product(),
                )
                .map(|(ins, outs)| {
                    let ampl: C64
                        = ins.iter().chain(outs.iter())
                        .any(|(_, s)| *s == Zero)
                        .then_some(c!(1.0))
                        .unwrap_or(a);
                    KetBra::new(ampl, outs, ins)
                })
                .collect()
            };
        Self { kind: ElementKind::H(a), terms }
    }

    /// Create a swap.
    ///
    /// ```math
    /// \text{SWAP}([i_0, i_1])
    ///     = \ket{0_{i_0} 0_{i_1}} \bra{0_{i_0} 0_{i_1}}
    ///     + \ket{1_{i_0} 0_{i_1}} \bra{0_{i_0} 1_{i_1}}
    ///     + \ket{0_{i_0} 1_{i_1}} \bra{1_{i_0} 0_{i_1}}
    ///     + \ket{1_{i_0} 1_{i_1}} \bra{1_{i_0} 1_{i_1}}
    /// ```
    ///
    /// If the two given wires are equal, return the output of [`Self::id`]
    /// instead.
    pub fn swap(wires: [usize; 2]) -> Self {
        let i0: usize = wires[0];
        let i1: usize = wires[1];
        if i0 == i1 {
            return Self::id(i0);
        }
        let terms: Vec<KetBra>
            = vec![
                KetBra::new(
                    c!(1.0),
                    [(i0, Zero), (i1, Zero)],
                    [(i0, Zero), (i1, Zero)],
                ),
                KetBra::new(
                    c!(1.0),
                    [(i0, One ), (i1, Zero)],
                    [(i0, Zero), (i1, One )],
                ),
                KetBra::new(
                    c!(1.0),
                    [(i0, Zero), (i1, One )],
                    [(i0, One ), (i1, Zero)],
                ),
                KetBra::new(
                    c!(1.0),
                    [(i0, One ), (i1, One )],
                    [(i0, One ), (i1, One )],
                ),
            ];
        Self { kind: ElementKind::Swap(wires[0], wires[1]), terms }
    }

    /// Create a cup (Bell state).
    ///
    /// ```math
    /// \text{cup}([j_0, j_1])
    ///     = Z([], [j_0, j_1], 0)
    ///     = X([], [j_0, j_1], 0)
    ///     = \ket{0_{i_0} 0_{i_1}} + \ket{1_{i_0} 1_{i_1}}
    /// ```
    ///
    /// Defaults to the `$\ket{00} + \ket{11}$` Bell state if no phase is given.
    ///
    /// *Panics if the two given wires are equal.*
    pub fn cup(outs: [usize; 2], phase: Option<Spider>) -> Self {
        if outs[0] == outs[1] {
            panic!("Element::cup: wires cannot be equal");
        }
        match phase {
            Some(Spider::Z(ph)) => Element::Z([], outs, Some(ph)),
            Some(Spider::X(ph)) => Element::X([], outs, Some(ph)),
            None => {
                let mut elem = Self::Z([], outs, None);
                elem.kind = ElementKind::Cup(outs[0], outs[1]);
                elem
            },
        }
    }

    /// Create a cap (Bell effect).
    ///
    /// ```math
    /// \text{cap}([i_0, i_1])
    ///     = Z([i_0, i_1], [], 0)
    ///     = X([i_0, i_1], [], 0)
    ///     = \bra{0_{i_0} 0_{i_1}} + \bra{1_{i_0} 1_{i_1}}
    /// ```
    ///
    /// Defaults to the `$\ket{00} + \ket{11}$` Bell effect if no phase is
    /// given.
    ///
    /// *Panics if the two given wires are equal.*
    pub fn cap(ins: [usize; 2], phase: Option<Spider>) -> Self {
        if ins[0] == ins[1] {
            panic!("Element::cap: wires cannot be equal");
        }
        match phase {
            Some(Spider::Z(ph)) => Element::Z(ins, [], Some(ph)),
            Some(Spider::X(ph)) => Element::X(ins, [], Some(ph)),
            None => {
                let mut elem = Self::Z(ins, [], None);
                elem.kind = ElementKind::Cap(ins[0], ins[1]);
                elem
            },
        }
    }

    /// Alias for [`Self::new`] for semantics purposes.
    pub fn from_sum<I>(terms: I) -> KBResult<Self>
    where I: IntoIterator<Item = KetBra>
    {
        Self::new(terms)
    }

    /// Compute the Kronecker product of `self` and `rhs`.
    ///
    /// No inner products are evaluated; the set of wire indices contained
    /// within `self` and `rhs`'s kets and bras must be disjoint, otherwise this
    /// operation fails.
    pub fn kron(&self, rhs: &Self) -> KBResult<Self> {
        let terms: Vec<KetBra>
            = Itertools::cartesian_product(
                self.terms.iter(),
                rhs.terms.iter(),
            )
            .map(|(l, r)| l.kron(r))
            .collect::<KBResult<Vec<KetBra>>>()?;
        Self::new(terms)
    }

    /// Compute the Kronecker product of multiple `Element`s.
    ///
    /// No inner products are evaluated; the set of wire indices contained
    /// within `self` and `rhs`'s kets and bras must be disjoint, otherwise this
    /// operation fails.
    pub fn multikron<'a, I>(elems: I) -> KBResult<Self>
    where I: IntoIterator<Item = &'a Element>
    {
        let terms: Vec<KetBra>
            = elems.into_iter()
            .cloned()
            .multi_cartesian_product()
            .map(KetBra::into_multikron)
            .collect::<KBResult<Vec<KetBra>>>()?;
        Self::new(terms)
    }

    /// Compute the Kronecker product of multiple `Element`s, consuming all
    /// inputs.
    ///
    /// No inner products are evaluated; the set of wire indices contained
    /// within `self` and `rhs`'s kets and bras must be disjoint, otherwise this
    /// operation fails.
    pub fn into_multikron<I>(elems: I) -> KBResult<Self>
    where I: IntoIterator<Item = Element>
    {
        let terms: Vec<KetBra>
            = elems.into_iter()
            .multi_cartesian_product()
            .map(KetBra::into_multikron)
            .collect::<KBResult<Vec<KetBra>>>()?;
        Self::new(terms)
    }

    /// Compute the dot-product of `self` and `rhs.
    ///
    /// Unmatched kets and bras are left untouched and "passed through" the
    /// product. See also [`KetBra::dot`]. The product must not leave multiple
    /// states on the same wire.
    pub fn dot(&self, rhs: &Self) -> KBResult<Self> {
        let terms: Vec<KetBra>
            = Itertools::cartesian_product(
                self.terms.iter(),
                rhs.terms.iter(),
            )
            .map(|(l, r)| l.dot(r))
            .collect::<KBResult<Vec<Element>>>()?
            .into_iter()
            .flatten()
            .collect();
        Self::new(terms)
    }

    /// Compute the dot-product of multiple `Element`s.
    ///
    /// Unmatched kets and bras are left untouched and "passed through" the
    /// product. See also [`KetBra::dot`]. The product must not leave multiple
    /// states on the same wire.
    pub fn multidot<'a, I>(elems: I) -> KBResult<Self>
    where I: IntoIterator<Item = &'a Element>
    {
        let mut acc = Self::one();
        for elem in elems.into_iter() {
            acc = acc.dot(elem)?;
        }
        Ok(acc)
    }

    /// Compute the dot-product of multiple `Element`s, consuming all inputs.
    ///
    /// Unmatched kets and bras are left untouched and "passed through" the
    /// product. See also [`KetBra::dot`]. The product must not leave multiple
    /// states on the same wire.
    pub fn into_multidot<I>(elems: I) -> KBResult<Self>
    where I: IntoIterator<Item = Element>
    {
        let mut acc = Self::one();
        for elem in elems.into_iter() {
            acc = acc.dot(&elem)?;
        }
        Ok(acc)
    }

    /// Normalize the amplitudes of all terms in `self` such that the majority
    /// will be `1`.
    ///
    /// Relative magnitudes and phases between amplitudes are preserved.
    pub fn normalize(&mut self) {
        let mut ampl_counts: Vec<(C64, usize)> = Vec::new();
        self.terms.iter()
            .for_each(|term| {
                if let Some(k)
                    = ampl_counts.iter()
                    .enumerate()
                    .find_map(|(k, (a, _))| (*a == term.ampl).then_some(k))
                {
                    let (a, count) = ampl_counts.swap_remove(k);
                    ampl_counts.push((a, count + 1));
                } else {
                    ampl_counts.push((term.ampl, 1));
                }
            });
        if let Some((norm_ampl, _))
            = ampl_counts.into_iter()
            .max_by(|(ampl_l, count_l), (ampl_r, count_r)| {
                count_l.cmp(count_r)
                    .then(ampl_r.im.total_cmp(&ampl_l.im))
            })
        {
            self.terms.iter_mut()
                .for_each(|term| term.ampl /= norm_ampl);
        }
    }

    /// Normalize the amplitudes of all terms in `self` such that the majority
    /// will be `1`, and return `self`.
    ///
    /// Ratios and relative phases between amplitudes are preserved.
    pub fn normalized(mut self) -> Self { self.normalize(); self }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut terms: Vec<&KetBra> = self.terms.iter().collect();
        terms.sort_by(|l, r| l.state_ord(r));
        let n: usize = terms.len();
        for (k, term) in terms.iter().enumerate() {
            term.fmt(f)?;
            if k < n - 1 { write!(f, " + ")?; }
        }
        Ok(())
    }
}

/// Represents a series of [`Element`]s as "slices" of a ZX-diagram, each of
/// whose inputs and outputs locally span the entire lateral space of wires,
/// (i.e. such that every element's output wires match the input wires of the
/// next one).
#[derive(Clone, Debug)]
pub struct Diagram {
    slices: Vec<Element>,
    scalar: C64,
}

impl Deref for Diagram {
    type Target = Vec<Element>;

    fn deref(&self) -> &Self::Target { &self.slices }
}

impl DerefMut for Diagram {
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.slices }
}

impl FromIterator<Element> for Diagram {
    fn from_iter<I>(iter: I) -> Self
    where I: IntoIterator<Item = Element>
    {
        Self {
            slices: iter.into_iter().collect(),
            scalar: 1.0_f64.into(),
        }
    }
}

/// Iterator type for [`Diagram`]s over each [`Element`] slice.
pub struct DiagramIntoIter {
    iter: std::vec::IntoIter<Element>
}

impl Iterator for DiagramIntoIter {
    type Item = Element;

    fn next(&mut self) -> Option<Self::Item> { self.iter.next() }
}

impl IntoIterator for Diagram {
    type Item = Element;
    type IntoIter = DiagramIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        DiagramIntoIter { iter: self.slices.into_iter() }
    }
}

impl From<Element> for Diagram {
    fn from(elem: Element) -> Self {
        Self {
            slices: vec![elem],
            scalar: 1.0_f64.into()
        }
    }
}

impl Diagram {
    /// Create a new diagram from a list of [`Element`]s.
    ///
    /// `Element`s should follow the order in which they would be drawn in a
    /// real ZX-diagram, i.e. the reverse of the order in which dot-products are
    /// taken.
    pub fn new<I>(elems: I) -> Self
    where I: IntoIterator<Item = Element>
    {
         Self {
            slices: elems.into_iter().collect(),
            scalar: 1.0_f64.into(),
        }
    }

    /// Set an overall scalar factor.
    pub fn with_scalar(mut self, z: C64) -> Self {
        self.scalar = z;
        self
    }

    /// Set an overall scalar factor.
    pub fn set_scalar(&mut self, z: C64) -> &mut Self {
        self.scalar = z;
        self
    }

    /// Fold all slices of `self` into a single [`Element`] via the
    /// [dot-product][Element::dot], in the diagram order.
    pub fn contract(&self) -> KBResult<Element> {
        Element::multidot(self.slices.iter().rev())
            .map(|elem| {
                Element::collect_unchecked(
                    elem.into_iter()
                        .map(|mut kb| {
                            kb.ampl *= self.scalar;
                            kb
                        })
                )
            })
    }

    /// Compose `self` with `rhs` by attaching the outputs of `rhs` with the
    /// inputs of `self`.
    ///
    /// This is equivalent to appending `self` to the the end of `rhs`,
    /// analogous to the usual `$\circ$` composition operation.
    pub fn compose(&self, rhs: &Self) -> Self {
        Self {
            slices:
                rhs.slices.iter()
                .chain(self.slices.iter())
                .cloned()
                .collect(),
            scalar: self.scalar * rhs.scalar,
        }
    }

    /// Join `rhs` to the end of `self`, attaching the inputs of `rhs` to the
    /// outputs of `self.
    ///
    /// This is the reflexive operation to [`Self::compose`].
    pub fn join(&self, rhs: &Self) -> Self {
        Self {
            slices:
                self.slices.iter()
                .chain(rhs.slices.iter())
                .cloned()
                .collect(),
            scalar: self.scalar * rhs.scalar,
        }
    }

    /// Return sets of all input and output wire indices.
    pub fn ins_outs(&self) -> (HashSet<usize>, HashSet<usize>) {
        let mut ins: HashSet<usize> = HashSet::default();
        let mut outs: HashSet<usize> = HashSet::default();
        for element in self.slices.iter() {
            element.ins()
                .into_iter()
                .for_each(|k| {
                    if outs.contains(&k) {
                        outs.remove(&k);
                    } else {
                        ins.insert(k);
                    }
                });
            element.outs()
                .into_iter()
                .for_each(|k| { outs.insert(k); });
        }
        (ins, outs)
    }

    /// Convert `self` to a [`graph::Diagram`] representation.
    ///
    /// Note that the indices of input and output wires are not preserved, only
    /// their order.
    pub fn as_graph(&self) -> KBResult<graph::Diagram> {
        let mut graph = graph::Diagram::new();
        let mut node_id: graph::NodeId;
        let (ins, outs) = self.ins_outs();
        let mut inputs: Vec<usize> = ins.into_iter().collect();
        inputs.sort();
        let mut outputs: Vec<usize> = outs.into_iter().collect();
        outputs.sort();
        let mut wires
            = HashMap::<usize, graph::NodeId>::default(); // wire idx -> node idx
        let mut to_remove: Vec<graph::NodeId> = Vec::new(); // placeholder nodes

        // add inputs
        for idx in inputs.into_iter() {
            node_id = graph.add_node(graph::NodeDef::Input);
            wires.insert(idx, node_id);
        }

        // add elements/wires
        for elem in self.slices.iter() {
            match elem.kind {
                ElementKind::Z(_) | ElementKind::X(_) | ElementKind::H(_) => {
                    let node_def: graph::NodeDef
                        = match elem.kind {
                            ElementKind::Z(ph) => graph::NodeDef::Z(ph),
                            ElementKind::X(ph) => graph::NodeDef::X(ph),
                            ElementKind::H(a) => graph::NodeDef::H(a),
                            _ => unreachable!(),
                        };
                    node_id = graph.add_node(node_def);
                    for in_wire in elem.ins().into_iter() {
                        if let Some(prev_id) = wires.remove(&in_wire) {
                            graph.add_wire(prev_id, node_id).ok();
                        } else {
                            unreachable!()
                        }
                    }
                    for out_wire in elem.outs().into_iter() {
                        if wires.insert(out_wire, node_id).is_some() {
                            unreachable!()
                        }
                    }
                },
                ElementKind::Swap(w1, w2) => {
                    let Some(n1) = wires.remove(&w1) else { unreachable!() };
                    let Some(n2) = wires.remove(&w2) else { unreachable!() };
                    wires.insert(w1, n2);
                    wires.insert(w2, n1);
                },
                ElementKind::Cup(w1, w2) => {
                    let n1 = graph.add_node(graph::NodeDef::Z(0.0));
                    to_remove.push(n1);
                    wires.insert(w1, n1);
                    let n2 = graph.add_node(graph::NodeDef::Z(0.0));
                    to_remove.push(n2);
                    wires.insert(w2, n2);
                    graph.add_wire(n1, n2).ok();
                },
                ElementKind::Cap(w1, w2) => {
                    let Some(n1) = wires.remove(&w1) else { unreachable!() };
                    let Some(n2) = wires.remove(&w2) else { unreachable!() };
                    graph.add_wire(n1, n2).ok();
                },
                ElementKind::Unknown => {
                    return Err(KBError::DiagramConversionNonGenerator);
                },
            }
        }

        // add outputs
        for idx in outputs.into_iter() {
            node_id = graph.add_node(graph::NodeDef::Output);
            if let Some(prev_id) = wires.remove(&idx) {
                graph.add_wire(prev_id, node_id).ok();
            } else {
                unreachable!()
            }
        }

        // all nodes to remove here are placeholders for cups; they have
        // exactly two neighbors
        let mut n1: graph::NodeId;
        let mut n2: graph::NodeId;
        for remove in to_remove.into_iter() {
            let mut neighbors = graph.neighbors_of(remove).unwrap();
            n1 = neighbors.next().unwrap().0;
            n2 = neighbors.next().unwrap().0;
            graph.remove_node(remove);
            graph.add_wire(n1, n2).ok();
        }

        Ok(graph)
    }

    /// Return an object containing an encoding of `self` in the [dot
    /// language][dot-lang].
    ///
    /// Rendering this object using the default formatter will result in a full
    /// dot string representation of the diagram.
    ///
    /// [dot-lang]: https://en.wikipedia.org/wiki/DOT_(graph_description_language)
    pub fn to_graphviz(&self, name: &str) -> KBResult<tabbycat::Graph> {
        use tabbycat::*;
        use tabbycat::attributes::*;
        use crate::vizdefs::*;
        let mut id_gen = 0..;
        let mut node_id: usize;

        // initial declarations
        let mut statements
            = StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rankdir(RankDir::LR)),
            )
            .add_attr(
                AttrType::Node,
                AttrList::new()
                    .add_pair(fontname(FONT))
                    .add_pair(fontsize(FONTSIZE))
                    .add_pair(margin(NODE_MARGIN))
                    ,
            );
        let (ins, outs) = self.ins_outs();
        let mut wires = HashMap::<usize, usize>::default(); // wire idx -> node id

        // ensure all inputs are in a subgraph at the same rank
        let mut inputs_subgraph_stmt
            = StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rank(RankType::Source)),
            );
        let mut inputs: Vec<usize> = ins.into_iter().collect();
        inputs.sort();
        let mut prev: Option<usize> = None;
        for idx in inputs.into_iter() {
            node_id = id_gen.next().unwrap();
            wires.insert(idx, node_id);
            inputs_subgraph_stmt
                = inputs_subgraph_stmt.add_node(
                    node_id.into(),
                    None,
                    Some(
                        AttrList::new()
                            .add_pair(label(format!("In {}", idx)))
                            .add_pair(shape(Shape::Plaintext))
                    ),
                );
            if let Some(prev_idx) = prev {
                inputs_subgraph_stmt
                    = inputs_subgraph_stmt.add_edge(
                        Edge::head_node(
                            prev_idx.into(),
                            Some(Port::compass(Compass::South)),
                        )
                        .line_to_node(
                            node_id.into(),
                            Some(Port::compass(Compass::North)),
                        )
                        .add_attrpair(style(Style::Invisible))
                    );
            }
            prev = Some(node_id);
        }
        statements
            = statements.add_subgraph(SubGraph::cluster(inputs_subgraph_stmt));

        // add elements/wires
        let mut u_idx = 0..;
        for elem in self.slices.iter() {
            match elem.kind {
                ElementKind::Z(_)
                | ElementKind::X(_)
                | ElementKind::H(_)
                | ElementKind::Unknown
                => {
                    node_id = id_gen.next().unwrap();
                    let attrs
                        = match elem.kind {
                            ElementKind::Z(ph) => {
                                AttrList::new()
                                    .add_pair(label(format_phase(ph)))
                                    .add_pair(shape(Shape::Circle))
                                    .add_pair(height(CIRCLE_HEIGHT))
                                    .add_pair(style(Style::Filled))
                                    .add_pair(fillcolor(Z_COLOR))
                            },
                            ElementKind::X(ph) => {
                                AttrList::new()
                                    .add_pair(label(format_phase(ph)))
                                    .add_pair(shape(Shape::Circle))
                                    .add_pair(height(CIRCLE_HEIGHT))
                                    .add_pair(style(Style::Filled))
                                    .add_pair(fillcolor(X_COLOR))
                            },
                            ElementKind::H(a) => {
                                let a_label
                                    = if a == c!(-1.0) {
                                        "".to_string()
                                    } else {
                                        format!("{}", a)
                                    };
                                AttrList::new()
                                    .add_pair(label(a_label))
                                    .add_pair(shape(Shape::Square))
                                    .add_pair(height(SQUARE_HEIGHT))
                                    .add_pair(style(Style::Filled))
                                    .add_pair(fillcolor(H_COLOR))
                            },
                            ElementKind::Unknown => {
                                let u_label
                                    = format!(
                                        "U{}",
                                        subscript_str(u_idx.next().unwrap())
                                    );
                                AttrList::new()
                                    .add_pair(label(u_label))
                                    .add_pair(shape(Shape::Rectangle))
                                    .add_pair(style(Style::Filled))
                                    .add_pair(fillcolor(Color::White))
                            },
                            _ => unreachable!(),
                        };
                    statements
                        = statements.add_node(
                            node_id.into(), None, Some(attrs));
                    for in_wire in elem.ins().into_iter() {
                        if let Some(prev_id) = wires.remove(&in_wire) {
                            statements
                                = statements.add_edge(
                                    Edge::head_node(prev_id.into(), None)
                                        .line_to_node(node_id.into(), None)
                                );
                        } else {
                            unreachable!()
                        }
                    }
                    for out_wire in elem.outs().into_iter() {
                        if wires.insert(out_wire, node_id).is_some() {
                            unreachable!()
                        }
                    }
                },
                ElementKind::Swap(w1, w2) => {
                    let Some(n1) = wires.remove(&w1) else { unreachable!() };
                    let Some(n2) = wires.remove(&w2) else { unreachable!() };
                    wires.insert(w1, n2);
                    wires.insert(w2, n1);
                },
                ElementKind::Cup(w1, w2) => {
                    let n1 = id_gen.next().unwrap();
                    let attrs1
                        = AttrList::new().add_pair(style(Style::Invisible));
                    statements
                        = statements.add_node(n1.into(), None, Some(attrs1));
                    wires.insert(w1, n1);
                    let n2 = id_gen.next().unwrap();
                    let attrs2
                        = AttrList::new().add_pair(style(Style::Invisible));
                    statements
                        = statements.add_node(n2.into(), None, Some(attrs2));
                    wires.insert(w2, n2);
                    statements
                        = statements.add_edge(
                            Edge::head_node(
                                n1.into(),
                                None,
                            )
                            .line_to_node(
                                n2.into(),
                                None,
                            )
                        );
                },
                ElementKind::Cap(w1, w2) => {
                    let Some(n1) = wires.remove(&w1) else { unreachable!() };
                    let Some(n2) = wires.remove(&w2) else { unreachable!() };
                    statements
                        = statements.add_edge(
                            Edge::head_node(
                                n1.into(),
                                None,
                            )
                            .line_to_node(
                                n2.into(),
                                None,
                            )
                        );
                },
            }
        }

        // ensure all outputs are in a subgraph at the same rank
        let mut outputs_subgraph_stmt
            = StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rank(RankType::Sink)),
            );
        let mut outputs: Vec<usize> = outs.into_iter().collect();
        outputs.sort();
        let mut prev: Option<usize> = None;
        for idx in outputs.into_iter() {
            node_id = id_gen.next().unwrap();
            outputs_subgraph_stmt
                = outputs_subgraph_stmt.add_node(
                    node_id.into(),
                    None,
                    Some(
                        AttrList::new()
                            .add_pair(label(format!("Out {}", idx)))
                            .add_pair(shape(Shape::Plaintext))
                    ),
                );
            if let Some(prev_idx) = prev {
                outputs_subgraph_stmt
                    = outputs_subgraph_stmt.add_edge(
                        Edge::head_node(
                            prev_idx.into(),
                            Some(Port::compass(Compass::South)),
                        )
                        .line_to_node(
                            node_id.into(),
                            Some(Port::compass(Compass::North)),
                        )
                        .add_attrpair(style(Style::Invisible))
                    );
            }
            prev = Some(node_id);
            if let Some(prev_id) = wires.remove(&idx) {
                statements
                    = statements.add_edge(
                        Edge::head_node(prev_id.into(), None)
                            .line_to_node(node_id.into(), None)
                    );
            } else {
                unreachable!()
            }
        }
        statements
            = statements.add_subgraph(SubGraph::cluster(outputs_subgraph_stmt));

        GraphBuilder::default()
            .graph_type(GraphType::Graph)
            .strict(false)
            .id(Identity::quoted(name))
            .stmts(statements)
            .build()
            .map_err(KBError::GraphVizError)
    }

    /// Like [`to_graphviz`][Self::to_graphviz], but render directly to a string
    /// and write it to `path`.
    pub fn save_graphviz<P>(&self, name: &str, path: P) -> KBResult<()>
    where P: AsRef<Path>
    {
        let graphviz = self.to_graphviz(name)?;
        fs::OpenOptions::new()
            .write(true)
            .append(false)
            .create(true)
            .truncate(true)
            .open(path)?
            .write_all(format!("{}", graphviz).as_bytes())?;
        Ok(())
    }
}

impl fmt::Display for Diagram {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (ins, outs) = self.ins_outs();
        if !ins.is_empty() {
            let mut ins: Vec<usize> = ins.into_iter().collect();
            ins.sort();
            writeln!(f, "ins: {:?}", ins)?;
        } else {
            writeln!(f, "ins: [ ]")?;
        }
        if !self.slices.is_empty() {
            writeln!(f, "slices: {{")?;
            for slice in self.slices.iter() {
                write!(f, "  ")?;
                slice.fmt(f)?;
                writeln!(f)?;
            }
            writeln!(f, "}}")?;
        } else {
            writeln!(f, "slices: {{ }}")?;
        }
        if !outs.is_empty() {
            let mut outs: Vec<usize> = outs.into_iter().collect();
            outs.sort();
            write!(f, "outs: {:?}", outs)?;
        } else {
            write!(f, "outs: [ ]")?;
        }
        Ok(())
    }
}

