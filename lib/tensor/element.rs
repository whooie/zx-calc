use num_complex::Complex64 as C64;
use crate::phase::Phase;
use super::dense::{ De, Q, Tensor };
use super::sparse::{ Sp, State, Basis, KetBra, Terms };

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

/// `ElementData` is a sealed traits to govern tensor storage types -- we want
/// to take advantage of the type system to make `Element` (and extensions
/// thereof) generic over storage types, but ensure strict relationships between
/// constructor methods and `Kind` without having to make `Kind` public.
pub(crate) mod private { pub trait ElementDataSeal { } }
pub(crate) use private::ElementDataSeal;

/// Trait for storage types to back [`Element`].
///
/// This trait is sealed to ensure strict relationships between [`Element`] and
/// its backing data.
pub trait ElementData: Sized + ElementDataSeal {
    /// Associated error type for general operations.
    type Error;

    /// Returned by [`Self::input_iter`], giving the indices of all input wires.
    type InputIter<'a>: Iterator<Item = usize> where Self: 'a;

    /// Returned by [`Self::output_iter`], giving the indices of all output
    /// wires.
    type OutputIter<'a>: Iterator<Item = usize> where Self: 'a;

    /// Return an iterator over all input (bra) wire indices.
    fn input_iter(&self) -> Self::InputIter<'_>;

    /// Return an iterator over all output (ket) wire indices.
    fn output_iter(&self) -> Self::OutputIter<'_>;

    /// Return `true` if the `k`-th wire index connects to `self` as an input
    /// (bra).
    fn has_input(&self, k: usize) -> bool {
        self.input_iter().any(|j| j == k)
    }

    /// Return `true` if the `k`-th wire index connects to `self` as an output
    /// (ket).
    fn has_output(&self, k: usize) -> bool {
        self.output_iter().any(|j| j == k)
    }

    /// Construct from an arbitrary scalar.
    fn scalar<T>(a: T) -> Self
    where T: Into<C64>;

    /// If `self` is a scalar, extract its value.
    fn as_scalar(&self) -> Option<C64>;

    /// Return `true` if `self` is a scalar (rank-0 tensor).
    fn is_scalar(&self) -> bool { self.as_scalar().is_some() }

    /// Return an identity element, corresponding to an empty wire.
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     id(<i>i</i>) =
    ///       ∣0<sub><i>i</i></sub>⟩⟨0<sub><i>i</i></sub>∣
    ///       + ∣1<sub><i>i</i></sub>⟩⟨1<sub><i>i</i></sub>∣
    ///   </p>
    /// </blockquote>
    fn id(i: usize) -> Self;

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
    fn z<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>;

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
    fn x<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>;

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
    /// If the total number of input and output wires is 2 and `a` is –1, an
    /// extra scalar factor of 1/√2 should be added so that the result will
    /// coincide with the usual definition of the Hadamard gate.
    fn h<I, J>(ins: I, outs: J, a: Option<C64>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>;

    /// Create a swap.
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     swap(<i>i</i><sub>0</sub>, <i>i</i><sub>1</sub>) =
    ///       ∣0<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>⟩⟨0<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>∣
    ///       + ∣0<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>⟩⟨1<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>∣
    ///       + ∣1<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>⟩⟨0<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>∣
    ///       + ∣1<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>⟩⟨1<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>∣
    ///   </p>
    /// </blockquote>
    fn swap(i0: usize, i1: usize) -> Self;

    /// Create a cup (Bell state).
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     cup(<i>i</i><sub>0</sub>, <i>i</i><sub>1</sub>) =
    ///       ∣0<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>⟩
    ///       + ∣1<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>⟩
    ///   </p>
    /// </blockquote>
    fn cup(i0: usize, i1: usize) -> Self;

    /// Create a cap (Bell effect).
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     cap(<i>i</i><sub>0</sub>, <i>i</i><sub>1</sub>) =
    ///       ⟨0<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>∣
    ///       + ⟨1<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>∣
    ///   </p>
    /// </blockquote>
    fn cap(i0: usize, i1: usize) -> Self;

    /// Compute the dot product `self · rhs`, consuming both.
    fn dot(self, rhs: Self) -> Result<Self, Self::Error>;

    /// Multiply by a scalar, modifying `self` in place.
    fn scalar_mul<C>(&mut self, scalar: C)
    where C: Into<C64>;

    /// Conjugate `self` in place, converting kets to bras and bras to kets, and
    /// taking all tensor elements to their complex conjugates.
    fn conjugate(&mut self);
}

/// A lateral slice of a diagram.
///
/// Strictly, an `Element` is a single tensor parameterized by a particular
/// storage type. This crate offers two: `De`, a dense array backed by a
/// [`Tensor`]; and `Sp`, a sparse array backed by a vector of [`KetBra`]s.
#[derive(Copy, Clone, Debug)]
pub struct Element<A> {
    pub(crate) kind: Kind, // used only for diagram-rendering purposes
    pub(crate) data: A,
}

/// An [`Element`] with dense tensor storage.
pub type ElementDe = Element<De>;

/// An [`Element`] with sparse tensor storage.
pub type ElementSp = Element<Sp>;

impl<A> From<C64> for Element<A>
where A: From<C64>
{
    fn from(a: C64) -> Self {
        let kind = Kind::Scalar(a);
        let data = A::from(a);
        Self { kind, data }
    }
}

impl<A> Element<A>
where A: ElementData
{
    /// Return an iterator over all input (bra) wire indices.
    pub fn input_iter(&self) -> A::InputIter<'_> { self.data.input_iter() }

    /// Return an iterator over all output (ket) wire indices.
    pub fn output_iter(&self) -> A::OutputIter<'_> { self.data.output_iter() }

    /// Return `true` if the `k`-th wire index connects to `self` as an input
    /// (bra).
    pub fn has_input(&self, k: usize) -> bool { self.data.has_input(k) }

    /// Return `true` if the `k`-th wire index connects to `self` as an output
    /// (ket).
    pub fn has_output(&self, k: usize) -> bool { self.data.has_output(k) }

    /// Construct from an arbitrary scalar.
    pub fn scalar<T>(a: T) -> Self
    where T: Into<C64>
    {
        let a = a.into();
        let kind = Kind::Scalar(a);
        let data = A::scalar(a);
        Self { kind, data }
    }

    /// If `self` is a scalar, extract its value.
    pub fn as_scalar(&self) -> Option<C64> { self.data.as_scalar() }

    /// Return `true` if `self` is a scalar (rank-0 tensor).
    pub fn is_scalar(&self) -> bool { self.data.is_scalar() }

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
        let data = A::id(i);
        Self { kind, data }
    }

    /// Create a Z-spider.
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     <i>Z</i>({<i>i</i><sub>0</sub>, ..., <i>i</i><sub><i>n</i></sub>}, {<i>j</i><sub>0</sub>, ..., <i>j</i><sub><i>m</i></sub>}, <i>α</i>) =
    ///       ∣0<sub><i>j</i><sub>0</sub></sub>, ..., 0<sub><i>j</i><sub><i>m</i></sub></sub>⟩⟨0<sub><i>i</i><sub>0</sub></sub>, ..., 0<sub><i>i</i><sub><i>n</i></sub></sub>∣
    ///       + <i>e</i><sup><i>iα</i></sup>
    ///         ∣1<sub><i>j</i><sub>0</sub></sub>, ..., 1<sub><i>j</i><sub><i>m</i></sub></sub>⟩⟨1<sub><i>i</i><sub>0</sub></sub>, ..., 1<sub><i>i</i><sub><i>n</i></sub></sub>∣
    ///   </p>
    /// </blockquote>
    ///
    /// `phase` defaults to 0 and duplicate wire indices are ignored.
    pub fn z<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let kind = Kind::Z(phase.unwrap_or_else(Phase::zero));
        let data = A::z(ins, outs, phase);
        Self { kind, data }
    }

    /// Create an X-spider.
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     <i>X</i>({<i>i</i><sub>0</sub>, ..., <i>i</i><sub><i>n</i></sub>}, {<i>j</i><sub>0</sub>, ..., <i>j</i><sub><i>m</i></sub>}, <i>α</i>) =
    ///       ∣+<sub><i>j</i><sub>0</sub></sub>, ..., +<sub><i>j</i><sub><i>m</i></sub></sub>⟩⟨+<sub><i>i</i><sub>0</sub></sub>, ..., +<sub><i>i</i><sub><i>n</i></sub></sub>∣
    ///       + <i>e</i><sup><i>iα</i></sup>
    ///         ∣–<sub><i>j</i><sub>0</sub></sub>, ..., –<sub><i>j</i><sub><i>m</i></sub></sub>⟩⟨–<sub><i>i</i><sub>0</sub></sub>, ..., –<sub><i>i</i><sub><i>n</i></sub></sub>∣
    ///   </p>
    /// </blockquote>
    ///
    /// `phase` defaults to 0 and duplicate wire indices are ignored.
    pub fn x<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let kind = Kind::X(phase.unwrap_or_else(Phase::zero));
        let data = A::x(ins, outs, phase);
        Self { kind, data }
    }

    /// Create an H-box.
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     <i>H</i>({<i>i</i><sub>0</sub>, ..., <i>i</i><sub><i>n</i></sub>}, {<i>j</i><sub>0</sub>, ..., <i>j</i><sub><i>m</i></sub>}, <i>a</i>) =
    ///       Σ<sub><i>s</i><sub><i>k</i></sub> ∊ {0, 1}</sub>
    ///         <i>a</i><sup>Π<sub><i>k</i></sub> <i>s</i><sub><i>k</i></sub></sup>
    ///           ∣<i>s</i><sub><i>j</i><sub>0</sub></sub>, ..., <i>s</i><sub><i>j</i><sub><i>m</i></sub></sub>⟩⟨<i>s</i><sub><i>i</i><sub>0</sub></sub>, ..., <i>s</i><sub><i>i</i><sub><i>n</i></sub></sub>∣
    ///   </p>
    /// </blockquote>
    ///
    /// If the total number of input and output wires is 2 and `a` is –1, an
    /// extra scalar factor of 1/√2 is added so that the result will coincide
    /// with the usual definition of the Hadamard gate. `a` defaults to –1 and
    /// duplicate wire indices are ignored.
    pub fn h<I, J>(ins: I, outs: J, a: Option<C64>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let kind = Kind::H(a.unwrap_or_else(|| -C64::from(1.0)));
        let data = A::h(ins, outs, a);
        Self { kind, data }
    }

    /// Create a swap.
    ///
    /// <blockquote>
    ///   <p style="font-size:20px">
    ///     swap(<i>i</i><sub>0</sub>, <i>i</i><sub>1</sub>) =
    ///       ∣0<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>⟩⟨0<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>∣
    ///       + ∣0<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>⟩⟨1<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>∣
    ///       + ∣1<sub><i>i</i><sub>0</sub></sub>, 0<sub><i>i</i><sub>1</sub></sub>⟩⟨0<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>∣
    ///       + ∣1<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>⟩⟨1<sub><i>i</i><sub>0</sub></sub>, 1<sub><i>i</i><sub>1</sub></sub>∣
    ///   </p>
    /// </blockquote>
    ///
    /// If `i0 == i1`, an [`id`][Self::id] is returned.
    pub fn swap(i0: usize, i1: usize) -> Self {
        if i0 == i1 { return Self::id(i0); }
        let kind = Kind::Swap(i0, i1);
        let data = A::swap(i0, i1);
        Self { kind, data }
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
    /// If `i0 == i1`, return a unary [`z`][Self::z] with phase 0.
    pub fn cup(i0: usize, i1: usize) -> Self {
        if i0 == i1 { return Self::z([], [i0], None); }
        let kind = Kind::Cup(i0, i1);
        let data = A::cup(i0, i1);
        Self { kind, data }
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
    /// If `i0 == i1`, return a unary [`z`][Self::z] with phase 0.
    pub fn cap(i0: usize, i1: usize) -> Self {
        if i0 == i1 { return Self::z([i0], [], None); }
        let kind = Kind::Cap(i0, i1);
        let data = A::cap(i0, i1);
        Self { kind, data }
    }

    /// Compute the dot product `self · rhs`, consuming both.
    pub fn dot(self, rhs: Self) -> Result<Self, A::Error> {
        let data = self.data.dot(rhs.data)?;
        let kind =
            if let Some(a) = data.as_scalar() {
                Kind::Scalar(a)
            } else {
                Kind::Unknown
            };
        Ok(Self { kind, data })
    }

    /// Thin wrapper around [`dot`][Self::dot] with the operands reversed (i.e.
    /// computing `rhs · self`) in order to conform to a left-to-right style of
    /// composition.
    pub fn then(self, rhs: Self) -> Result<Self, A::Error> { rhs.dot(self) }

    /// Multiply by a scalar, modifying `self` in place.
    pub fn scalar_mul_mut<C>(&mut self, scalar: C)
    where C: Into<C64>
    {
        const EPSILON: f64 = 1e-12;
        let scalar = scalar.into();
        if (scalar - 1.0).norm() < EPSILON { return; }
        self.data.scalar_mul(scalar);
        self.kind = Kind::Unknown;
    }

    /// Multiply by a scalar, consuming `self` and returning the result.
    pub fn scalar_mul<C>(mut self, scalar: C) -> Self
    where C: Into<C64>
    {
        self.scalar_mul_mut(scalar);
        self
    }

    /// Conjugate `self` in place, converting kets to bras and bras to kets, and
    /// taking all tensor elements to their complex conjugates.
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
        self.data.conjugate();
    }

    /// Return the conjugate of `self`, converting kets to bras and bras to
    /// kets, and taking all tensor elements to their complex conjugates.
    pub fn adjoint(mut self) -> Self {
        self.adjoint_mut();
        self
    }
}

impl<A> std::fmt::Display for Element<A>
where A: std::fmt::Display
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f)
    }
}

impl From<Tensor> for Element<De> {
    fn from(data: Tensor) -> Self {
        let kind =
            if let Some(a) = data.as_scalar() {
                Kind::Scalar(a)
            } else {
                Kind::Unknown
            };
        let data = De::from(data);
        Self { kind, data }
    }
}

impl Element<De> {
    /// Create a new `Element` from a function over index values.
    ///
    /// This function is meant to generate tensor elements from an iterator over
    /// all possible index values of all (unique) indices passed to this
    /// function. The positions of the index values passed to the generating
    /// function correspond to the *unique* indices passed to this function *in
    /// the order in which they're passed*.
    ///
    /// Each index corresponds to a qubit wire; hence every index ranges over
    /// values {0, 1}.
    ///
    /// # Example
    /// TODO
    pub fn new_dense<I, F>(indices: I, elems: F) -> Self
    where
        I: IntoIterator<Item = Q>,
        F: FnMut(&[usize]) -> C64,
    {
        let kind = Kind::Unknown;
        let data = De::new_elems(indices, elems);
        Self { kind, data }
    }

    /// Sort axes in place so that wire indices are placed in ascending order,
    /// with all kets to the left of all bras.
    pub fn sort_indices(&mut self) { self.data.sort_indices() }

    /// Sort axes so that wire indices are placed in ascending order, with all
    /// kets to the left of all bras, consuming `self`.
    pub fn sorted_indices(mut self) -> Self {
        self.sort_indices();
        self
    }

    /// Return `true` if `self` and `other` denote the same tensor.
    ///
    /// Specifically, `self` and `other` must have identical indices and the
    /// modulus of the difference between any two corresponding tensor elements
    /// must be less than `1e-12`. Mutable references are required to sort the
    /// indices of both before the comparison is made.
    pub fn approx_eq(&mut self, other: &mut Self) -> bool {
        self.data.approx_eq(&mut other.data, Some(1e-12))
    }

    /// Convert to a [sparse][Sp] representation.
    pub fn to_sparse(&self) -> Element<Sp> {
        match self.kind {
            Kind::Scalar(a) => Element::scalar(a),
            Kind::Id(i) => Element::id(i),
            Kind::Z(ph) =>
                Element::z(self.input_iter(), self.output_iter(), Some(ph)),
            Kind::X(ph) =>
                Element::x(self.input_iter(), self.output_iter(), Some(ph)),
            Kind::H(a) =>
                Element::h(self.input_iter(), self.output_iter(), Some(a)),
            Kind::Swap(a, b) => Element::swap(a, b),
            Kind::Cup(a, b) => Element::cup(a, b),
            Kind::Cap(a, b) => Element::cap(a, b),
            Kind::Unknown => {
                let data = self.data.to_sparse();
                if let Some(a) = data.as_scalar() {
                    Element { kind: Kind::Scalar(a), data }
                } else {
                    Element { kind: Kind::Unknown, data }
                }
            }
        }
    }
}

impl From<KetBra> for Element<Sp> {
    fn from(data: KetBra) -> Self {
        let kind =
            if let Some(a) = data.as_scalar() {
                Kind::Scalar(a)
            } else {
                Kind::Unknown
            };
        let data = Sp::from(data);
        Self { kind, data }
    }
}

impl Element<Sp> {
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
    pub fn new_sparse<I, J, F>(ins: I, outs: J, ampl: F, basis: Basis) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
        F: FnMut(&[State], &[State]) -> Option<C64>,
    {
        let kind = Kind::Unknown;
        let data = Sp::new_terms(ins, outs, ampl, basis);
        Self { kind, data }
    }

    /// Return an iterator over all ket-bra terms.
    ///
    /// All terms are guaranteed to have states on the same set of input and
    /// output wire indices.
    ///
    /// The iterator item type is `&`[`KetBra`].
    pub fn terms(&self) -> Terms<'_> { self.data.terms() }

    /// Sort terms in place by their indices and states lexicographically.
    pub fn sort_terms(&mut self) { self.data.sort_terms(); }

    /// Sort terms by their indices and states lexicographically, consuming
    /// `self`.
    pub fn sorted_terms(mut self) -> Self {
        self.sort_terms();
        self
    }

    /// Combine like terms in place, ensuring minimal storage without changing
    /// basis.
    pub fn simplify(&mut self) { self.data.simplify(); }

    /// Combine like terms, ensuring minimal storage without changing basis,
    /// consuming `self`.
    pub fn simplified(mut self) -> Self {
        self.simplify();
        self
    }

    /// Re-express all terms in a single basis, in place.
    ///
    /// Terms are [simplified][Self::simplify] afterward.
    pub fn make_basis(&mut self, basis: Basis) { self.data.make_basis(basis); }

    /// Re-express all terms in a single basis, consuming `self`.
    ///
    /// Terms are [simplified][Self::simplify] afterward.
    pub fn into_basis(mut self, basis: Basis) -> Self {
        self.make_basis(basis);
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
        self.data.approx_eq(&mut other.data, thresh)
    }

    /// Convert to a [dense][`De`] representation.
    pub fn to_dense(&self) -> Element<De> {
        match self.kind {
            Kind::Scalar(a) => Element::scalar(a),
            Kind::Id(i) => Element::id(i),
            Kind::Z(ph) =>
                Element::z(self.input_iter(), self.output_iter(), Some(ph)),
            Kind::X(ph) =>
                Element::x(self.input_iter(), self.output_iter(), Some(ph)),
            Kind::H(a) =>
                Element::h(self.input_iter(), self.output_iter(), Some(a)),
            Kind::Swap(a, b) => Element::swap(a, b),
            Kind::Cup(a, b) => Element::cup(a, b),
            Kind::Cap(a, b) => Element::cap(a, b),
            Kind::Unknown => {
                let data = self.data.to_dense();
                if let Some(a) = data.as_scalar() {
                    Element { kind: Kind::Scalar(a), data }
                } else {
                    Element { kind: Kind::Unknown, data }
                }
            }
        }
    }
}

