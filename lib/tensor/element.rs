use std::fmt;
use ndarray::{ self as nd, Dimension };
use num_complex::Complex64 as C64;
use crate::{
    phase::Phase,
    tensor::{ Q, Tensor, TensorResult, KetIndices, BraIndices },
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

/// Represents a single element of a diagram as a tensor.
///
/// Note what while every named object in the ZX(H)-calculus (i.e. spiders,
/// cups, caps, and H-boxes) is an `Element` and can be created as such, not
/// every `Element` will be a named object.
#[derive(Clone, Debug)]
pub struct Element {
    pub(crate) kind: Kind, // used only for rendering purposes
    pub(crate) data: Tensor,
}

impl From<Tensor> for Element {
    fn from(data: Tensor) -> Self {
        if let Some(a) = data.as_scalar() {
            Self { kind: Kind::Scalar(a), data }
        } else {
            Self { kind: Kind::Unknown, data }
        }
    }
}

impl From<C64> for Element {
    fn from(val: C64) -> Self {
        Self { kind: Kind::Scalar(val), data: val.into() }
    }
}

impl Element {
    /// Create a new `Element` from an arbitrary tensor.
    pub fn new(data: Tensor) -> Self { Self::from(data) }

    /// Create a new `Element` from an arbitrary scalar.
    pub fn new_scalar<T>(val: T) -> Self
    where T: Into<C64>
    {
        Self::from(val.into())
    }

    /// Return `true` if `self` is a 0-ary generator or rank-0 tensor.
    pub fn is_scalar(&self) -> bool { self.data.is_scalar() }

    /// If `self` is a 0-ary generator or rank-0 tensor, extract its value as a
    /// single scalar.
    pub fn as_scalar(&self) -> Option<C64> { self.data.as_scalar() }

    /// Return an iterator over all input wire indices.
    pub fn ins(&self) -> BraIndices<'_> { self.data.bra_indices() }

    /// Return an iterator over all output wire indices.
    pub fn outs(&self) -> KetIndices<'_> { self.data.ket_indices() }

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
        let indices: Vec<Q> = vec![Q::Ket(i), Q::Bra(i)];
        let array: nd::Array2<C64> = nd::Array2::eye(2);
        let data = Tensor::from_array_unchecked(indices, array);
        let kind = Kind::Id(i);
        Self { kind, data }
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
        let mut indices: Vec<Q> = Vec::new();
        let phase = phase.unwrap_or_else(Phase::zero);
        let kind = Kind::Z(phase);
        outs.into_iter()
            .for_each(|idx| {
                let q = Q::Ket(idx);
                if !indices.contains(&q) { indices.push(q); }
            });
        ins.into_iter()
            .for_each(|idx| {
                let q = Q::Bra(idx);
                if !indices.contains(&q) { indices.push(q); }
            });
        if indices.is_empty() {
            let data = Tensor::new_scalar(1.0 + phase.cis());
            Self { kind, data }
        } else {
            let shape: Vec<usize> = vec![2; indices.len()];
            let mut array: nd::ArrayD<C64> = nd::ArrayD::zeros(shape);
            let elem = array.first_mut().unwrap();
            *elem = 1.0.into();
            let elem = array.last_mut().unwrap();
            *elem = phase.cis();
            let data = Tensor::from_array_unchecked(indices, array);
            Self { kind, data }
        }
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
        let mut indices: Vec<Q> = Vec::new();
        let phase = phase.unwrap_or_else(Phase::zero);
        let kind = Kind::X(phase);
        outs.into_iter()
            .for_each(|idx| {
                let q = Q::Ket(idx);
                if !indices.contains(&q) { indices.push(q); }
            });
        ins.into_iter()
            .for_each(|idx| {
                let q = Q::Bra(idx);
                if !indices.contains(&q) { indices.push(q); }
            });
        if indices.is_empty() {
            let data = Tensor::new_scalar(1.0 + phase.cis());
            Self { kind, data }
        } else {
            let sqrt_2_pow = indices.len() as i32;
            let norm = std::f64::consts::FRAC_1_SQRT_2.powi(sqrt_2_pow);
            let ph = phase.cis();
            let data =
                Tensor::new_unchecked(
                    indices,
                    |idx| {
                        let num_ones: usize = idx.iter().copied().sum();
                        let sign: f64 =
                            if num_ones % 2 == 1 { -1.0 } else { 1.0 };
                        (1.0 + sign * ph) * norm
                    },
                );
            Self { kind, data }
        }
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
        let mut indices: Vec<Q> = Vec::new();
        let arg = a.unwrap_or_else(|| -C64::from(1.0));
        let kind = Kind::H(arg);
        outs.into_iter()
            .for_each(|idx| {
                let q = Q::Ket(idx);
                if !indices.contains(&q) { indices.push(q); }
            });
        ins.into_iter()
            .for_each(|idx| {
                let q = Q::Bra(idx);
                if !indices.contains(&q) { indices.push(q); }
            });
        if indices.is_empty() {
            let data = Tensor::new_scalar(arg);
            Self { kind, data }
        } else {
            let norm =
                if indices.len() == 2 && (arg + 1.0).norm() < EPSILON {
                    std::f64::consts::FRAC_1_SQRT_2
                } else {
                    1.0
                };
            let data =
                Tensor::new_unchecked(
                    indices,
                    |idx| {
                        if idx.iter().all(|ix| *ix == 1) {
                            arg * norm
                        } else {
                            norm.into()
                        }
                    },
                );
            Self { kind, data }
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
        if i0 == i1 {
            Self::id(i0)
        } else {
            let kind = Kind::Swap(i0, i1);
            let indices: Vec<Q> =
                vec![Q::Ket(i0), Q::Ket(i1), Q::Bra(i0), Q::Bra(i1)];
            let mut array: nd::Array4<C64> = nd::Array4::zeros((2, 2, 2, 2));
            array[[0, 0, 0, 0]] = 1.0.into();
            array[[0, 1, 1, 0]] = 1.0.into();
            array[[1, 0, 0, 1]] = 1.0.into();
            array[[1, 1, 1, 1]] = 1.0.into();
            let data = Tensor::from_array_unchecked(indices, array);
            Self { kind, data }
        }
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
        if i0 == i1 { panic!("Element::cup: wires cannot be equal"); }
        let kind = Kind::Cup(i0, i1);
        let indices: Vec<Q> = vec![Q::Ket(i0), Q::Ket(i1)];
        let mut array: nd::Array2<C64> = nd::Array2::zeros((2, 2));
        array[[0, 0]] = 1.0.into();
        array[[1, 1]] = 1.0.into();
        let data = Tensor::from_array_unchecked(indices, array);
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
    /// *Panics if `i0 == i1`.*
    pub fn cap(i0: usize, i1: usize) -> Self {
        if i0 == i1 { panic!("Element::cap: wires cannot be equal"); }
        let kind = Kind::Cup(i0, i1);
        let indices: Vec<Q> = vec![Q::Bra(i0), Q::Bra(i1)];
        let mut array: nd::Array2<C64> = nd::Array2::zeros((2, 2));
        array[[0, 0]] = 1.0.into();
        array[[1, 1]] = 1.0.into();
        let data = Tensor::from_array_unchecked(indices, array);
        Self { kind, data }
    }

    /// Compute the dot product `self · rhs`, consuming both.
    pub fn dot(self, rhs: Self) -> TensorResult<Self> {
        let data = self.data.contract(rhs.data)?;
        let kind =
            if let Some(scalar) = data.as_scalar() {
                Kind::Scalar(scalar)
            } else {
                Kind::Unknown
            };
        Ok(Self { kind, data })
    }

    /// Thin wrapper around [`dot`][Self::dot] but with the operands reversed
    /// (i.e. computing `rhs · self`) in order to conform to a left-to-right
    /// style of composition.
    pub fn then(self, rhs: Self) -> TensorResult<Self> {
        let data = rhs.data.contract(self.data)?;
        let kind =
            if let Some(scalar) = data.as_scalar() {
                Kind::Scalar(scalar)
            } else {
                Kind::Unknown
            };
        Ok(Self { kind, data })
    }

    /// Multiply by a scalar, modifying `self` in place.
    pub fn scalar_mul_mut(&mut self, scalar: C64) {
        const EPSILON: f64 = 1e-12;
        if (scalar + 1.0).norm() < EPSILON { return; }
        self.data.scalar_mul_inplace(scalar);
        self.kind = Kind::Unknown;
    }

    /// Multiply by a scalar, consuming `self`.
    pub fn scalar_mul(mut self, scalar: C64) -> Self {
        self.data.scalar_mul_inplace(scalar);
        self
    }

    /// Conjugate `self` in place, converting kets to bras and bras to kets, and
    /// replacing each element of the inner tensor with its complex conjugate.
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
        self.data.conj_inplace();
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
    /// Specifically, `self` and `other` must have identical indices and the
    /// modulus of the difference between any two corresponding tensor elements
    /// must be less than `thresh`. Mutable references are required to sort the
    /// indices of both before the comparison is made.
    pub fn approx_eq(&mut self, other: &mut Self, thresh: Option<f64>) -> bool {
        self.data.sort_indices();
        other.data.sort_indices();
        self.data.approx_eq(&other.data, thresh)
    }

    /// Convert to a [`ketbra::Element`][crate::ketbra2::Element`].
    pub fn as_ketbras(&self) -> crate::ketbra2::Element {
        use crate::ketbra2 as kb;
        match self.kind {
            Kind::Scalar(a) =>
                kb::Element::from_terms([kb::KetBra::new(a, [], [])]),
            Kind::Id(i) => kb::Element::z([i], [i], None),
            Kind::Z(ph) => kb::Element::z(self.ins(), self.outs(), Some(ph)),
            Kind::X(ph) => kb::Element::x(self.ins(), self.outs(), Some(ph)),
            Kind::H(a) => kb::Element::h(self.ins(), self.outs(), Some(a)),
            Kind::Swap(a, b) => kb::Element::swap(a, b),
            Kind::Cup(a, b) => kb::Element::cup(a, b),
            Kind::Cap(a, b) => kb::Element::cap(a, b),
            Kind::Unknown => {
                fn ixstate(b: usize) -> kb::State {
                    if b == 0 { kb::State::Zero } else { kb::State::One }
                }

                if let Some(a) = self.data.as_scalar() {
                    kb::Element::from_terms([kb::KetBra::new(a, [], [])])
                } else {
                    let arr = self.data.as_array().unwrap();
                    let idxs = self.data.indices().unwrap();
                    let terms =
                        arr.indexed_iter()
                        .map(|(ix, a)| {
                            let ix_arr = ix.as_array_view();
                            let ix_slice = ix_arr.as_slice().unwrap();
                            let kets =
                                idxs.iter().zip(ix_slice)
                                .filter(|(idx, _)| idx.is_ket())
                                .map(|(idx, ix)| (idx.wire_index(), ixstate(*ix)));
                            let bras =
                                idxs.iter().zip(ix_slice)
                                .filter(|(idx, _)| idx.is_bra())
                                .map(|(idx, ix)| (idx.wire_index(), ixstate(*ix)));
                            kb::KetBra::new(*a, kets, bras)
                        });
                    kb::Element::from_terms(terms)
                }
            },
        }
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.data.fmt(f)
    }
}

