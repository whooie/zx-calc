//! Dense tensor support for [`Element`][super::Element].

use ndarray::{ self as nd, Dimension };
use num_complex::Complex64 as C64;
use thiserror::Error;
use crate::phase::Phase;
use super::{ ElementDataSeal, ElementData };

#[derive(Debug, Error)]
pub enum DeError {
    #[error("duplicate index {0:?}")]
    DuplicateIndex(Q),

    #[error("non-matching indices {0:?} and shape {1:?}")]
    IncompatibleShape(Box<[Q]>, Box<[usize]>),

    #[error("un-matched duplicate index in contraction {0:?}")]
    ContractDuplicateIndex(Q),
}
pub type DeResult<T> = Result<T, DeError>;

pub(crate) mod tensor;
pub use tensor::*;

/// Dense storage type for use with [`Element`][super::Element].
#[derive(Clone, Debug)]
pub struct De(pub(crate) Tensor);

impl From<Tensor> for De {
    fn from(data: Tensor) -> Self { Self(data) }
}

impl From<C64> for De {
    fn from(a: C64) -> Self { Self(Tensor::new_scalar(a)) }
}

impl ElementDataSeal for De { }
impl ElementData for De {
    type Error = DeError;
    type InputIter<'a> = BraIndices<'a>;
    type OutputIter<'a> = KetIndices<'a>;

    fn input_iter(&self) -> BraIndices<'_> { self.0.bra_indices() }

    fn output_iter(&self) -> KetIndices<'_> { self.0.ket_indices() }

    fn has_input(&self, k: usize) -> bool { self.0.has_index(&Q::Bra(k)) }

    fn has_output(&self, k: usize) -> bool { self.0.has_index(&Q::Ket(k)) }

    fn scalar<T>(a: T) -> Self
    where T: Into<C64>
    {
        Self(Tensor::new_scalar(a.into()))
    }

    fn as_scalar(&self) -> Option<C64> { self.0.as_scalar() }

    fn is_scalar(&self) -> bool { self.0.is_scalar() }

    fn id(i: usize) -> Self {
        let indices: Vec<Q> = vec![Q::Ket(i), Q::Bra(i)];
        let array: nd::Array2<C64> = nd::Array2::eye(2);
        let data = Tensor::from_array_unchecked(indices, array);
        Self(data)
    }

    fn z<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let mut indices: Vec<Q> = Vec::new();
        let phase = phase.unwrap_or_else(Phase::zero);
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
            Self(data)
        } else {
            let shape: Vec<usize> = vec![2; indices.len()];
            let mut array: nd::ArrayD<C64> = nd::ArrayD::zeros(shape);
            let elem = array.first_mut().unwrap();
            *elem = 1.0.into();
            let elem = array.last_mut().unwrap();
            *elem = phase.cis();
            let data = Tensor::from_array_unchecked(indices, array);
            Self(data)
        }
    }

    fn x<I, J>(ins: I, outs: J, phase: Option<Phase>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        let mut indices: Vec<Q> = Vec::new();
        let phase = phase.unwrap_or_else(Phase::zero);
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
            Self(data)
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
            Self(data)
        }
    }

    fn h<I, J>(ins: I, outs: J, a: Option<C64>) -> Self
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        const EPSILON: f64 = 1e-12;
        let mut indices: Vec<Q> = Vec::new();
        let arg = a.unwrap_or_else(|| -C64::from(1.0));
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
            Self(data)
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
            Self(data)
        }
    }

    fn swap(i0: usize, i1: usize) -> Self {
        if i0 == i1 {
            Self::id(i0)
        } else {
            let indices: Vec<Q> =
                vec![Q::Ket(i0), Q::Ket(i1), Q::Bra(i0), Q::Bra(i1)];
            let mut array: nd::Array4<C64> = nd::Array4::zeros((2, 2, 2, 2));
            array[[0, 0, 0, 0]] = 1.0.into();
            array[[0, 1, 1, 0]] = 1.0.into();
            array[[1, 0, 0, 1]] = 1.0.into();
            array[[1, 1, 1, 1]] = 1.0.into();
            let data = Tensor::from_array_unchecked(indices, array);
            Self(data)
        }
    }

    fn cup(i0: usize, i1: usize) -> Self {
        if i0 == i1 { return Self::z([], [i0], None); }
        let indices: Vec<Q> = vec![Q::Ket(i0), Q::Ket(i1)];
        let mut array: nd::Array2<C64> = nd::Array2::zeros((2, 2));
        array[[0, 0]] = 1.0.into();
        array[[1, 1]] = 1.0.into();
        let data = Tensor::from_array_unchecked(indices, array);
        Self(data)
    }

    fn cap(i0: usize, i1: usize) -> Self {
        if i0 == i1 { return Self::z([i0], [], None); }
        let indices: Vec<Q> = vec![Q::Bra(i0), Q::Bra(i1)];
        let mut array: nd::Array2<C64> = nd::Array2::zeros((2, 2));
        array[[0, 0]] = 1.0.into();
        array[[1, 1]] = 1.0.into();
        let data = Tensor::from_array_unchecked(indices, array);
        Self(data)
    }

    fn dot(self, rhs: Self) -> Result<Self, Self::Error> {
        let data = self.0.contract(rhs.0)?;
        Ok(Self(data))
    }

    fn scalar_mul<C>(&mut self, scalar: C)
    where C: Into<C64>
    {
        const EPSILON: f64 = 1e-12;
        let scalar = scalar.into();
        if (scalar - 1.0).norm() < EPSILON { return; }
        self.0.scalar_mul_inplace(scalar);
    }

    fn conjugate(&mut self) { self.0.conj_inplace(); }
}

impl De {
    /// Create a new `De` using a function over index values.
    ///
    /// This function is meant to generate tensor elements from an iterator over
    /// all possible index values of all (unique) indices passed to this
    /// function. The positions of the index values passed to the generating
    /// function correspond to the *unique* indices passed to this function *in
    /// the order in which they're passed*.
    ///
    /// Each index corresponds to a qubit wire; hence indices all range over
    /// values {0, 1}.
    pub fn new_elems<I, F>(indices: I, elems: F) -> Self
    where
        I: IntoIterator<Item = Q>,
        F: FnMut(&[usize]) -> C64,
    {
        Self(Tensor::new(indices, elems))
    }

    /// Sort axes in place so that wire indices are placed in ascending order,
    /// with all kets to the left of all bras.
    pub fn sort_indices(&mut self) { self.0.sort_indices() }

    /// Return `true` if `self` and `other` denote the same tensor to within an
    /// optional threshold, which defaults to `1e-12`.
    ///
    /// Specifically, `self` and `other` must have identical indices and the
    /// modulus of the difference between any two corresponding tensor elements
    /// must be less than `thresh`. Mutable references are required to sort the
    /// indices of both before the comparison is made.
    pub fn approx_eq(&mut self, other: &mut Self, thresh: Option<f64>)
        -> bool
    {
        self.0.sort_indices();
        other.0.sort_indices();
        self.0.approx_eq(&other.0, thresh)
    }

    /// Convert to a corresponding sparse representation.
    pub fn to_sparse(&self) -> super::sparse::Sp {
        use super::sparse as sp;

        fn ixstate(b: usize) -> sp::State {
            if b == 0 { sp::State::Zero } else { sp::State::One }
        }

        if let Some(a) = self.0.as_scalar() {
            sp::Sp::scalar(a)
        } else {
            let arr = self.0.as_array().unwrap();
            let idxs = self.0.indices().unwrap();
            let terms: Vec<sp::KetBra> =
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
                    sp::KetBra::new(*a, kets, bras)
                })
                .collect();
            sp::Sp(terms)
        }
    }
}

impl std::fmt::Display for De {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

