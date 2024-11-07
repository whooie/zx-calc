use std::fmt;
use ndarray::{ self as nd, Dimension };
use num_complex::Complex64 as C64;
use super::{ TensorError, TensorResult };

use TensorError::*;

/// An index for a qubit wire on a specific side of a ket-bra.
///
/// Each index of this type is identified with a specific "half" of a ket-bra,
/// meaning that it names a particular qubit wire either on the output (ket) or
/// input (bra) side of the linear transformation as a whole.
///
/// [`Ord`] is implemented to sort all ket indices before all bra indices,
/// deferring to the index values to determined orderings within those
/// subgroups.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum Q {
    /// A qubit index on the ket side of a ket-bra.
    Ket(usize),

    /// A qubit index on the bra side of a ket-bra.
    Bra(usize),
}

impl Q {
    /// Return `true` if `self` is `Ket`.
    pub fn is_ket(&self) -> bool { matches!(self, Self::Ket(_)) }

    /// Return `true` if `self` is `Ket` and the wire index satisfies a
    /// predicate.
    pub fn is_ket_and<F>(&self, pred: F) -> bool
    where F: FnOnce(usize) -> bool
    {
        if let Self::Ket(k) = self {
            pred(*k)
        } else {
            false
        }
    }

    /// Return `true` if `self` is `Bra`.
    pub fn is_bra(&self) -> bool { matches!(self, Self::Bra(_)) }

    /// Return `true` if `self` is `Bra` and the wire index satisfies a
    /// predicate.
    pub fn is_bra_and<F>(&self, pred: F) -> bool
    where F: FnOnce(usize) -> bool
    {
        if let Self::Bra(k) = self {
            pred(*k)
        } else {
            false
        }
    }

    /// Conjugate `self`, changing `Ket`s to `Bra`s and `Bra`s to `Ket`s.
    pub fn conj(self) -> Self {
        match self {
            Self::Ket(k) => Self::Bra(k),
            Self::Bra(k) => Self::Ket(k),
        }
    }

    /// Conjugate `self` in place, changing `Ket`s to `Bra`s and `Bra`s to
    /// `Ket`s.
    pub fn conj_inplace(&mut self) {
        match self {
            Self::Ket(k) => { *self = Self::Bra(*k); },
            Self::Bra(k) => { *self = Self::Ket(*k); },
        }
    }

    /// Return `true` if `self` is a `Bra` and `other` is a `Ket` with the same
    /// wire index.
    pub fn matches_with(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (Self::Bra(l), Self::Ket(r)) if l == r,
        )
    }

    /// Return `self` as a `Ket`, maintaining wire index.
    pub fn as_ket(self) -> Self {
        match self {
            Self::Ket(k) => Self::Ket(k),
            Self::Bra(k) => Self::Ket(k),
        }
    }

    /// Return `self` as a `Bra`, maintaining wire index.
    pub fn as_bra(self) -> Self {
        match self {
            Self::Ket(k) => Self::Bra(k),
            Self::Bra(k) => Self::Bra(k),
        }
    }

    /// Return the wire index.
    pub fn wire_index(&self) -> usize {
        match self {
            Self::Ket(k) => *k,
            Self::Bra(k) => *k,
        }
    }

    /// Return `true` if `self` and `other` have the same wire index, regardless
    /// of ket-bra-ness.
    pub fn same_wire(&self, other: &Self) -> bool {
        self.wire_index() == other.wire_index()
    }
}

impl PartialOrd for Q {
    fn partial_cmp(&self, rhs: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(rhs))
    }
}

impl Ord for Q {
    fn cmp(&self, rhs: &Self) -> std::cmp::Ordering {
        match (self, rhs) {
            (Self::Ket(_), Self::Bra(_)) => std::cmp::Ordering::Less,
            (Self::Bra(_), Self::Ket(_)) => std::cmp::Ordering::Greater,
            (Self::Ket(l), Self::Ket(r)) => l.cmp(r),
            (Self::Bra(l), Self::Bra(r)) => l.cmp(r),
        }
    }
}

impl fmt::Display for Q {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

/// A dynamically dimensioned array, with data stored behind an atomic reference
/// counter.
pub type ArcArrayD<T> = nd::ArcArray<T, nd::IxDyn>;

#[derive(Clone, PartialEq, Debug)]
enum TensorData {
    Scalar(C64),
    Tensor(Vec<Q>, ArcArrayD<C64>),
}

impl From<C64> for TensorData {
    fn from(val: C64) -> Self { Self::Scalar(val) }
}

impl TensorData {
    fn new<I, F>(indices: I, mut elems: F) -> Self
    where
        I: IntoIterator<Item = Q>,
        F: FnMut(&[usize]) -> C64,
    {
        let mut idxs: Vec<Q> = Vec::new();
        indices.into_iter()
            .for_each(|idx| { if !idxs.contains(&idx) { idxs.push(idx); } });
        if idxs.is_empty() {
            let scalar: C64 = elems(&[]);
            Self::Scalar(scalar)
        } else {
            let shape: Vec<usize> = vec![2; idxs.len()];
            let data: nd::ArrayD<C64> =
                nd::ArrayD::from_shape_fn(
                    shape,
                    |ix| elems(ix.as_array_view().as_slice().unwrap())
                );
            Self::Tensor(idxs, data.into_shared())
        }
    }

    pub(crate) fn new_unchecked<F>(indices: Vec<Q>, mut elems: F) -> Self
    where F: FnMut(&[usize]) -> C64
    {
        if indices.is_empty() {
            Self::Scalar(elems(&[]))
        } else {
            let shape: Vec<usize> = vec![2; indices.len()];
            let data: nd::ArrayD<C64> =
                nd::ArrayD::from_shape_fn(
                    shape,
                    |idxs| elems(idxs.as_array_view().as_slice().unwrap())
                );
            Self::Tensor(indices, data.into_shared())
        }
    }

    fn new_scalar(val: C64) -> Self { Self::Scalar(val) }

    fn from_array<I, S, D>(indices: I, array: nd::ArrayBase<S, D>)
        -> TensorResult<Self>
    where
        I: IntoIterator<Item = Q>,
        S: nd::DataOwned<Elem = C64>,
        D: nd::Dimension,
    {
        let indices: Vec<Q> = indices.into_iter().collect();
        let mb_dup =
            indices.iter().enumerate()
            .find_map(|(k, idx)| {
                indices.iter().skip(k + 1).find(|idx2| *idx2 == idx)
            });
        if let Some(dup) = mb_dup {
            return Err(DuplicateIndex(*dup));
        }
        let shape = array.shape();
        if indices.len() == shape.len() && shape.iter().all(|dim| *dim == 2) {
            if indices.is_empty() {
                Ok(Self::Scalar(*array.into_iter().next().unwrap()))
            } else {
                Ok(Self::Tensor(indices, array.into_dyn().into_shared()))
            }
        } else if shape.len() == 1 && shape[0] == 1_usize << indices.len() {
            let tensor_shape: Vec<usize> = vec![2; indices.len()];
            let data = array.into_dyn().into_shared();
            data.reshape(tensor_shape);
            Ok(Self::Tensor(indices, data))
        } else {
            let idx_strs: Box<[Q]> = indices.into();
            let array_shape: Box<[usize]> = shape.iter().copied().collect();
            Err(IncompatibleShape(idx_strs, array_shape))
        }
    }

    fn from_array_unchecked<S, D>(
        indices: Vec<Q>,
        array: nd::ArrayBase<S, D>,
    ) -> Self
    where
        S: nd::DataOwned<Elem = C64>,
        D: nd::Dimension,
    {
        if indices.is_empty() {
            Self::Scalar(*array.into_iter().next().unwrap())
        } else {
            Self::Tensor(indices, array.into_dyn().into_shared())
        }
    }

    fn is_scalar(&self) -> bool { matches!(self, Self::Scalar(_)) }

    fn as_scalar(&self) -> Option<C64> {
        match self {
            Self::Scalar(a) => Some(*a),
            _ => None,
        }
    }

    fn is_tensor(&self) -> bool { matches!(self, Self::Tensor(..)) }

    fn as_array(&self) -> Option<&ArcArrayD<C64>> {
        match self {
            Self::Tensor(_, arr) => Some(arr),
            _ => None
        }
    }

    fn has_index(&self, index: &Q) -> bool {
        match self {
            Self::Scalar(_) => false,
            Self::Tensor(idxs, _) => idxs.contains(index),
        }
    }

    fn index_pos(&self, index: &Q) -> Option<usize> {
        match self {
            Self::Scalar(_) => None,
            Self::Tensor(idxs, _) => {
                idxs.iter().enumerate()
                    .find_map(|(k, idx)| (idx == index).then_some(k))
            },
        }
    }

    fn get_index(&self, k: usize) -> Option<&Q> {
        match self {
            Self::Scalar(_) => None,
            Self::Tensor(idxs, _) => idxs.get(k),
        }
    }

    unsafe fn get_index_mut(&mut self, k: usize) -> Option<&mut Q> {
        match self {
            Self::Scalar(_) => None,
            Self::Tensor(idxs, _) => idxs.get_mut(k),
        }
    }

    fn indices(&self) -> Option<&Vec<Q>> {
        match self {
            Self::Scalar(_) => None,
            Self::Tensor(idxs, _) => Some(idxs),
        }
    }

    unsafe fn indices_mut(&mut self) -> Option<&mut Vec<Q>> {
        match self {
            Self::Scalar(_) => None,
            Self::Tensor(idxs, _) => Some(idxs),
        }
    }

    fn ket_indices(&self) -> KetIndices<'_> {
        match self {
            Self::Scalar(_) => KetIndicesData::Scalar.into(),
            Self::Tensor(idxs, _) => KetIndicesData::Tensor(idxs.iter()).into(),
        }
    }

    fn bra_indices(&self) -> BraIndices<'_> {
        match self {
            Self::Scalar(_) => BraIndicesData::Scalar.into(),
            Self::Tensor(idxs, _) => BraIndicesData::Tensor(idxs.iter()).into(),
        }
    }

    fn rank(&self) -> usize {
        match self {
            Self::Scalar(_) => 0,
            Self::Tensor(idxs, _) => idxs.len(),
        }
    }

    unsafe fn map_indices<F>(self, map: F) -> TensorData
    where F: FnMut(Q) -> Q
    {
        match self {
            Self::Scalar(a) => TensorData::Scalar(a),
            Self::Tensor(idxs, data) => {
                let idxs = idxs.into_iter().map(map).collect();
                TensorData::Tensor(idxs, data)
            },
        }
    }

    fn map<F>(&self, mut f: F) -> Self
    where F: FnMut(&[usize], C64) -> C64,
    {
        match self {
            Self::Scalar(a) => Self::Scalar(f(&[], *a)),
            Self::Tensor(idxs, data) => {
                let new_data =
                    data.indexed_iter()
                    .map(|(ix, a)| f(ix.as_array_view().as_slice().unwrap(), *a))
                    .collect::<nd::Array1<C64>>()
                    .into_dyn()
                    .into_shape(data.raw_dim())
                    .unwrap();
                Self::Tensor(idxs.clone(), new_data.into_shared())
            },
        }
    }

    fn map_inplace<F>(&mut self, mut f: F)
    where F: FnMut(&[usize], C64) -> C64
    {
        match self {
            Self::Scalar(a) => { *a = f(&[], *a); },
            Self::Tensor(_, data) => {
                data.indexed_iter_mut()
                    .for_each(|(ix, a)| {
                        *a = f(ix.as_array_view().as_slice().unwrap(), *a);
                    });
            },
        }
    }

    fn conj(&self) -> Self {
        match self {
            Self::Scalar(a) => Self::Scalar(a.conj()),
            Self::Tensor(idxs, a) => {
                let idxs_new: Vec<Q> =
                    idxs.iter()
                    .map(|idx| idx.conj())
                    .collect();
                let data_new: ArcArrayD<C64> = a.mapv(|ak| ak.conj()).into();
                Self::Tensor(idxs_new, data_new)
            },
        }
    }

    fn conj_inplace(&mut self) {
        match self {
            Self::Scalar(a) => { *a = a.conj(); },
            Self::Tensor(idxs, a) => {
                idxs.iter_mut().for_each(|idx| { idx.conj_inplace(); });
                a.map_inplace(|ak| { *ak = ak.conj(); });
            },
        }
    }

    fn swap_indices(&mut self, a: &Q, b: &Q) {
        if let Self::Tensor(idxs, data) = self {
            let pos_a =
                idxs.iter().enumerate()
                .find_map(|(k, idx)| (idx == a).then_some(k));
            let pos_b =
                idxs.iter().enumerate()
                .find_map(|(k, idx)| (idx == b).then_some(k));
            if let Some((k_a, k_b)) = pos_a.zip(pos_b) {
                idxs.swap(k_a, k_b);
                data.swap_axes(k_a, k_b);
            }
        }
    }

    fn swap_indices_pos(&mut self, a: usize, b: usize) {
        if let Self::Tensor(idxs, data) = self {
            if a >= idxs.len() || b >= idxs.len() { return; }
            idxs.swap(a, b);
            data.swap_axes(a, b);
        }
    }

    fn sort_indices(&mut self) {
        // have to use bubble sort because we can only swap axes on the arrays
        if let Self::Tensor(idxs, data) = self {
            let mut n = idxs.len();
            let mut swapped = true;
            while swapped {
                swapped = false;
                for i in 1..=n - 1 {
                    if idxs[i - 1] > idxs[i] {
                        idxs.swap(i - 1, i);
                        data.swap_axes(i - 1, i);
                        swapped = true;
                    }
                }
                n -= 1;
            }
        }
    }

    fn sorted_indices(mut self) -> Self {
        self.sort_indices();
        self
    }

    fn do_contract(
        idx_common: Vec<usize>,
        mut idxs_a: Vec<Q>,
        mut a: ArcArrayD<C64>,
        mut idxs_b: Vec<Q>,
        mut b: ArcArrayD<C64>,
    ) -> Self
    {
        // swap common indices and corresponding axes to the rightmost positions
        // in a and the leftmost in b
        let n_idx_a = idxs_a.len();
        let n_common = idx_common.len();
        let n_idx_b = idxs_b.len();
        let mut k_src: usize;
        for (k_targ, idx) in idx_common.iter().enumerate() {
            k_src =
                idxs_a.iter().enumerate()
                .find_map(|(k_src, idx_src)| {
                    idx_src.is_bra_and(|ix| ix == *idx).then_some(k_src)
                })
                .unwrap();
            idxs_a.swap(k_src, n_idx_a - n_common + k_targ);
            a.swap_axes(k_src, n_idx_a - n_common + k_targ);

            k_src =
                idxs_b.iter().enumerate()
                .find_map(|(k_src, idx_src)| {
                    idx_src.is_ket_and(|ix| ix == *idx).then_some(k_src)
                })
                .unwrap();
            idxs_b.swap(k_src, k_targ);
            b.swap_axes(k_src, k_targ);
        }

        // reshape to fuse all common and non-common indices together so we can
        // use matmul and then reshape the result to un-fuse
        let dim_noncomm_a = 1_usize << (n_idx_a - n_common);
        let dim_comm = 1_usize << n_common;
        let dim_noncomm_b = 1_usize << (n_idx_b - n_common);
        let a: nd::CowArray<C64, nd::Ix2> =
            a.as_standard_layout()
            .into_shape((dim_noncomm_a, dim_comm))
            .unwrap();
        let b: nd::CowArray<C64, nd::Ix2> =
            b.as_standard_layout()
            .into_shape((dim_comm, dim_noncomm_b))
            .unwrap();
        let c: nd::Array2<C64> = a.dot(&b);
        let new_shape: Vec<usize> = vec![2; n_idx_a + n_idx_b - 2 * n_common];
        if new_shape.is_empty() {
            let c_val = c.into_iter().next().unwrap();
            Self::Scalar(c_val)
        } else {
            let new_idxs: Vec<Q> =
                idxs_a.into_iter().take(n_idx_a - n_common)
                .chain(idxs_b.into_iter().skip(n_common))
                .collect();
            let c = c.into_shape(new_shape).unwrap();
            Self::Tensor(new_idxs, c.into_shared())
        }
    }

    fn contract(self, rhs: Self) -> TensorResult<Self> {
        match (self, rhs) {
            (Self::Scalar(a), Self::Scalar(b)) => Ok(Self::Scalar(a * b)),
            (Self::Scalar(a), Self::Tensor(idxs, mut b)) => {
                b.mapv_inplace(|bk| a * bk);
                Ok(Self::Tensor(idxs, b))
            },
            (Self::Tensor(idxs, mut a), Self::Scalar(b)) => {
                a.mapv_inplace(|ak| ak * b);
                Ok(Self::Tensor(idxs, a))
            },
            (Self::Tensor(idxs_a, a), Self::Tensor(idxs_b, b)) => {
                let mut idx_common: Vec<usize> =
                    Vec::with_capacity(idxs_a.len().max(idxs_b.len()));
                // check for duplicate bras
                for idx_a in idxs_a.iter() {
                    if idx_a.is_ket() { continue; }
                    let has_match =
                        idxs_b.iter().any(|idx_b| idx_a.matches_with(idx_b));
                    let has_dup = idxs_b.contains(idx_a);
                    if has_match {
                        idx_common.push(idx_a.wire_index());
                    } else if has_dup {
                        return Err(ContractDuplicateIndex(*idx_a));
                    }
                }
                // check for duplicate kets
                for idx_b in idxs_b.iter() {
                    if idx_b.is_bra() { continue; }
                    let has_match = idx_common.contains(&idx_b.wire_index());
                    let has_dup = idxs_a.contains(idx_b);
                    if !has_match && has_dup {
                        return Err(ContractDuplicateIndex(*idx_b));
                    }
                }
                Ok(Self::do_contract(idx_common, idxs_a, a, idxs_b, b))
            },
        }
    }

    fn scalar_mul_inplace(&mut self, scalar: C64) {
        match self {
            Self::Scalar(a) => { *a *= scalar; },
            Self::Tensor(_, data) => { *data *= scalar; },
        }
    }

    fn scalar_mul(mut self, scalar: C64) -> Self {
        self.scalar_mul_inplace(scalar);
        self
    }

    fn into_flat(self) -> (Vec<Q>, nd::Array1<C64>) {
        match self {
            Self::Scalar(a) => (Vec::new(), nd::array![a]),
            Self::Tensor(idxs, data) => (idxs, data.into_iter().collect()),
        }
    }

    fn into_raw(self) -> (Vec<Q>, nd::ArrayD<C64>) {
        match self {
            Self::Scalar(a) => (Vec::new(), nd::array![a].into_dyn()),
            Self::Tensor(idxs, data) => (idxs, data.into_owned()),
        }
    }

    fn approx_eq(&self, other: &Self, thresh: Option<f64>) -> bool {
        let eps = thresh.unwrap_or(1e-12);
        match (self, other) {
            (Self::Scalar(l), Self::Scalar(r)) => (*l - *r).norm() < eps,
            (Self::Tensor(idxs_l, data_l), Self::Tensor(idxs_r, data_r)) => {
                idxs_l == idxs_r
                    && data_l.iter().zip(data_r)
                        .all(|(l, r)| (*l - *r).norm() < eps)
            },
            _ => false
        }
    }
}

impl fmt::Display for TensorData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Scalar(a) => {
                a.fmt(f)?;
                write!(f, " {{ }}")?;
            },
            Self::Tensor(idxs, a) => {
                a.fmt(f)?;
                write!(f, "\n{{ ")?;
                let n_idxs = idxs.len();
                for (k, idx) in idxs.iter().enumerate() {
                    idx.fmt(f)?;
                    if k < n_idxs - 1 { write!(f, ", ")?; }
                }
                write!(f, " }}")?;
            },
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
enum KetIndicesData<'a> {
    Scalar,
    Tensor(std::slice::Iter<'a, Q>),
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
            KetIndicesData::Tensor(ref mut iter) => {
                iter.find_map(|idx| idx.is_ket().then(|| idx.wire_index()))
            },
        }
    }
}

impl<'a> DoubleEndedIterator for KetIndices<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        match &mut self.0 {
            KetIndicesData::Scalar => None,
            KetIndicesData::Tensor(ref mut iter) => {
                iter.rfind(|idx| idx.is_ket()).map(|idx| idx.wire_index())
            },
        }
    }
}

impl<'a> std::iter::FusedIterator for KetIndices<'a> { }

#[derive(Clone, Debug)]
enum BraIndicesData<'a> {
    Scalar,
    Tensor(std::slice::Iter<'a, Q>),
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
            BraIndicesData::Tensor(ref mut iter) => {
                iter.find_map(|idx| idx.is_bra().then(|| idx.wire_index()))
            },
        }
    }
}

impl<'a> DoubleEndedIterator for BraIndices<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        match &mut self.0 {
            BraIndicesData::Scalar => None,
            BraIndicesData::Tensor(ref mut iter) => {
                iter.rfind(|idx| idx.is_bra()).map(|idx| idx.wire_index())
            },
        }
    }
}

impl<'a> std::iter::FusedIterator for BraIndices<'a> { }

/// Basic implementation of an abstract tensor object.
///
/// A `Tensor` consists of some number of complex numbers and a series of [`Q`]
/// wire indices.
///
/// This implementation distinguishes between rank-0 (scalar) and rank > 0
/// (array) quantities for some small operational benefits, but otherwise lived
/// in an abstraction layer one step above concrete array representations.
/// Tensors are thread-safe and relatively cheap to clone. No parallelism or GPU
/// computation is offered.
#[derive(Clone, PartialEq, Debug)]
pub struct Tensor(TensorData);

impl From<TensorData> for Tensor {
    fn from(data: TensorData) -> Self { Self(data) }
}

impl From<C64> for Tensor {
    fn from(val: C64) -> Self { Self(val.into()) }
}

impl Tensor {
    /// Create a new tensor using a function over given indices.
    ///
    /// This function is meant to generate tensor elements from an iterator over
    /// all possible index values of all (unique) indices passed to this
    /// function. The positions of the index values passed to the generating
    /// function correspond to the *unique* indices passed to this function *in
    /// the order in which they're passed*.
    ///
    /// Each index corresponds to a qubit wire; hence indices all range over
    /// values {0, 1}.
    pub fn new<I, F>(indices: I, elems: F) -> Self
    where
        I: IntoIterator<Item = Q>,
        F: FnMut(&[usize]) -> C64,
    {
        TensorData::new(indices, elems).into()
    }

    pub(crate) fn new_unchecked<F>(indices: Vec<Q>, elems: F) -> Self
    where F: FnMut(&[usize]) -> C64
    {
        TensorData::new_unchecked(indices, elems).into()
    }

    /// Create a new rank-0 (scalar) tensor.
    pub fn new_scalar(val: C64) -> Self { TensorData::new_scalar(val).into() }

    /// Create a new tensor from an n-dimensional array.
    ///
    /// Fails if duplicate indices are provided or the dimensions of the array
    /// do not match those of the indices. The array's dimensions match the
    /// indices if the array has an axis of length 2 for every index or is
    /// one-dimensional with length 2<sup><i>n</i></sup>, where *n* is the
    /// number of indices.
    pub fn from_array<I, S, D>(indices: I, array: nd::ArrayBase<S, D>)
        -> TensorResult<Self>
    where
        I: IntoIterator<Item = Q>,
        S: nd::DataOwned<Elem = C64>,
        D: nd::Dimension,
    {
        TensorData::from_array(indices, array).map(Self::from)
    }

    pub(crate) fn from_array_unchecked<S, D>(
        indices: Vec<Q>,
        array: nd::ArrayBase<S, D>,
    ) -> Self
    where
        S: nd::DataOwned<Elem = C64>,
        D: nd::Dimension,
    {
        TensorData::from_array_unchecked(indices, array).into()
    }

    /// Return `true` if `self` has rank 0.
    pub fn is_scalar(&self) -> bool { self.0.is_scalar() }

    /// If `self` has rank 0, return its value as a single scalar.
    pub fn as_scalar(&self) -> Option<C64> { self.0.as_scalar() }

    /// Return `true` if `self` has rank > 0.
    pub fn is_tensor(&self) -> bool { self.0.is_tensor() }

    /// If `self` has rank > 0, return its inner data as an ordinary array.
    pub fn as_array(&self) -> Option<&ArcArrayD<C64>> { self.0.as_array() }

    /// Return `true` if `self` has the given index.
    pub fn has_index(&self, index: &Q) -> bool { self.0.has_index(index) }

    /// Return the position of an index if it exists.
    pub fn index_pos(&self, index: &Q) -> Option<usize> {
        self.0.index_pos(index)
    }

    /// Return a reference to the index at axis position `k`, if it exists.
    pub fn get_index(&self, k: usize) -> Option<&Q> { self.0.get_index(k) }

    pub(crate) unsafe fn get_index_mut(&mut self, k: usize) -> Option<&mut Q> {
        self.0.get_index_mut(k)
    }

    /// Return a reference to all indices, if `self` has rank > 0.
    pub fn indices(&self) -> Option<&Vec<Q>> { self.0.indices() }

    pub(crate) unsafe fn indices_mut(&mut self) -> Option<&mut Vec<Q>> {
        self.0.indices_mut()
    }

    /// Return an iterator over all ket wire indices.
    ///
    /// The iterator item type is `usize`.
    pub fn ket_indices(&self) -> KetIndices<'_> { self.0.ket_indices() }

    /// Return an iterator over all bra wire indices.
    ///
    /// The iterator item type is `usize`.
    pub fn bra_indices(&self) -> BraIndices<'_> { self.0.bra_indices() }

    /// Return the rank of `self`.
    pub fn rank(&self) -> usize { self.0.rank() }


    pub(crate) unsafe fn map_indices<F>(self, f: F) -> Tensor
    where F: FnMut(Q) -> Q
    {
        self.0.map_indices(f).into()
    }

    /// Apply a mapping function to the (indexed) elements of `self`, returning
    /// a new `Tensor` with the same indices.
    pub fn map<F>(&self, f: F) -> Self
    where F: FnMut(&[usize], C64) -> C64,
    {
        self.0.map(f).into()
    }

    /// Apply a mapping function to the (indexed) elements of `self` in place,
    /// reassigning to the output of the function.
    pub fn map_inplace<F>(&mut self, f: F)
    where F: FnMut(&[usize], C64) -> C64
    {
        self.0.map_inplace(f);
    }

    /// Return the element-wise conjugate of `self` as a new `Tensor`.
    pub fn conj(&self) -> Self { self.0.conj().into() }

    /// Conjugate the elements of `self` in place.
    pub fn conj_inplace(&mut self) { self.0.conj_inplace(); }

    /// Swap two indices in the underlying representation of `self`.
    ///
    /// Note that this operation has no bearing over the computational
    /// complexity of tensor contractions. If either given index does not exist,
    /// no swap is performed.
    pub fn swap_indices(&mut self, a: &Q, b: &Q) { self.0.swap_indices(a, b); }

    /// Swap two index positions in the underlying representation of `self`.
    ///
    /// Note that this operation has no bearing over the computational
    /// complexity of tensor contractions. If either axis position does not
    /// exist, no swap is performed.
    pub fn swap_indices_pos(&mut self, a: usize, b: usize) {
        self.0.swap_indices_pos(a, b);
    }

    /// Apply an ordering to the set of indices in the underlying representation
    /// of `self`.
    ///
    /// Note that this operation has no bearing over the computational
    /// complexity of tensor contractions.
    pub fn sort_indices(&mut self) { self.0.sort_indices() }

    /// Apply an ordering to the set of indices in the underlying representation
    /// of `self`.
    ///
    /// Note that this operation has no bearing over the computational
    /// complexity of tensor contractions.
    pub fn sorted_indices(self) -> Self { self.0.sorted_indices().into() }

    /// Contract `self` with `rhs` over all common indices, consuming both.
    ///
    /// this operation results in a new `Tensor` whose indices comprise all
    /// non-common indices belonging to `self` and `rhs`. All indices originally
    /// belonging to `self` are placed before those from `rhs`, but the ordering
    /// of indices within these groups is not preserved.
    pub fn contract(self, rhs: Self) -> TensorResult<Self> {
        self.0.contract(rhs.0).map(Self::from)
    }

    /// Multiply by a scalar, modifying `self` in place.
    pub fn scalar_mul_inplace(&mut self, scalar: C64) {
        self.0.scalar_mul_inplace(scalar);
    }

    /// Multiply by a scalar, consuming `self`.
    pub fn scalar_mul(self, scalar: C64) -> Self {
        self.0.scalar_mul(scalar).into()
    }

    /// Flatten `self` into a list of indices and an ordinary 1D array.
    ///
    /// If `self` is a scalar, the index list is empty.
    pub fn into_flat(self) -> (Vec<Q>, nd::Array1<C64>) { self.0.into_flat() }

    /// Unwrap `self` into a list of indices and an ordinary N-dimensional
    /// array.
    ///
    /// If `self` is a scalar, the list is empty.
    pub fn into_raw(self) -> (Vec<Q>, nd::ArrayD<C64>) { self.0.into_raw() }

    /// Return `true` if `self` and `other` denote the same tensor, up to an
    /// optional threshold.
    ///
    /// Specifically, `self` and `other` must have identical indices and the
    /// modulus of the difference between any two corresponding tensor elements
    /// must be less than `thresh`, which defaults to `1e-12`.
    pub fn approx_eq(&self, other: &Self, thresh: Option<f64>) -> bool {
        self.0.approx_eq(&other.0, thresh)
    }
}

impl fmt::Display for Tensor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { self.0.fmt(f) }
}

