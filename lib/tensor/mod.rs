//! Tensor-based tools to compute end producs of diagrams in the ZX(H)-calculus.
//!
//! # TODO

use thiserror::Error;

#[derive(Debug, Error)]
pub enum TensorError {
    #[error("duplicate index {0:?}")]
    DuplicateIndex(Q),

    #[error("non-matching indices {0:?} and shape {1:?}")]
    IncompatibleShape(Box<[Q]>, Box<[usize]>),

    #[error("un-matched duplicate index in contraction {0:?}")]
    ContractDuplicateIndex(Q),

    #[error("encountered a non-generator element in graph conversion")]
    GraphConvNonGenerator,

    #[error("I/O error: {0}")]
    IOError(#[from] std::io::Error),
}
pub type TensorResult<T> = Result<T, TensorError>;

#[allow(clippy::module_inception)]
pub(crate) mod tensor;
pub use tensor::*;

pub(crate) mod element;
pub use element::*;

pub(crate) mod diagram;
#[allow(unused_imports)]
pub use diagram::*;

pub mod circuit;

