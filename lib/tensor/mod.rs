//! Tensor-based tools to create and compute end products of diagrams in the
//! ZX(H)-calculus.
//!
//! Diagrams are represented as "sliced" series of diagram elements. Each
//! element has a tensor representation, and the overall linear map denoted by
//! a diagram is simply the contraction of all the tensors it comprises, with
//! diagram inputs and outputs inferred from leftover tensor indices.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum TensorError {
    #[error("encountered non-matched duplicate wire index {0}")]
    DuplicateWireIndex(usize),

    #[error("encountered a non-generator element in graph conversion")]
    GraphConvNonGenerator,

    #[error("dense error: {0}")]
    DeError(#[from] dense::DeError),

    #[error("sparse error: {0}")]
    SpError(#[from] sparse::SpError),

    #[error("I/O error: {0}")]
    IOError(#[from] std::io::Error),
}
pub type TensorResult<T> = Result<T, TensorError>;

pub mod dense;
pub use dense::{ De, Q, Tensor };

pub mod sparse;
pub use sparse::{ Sp, Basis, State, KetBra };

pub(crate) mod element;
pub use element::*;

pub(crate) mod diagram;
pub use diagram::*;

pub mod circuit;
