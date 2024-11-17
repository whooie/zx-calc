//! Tensor-based tools to create and compute end products of diagrams in the
//! ZX(H)-calculus.
//!
//! Diagrams are represented as "sliced" series of diagram elements. Each
//! element has a tensor representation, and the overall linear map denoted by
//! a diagram is simply the contraction of all the tensors it comprises, with
//! diagram inputs and outputs inferred from leftover tensor indices.

use thiserror::Error;

/// Errors for fallible operations on tensor-backed diagrams (sparse and dense).
#[derive(Debug, Error)]
pub enum TensorError {
    /// Returned when converting a tensor-backed diagram to a graph-backed
    /// diagram and an unmatched duplicate wire index is encountered.
    #[error("encountered unmatched duplicate wire index {0}")]
    DuplicateWire(usize),

    /// Returned when converting a tensor-backed diagram to a graph-backed
    /// diagram and a diagram element is encountered for which no single
    /// generator equivalent can be chosen.
    #[error("encountered a non-generator element in graph conversion")]
    NonGenerator,

    /// Returned by operations on an underlying [dense][De] tensor
    /// representation.
    #[error("{0}")]
    DeError(#[from] dense::DeError),

    /// Returned by operations on an underlying [sparse][Sp] tensor
    /// representation.
    #[error("{0}")]
    SpError(#[from] sparse::SpError),

    /// I/O error when writing a diagram to a file.
    #[error("{0}")]
    IOError(#[from] std::io::Error),
}
/// Results with [`TensorError`] errors.
pub type TensorResult<T> = Result<T, TensorError>;

pub mod dense;
pub use dense::{ De, Q, Tensor };

pub mod sparse;
pub use sparse::{ Sp, Basis, State, KetBra };

pub(crate) mod element;
pub use element::*;

pub(crate) mod diagram;
pub use diagram::*;

