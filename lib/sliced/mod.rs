use thiserror::Error;

#[derive(Debug, Error)]
pub enum SlicedError {
    #[error("encountered non-matched duplicate wire index {0}")]
    DuplicateWireIndex(usize),

    #[error("encountered a non-generator element in graph conversion")]
    GraphConvNonGenerator,

    #[error("tensor error: {0}")]
    TensorError(#[from] dense::TensorError),

    #[error("ket-bra error: {0}")]
    KBError(#[from] sparse::KBError),

    #[error("I/O error: {0}")]
    IOError(#[from] std::io::Error),
}
pub type SlicedResult<T> = Result<T, SlicedError>;

pub mod dense;
pub use dense::{ De, Q, Tensor };

pub mod sparse;
pub use sparse::{ Sp, Basis, State, KetBra };

pub(crate) mod element;
pub use element::*;

pub(crate) mod diagram;
pub use diagram::*;

pub mod circuit;
