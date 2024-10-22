//! Tools to create and compute end product of diagrams in the ZX(H)-calculus
//! based on naive ketbra multiplication.
//!
//! You will most likely only want to deal directly with [`Element`] and
//! [`Diagram`].

use thiserror::Error;

#[derive(Debug, Error)]
pub enum KBError {
    #[error("duplicate ket key in dot result: {0}")]
    DotDuplicateKetKey(usize),

    #[error("duplicate bra key in dot result: {0}")]
    DotDuplicateBraKey(usize),

    #[error("diagram: encountered a non-generator element")]
    DiagramConversionNonGenerator,

    #[error("error constructing GraphViz representation: {0}")]
    GraphVizError(String),

    #[error("I/O error: {0}")]
    IOError(#[from] std::io::Error),
}
pub type KBResult<T> = Result<T, KBError>;

pub(crate) mod state;
pub use state::*;

#[allow(clippy::module_inception)]
pub(crate) mod ketbra;
pub use ketbra::*;

pub(crate) mod element;
pub use element::*;

pub(crate) mod diagram;
pub use diagram::*;

pub mod circuit;

