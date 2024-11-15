//! Graph-based tools to create and compute end products of diagrams in the
//! ZX(H)-calculus.
//!
//! Diagrams are represented as an undirected, unweighted graph with parallel
//! edges and data attached to its nodes. See
//! [`graph_diagram!`][crate::graph_diagram] for example usage and abbreviated
//! syntax.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum GraphError {
    #[error("error removing node: missing node {0}")]
    RemoveNodeMissingNode(NodeId),

    #[error("error adding wire: missing node {0}")]
    AddWireMissingNode(NodeId),

    #[error("error adding wire: input/output node {0} is already connected")]
    AddWireConnectedIO(NodeId),

    #[error("error removing wire(s): missing node {0}")]
    RemoveWireMissingNode(NodeId),

    #[error("error applying state/effect: node {0} is not an input/output")]
    ApplyStateNotIO(NodeId),

    #[error("error applying state/effect: missing qubit {0}")]
    ApplyStateMissingQubit(QubitId),

    #[error("error applying bell state/effect: nodes {0}, {1} are not both inputs/outputs")]
    ApplyBellNotIO(NodeId, NodeId),

    #[error("error applying bell state/effect: missing qubit {0}")]
    ApplyBellMissingQubit(QubitId),

    #[error("error applying bell state/effect: cannot apply to a single qubit")]
    ApplyBellSameQubit,

    #[error("error applying bell state/effect: input/output node {0} is disconnected")]
    ApplyBellDisconnected(NodeId),

    #[error("error in composition: cannot match {0} free output(s) with {1} free input(s)")]
    NonMatchingIO(usize, usize),

    #[error("I/O error: {0}")]
    IOError(#[from] std::io::Error),
}
pub type GraphResult<T> = Result<T, GraphError>;

pub(crate) mod spider;
pub use spider::*;

pub mod zx;
pub use zx::{ ZX, ZXNode, ZXWire };

pub mod zh;
pub use zh::{ ZH, ZHNode };

pub mod clifft;
pub use clifft::{ CT, CTNode, CTWire, TPhase, Complex };

pub mod rules;

pub(crate) mod diagram;
pub use diagram::*;

/// Identifies a node in a diagram.
pub type NodeId = usize;

/// Identifies an input or output qubit in a diagram.
pub type QubitId = usize;

