//! Graph-based tools to create and compute end products of diagrams in the
//! ZX(H)-calculus.
//!
//! Diagrams are represented as an undirected, unweighted graph with parallel
//! edges and data attached to its nodes. See
//! [`graph_diagram!`][crate::graph_diagram] for example usage and abbreviated
//! syntax.

use thiserror::Error;

/// Errors for fallible operations on graph-backed diagrams.
#[derive(Debug, Error)]
pub enum GraphError {
    /// Returned when a node does not exist under a given ID.
    #[error("missing node {0}")]
    MissingNode(NodeId),

    /// Returned when a diagram input or output does not exist under a given
    /// index.
    #[error("missing qubit {0}")]
    MissingQubit(QubitId),

    /// Returned when attempting to attach a new wire to an `Input` or `Output`
    /// node that already has a wire attached.
    #[error("input/output node {0} is already connected")]
    ConnectedIO(NodeId),

    /// Returned when attempting to assign a state to a diagram input or output
    /// that has no neighbors.
    #[error("cannot apply state to a disconnected input/output node")]
    DisconnectedIO,

    /// Returned when attempting to assign a state to a diagram input or output
    /// that does not exist or already has a state assigned.
    #[error("node {0} is not an input/output")]
    NotIO(NodeId),

    /// Returned when attempting to remove an input or output state node that
    /// does not exist or is not an input or output state to a diagram.
    #[error("node {0} is not an input/output state")]
    NotIOState(NodeId),

    /// Returned when attempting to attach a Bell state/effect wire to two nodes
    /// that do not both exist, or are not both inputs or outputs.
    #[error("nodes {0}, {1} are not both inputs/outputs")]
    BellNotIO(NodeId, NodeId),

    /// Returned when attempting to attach a Bell state/effect wire to the same
    /// input or output node.
    #[error("cannot apply Bell state/effect to a single qubit")]
    SameQubit,

    /// Returned when composition is attempted between two diagrams that have
    /// non-matching inputs and outputs.
    #[error("cannot match {0} free output(s) with {1} free input(s)")]
    ComposeIO(usize, usize),

    /// I/O error when writing a diagram to a file.
    #[error("{0}")]
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

