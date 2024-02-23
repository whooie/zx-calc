//! Graph-based tools to create and compute end products of diagrams in the
//! ZX(H)-calculus.
//!
//! Diagrams in the calculus are represented as an undirected, unweighted graph
//! with data attached to its nodes.
//!
//! See \[[1][pennylane]\] and \[[2][arxiv]\] for more info.
//!
//! \[1\] <https://pennylane.ai/qml/demos/tutorial_zx_calculus>
//!
//! \[2\] <https://arxiv.org/abs/2012.13966>
//!
//! [pennylane]: https://pennylane.ai/qml/demos/tutorial_zx_calculus
//! [arxiv]: https://arxiv.org/abs/2012.13966

use std::{
    collections::{ HashMap, HashSet },
    f64::consts::{ PI, TAU },
    fs,
    io::Write,
    ops::{ Deref, DerefMut },
    path::Path,
};
use num_complex::Complex64 as C64;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum GraphError {
    #[error("error adding wire: missing node {0}")]
    AddWireMissingNode(usize),

    #[error("error adding wire: input/output node {0} is already connected")]
    AddWireConnectedInputOutput(usize),

    #[error("error constructing GraphViz representation: {0}")]
    GraphVizError(String),

    #[error("I/O error: {0}")]
    IOError(#[from] std::io::Error),
}
pub type GraphResult<T> = Result<T, GraphError>;

/// A definition for a single node in a diagram.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum NodeDef {
    /// A Z-spider, parameterized by a real phase.
    Z(f64),
    /// An X-spider, parameterized by a real phase.
    X(f64),
    /// An H-box, parameterized by a general complex number.
    H(C64),
    /// Termination of a wire as an input to the diagram.
    Input,
    /// Termination of a wire as an output of the diagram.
    Output,
}

impl NodeDef {
    pub fn z() -> Self { Self::Z(0.0) }

    pub fn x() -> Self { Self::X(0.0) }

    pub fn h() -> Self { Self::H((-1.0).into()) }

    pub fn input() -> Self { Self::Input }

    pub fn output() -> Self { Self::Output }
}

/// The type of a single spider.
#[derive(Copy, Clone, Debug)]
pub enum Spider {
    Z,
    X,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Node {
    Z(f64),
    X(f64),
    H(C64),
    Input(usize),
    Output(usize),
}

impl From<Node> for NodeDef {
    fn from(node: Node) -> Self {
        match node {
            Node::Z(ph) => Self::Z(ph),
            Node::X(ph) => Self::X(ph),
            Node::H(a) => Self::H(a),
            Node::Input(_) => Self::Input,
            Node::Output(_) => Self::Output,
        }
    }
}

impl Node {
    /// Create a new spider of a certain kind with a certain phase.
    ///
    /// This function is guaranteed to return either [`Self::Z`] or [`Self::X`].
    pub fn new_spider(color: Spider, phase: f64) -> Self {
        match color {
            Spider::Z => Self::Z(phase),
            Spider::X => Self::X(phase),
        }
    }

    /// Return `true` if `self` is `Z`.
    pub fn is_z(&self) -> bool { matches!(self, Self::Z(_)) }

    /// Return `true` if `self` is `X`.
    pub fn is_x(&self) -> bool { matches!(self, Self::X(_)) }

    /// Return `true` if `self` is `H`.
    pub fn is_h(&self) -> bool { matches!(self, Self::H(_)) }

    /// Return `true` if `self` is `Input`.
    pub fn is_input(&self) -> bool { matches!(self, Self::Input(_)) }

    /// Return `true` if `self` is `Output`.
    pub fn is_output(&self) -> bool { matches!(self, Self::Output(_)) }

    /// Return `true` if `self` is a generator (`Z`, `X`, or `H`).
    pub fn is_generator(&self) -> bool {
        self.is_z() || self.is_x() || self.is_h()
    }

    /// Return `true` if `self` is `Z`, `X`, or `H` and has the default argument
    /// value. `Input` and `Output` return `false`.
    pub fn has_defarg(&self) -> bool {
        match *self {
            Self::Z(ph) => (ph % TAU).abs() < 1e-15,
            Self::X(ph) => (ph % TAU).abs() < 1e-15,
            Self::H(a) => (a + 1.0).norm() < 1e-15,
            Self::Input(_) => false,
            Self::Output(_) => false,
        }
    }
}

macro_rules! isomorphism {
    (
        $docstring:literal,
        $name:ident ($iso_to:ident),
        derive: { $($derive:ident),* $(,)? },
        from: { $($from:ident),* $(,)? } $(,)?
    ) => {
        #[doc = $docstring]
        #[derive($($derive),*)]
        pub struct $name(pub $iso_to);

        impl From<$iso_to> for $name {
            fn from(x: $iso_to) -> Self { Self(x) }
        }

        impl From<$name> for $iso_to {
            fn from(x: $name) -> Self { x.0 }
        }

        impl AsRef<$iso_to> for $name {
            fn as_ref(&self) -> &$iso_to { &self.0 }
        }

        impl AsMut<$iso_to> for $name {
            fn as_mut(&mut self) -> &mut $iso_to { &mut self.0 }
        }

        impl Deref for $name {
            type Target = $iso_to;

            fn deref(&self) -> &Self::Target { &self.0 }
        }

        impl DerefMut for $name {
            fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
        }

        $(
            impl From<$from> for $name {
                fn from(x: $from) -> Self { Self(x.into()) }
            }
        )*
    }
}

isomorphism!(
    "Sugared `usize` representing the ID of a single node.",
    NodeId (usize),
    derive: { Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash },
    from: { u8, u16 },
);

isomorphism!(
    "Sugared `usize` representing the ID of a single wire.",
    WireId (usize),
    derive: { Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash },
    from: { u8, u16 },
);

isomorphism!(
    "Sugared `usize` representing the ID of a single qubit.",
    QubitId (usize),
    derive: { Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash },
    from: { u8, u16 },
);

/// Represents the connection of one node to another.
///
/// Note that multiple wires are undirected and, since multiple edges are
/// allowed between nodes, a particular `Wire` value may not be unique in a
/// diagram even though node IDs are guaranteed to be unique.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct Wire(pub NodeId, pub NodeId);

impl<A, B> From<(A, B)> for Wire
where
    A: Into<NodeId>,
    B: Into<NodeId>,
{
    fn from(ab: (A, B)) -> Self {
        let (a, b) = ab;
        Self(a.into(), b.into())
    }
}

impl<A, B> From<Wire> for (A, B)
where
    NodeId: Into<A>,
    NodeId: Into<B>,
{
    fn from(wire: Wire) -> Self {
        let Wire(a, b) = wire;
        (a.into(), b.into())
    }
}

impl Wire {
    pub fn has_node<Id>(&self, id: Id) -> bool
    where Id: Into<NodeId>
    {
        let id = id.into();
        self.0 == id || self.1 == id
    }

    pub fn other_end_of<Id>(&self, start: Id) -> Option<NodeId>
    where Id: Into<NodeId>
    {
        let s = start.into();
        match *self {
            Wire(a, b) if a == s => Some(b),
            Wire(a, b) if b == s => Some(a),
            _ => None,
        }
    }
}

/// Like [`Wire`], but holding the data associated with the wire's endpoints
/// instead of their IDs.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct WireEnds(pub Node, pub Node);

/// Either a phaseless spider with two wires or two adjacent H-boxes of arity 2
/// and default argument.
#[derive(Copy, Clone, Debug)]
enum Identity {
    Spider(NodeId),
//  Spider(spider)
    HBox(NodeId, NodeId),
//  HBox(hbox1, hbox2)
}

/// Two spiders of the same color and arbitrary arity and phase with any number
/// of adjoining wires.
#[derive(Copy, Clone, Debug)]
struct Fuse(NodeId, NodeId);
//     Fuse(spider1, spider2)

/// Two spiders of opposite color and arbitrary arity and phase with exactly two
/// adjoining wires.
#[derive(Copy, Clone, Debug)]
struct Hopf(NodeId, NodeId);
//     Hopf(spider1, spider2)

/// Two spiders of the same color and arbitrary arity and phase with exactly two
/// adjoining wires, each with a single H-box with default argument.
#[derive(Copy, Clone, Debug)]
struct H2Hopf(NodeId, NodeId);
//     H2Hopf(hbox1, hbox2)

/// A spider of arbitrary arity and phase with a single self-loop containing a
/// single H-box with default argument.
#[derive(Copy, Clone, Debug)]
struct HLoop(NodeId, NodeId);
//     HLoop(spider, hbox)

/// Two spiders of arbitrarity arity and phase π/2 with the same color
/// sandwiching an H-box with default argument on a single wire.
#[derive(Copy, Clone, Debug)]
struct IHSandwich(NodeId, NodeId, NodeId);
//     IHSandwich(spider1, hbox, spider2)

/// Two spiders of arbitrary arity and phase with opposite colors sandwiching an
/// H-box with default argument on a single wire.
#[derive(Copy, Clone, Debug)]
struct HMove(NodeId, NodeId, NodeId);
//     HMove(spider1, hbox, spider2)

/// A spider of arbitrary arity and phase surrounded by H-boxes of default
/// argument.
#[derive(Clone, Debug)]
struct ColorSwap(NodeId, Vec<NodeId>);

/// A spider of arbitrary arity and phase with oppositely colored spiders of
/// arity 2 and phase π on all but one of its wires.
#[derive(Clone, Debug)]
struct PiCommute(NodeId, Vec<NodeId>, WireId);

/// `n` phaseless Z-spiders with all-to-all connectivity to `m` phaseless
/// X-spiders, where each Z-spider has arity `m + 1` and each X-spider has arity
/// `n + 1`.
#[derive(Clone, Debug)]
struct BiAlgebra(Vec<NodeId>, Vec<NodeId>);

/// Four spiders in a checkerboard square, each with arity 3, such that the
/// phases of the spiders of each color all have idential phases that are an
/// integer multiple of π.
#[derive(Copy, Clone, Debug)]
struct BitBiAlgebra((NodeId, NodeId), (NodeId, NodeId));

/// A spider of arbitrary arity and phase connected to a single oppositely
/// colored state/effect with phase equal to an integer multiple of π.
#[derive(Clone, Debug)]
struct StateCopy(NodeId, NodeId, Vec<WireId>);

/// H-box of arbitrary arity and argument connected to an H-box of aribtrary
/// arity and default argument via a single wire carrying a single H-box of
/// default argument.
#[derive(Copy, Clone, Debug)]
struct HFuse(NodeId, NodeId, NodeId);

/// An H-box of arbitrary arity connected to any number of X-states/effects,
/// each with phase π.
#[derive(Clone, Debug)]
struct HAbsorb(NodeId, Vec<NodeId>);

/// An H-box of arbitrary arity and argument connected to a single phaseless
/// X-state/effect.
#[derive(Clone, Debug)]
struct HExplode(NodeId, NodeId, Vec<WireId>);

/// An H-box of arity 1 with argument ±1.
#[derive(Copy, Clone, Debug)]
struct HState(NodeId);

/// An H-box of arbitrary arity with argument 1.
#[derive(Copy, Clone, Debug)]
struct HMultiState(NodeId);

/// An H-box of arity 1 and argument -1 connected to an H-box of arbitrary arity
/// and argument -1.
#[derive(Copy, Clone, Debug)]
struct HStateCopy(NodeId);

/// A phaseless Z-spider and an H-box of default argument, each of arity 3,
/// connected by exactly two adjoining wires.
#[derive(Copy, Clone, Debug)]
struct HHopf(NodeId, NodeId);

/// `n` phaseless Z-spiders with all-to-all connectivity to `m` default-argument
/// H-boxes, where each Z-spider has arity `m + 1` and each H-box has arity `n +
/// 1`.
#[derive(Clone, Debug)]
struct HBiAlgebra(Vec<NodeId>, Vec<NodeId>);

/// An X-spider of arity 2 and phase π connected to a phaseless Z-spider of
/// arity 3 over two wires, each with an H-box of arbitrary argument.
#[derive(Copy, Clone, Debug)]
struct HStateAvg((NodeId, NodeId), (NodeId, NodeId));

/// Two phaseless Z-spiders of arbitrary arities connected to each other by an
/// arbitrary number of wires, each with an H-box of arbitrary argument.
#[derive(Clone, Debug)]
struct HMul(NodeId, Vec<NodeId>, NodeId);

/// An arbitrary of H-boxes, each with arity 1, connected to a phaseless
/// Z-spider.
#[derive(Clone, Debug)]
struct HStateMul(Vec<NodeId>, NodeId);

type NodePath = [(NodeId, Node)];
type PathNodeSpec = fn(&Diagram, &NodePath, &NodeId, &Node) -> bool;

/// Represents a diagram in the calculus.
///
/// Every node and edge is given a unique index necessary to identify it when
/// adding or removing components of the graph. Note that multiple edges may
/// exist between two nodes.
///
/// All accessors are O(1).
#[derive(Clone, Debug)]
pub struct Diagram {
    nodes: HashMap<NodeId, Node>,
    wires: HashMap<WireId, Wire>,
    wires_from: HashMap<NodeId, HashSet<WireId>>,
    wires_between: HashMap<Wire, HashSet<WireId>>,
    node_id: usize,
    input_counter: usize,
    output_counter: usize,
    edge_id: usize,
}

impl Default for Diagram {
    fn default() -> Self { Self::new() }
}

impl Diagram {
    /// Create a new, empty diagram.
    pub fn new() -> Self {
        Self {
            nodes: HashMap::new(),
            wires: HashMap::new(),
            wires_from: HashMap::new(),
            wires_between: HashMap::new(),
            node_id: 0,
            input_counter: 0,
            output_counter: 0,
            edge_id: 0,
        }
    }

    /// Create a new, disconnected diagram given a set of nodes.
    ///
    /// All nodes passed to this function are given unique IDs starting from 0,
    /// i.e. if `n` nodes are passed to this function then the first node seen
    /// will have ID `0`, and the last will have ID `n - 1`. Input and output
    /// nodes will be given associated qubit numbers in the order they are seen.
    pub fn from_nodes<I>(nodes: I) -> Self
    where I: IntoIterator<Item = NodeDef>
    {
        let mut input_counter: usize = 0;
        let mut output_counter: usize = 0;
        let nodes: HashMap<NodeId, Node>
            = nodes.into_iter()
            .enumerate()
            .map(|(k, node_def)| {
                let node
                    = match node_def {
                        NodeDef::Z(ph) => Node::Z(ph),
                        NodeDef::X(ph) => Node::X(ph),
                        NodeDef::H(a) => Node::H(a),
                        NodeDef::Input => {
                            input_counter += 1;
                            Node::Input(input_counter - 1)
                        },
                        NodeDef::Output => {
                            output_counter += 1;
                            Node::Output(output_counter - 1)
                        },
                    };
                (k.into(), node)
            })
            .collect();
        let node_id = nodes.len();
        Self {
            nodes,
            wires: HashMap::new(),
            wires_from: HashMap::new(),
            wires_between: HashMap::new(),
            node_id,
            input_counter,
            output_counter,
            edge_id: 0,
        }
    }

    /// Return the number of nodes.
    pub fn count_nodes(&self) -> usize { self.nodes.len() }

    /// Return the number of wires.
    pub fn count_wires(&self) -> usize { self.wires.len() }

    /// Get the node associated with a particular ID if it exists.
    pub fn get_node<Id>(&self, node_id: Id) -> Option<&Node>
    where Id: Into<NodeId>
    {
        self.nodes.get(&node_id.into())
    }

    /// Get the number of wires attached to a node, if it exists.
    pub fn degree<Id>(&self, node_id: Id) -> Option<usize>
    where Id: Into<NodeId>
    {
        self.wires_from.get(&node_id.into()).map(|wires| wires.len())
    }

    /// Return `true` if a node has neighbors, if it exists.
    pub fn is_connected<Id>(&self, node_id: Id) -> Option<bool>
    where Id: Into<NodeId>
    {
        self.wires_from.get(&node_id.into()).map(|wires| !wires.is_empty())
    }

    /// Get the number of wires connecting two nodes, if they both exist.
    pub fn mutual_degree<A, B>(&self, a: A, b: B) -> Option<usize>
    where
        A: Into<NodeId>,
        B: Into<NodeId>,
    {
        let a = a.into();
        let b = b.into();
        if self.nodes.contains_key(&a) && self.nodes.contains_key(&b) {
            if let Some(wires) = self.wires_between.get(&Wire(a, b)) {
                Some(wires.len())
            } else {
                Some(0)
            }
        } else {
            None
        }
    }

    /// Add a node to the diagram and return its ID.
    pub fn add_node(&mut self, node_def: NodeDef) -> NodeId {
        let id = self.node_id;
        let node
            = match node_def {
                NodeDef::Z(ph) => Node::Z(ph),
                NodeDef::X(ph) => Node::X(ph),
                NodeDef::H(a) => Node::H(a),
                NodeDef::Input => {
                    self.input_counter += 1;
                    Node::Input(self.input_counter - 1)
                },
                NodeDef::Output => {
                    self.output_counter += 1;
                    Node::Output(self.output_counter - 1)
                },
            };
        self.nodes.insert(id.into(), node);
        self.wires_from.insert(id.into(), HashSet::new());
        self.node_id += 1;
        id.into()
    }

    /// Remove the node associated with a particular ID and return its data if
    /// it exists.
    ///
    /// This method also removes all wires with an endpoint at the node.
    pub fn remove_node<Id>(&mut self, node_id: Id) -> Option<Node>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        if let Some(wires_from) = self.wires_from(id) {
            let wire_ids: Vec<WireId> = wires_from.map(|(id, _)| id).collect();
            wire_ids.iter().for_each(|id| { self.wires.remove(id); });
            self.wires_from.values_mut()
                .for_each(|wires_from| {
                    wire_ids.iter().for_each(|id| { wires_from.remove(id); });
                });
            self.wires_between.values_mut()
                .for_each(|wires_btw| {
                    wire_ids.iter().for_each(|id| { wires_btw.remove(id); });
                });
            self.nodes.remove(&id)
        } else {
            None
        }
    }

    /// Remove all nodes that have no wires going to any other node in the
    /// diagram and return their IDs and data.
    pub fn remove_disconnected(&mut self) -> Vec<(NodeId, Node)> {
        let remove_ids: Vec<NodeId>
            = self.wires_from.iter()
            .filter_map(|(node_id, wires_from)| {
                wires_from.is_empty().then_some(*node_id)
            })
            .collect();
        remove_ids.into_iter()
            .filter_map(|id| self.remove_node(id).map(|node| (id, node)))
            .collect()
    }

    /// Get the wire associated with a particular ID if it exists.
    pub fn get_wire<Id>(&self, edge_id: Id) -> Option<&Wire>
    where Id: Into<WireId>
    {
        self.wires.get(&edge_id.into())
    }

    /// Add a wire between two nodes and return its ID.
    ///
    /// Fails if one or both of the nodes don't exist.
    pub fn add_wire<A, B>(&mut self, a: A, b: B) -> GraphResult<WireId>
    where
        A: Into<NodeId>,
        B: Into<NodeId>,
    {
        let a = a.into();
        let b = b.into();
        self.nodes.get(&a)
            .ok_or(GraphError::AddWireMissingNode(a.0))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(a).unwrap())
                    .then_some(())
                    .ok_or(GraphError::AddWireConnectedInputOutput(a.0))
            })?;
        self.nodes.get(&b)
            .ok_or(GraphError::AddWireMissingNode(b.0))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(b).unwrap())
                    .then_some(())
                    .ok_or(GraphError::AddWireConnectedInputOutput(b.0))
            })?;
        let id = self.edge_id;
        self.wires.insert(id.into(), Wire(a, b));
        if let Some(wires_from_a) = self.wires_from.get_mut(&a) {
            wires_from_a.insert(id.into());
        } else {
            self.wires_from.insert(a, [id.into()].into());
        }
        if let Some(wires_from_b) = self.wires_from.get_mut(&b) {
            wires_from_b.insert(id.into());
        } else {
            self.wires_from.insert(b, [id.into()].into());
        }
        if let Some(wires_btw) = self.wires_between.get_mut(&Wire(a, b)) {
            wires_btw.insert(id.into());
        } else {
            self.wires_between.insert(Wire(a, b), [id.into()].into());
            self.wires_between.insert(Wire(b, a), [id.into()].into());
        }
        self.edge_id += 1;
        Ok(id.into())
    }

    /// Add a wire with an attached [`Input`][Node::Input] to a pre-existing
    /// node and return the IDs of both the input and the wire as well as the
    /// index of the input qubit.
    ///
    /// The pre-existing node cannot already be an input or output.
    pub fn add_input_wire<Id>(&mut self, node_id: Id)
        -> GraphResult<(NodeId, WireId, QubitId)>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        self.nodes.get(&id)
            .ok_or(GraphError::AddWireMissingNode(id.0))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(GraphError::AddWireConnectedInputOutput(id.0))
            })?;
        let input_id = self.add_node(NodeDef::Input);
        let wire_id = self.add_wire(id, input_id)?;
        let qubit_id = self.input_counter - 1;
        Ok((input_id, wire_id, qubit_id.into()))
    }

    /// Add a wire with an attached [`Output`][Node::Output] to a pre-existing
    /// node and return the IDs of both the output and the wire as well as the
    /// index of the output qubit.
    ///
    /// The pre-existing node cannot already be an input or output.
    pub fn add_output_wire<Id>(&mut self, node_id: Id)
        -> GraphResult<(NodeId, WireId, QubitId)>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        self.nodes.get(&id)
            .ok_or(GraphError::AddWireMissingNode(id.0))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(GraphError::AddWireConnectedInputOutput(id.0))
            })?;
        let output_id = self.add_node(NodeDef::Output);
        let wire_id = self.add_wire(id, output_id)?;
        let qubit_id = self.output_counter - 1;
        Ok((output_id, wire_id, qubit_id.into()))
    }

    /// Remove the wire associated with a particular ID and return its endpoints
    /// if it exists.
    ///
    /// Note that this will not remove any nodes that become disconnected as a
    /// result of this operation. To remove an input/output node (which have
    /// enforced degree 1), use [`Self::remove_node`] instead.
    pub fn remove_wire<Id>(&mut self, wire_id: Id) -> Option<Wire>
    where Id: Into<WireId>
    {
        let id = wire_id.into();
        if self.wires.contains_key(&id) {
            self.wires_from.values_mut()
                .for_each(|wires_from| { wires_from.remove(&id); });
            self.wires_between.values_mut()
                .for_each(|wires_between| { wires_between.remove(&id); });
            self.wires.remove(&id)
        } else {
            None
        }
    }

    /// Return an iterator over all nodes in the diagram.
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
    pub fn nodes(&self) -> Nodes<'_> {
        Nodes { iter: self.nodes.iter() }
    }

    /// Return an iterator over all interior (non-input/output) nodes in the
    /// diagram.
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
    pub fn nodes_inner(&self) -> NodesInner<'_> {
        NodesInner { iter: self.nodes() }
    }

    /// Return an iterator over all diagram inputs.
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
    pub fn inputs(&self) -> Inputs<'_> {
        Inputs { iter: self.nodes.iter() }
    }

    /// Return an iterator over all diagram outputs.
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
    pub fn outputs(&self) -> Outputs<'_> {
        Outputs { iter: self.nodes.iter() }
    }

    /// Return an iterator over all wires in the diagram.
    ///
    /// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
    pub fn wires(&self) -> Wires<'_> {
        Wires { iter: self.wires.iter() }
    }

    /// Return an iterator over all interior (not connected to an input/output
    /// node) wires in the diagram.
    ///
    /// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
    pub fn wires_inner(&self) -> WiresInner<'_> {
        WiresInner { diagram: self, iter: self.wires() }
    }

    /// Return an iterator over all wires in the diagram with their endpoints'
    /// associated data.
    ///
    /// The iterator item type is `(`[`WireId`]`, `[`WireEnds`]`)`.
    pub fn wires_ends(&self) -> WiresEnds<'_> {
        WiresEnds { diagram: self, iter: self.wires.iter() }
    }

    /// Return an iterator over all interior (not connected to an input/output
    /// node) wires in the diagram with their associated data.
    ///
    /// The iterator item type is `(`[`WireId`]`, `[`WireEnds`]`)`.
    pub fn wires_ends_inner(&self) -> WiresEndsInner<'_> {
        WiresEndsInner { iter: self.wires_ends() }
    }

    /// Return an iterator over all wires with an endpoint at `node_id`, if the
    /// node exists.
    ///
    /// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
    pub fn wires_from<Id>(&self, node_id: Id) -> Option<WiresFrom<'_>>
    where Id: Into<NodeId>
    {
        self.wires_from.get(&node_id.into())
            .map(|wire_ids| {
                WiresFrom { diagram: self, iter: wire_ids.iter() }
            })
    }

    /// Return an iterator over all interior (not connected to an input/output
    /// node) wires with an endpoint at `node_id`, if the node exists and is not
    /// an input or output node.
    ///
    /// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
    pub fn wires_from_inner<Id>(&self, node_id: Id)
        -> Option<WiresFromInner<'_>>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        if let Some(node) = self.nodes.get(&id) {
            node.is_generator()
                .then(|| WiresFromInner { iter: self.wires_from(id).unwrap() })
        } else {
            None
        }
    }

    /// Return an iterator over all wires with endpoints at `a` and `b`, if both
    /// nodes exist and there is at least one wire between them.
    ///
    /// The iterator item type is [`WireId`].
    pub fn wires_between<A, B>(&self, a: A, b: B) -> Option<WiresBetween<'_>>
    where
        A: Into<NodeId>,
        B: Into<NodeId>,
    {
        self.wires_between.get(&(a, b).into())
            .and_then(|wire_ids| {
                (!wire_ids.is_empty())
                    .then_some(WiresBetween { iter: wire_ids.iter() })
            })
    }

    /// Return an iterator over all interior (not connected to an input/output
    /// node) wires with endpoints at `a` and `b`, if both nodes exist,
    /// neither are input or output nodes, and there is at least one wire
    /// between them.
    ///
    /// The iterator item type is [`WireId`].
    pub fn wires_between_inner<A, B>(&self, a: A, b: B)
        -> Option<WiresBetweenInner<'_>>
    where
        A: Into<NodeId>,
        B: Into<NodeId>,
    {
        let a = a.into();
        let b = b.into();
        (
            self.nodes.get(&a).is_some_and(|node| node.is_generator())
            && self.nodes.get(&b).is_some_and(|node| node.is_generator())
        ).then(|| {
            self.wires_between(a, b)
                .map(|iter| WiresBetweenInner { diagram: self, iter })
        })
        .flatten()
    }

    /// Return an iterator over all nodes neighboring the node at `node_id`, if
    /// it exists.
    pub fn neighbors_of<Id>(&self, node_id: Id) -> Option<Neighbors<'_>>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        self.wires_from.get(&id)
            .map(|wire_ids| {
                Neighbors { diagram: self, source: id, iter: wire_ids.iter() }
            })
    }

    /// Return an iterator over all interior (non-input/output) nodes
    /// neighboring the node at `node_id`, if it exists and is not an input or
    /// output node.
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
    pub fn neighbors_of_inner<Id>(&self, node_id: Id)
        -> Option<NeighborsInner<'_>>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        if self.nodes.get(&id).is_some_and(|node| node.is_generator()) {
            Some(NeighborsInner { iter: self.neighbors_of(id).unwrap() })
        } else {
            None
        }
    }

    /// Recursively find a finite sequence of nodes in the diagram satisfying a
    /// sequence of predicates, given a locally defined set of nodes to search
    /// over at each point in the path. The initial call to this function can
    /// pass any iterator as the set of nodes to search, but recursive calls
    /// will only use the set of neighbors of the last matching path node. Nodes
    /// with IDs matching any of those already in the path are automatically
    /// skipped. Any matching path can only have as many nodes as the number of
    /// predicates provided. Nodes matching a predicate are pushed onto the
    /// accumulator, which is returned as `Some` if all predicates are matched,
    /// otherwise `None` is returned instead.
    fn find_path<I>(
        &self,
        nodes: I,
        specs: &[PathNodeSpec],
        mut acc: Vec<(NodeId, Node)>,
    ) -> Option<Vec<(NodeId, Node)>>
    where I: Iterator<Item = (NodeId, Node)>
    {
        if let Some(spec) = specs.first() {
            for (id, node) in nodes {
                if acc.iter().any(|(path_id, _)| *path_id == id) {
                    continue;
                }
                if spec(self, &acc, &id, &node) {
                    acc.push((id, node));
                    return self.find_path(
                        self.neighbors_of(id).unwrap(),
                        &specs[1..],
                        acc,
                    );
                }
            }
            None
        } else {
            Some(acc)
        }
    }

    fn find_identity(&self) -> Option<Identity> {
        for (id, node) in self.nodes_inner() {
            if node.has_defarg() && self.degree(id).unwrap() == 2 {
                if let Node::Z(_) = node {
                    return Some(Identity::Spider(id));
                } else if let Node::X(_) = node {
                    return Some(Identity::Spider(id));
                } else if let Node::H(_) = node {
                    for (id2, node2) in self.neighbors_of(id).unwrap() {
                        if
                            node2.is_h()
                            && node2.has_defarg()
                            && self.degree(id).unwrap() == 2
                        {
                            return Some(Identity::HBox(id, id2));
                        }
                    }
                }
            }
        }
        None
    }

    /// Remove all identities from the diagram.
    ///
    /// The following structures are counted as identities:
    /// - a Z- or X-spider with two wires and zero phase
    /// - adjacent H-boxes, each with two wires and argument -1
    pub fn simplify_identity(&mut self) -> bool {
        if let Some(identity) = self.find_identity() {
            match identity {
                Identity::Spider(s) => {
                    let neighbors: Vec<NodeId>
                        = self.neighbors_of(s).unwrap()
                        .map(|(id, _)| id)
                        .collect();
                    self.remove_node(s);
                    self.add_wire(neighbors[0], neighbors[1]).ok();
                },
                Identity::HBox(h1, h2) => {
                    let neighbors: Vec<NodeId>
                        = self.neighbors_of(h1).unwrap()
                        .chain(self.neighbors_of(h2).unwrap())
                        .filter_map(|(id, _)| {
                            (id != h1 && id != h2).then_some(id)
                        })
                        .collect();
                    self.remove_node(h1);
                    self.remove_node(h2);
                    self.add_wire(neighbors[0], neighbors[1]).ok();
                },
            }
            true
        } else {
            false
        }
    }

    fn find_fuse(&self) -> Option<Fuse> {
        for (id, node) in self.nodes_inner() {
            if node.is_z() {
                for (id2, node2) in self.neighbors_of_inner(id).unwrap() {
                    if id2 == id {
                        continue;
                    }
                    if node2.is_z() {
                        return Some(Fuse(id, id2));
                    }
                }
            } else if node.is_x() {
                for (id2, node2) in self.neighbors_of_inner(id).unwrap() {
                    if id2 == id {
                        continue;
                    }
                    if node2.is_x() {
                        return Some(Fuse(id, id2));
                    }
                }
            }
        }
        None
    }

    /// Fuse adjacent spiders of the same color.
    pub fn simplify_fuse(&mut self) -> bool {
        if let Some(Fuse(s1, s2)) = self.find_fuse() {
            let neighbors2: Vec<NodeId>
                = self.neighbors_of(s2).unwrap()
                .filter_map(|(id, _)| (id != s1).then_some(id))
                .collect();
            let self_loops: usize
                = neighbors2.iter()
                .filter(|id| **id == s2)
                .count();
            let new
                = match (self.nodes[&s1], self.nodes[&s2]) {
                    (Node::Z(ph1), Node::Z(ph2)) => Node::Z((ph1 + ph2) % TAU),
                    (Node::X(ph1), Node::X(ph2)) => Node::X((ph1 + ph2) % TAU),
                    _ => unreachable!(),
                };
            self.remove_node(s2);
            if let Some(node) = self.nodes.get_mut(&s1) { *node = new; }
            neighbors2.into_iter()
                .for_each(|id| { self.add_wire(s1, id).ok(); });
            (0..self_loops)
                .for_each(|_| { self.add_wire(s1, s1).ok(); });
            true
        } else {
            false
        }
    }

    fn find_hopf(&self) -> Option<Hopf> {
        for (id, node) in self.nodes_inner() {
            if node.is_z() {
                for (id2, node2) in self.neighbors_of_inner(id).unwrap() {
                    if
                        node2.is_x()
                        && self.mutual_degree(id, id2).unwrap() == 2
                    {
                        return Some(Hopf(id, id2));
                    }
                }
            } else if node.is_x() {
                for (id2, node2) in self.neighbors_of_inner(id).unwrap() {
                    if
                        node2.is_z()
                        && self.mutual_degree(id, id2).unwrap() == 2
                    {
                        return Some(Hopf(id, id2));
                    }
                }
            }
        }
        None
    }

    /// Apply the Hopf rule, disconnecting two spiders of opposite color
    /// joined by two wires.
    pub fn simplify_hopf(&mut self) -> bool {
        if let Some(Hopf(s1, s2)) = self.find_hopf() {
            let wires: Vec<WireId>
                = self.wires_between(s1, s2).unwrap()
                .collect();
            wires.into_iter()
                .for_each(|id| { self.remove_wire(id); });
            true
        } else {
            false
        }
    }

    fn find_h2hopf(&self) -> Option<H2Hopf> {
        let n0 = |_: &Diagram, _: &NodePath, _: &NodeId, n: &Node| {
            n.is_z() || n.is_x()
        };
        let n1 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_h() && n.has_defarg() && dg.degree(*id).unwrap() == 2
        };
        let n2 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            (p[0].1.is_z() && n.is_z()) || (p[0].1.is_x() && n.is_x())
                && dg.mutual_degree(*id, p[0].0).is_some_and(|deg| deg == 0)
        };
        let n3 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            n.is_h() && n.has_defarg() && dg.degree(*id).unwrap() == 2
                && dg.mutual_degree(*id, p[0].0).is_some_and(|deg| deg == 1)
        };
        self.find_path(self.nodes_inner(), &[n0, n1, n2, n3], Vec::new())
            .map(|path| H2Hopf(path[1].0, path[3].0))
    }

    /// Apply a variant of the Hopf rule, disconnecting two spiders of arbitrary
    /// arity and phase joined by two wires, each with a single H-box of default
    /// argument.
    pub fn simplify_h2hopf(&mut self) -> bool {
        if let Some(H2Hopf(h1, h2)) = self.find_h2hopf() {
            self.remove_node(h1);
            self.remove_node(h2);
            true
        } else {
            false
        }
    }

    fn find_hloop(&self) -> Option<HLoop> {
        for (id, node) in self.nodes_inner() {
            if node.is_z() {
                for (id2, node2) in self.neighbors_of_inner(id).unwrap() {
                    if
                        node2.is_h()
                        && node2.has_defarg()
                        && self.degree(id2).unwrap() == 2
                        && self.mutual_degree(id, id2).unwrap() == 2
                    {
                        return Some(HLoop(id, id2));
                    }
                }
            }
        }
        None
    }

    /// Simplify a spider connected to an H-box in a single loop by removing the
    /// loop and adding π to the spider's phase.
    pub fn simplify_hloop(&mut self) -> bool {
        if let Some(HLoop(s, h)) = self.find_hloop() {
            let new
                = if let Node::Z(ph) = self.nodes[&s] {
                    Node::Z((ph + PI) % TAU)
                } else {
                    unreachable!()
                };
            self.remove_node(h);
            if let Some(node) = self.nodes.get_mut(&s) { *node = new; }
            true
        } else {
            false
        }
    }

    fn find_ih_sandwich(&self) -> Option<IHSandwich> {
        let n0 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            match n {
                Node::Z(ph) | Node::X(ph) => {
                    (ph.rem_euclid(TAU) / PI - 0.5).abs() < 1e-15
                        && dg.degree(*id).unwrap() == 2
                },
                _ => false,
            }
        };
        let n1 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_h() && n.has_defarg() && dg.degree(*id).unwrap() == 2
        };
        let n2 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            match (&p[0].1, n) {
                (Node::Z(_), Node::Z(ph)) | (Node::X(_), Node::X(ph)) => {
                    (ph.rem_euclid(TAU) / PI - 0.5).abs() < 1e-15
                        && dg.degree(*id).unwrap() == 2
                },
                _ => false,
            }
        };
        self.find_path(self.nodes_inner(), &[n0, n1, n2], Vec::new())
            .map(|path| IHSandwich(path[0].0, path[1].0, path[2].0))
    }

    /// Simplify two spiders of the same color with π/2 phase sandwiching an
    /// H-box on a single wire by rewriting the H-box as its Euler angles and
    /// fusing spiders appropriately.
    pub fn simplify_ih_sandwich(&mut self) -> bool {
        if let Some(IHSandwich(s1, h, s2)) = self.find_ih_sandwich() {
            let neighbors: Vec<NodeId>
                = self.neighbors_of(s1).unwrap()
                .chain(self.neighbors_of(s2).unwrap())
                .filter_map(|(id, _)| (id != h).then_some(id))
                .collect();
            let new
                = if self.nodes[&s1].is_z() {
                    Node::X(-PI / 2.0)
                } else {
                    Node::Z(-PI / 2.0)
                };
            self.remove_node(s1);
            self.remove_node(s2);
            if let Some(node) = self.nodes.get_mut(&h) { *node = new; }
            self.add_wire(neighbors[0], h).ok();
            self.add_wire(neighbors[1], h).ok();
            true
        } else {
            false
        }
    }

    fn find_hmove(&self) -> Option<HMove> {
        let n0 = |_: &Diagram, _: &NodePath, _: &NodeId, n: &Node| {
            n.is_z() || n.is_x()
        };
        let n1 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            n.is_h() && n.has_defarg()
                && dg.mutual_degree(*id, p[0].0).unwrap() == 1
        };
        let n2 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            ((p[0].1.is_z() && n.is_x()) || (p[0].1.is_x() && n.is_z()))
                && dg.mutual_degree(*id, p[1].0).unwrap() == 1
        };
        self.find_path(self.nodes_inner(), &[n0, n1, n2], Vec::new())
            .map(|path| HMove(path[0].0, path[1].0, path[2].0))
    }

    /// Simplify two spiders of opposite color sandwiching an H-box on a single
    /// wire by moving the H-box through the spider with fewer wires.
    #[allow(unused_variables)]
    pub fn simplify_hmove(&mut self) -> bool {
        if let Some(hmove) = self.find_hmove() {
            todo!()
        } else {
            false
        }
    }

    fn find_color_swap(&self) -> Option<ColorSwap> {
        todo!()
    }

    /// Simplify a spider surrounded by H-boxes by swapping the color of the
    /// spider and removing all H-boxes.
    #[allow(unused_variables)]
    pub fn simplify_color_swap(&mut self) -> bool {
        if let Some(color_swap) = self.find_color_swap() {
            todo!()
        } else {
            false
        }
    }

    fn find_pi_commute(&self) -> Option<PiCommute> {
        todo!()
    }

    /// Simplify by applying the π-commute rule.
    #[allow(unused_variables)]
    pub fn simplify_pi_commute(&mut self) -> bool {
        if let Some(pi_commute) = self.find_pi_commute() {
            todo!()
        } else {
            false
        }
    }

    fn find_bialgebra(&self) -> Option<BiAlgebra> {
        todo!()
    }

    /// Simplify by applying the bialgebra rule.
    #[allow(unused_variables)]
    pub fn simplify_bialgebra(&mut self) -> bool {
        if let Some(bialgebra) = self.find_bialgebra() {
            todo!()
        } else {
            false
        }
    }

    fn find_bit_bialgebra(&self) -> Option<BitBiAlgebra> {
        todo!()
    }

    /// Simplify by applying the 2-2 special case of the biagebra rule where
    /// the spiders have phases equal to integer multiples of π.
    #[allow(unused_variables)]
    pub fn simplify_bit_bialgebra(&mut self) -> bool {
        if let Some(bit_bialgebra) = self.find_bit_bialgebra() {
            todo!()
        } else {
            false
        }
    }

    fn find_state_copy(&self) -> Option<StateCopy> {
        todo!()
    }

    /// Simplify by copying a state or effect through an oppositely colored
    /// spider.
    #[allow(unused_variables)]
    pub fn simplify_state_copy(&mut self) -> bool {
        if let Some(state_copy) = self.find_state_copy() {
            todo!()
        } else {
            false
        }
    }

    fn find_hfuse(&self) -> Option<HFuse> {
        todo!()
    }

    /// Simplify by applying the fuse rule for H-boxes.
    #[allow(unused_variables)]
    pub fn simplify_hfuse(&mut self) -> bool {
        if let Some(hfuse) = self.find_hfuse() {
            todo!()
        } else {
            false
        }
    }

    fn find_habsorb(&self) -> Option<HAbsorb> {
        todo!()
    }

    /// Simplify by absorbing all X-states/effects of phase π into a neighboring
    /// H-box with arbitrary arity and argument.
    #[allow(unused_variables)]
    pub fn simplify_habsorb(&mut self) -> bool {
        if let Some(habsorb) = self.find_habsorb() {
            todo!()
        } else {
            false
        }
    }

    fn find_hexplode(&self) -> Option<HExplode> {
        todo!()
    }

    /// Simplify by exploding a single phaseless X-state/effect through a
    /// neighboring H-box with arbitrary arity and argument.
    #[allow(unused_variables)]
    pub fn simplify_hexplode(&mut self) -> bool {
        if let Some(hexplode) = self.find_hexplode() {
            todo!()
        } else {
            false
        }
    }

    fn find_hstate(&self) -> Option<HState> {
        todo!()
    }

    /// Simplify by converting an H-box with argument ±1 and one wire to a
    /// state represented by a Z-spider.
    #[allow(unused_variables)]
    pub fn simplify_hstate(&mut self) -> bool {
        if let Some(hstate) = self.find_hstate() {
            todo!()
        } else {
            false
        }
    }

    fn find_hmultistate(&self) -> Option<HMultiState> {
        todo!()
    }

    /// Simplify by expanding an H-box of arbitrary arity and argument 1 into a
    /// number of phaseless Z-states/effects equal to the original arity of the
    /// H-box.
    #[allow(unused_variables)]
    pub fn simplify_hmultistate(&mut self) -> bool {
        if let Some(hmultistate) = self.find_hmultistate() {
            todo!()
        } else {
            false
        }
    }

    fn find_hstate_copy(&self) -> Option<HStateCopy> {
        todo!()
    }

    /// Simplify by applying the all-H-box version of the state copy rule.
    #[allow(unused_variables)]
    pub fn simplify_hstate_copy(&mut self) -> bool {
        if let Some(hstate_copy) = self.find_hstate_copy() {
            todo!()
        } else {
            false
        }
    }

    fn find_hhopf(&self) -> Option<HHopf> {
        todo!()
    }

    /// Simplify by applying the H-box version of the Hopf rule.
    #[allow(unused_variables)]
    pub fn simplify_hhopf(&mut self) -> bool {
        if let Some(hhopf) = self.find_hhopf() {
            todo!()
        } else {
            false
        }
    }

    fn find_hbialgebra(&self) -> Option<HBiAlgebra> {
        todo!()
    }

    /// Simplify by applying the H-box version of the bialgebra rule.
    #[allow(unused_variables)]
    pub fn simplify_hbialgebra(&mut self) -> bool {
        if let Some(hbialgebra) = self.find_hbialgebra() {
            todo!()
        } else {
            false
        }
    }

    fn find_hstate_avg(&self) -> Option<HStateAvg> {
        todo!()
    }

    /// Simplify by applying the averaging rule.
    #[allow(unused_variables)]
    pub fn simplify_hstate_avg(&mut self) -> bool {
        if let Some(hstate_avg) = self.find_hstate_avg() {
            todo!()
        } else {
            false
        }
    }

    fn find_hmul(&self) -> Option<HMul> {
        todo!()
    }

    /// Simplify by applying the generalized multiplying rule.
    #[allow(unused_variables)]
    pub fn simplify_hmul(&mut self) -> bool {
        if let Some(hmul) = self.find_hmul() {
            todo!()
        } else {
            false
        }
    }

    fn find_hstate_mul(&self) -> Option<HStateMul> {
        todo!()
    }

    /// Simplify by applying the multiplying rule for H-box states.
    #[allow(unused_variables)]
    pub fn simplify_hstate_mul(&mut self) -> bool {
        if let Some(hstate_mul) = self.find_hstate_mul() {
            todo!()
        } else {
            false
        }
    }

    /// Simplify a diagram by applying all simplification rules until no further
    /// simplifications can be made.
    pub fn simplify(&mut self) { todo!() }

    /// Return an object containing an encoding of `self` in the [dot
    /// language][dot-lang].
    ///
    /// Rendering this object using the default formatter will result in a full
    /// dot string representation of the diagram.
    ///
    /// [dot-lang]: https://en.wikipedia.org/wiki/DOT_(graph_description_language)
    pub fn to_graphviz(&self, name: &str) -> GraphResult<tabbycat::Graph> {
        use tabbycat::*;
        use tabbycat::attributes::*;
        const SQUARE_HEIGHT: f64 = 0.15;
        const CIRCLE_HEIGHT: f64 = 0.20;
        // const Z_COLOR: Color = Color::White;
        // const X_COLOR: Color = Color::Gray50;
        // const H_COLOR: Color = Color::White;
        const Z_COLOR: Color = Color::Rgb(115, 150, 250);
        const X_COLOR: Color = Color::Rgb(230, 115, 125);
        const H_COLOR: Color = Color::Rgb(250, 205, 115);
        // initial declarations
        let mut statements
            = StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rankdir(RankDir::LR)),
            )
            .add_attr(
                AttrType::Node,
                AttrList::new()
                    .add_pair(fontname("DejaVu Sans"))
                    .add_pair(fontsize(10.0))
                    ,
            );
        // ensure all inputs are in a subgraph at the same rank
        let mut inputs_subgraph_stmt
            = StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rank(RankType::Min)),
            );
        let mut inputs: Vec<(NodeId, Node)> = self.inputs().collect();
        inputs.sort_by(|(l, _), (r, _)| l.cmp(r));
        for (NodeId(id), node) in inputs.into_iter() {
            let Node::Input(index) = node else { unreachable!() };
            inputs_subgraph_stmt
                = inputs_subgraph_stmt.add_node(
                    id.into(),
                    None,
                    Some(
                        AttrList::new()
                            .add_pair(label(format!("In {}", index)))
                            .add_pair(shape(Shape::Plaintext))
                    ),
                );
        }
        statements
            = statements.add_subgraph(SubGraph::cluster(inputs_subgraph_stmt));
        // ensure all outputs are in a subgraph at the same rank
        let mut outputs_subgraph_stmt
            = StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rank(RankType::Max)),
            );
        let mut outputs: Vec<(NodeId, Node)> = self.outputs().collect();
        outputs.sort_by(|(l, _), (r, _)| l.cmp(r));
        for (NodeId(id), node) in outputs.into_iter() {
            let Node::Output(index) = node else { unreachable!() };
            outputs_subgraph_stmt
                = outputs_subgraph_stmt.add_node(
                    id.into(),
                    None,
                    Some(
                        AttrList::new()
                            .add_pair(label(format!("Out {}", index)))
                            .add_pair(shape(Shape::Plaintext))
                    ),
                );
        }
        statements
            = statements.add_subgraph(SubGraph::cluster(outputs_subgraph_stmt));
        // add interior nodes
        for (NodeId(id), node) in self.nodes_inner() {
            let attrs: AttrList
                = match node {
                    Node::Z(ph) => {
                        let ph = ph.rem_euclid(TAU);
                        let ph_label
                            = if node.has_defarg() {
                                "".to_string()
                            } else if ph.rem_euclid(PI) < 1e-15 {
                                if (ph / PI - 1.0).abs() < 1e-15 {
                                    "π".to_string()
                                } else {
                                    format!("{:.0}π", ph / PI)
                                }
                            } else {
                                format!("{}π", ph / PI)
                            };
                        AttrList::new()
                            .add_pair(label(ph_label))
                            .add_pair(shape(Shape::Circle))
                            .add_pair(height(CIRCLE_HEIGHT))
                            .add_pair(style(Style::Filled))
                            .add_pair(fillcolor(Z_COLOR))
                    },
                    Node::X(ph) => {
                        let ph = ph.rem_euclid(TAU);
                        let ph_label
                            = if node.has_defarg() {
                                "".to_string()
                            } else if ph.rem_euclid(PI) < 1e-15 {
                                if (ph / PI - 1.0).abs() < 1e-15 {
                                    "π".to_string()
                                } else {
                                    format!("{:.0}π", ph / PI)
                                }
                            } else {
                                format!("{}π", ph / PI)
                            };
                        AttrList::new()
                            .add_pair(label(ph_label))
                            .add_pair(shape(Shape::Circle))
                            .add_pair(height(CIRCLE_HEIGHT))
                            .add_pair(style(Style::Filled))
                            .add_pair(fillcolor(X_COLOR))
                    },
                    Node::H(a) => {
                        let a_label
                            = if node.has_defarg() {
                                "".to_string()
                            } else {
                                format!("{}", a)
                            };
                        AttrList::new()
                            .add_pair(label(a_label))
                            .add_pair(shape(Shape::Square))
                            .add_pair(height(SQUARE_HEIGHT))
                            .add_pair(style(Style::Filled))
                            .add_pair(fillcolor(H_COLOR))
                    },
                    Node::Input(_) | Node::Output(_) => unreachable!(),
                };
            statements = statements.add_node(id.into(), None, Some(attrs));
        }
        // add wires
        for (_, Wire(NodeId(a), NodeId(b))) in self.wires() {
            statements
                = statements.add_edge(
                    Edge::head_node(a.into(), None)
                        .line_to_node(b.into(), None)
                );
        }
        GraphBuilder::default()
            .graph_type(GraphType::Graph)
            .strict(false)
            .id(Identity::quoted(name))
            .stmts(statements)
            .build()
            .map_err(GraphError::GraphVizError)
    }

    /// Like [`to_graphviz`][Self::to_graphviz], but render directly to a string
    /// and write it to `path`.
    pub fn save_graphviz<P>(&self, name: &str, path: P) -> GraphResult<()>
    where P: AsRef<Path>
    {
        let graphviz = self.to_graphviz(name)?;
        fs::OpenOptions::new()
            .write(true)
            .append(false)
            .create(true)
            .truncate(true)
            .open(path)?
            .write_all(format!("{}", graphviz).as_bytes())?;
        Ok(())
    }
}

/// Iterator type over the nodes in a diagram.
///
/// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
pub struct Nodes<'a> {
    iter: std::collections::hash_map::Iter<'a, NodeId, Node>
}

impl<'a> Iterator for Nodes<'a> {
    type Item = (NodeId, Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|(id, node)| (*id, *node))
    }
}

/// Iterator type over the interior (non-input/output) nodes in a diagram.
///
/// The iterator item type is `(`[`Node`]`, `[`Node`]`)`.
pub struct NodesInner<'a> {
    iter: Nodes<'a>
}

impl<'a> Iterator for NodesInner<'a> {
    type Item = (NodeId, Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, node)| node.is_generator())
    }
}

/// Iterator type over the inputs to a diagram.
///
/// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
pub struct Inputs<'a> {
    iter: std::collections::hash_map::Iter<'a, NodeId, Node>
}

impl<'a> Iterator for Inputs<'a> {
    type Item = (NodeId, Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, node)| node.is_input())
            .map(|(id, node)| (*id, *node))
    }
}

/// Iterator type over the outputs of a diagram.
///
/// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
pub struct Outputs<'a> {
    iter: std::collections::hash_map::Iter<'a, NodeId, Node>
}

impl<'a> Iterator for Outputs<'a> {
    type Item = (NodeId, Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, node)| node.is_output())
            .map(|(id, node)| (*id, *node))
    }
}

/// Iterator type over the wires in a diagram.
///
/// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
pub struct Wires<'a> {
    iter: std::collections::hash_map::Iter<'a, WireId, Wire>
}

impl<'a> Iterator for Wires<'a> {
    type Item = (WireId, Wire);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|(id, wire)| (*id, *wire))
    }
}

/// Iterator type over the interior (not connected to an input/output node)
/// wires in a diagram.
///
/// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
pub struct WiresInner<'a> {
    diagram: &'a Diagram,
    iter: Wires<'a>,
}

impl<'a> Iterator for WiresInner<'a> {
    type Item = (WireId, Wire);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, Wire(a, b))| {
            self.diagram.nodes[a].is_generator()
                && self.diagram.nodes[b].is_generator()
        })
    }
}


/// Iterator type over the endpoint data associated with the wires in a diagram.
///
/// The iterator item type is `(`[`WireId`]`, `[`WireEnds`]`)`.
pub struct WiresEnds<'a> {
    diagram: &'a Diagram,
    iter: std::collections::hash_map::Iter<'a, WireId, Wire>,
}

impl<'a> Iterator for WiresEnds<'a> {
    type Item = (WireId, WireEnds);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
            .map(|(id, Wire(a, b))| {
                (*id, WireEnds(self.diagram.nodes[a], self.diagram.nodes[b]))
            })
    }
}

/// Iterator type over the endpoint data associated with the interior (not
/// connected to an input/output node) wires in a diagram.
///
/// The iterator item type is `(`[`WireId`]`, `[`WireEnds`]`)`.
pub struct WiresEndsInner<'a> {
    iter: WiresEnds<'a>
}

impl<'a> Iterator for WiresEndsInner<'a> {
    type Item = (WireId, WireEnds);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, WireEnds(a, b))| {
            a.is_generator() && b.is_generator()
        })
    }
}

/// Iterator type over the wires attached to a particular node in a diagram.
///
/// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
pub struct WiresFrom<'a> {
    diagram: &'a Diagram,
    iter: std::collections::hash_set::Iter<'a, WireId>,
}

impl<'a> Iterator for WiresFrom<'a> {
    type Item = (WireId, Wire);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
            .map(|wire_id| (*wire_id, self.diagram.wires[wire_id]))
    }
}

/// Iterator type over the interior (not connected to an input/output node)
/// wires attached to a particular node in a diagram.
///
/// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
pub struct WiresFromInner<'a> {
    iter: WiresFrom<'a>
}

impl<'a> Iterator for WiresFromInner<'a> {
    type Item = (WireId, Wire);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, Wire(a, b))| {
            self.iter.diagram.nodes[a].is_generator()
                && self.iter.diagram.nodes[b].is_generator()
        })
    }
}

/// Iterator type over the wires attached to two particular nodes in a diagram.
///
/// The iterator item type is [`WireId`].
pub struct WiresBetween<'a> {
    iter: std::collections::hash_set::Iter<'a, WireId>,
}

impl<'a> Iterator for WiresBetween<'a> {
    type Item = WireId;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().copied()
    }
}

/// Iterator type over the interior (not connected to an input/output node)
/// wires attached to two particular nodes in a diagram.
///
/// The iterator item type is [`WireId`].
pub struct WiresBetweenInner<'a> {
    diagram: &'a Diagram,
    iter: WiresBetween<'a>,
}

impl<'a> Iterator for WiresBetweenInner<'a> {
    type Item = WireId;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|id| {
            let Wire(a, b) = &self.diagram.wires[id];
            self.diagram.nodes[a].is_generator()
                && self.diagram.nodes[b].is_generator()
        })
    }
}

/// Iterator type over the nodes neighboring a particular node in a diagram.
///
/// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
pub struct Neighbors<'a> {
    diagram: &'a Diagram,
    source: NodeId,
    iter: std::collections::hash_set::Iter<'a, WireId>,
}

impl<'a> Iterator for Neighbors<'a> {
    type Item = (NodeId, Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
            .map(|wire_id| {
                let other
                    = self.diagram.wires[wire_id].other_end_of(self.source)
                    .unwrap();
                (other, self.diagram.nodes[&other])
            })
    }
}

/// Iterator type over the interior (non-input/output) nodes neighboring a
/// particular node in a diagram.
///
/// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
pub struct NeighborsInner<'a> {
    iter: Neighbors<'a>
}

impl<'a> Iterator for NeighborsInner<'a> {
    type Item = (NodeId, Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, node)| node.is_generator())
    }
}

/// Create a diagram using abbreviated syntax.
///
/// The first block defines the nodes with the syntax
/// ```text
/// <node_label> = <node_type> [ [arg] ]
/// ```
/// where `<node_type>` is a variant or method of [`NodeDef`]. Labels given to
/// nodes in this block are then used to define wire connections in the second
/// block with the syntax
/// ```text
/// <node1> -- <node2>
/// ```
/// All nodes' labels and corresponding IDs are returned in a
/// [`HashMap`]`<&'static `[`str`]`, `[`NodeId`]`>` and all wires are given a
/// label `"<node1> - <node2>"` with labels and IDs returned in a
/// [`HashMap`]`<&'static `[`str`]`, `[`WireId`]`>`.
///
/// The total return type is [`Result`]`<(`[`Diagram`]`, `[`HashMap`]`<&'static
/// `[`str`]`, `[`NodeId`]`>, `[`HashMap`]`<&'static `[`str`]`, `[`WireId`]`>),
/// `[`GraphError`]`>`.
///
/// Up to hash maps and variable bindings, the usage
/// ```ignore
/// use std::f64::consts::PI;
///
/// diagram!(
///     nodes: {
///         i = input [],
///         o = output [],
///         n1 = Z [PI],
///         n2 = x [],
///         n3 = H [-1.0],
///     },
///     wires: {
///         i -- n1,
///         o -- n3,
///         n1 -- n2,
///         n1 -- n2,
///         n2 -- n3,
///     },
/// );
/// ```
/// is equivalent to
/// ```ignore
/// use std::f64::consts::PI;
///
/// let mut diagram = Diagram::new();
///
/// let i = diagram.add_node(NodeDef::Input);
/// let o = diagram.add_node(NodeDef::Output);
/// let n1 = diagram.add_node(NodeDef::Z(PI));
/// let n2 = diagram.add_node(NodeDef::x());
/// let n3 = diagram.add_node(NodeDef::H((-1.0).into()));
/// diagram.add_wire(i, n1)?;
/// diagram.add_wire(o, n3)?;
/// diagram.add_wire(n1, n2)?;
/// diagram.add_wire(n1, n2)?;
/// diagram.add_wire(n2, n3)?;
/// ```
#[macro_export]
macro_rules! diagram {
    (
        nodes: {
            $( $node_name:ident = $node_def:ident [$( $arg:expr ),* $(,)?] ),*
            $(,)?
        },
        wires: {
            $( $node1_name:ident -- $node2_name:ident ),* $(,)?
        } $(,)?
    ) => {
        {
            let mut diagram = $crate::graph::Diagram::new();
            $(
                let $node_name
                    = diagram.add_node(
                        NodeDef::$node_def($( ($arg).into() ),*)
                    );
            )*
            [$(
                diagram.add_wire($node1_name, $node2_name)
                    .map(|id| (stringify!($node1_name - $node2_name), id)),
            )*]
            .into_iter()
            .collect::<
                Result<
                    std::collections::HashMap<
                        &'static str,
                        $crate::graph::WireId
                    >,
                    $crate::graph::GraphError
                >
            >()
            .map(|edges| {
                let nodes:
                    std::collections::HashMap<
                        &'static str,
                        $crate::graph::NodeId
                    >
                    = [$( (stringify!($node_name), $node_name) ),*]
                    .into_iter()
                    .collect();
                (diagram, nodes, edges)
            })
        }
    }
    // (
    //     $( $node_name:ident = $node_def:ident [$( $arg:expr ),* $(,)?] ),*
    //     $( $node1_name:ident -- $node2_name:ident ),*
    //     $(,)?
    // ) => {
    //     diagram!(
    //         nodes: { $( $node_name = $node_def [$( $arg ),*] ),* }
    //         wires: { $( $node1_name -- $node2_name ),* }
    //     )
    // };
    // (
    //     $( $node_name:ident [$node_def:ident $(; $arg:expr )* $(;)?] );*
    //     $( $node1_name:ident -- $node2_name:ident );*
    //     $(;)?
    // ) => {
    //     diagram!(
    //         nodes: { $( $node_name = $node_def [$( $arg ),*] ),* }
    //         wires: { $( $node1_name -- $node2_name ),* }
    //     )
    // };
}

