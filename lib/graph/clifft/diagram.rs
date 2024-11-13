#![allow(unused_imports)]

use std::{ collections::VecDeque, fs, io::Write, path::Path };
use rustc_hash::FxHashMap;
use crate::graph::{
    GraphError,
    GraphResult,
    IONodeId,
    NodeId,
    QubitId,
    clifft::{ TNode, tphase::TPhase, complex::Complex },
};

use GraphError::*;

/// A wire in a [`CTDiagram`] connecting to a certain node ID.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Wire {
    /// A normal, empty wire.
    E(NodeId),
    /// A Hadamard wire.
    H(NodeId),
    /// A star wire.
    S(NodeId),
}

impl Wire {
    /// Return `true` if `self` is `E`.
    pub fn is_e(&self) -> bool { matches!(self, Self::E(_)) }

    /// Return `true` if `self` is `E` and the underlying node ID satisfies some
    /// predicate.
    pub fn is_e_and<F>(&self, pred: F) -> bool
    where F: FnOnce(NodeId) -> bool
    {
        match self {
            Self::E(id) => pred(*id),
            _ => false,
        }
    }

    /// Return `true` if `self` is `E` with an underlying node ID equal to `id`.
    pub fn has_e_id(&self, id: NodeId) -> bool {
        match self {
            Self::E(k) => *k == id,
            _ => false,
        }
    }

    pub(crate) fn make_e(&mut self) {
        match self {
            Self::E(_) => { },
            Self::H(id) | Self::S(id) => { *self = Self::E(*id); },
        }
    }

    /// Return `true` if `self` is `H`.
    pub fn is_h(&self) -> bool { matches!(self, Self::H(_)) }

    /// Return `true` if `self` is `H` and the underlying node ID satisfies some
    /// predicate.
    pub fn is_h_and<F>(&self, pred: F) -> bool
    where F: FnOnce(NodeId) -> bool
    {
        match self {
            Self::H(id) => pred(*id),
            _ => false,
        }
    }

    /// Return `true` if `self` is `H` with an underlying node ID equal to `id`.
    pub fn has_h_id(&self, id: NodeId) -> bool {
        match self {
            Self::H(k) => *k == id,
            _ => false,
        }
    }

    pub(crate) fn make_h(&mut self) {
        match self {
            Self::H(_) => { },
            Self::E(id) | Self::S(id) => { *self = Self::H(*id); },
        }
    }

    /// Return `true` if `self` is `S`.
    pub fn is_s(&self) -> bool { matches!(self, Self::S(_)) }

    /// Return `true` if `self` is `S` and the underlying node ID satisfies some
    /// predicate.
    pub fn is_s_and<F>(&self, pred: F) -> bool
    where F: FnOnce(NodeId) -> bool
    {
        match self {
            Self::S(id) => pred(*id),
            _ => false,
        }
    }

    /// Return `true` if `self` is `S` with an underlying node ID equal to `id`.
    pub fn has_s_id(&self, id: NodeId) -> bool {
        match self {
            Self::S(k) => *k == id,
            _ => false,
        }
    }

    pub(crate) fn make_s(&mut self) {
        match self {
            Self::S(_) => { },
            Self::E(id) | Self::H(id) => { *self = Self::S(*id); },
        }
    }

    /// Return `true` if `self` has an underlying node ID equal to `id`.
    pub fn has_id(&self, id: NodeId) -> bool {
        match self {
            Self::E(k) | Self::H(k) | Self::S(k) => *k == id
        }
    }

    /// Return the underlying node ID if `self` is `E`.
    pub fn e_id(&self) -> Option<NodeId> {
        match self {
            Self::E(id) => Some(*id),
            _ => None,
        }
    }

    /// Return the underlying node ID if `self` is `H`.
    pub fn h_id(&self) -> Option<NodeId> {
        match self {
            Self::H(id) => Some(*id),
            _ => None,
        }
    }

    /// Return the underlying node ID if `self` is `S`.
    pub fn s_id(&self) -> Option<NodeId> {
        match self {
            Self::S(id) => Some(*id),
            _ => None,
        }
    }

    /// Return the underlying node ID.
    pub fn id(&self) -> NodeId {
        match self {
            Self::E(id) | Self::H(id) | Self::S(id) => *id
        }
    }

    pub(crate) fn shift_id(&mut self, sh: usize) {
        match self {
            Self::E(id) | Self::H(id) | Self::S(id) => { *id += sh; }
        }
    }
}

/// Represents a diagram in the ZX-calculus, restricted to Clifford+*T*
/// "graph-like" forms, plus star wires.
///
/// A ZX-diagram is "graph-like" if it contains only Z-spiders connected to each
/// other via Hadamard wires. Here, we also allow star wires and limit all
/// phases to integer multiples of *π*/4.
#[derive(Clone, Debug)]
pub struct CTDiagram {
    pub(crate) nodes: Vec<Option<TNode>>,
    pub(crate) node_count: usize,
    pub(crate) wires: Vec<Option<Vec<Wire>>>,
    pub(crate) wire_count: usize,
    pub(crate) inputs: Vec<IONodeId>,
    pub(crate) outputs: Vec<IONodeId>,
    pub(crate) free: Vec<NodeId>,
    pub(crate) scalar: Complex,
}

impl Default for CTDiagram {
    fn default() -> Self { Self::new() }
}

impl CTDiagram {
    /// Create a new, empty diagram.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            node_count: 0,
            wires: Vec::new(),
            wire_count: 0,
            inputs: Vec::new(),
            outputs: Vec::new(),
            free: Vec::new(),
            scalar: Complex::ONE,
        }
    }

    /// Create a new, disconnected diagram from a set of nodes.
    ///
    /// All nodes passed to this function are given unique IDs starting from 0,
    /// i.e. if `n` nodes are passed to this function then the first node seen
    /// will have ID `0`, and the last will have ID `n - 1`. Input and output
    /// nodes will be given associated qubit numbers in the order they are seen,
    /// starting from zero.
    pub fn from_nodes<I>(nodes: I) -> Self
    where I: IntoIterator<Item = TNode>
    {
        let mut dg = Self::new();
        nodes.into_iter()
            .for_each(|def| { dg.add_node(def); });
        dg
    }

    /// Return the global scalar on the diagram.
    pub fn scalar(&self) -> Complex { self.scalar }

    /// Apply a mapping function to the global scalar.
    pub fn map_scalar<F>(&mut self, map: F)
    where F: FnOnce(Complex) -> Complex
    {
        self.scalar = map(self.scalar);
    }

    /// Return `true` if the gobal scalar is zero.
    pub fn scalar_is_zero(&self) -> bool { self.scalar == Complex::ZERO }

    /// Return the number of nodes.
    pub fn count_nodes(&self) -> usize { self.node_count }

    /// Return the number of Z-spiders.
    pub fn count_z(&self) -> usize {
        self.nodes.iter().filter(|mb_n| mb_n.is_some_and(|n| n.is_z())).count()
    }

    /// Return the number of inputs.
    pub fn count_inputs(&self) -> usize { self.inputs.len() }

    /// Return the number of free input wires.
    pub fn count_free_inputs(&self) -> usize {
        self.inputs.iter().filter(|input| input.is_free()).count()
    }

    /// Return the number of input wires with attached states.
    pub fn count_state_inputs(&self) -> usize {
        self.inputs.iter().filter(|input| input.is_state()).count()
    }

    /// Return the number of outputs.
    pub fn count_outputs(&self) -> usize { self.outputs.len() }

    /// Return the number of free output wires.
    pub fn count_free_outputs(&self) -> usize {
        self.outputs.iter().filter(|output| output.is_free()).count()
    }

    /// Return the number of output wires with attached effects.
    pub fn count_effect_outputs(&self) -> usize {
        self.outputs.iter().filter(|output| output.is_state()).count()
    }

    /// Get the node associated with a particular ID if it exists.
    pub fn get_node(&self, id: NodeId) -> Option<&TNode> {
        self.nodes.get(id).and_then(|mb_n| mb_n.as_ref())
    }

    pub(crate) fn get_node_mut(&mut self, id: NodeId) -> Option<&mut TNode> {
        self.nodes.get_mut(id).and_then(|mb_n| mb_n.as_mut())
    }

    /// Return `true` if a node exists with the given ID.
    pub fn has_node(&self, id: NodeId) -> bool {
        self.nodes.get(id).is_some_and(|mb_n| mb_n.is_some())
    }

    /// Return the qubit index corresponding to a node ID, if the node exists
    /// and is an input.
    pub fn get_input_index(&self, id: NodeId) -> Option<QubitId> {
        self.get_node(id)
            .and_then(|n| {
                (!n.is_output()).then(|| {
                    self.inputs.iter().enumerate()
                        .find_map(|(qid, ioid)| ioid.has_id(id).then_some(qid))
                        .expect("bad book-keeping: unrecorded input")
                })
            })
    }

    // like `get_input_index`, but for use after the node has already been
    // removed -- this is purely a search through `self.inputs`
    pub(crate) fn find_input_index(&self, id: NodeId) -> Option<QubitId> {
        self.inputs.iter().enumerate()
            .find_map(|(qid, ioid)| ioid.has_id(id).then_some(qid))
    }

    /// Return the qubit index corresponding to a node ID, if the node exists
    /// and is an output.
    pub fn get_output_index(&self, id: NodeId) -> Option<QubitId> {
        self.get_node(id)
            .and_then(|n| {
                (!n.is_input()).then(|| {
                    self.outputs.iter().enumerate()
                        .find_map(|(qid, ioid)| ioid.has_id(id).then_some(qid))
                        .expect("bad book-keeping: unrecorded output")
                })
            })
    }

    // like `get_output_index`, but for use after the node has already been
    // removed -- this is purely a search through `self.outputs`
    pub(crate) fn find_output_index(&self, id: NodeId) -> Option<QubitId> {
        self.outputs.iter().enumerate()
            .find_map(|(qid, ioid)| ioid.has_id(id).then_some(qid))
    }

    // get the raw neighbor IDs of a node, if it exists
    pub(crate) fn get_neighbors_mut(&mut self, id: NodeId)
        -> Option<&mut Vec<Wire>>
    {
        self.wires.get_mut(id).and_then(|mb_nnb| mb_nnb.as_mut())
    }

    /// Get the number of wires attached to a node, if it exists.
    pub fn arity(&self, id: NodeId) -> Option<usize> {
        self.wires.get(id)
            .and_then(|mb_nnb| mb_nnb.as_ref().map(|nnb| nnb.len()))
    }

    /// Return `true` if a node has neighbors, if it exists.
    pub fn is_connected(&self, id: NodeId) -> Option<bool> {
        self.wires.get(id)
            .and_then(|mb_nnb| mb_nnb.as_ref().map(|nnb| !nnb.is_empty()))
    }

    /// Get the number of wires connecting two nodes, if they both exist.
    pub fn mutual_arity(&self, a: NodeId, b: NodeId) -> Option<usize> {
        self.has_node(b).then_some(())?;
        self.wires.get(a)
            .and_then(|mb_nnb| {
                mb_nnb.as_ref()
                    .map(|nnb| {
                        let ma =
                            nnb.iter().filter(|wire| wire.has_id(b)).count();
                        if a == b { ma / 2 } else { ma }
                    })
            })
    }

    // pop a freed ID off the list or allocate a new one
    pub(crate) fn fresh_node_id(&mut self) -> NodeId {
        if let Some(id) = self.free.pop() {
            self.node_count += 1;
            id
        } else {
            let id = self.nodes.len();
            self.nodes.push(None);
            self.wires.push(None);
            self.node_count += 1;
            id
        }
    }

    /// Add a node to the diagram and return its ID.
    pub fn add_node(&mut self, node: TNode) -> NodeId {
        let id = self.fresh_node_id();
        if node.is_input() {
            self.inputs.push(IONodeId::Free(id));
        } else if node.is_output() {
            self.outputs.push(IONodeId::Free(id));
        }
        self.nodes[id] = Some(node);
        self.wires[id] = Some(Vec::new());
        id
    }

    /// Add a Z-spider to the diagram and return its ID.
    ///
    /// `phase` defaults to zero.
    pub fn add_z(&mut self, phase: Option<TPhase>) -> NodeId {
        self.add_node(TNode::Z(phase.unwrap_or(TPhase::T0)))
    }

    /// Add an input node to the diagram and return its ID.
    pub fn add_input(&mut self) -> NodeId {
        self.add_node(TNode::Input)
    }

    /// Add an output node to the diagram and return its ID.
    pub fn add_output(&mut self) -> NodeId {
        self.add_node(TNode::Output)
    }

    /// Remove the node associated with a particular ID and return its data if
    /// it exists.
    ///
    /// This method also removes all wires with an endpoint at the node.
    pub fn remove_node(&mut self, id: NodeId) -> GraphResult<TNode> {
        let node =
            self.nodes.get_mut(id)
            .ok_or(RemoveNodeMissingNode(id))?
            .take()
            .ok_or(RemoveNodeMissingNode(id))?;
        self.node_count -= 1;
        self.free.push(id);
        let nnb_of = self.wires[id].take().unwrap();
        let self_loops = nnb_of.iter().filter(|nb| nb.has_id(id)).count() / 2;
        self.wire_count -= nnb_of.len() - self_loops;
        for nb_of in nnb_of.into_iter() {
            while let Some(k) =
                self.wires[nb_of.id()].as_ref()
                    .and_then(|nnb| {
                        nnb.iter().enumerate()
                            .find_map(|(k, nb)| nb.has_id(id).then_some(k))
                    })
            {
                self.wires[nb_of.id()].as_mut().unwrap()
                    .swap_remove(k);
            }
        }
        if node.is_input() || self.inputs.iter().any(|ioid| ioid.has_id(id)) {
            let k = self.find_input_index(id).unwrap();
            self.inputs.remove(k);
        }
        if node.is_output() || self.outputs.iter().any(|ioid| ioid.has_id(id)) {
            let k = self.find_output_index(id).unwrap();
            self.outputs.remove(k);
        }
        Ok(node)
    }

    /// Like [`remove_node`][Self::remove_node], but returning a list of the
    /// node's neighbors alongside its data. Note that any H-boxes or stars
    /// living on wires connected to the node are removed as well.
    pub fn remove_node_nb(&mut self, id: NodeId)
        -> GraphResult<(TNode, Vec<Wire>)>
    {
        let node =
            self.nodes.get_mut(id)
            .ok_or(RemoveNodeMissingNode(id))?
            .take()
            .ok_or(RemoveNodeMissingNode(id))?;
        self.node_count -= 1;
        self.free.push(id);
        let nnb_of = self.wires[id].take().unwrap();
        let self_loops = nnb_of.iter().filter(|nb| nb.has_id(id)).count() / 2;
        self.wire_count -= nnb_of.len() - self_loops;
        for nb_of in nnb_of.iter() {
            while let Some(k) =
                self.wires[nb_of.id()].as_ref()
                    .and_then(|nnb| {
                        nnb.iter().enumerate()
                            .find_map(|(k, nb)| nb.has_id(id).then_some(k))
                    })
            {
                self.wires[nb_of.id()].as_mut().unwrap()
                    .swap_remove(k);
            }
        }
        if node.is_input() || self.inputs.iter().any(|ioid| ioid.has_id(id)) {
            let k = self.find_input_index(id).unwrap();
            self.inputs.remove(k);
        }
        if node.is_output() || self.outputs.iter().any(|ioid| ioid.has_id(id)) {
            let k = self.find_output_index(id).unwrap();
            self.outputs.remove(k);
        }
        Ok((node, nnb_of))
    }
}

