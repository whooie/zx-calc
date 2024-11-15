use num_complex::Complex64 as C64;
use num_traits::One;
use rustc_hash::FxHashMap;
use crate::graph::{ GraphError, GraphResult, NodeId, QubitId, ZX, ZH, CT };
use GraphError::*;

/// Associated types for a graph-backed diagram implementing some variant of the
/// ZX-calculus.
pub trait DiagramData {
    /// Node type.
    type Node: NodeData;

    /// Wire type.
    type Wire: WireData;

    /// Complex-valued diagram scalar type.
    type Scalar: ComplexRing;
}

/// Basic requirements for a node type in a diagram.
///
/// Node types are allowed to implement any set of generators for a diagram, but
/// must have variants for simple input/output nodes that indicate empty wires
/// at the boundaries of diagrams. Nodes must also have an "identity" that, when
/// placed on a single empty wire (i.e. the generator has arity 2), corresponds
/// to the identity on the wire.
pub trait NodeData: Clone + std::fmt::Debug {
    /// Create a new identity.
    fn new_id() -> Self;

    /// Create a new input node.
    fn new_input() -> Self;

    /// Create a new output node.
    fn new_output() -> Self;

    /// Return `true` if `self` is an input node.
    ///
    /// This must give consistent results with [`new_input`][Self::new_input].
    fn is_input(&self) -> bool;

    /// Return `true` if `self` is an output node.
    ///
    /// This must give consistent results with [`new_output`][Self::new_output].
    fn is_output(&self) -> bool;

    /// Return `true` if `self` is neither an [input][Self::new_input] nor an
    /// [output][Self::new_output].
    fn is_generator(&self) -> bool { !self.is_input() && !self.is_output() }

    /// Conjugate `self` in place, replacing it with its adjoint.
    ///
    /// This must take inputs to outputs and vice-versa.
    fn conjugate(&mut self);
}

/// Basic requirements for a wire type in a diagram.
///
/// Each wire is taken relative to a fixed central node, and is identified with
/// a single neighbor's ID, as used in e.g. [`Diagram::neighbors`]. This ID
/// must be accessible from any value belonging to the wire type, and the type's
/// implementation must allow the ID to be dynamically adjusted in order to
/// facilitate diagram composition and tensoring.
///
/// Additionally, wire types must have a variant associated with an ordinary,
/// empty wire in the base ZX-calculus.
pub trait WireData: Clone + std::fmt::Debug {
    /// Return the inner neighbor node ID.
    fn id(&self) -> NodeId;

    /// Return `true` if the inner neighbor node ID is equal to `id`.
    ///
    /// This must give consistent results with [`id`][Self::id].
    fn has_id(&self, id: NodeId) -> bool { self.id() == id }

    /// Return the inner neighbor node ID if it satisfies some predicate.
    fn id_if<F>(&self, pred: F) -> Option<NodeId>
    where F: FnOnce(NodeId) -> bool
    {
        let id = self.id();
        pred(id).then_some(id)
    }

    /// Apply a mapping function to the inner neighbor node ID in place.
    fn map_id<F>(&mut self, map: F)
    where F: FnOnce(NodeId) -> NodeId;

    /// Increase value of the inner neighbor node ID by `sh`.
    fn shift_id(&mut self, sh: usize) { self.map_id(|id| id + sh); }

    /// Create a new, empty wire variant with neighbor node ID `id`.
    fn new_empty(id: NodeId) -> Self;

    /// Return `true` if `self` corresponds to an empty wire.
    ///
    /// This must give consistent results with [`new_empty`][Self::new_empty].
    fn is_empty(&self) -> bool;
}

impl WireData for NodeId {
    fn id(&self) -> NodeId { *self }

    fn map_id<F>(&mut self, map: F)
    where F: FnOnce(NodeId) -> NodeId
    {
        *self = map(*self);
    }

    fn new_empty(id: NodeId) -> Self { id }

    fn is_empty(&self) -> bool { true }
}

/// A description of numerical values belonging to a complex-valued ring.
///
/// Specifically, this type must be closed under addition and multiplication
/// with identity elements for both. The addition must have inverse elements
/// (i.e. subtraction must be defined), but multiplication need not have them
/// (division need not exist).
pub trait ComplexRing:
    Copy
    + Clone
    + PartialEq<Self>
    + num_traits::Zero
    + num_traits::One
    + std::ops::Add<Self, Output = Self>
    + std::ops::AddAssign<Self>
    + std::ops::Sub<Self, Output = Self>
    + std::ops::SubAssign<Self>
    + std::ops::Mul<Self, Output = Self>
    + std::ops::MulAssign<Self>
    + std::fmt::Debug
{
    /// Return the complex conjugate of `self`.
    fn conj(self) -> Self;

    /// Return the square of the norm of `self`.
    fn norm2(self) -> Self { self * self.conj() }
}

impl ComplexRing for C64 {
    fn conj(self) -> Self { C64::conj(&self) }
}

/// The ID associated with a particular input/output node of a diagram.
///
/// An `IONodeId` can be either `Free`, where the node itself is properly an
/// [input][NodeData::new_input] or [output][NodeData::new_output] to the
/// diagram, or a `State`, where a particular unary generator has been assigned
/// to a previously free input/output.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum IONodeId {
    /// A proper free input/output to a diagram.
    Free(NodeId),
    /// A previously free input/output that now has a specified state.
    State(NodeId),
}

impl IONodeId {
    /// Return `true` if `self` is `Free`.
    pub fn is_free(&self) -> bool { matches!(self, Self::Free(_)) }

    /// Return `true` if `self` is `Free` and the underlying node ID satisfies
    /// some predicate.
    pub fn is_free_and<F>(&self, pred: F) -> bool
    where F: FnOnce(NodeId) -> bool
    {
        match self {
            Self::Free(id) => pred(*id),
            Self::State(_) => false,
        }
    }

    pub(crate) fn make_free(&mut self) {
        match self {
            Self::Free(_) => { },
            Self::State(id) => { *self = Self::Free(*id); },
        }
    }

    /// Return `true` if `self` is `State`.
    pub fn is_state(&self) -> bool { matches!(self, Self::State(_)) }

    /// Return `true` if `self` is `State` and the underlying node ID satisfies
    /// some predicate.
    pub fn is_state_and<F>(&self, pred: F) -> bool
    where F: FnOnce(NodeId) -> bool
    {
        match self {
            Self::Free(_) => false,
            Self::State(id) => pred(*id),
        }
    }

    pub(crate) fn make_state(&mut self) {
        match self {
            Self::Free(id) => { *self = Self::State(*id); },
            Self::State(_) => { },
        }
    }

    /// Return the underlying node ID.
    pub fn id(&self) -> NodeId {
        match self {
            Self::Free(id) => *id,
            Self::State(id) => *id,
        }
    }

    /// Return `true` if `self` has an underlying node ID equal to `id`.
    pub fn has_id(&self, id: NodeId) -> bool {
        match self {
            Self::Free(k) => *k == id,
            Self::State(k) => *k == id,
        }
    }

    /// Return `true` if `self` is `Free` and has an underlying node ID equal to
    /// `id`.
    pub fn has_free_id(&self, id: NodeId) -> bool {
        match self {
            Self::Free(k) => *k == id,
            Self::State(_) => false,
        }
    }

    /// Return `true` if `self` is `State` and has an underlying node ID equal
    /// to `id`.
    pub fn has_state_id(&self, id: NodeId) -> bool {
        match self {
            Self::Free(_) => false,
            Self::State(k) => *k == id,
        }
    }

    /// Return the underlying node ID if `self` is `Free`.
    pub fn free_id(&self) -> Option<NodeId> {
        match self {
            Self::Free(id) => Some(*id),
            Self::State(_) => None,
        }
    }

    /// Return the underlying node ID if `self` is `State`.
    pub fn state_id(&self) -> Option<NodeId> {
        match self {
            Self::Free(_) => None,
            Self::State(id) => Some(*id),
        }
    }

    pub(crate) fn shift_id(&mut self, sh: usize) {
        match self {
            Self::Free(id) => { *id += sh; },
            Self::State(id) => { *id += sh; },
        }
    }
}

/// Represents a diagram in a ZX-calculus(-like) graphical language.
///
/// Every node is given a unique [index][NodeId] for identification purposes and
/// parallel edges are, in general, allowed. Input and output qubits are
/// numbered dynamically, but the order in which they are added to the diagram
/// is always preserved.
///
/// This diagram type is parameterized by a type representing a particular
/// variant of the ZX-calculus, which specifies a node, wire, and scalar type
/// through the [`DiagramData`] trait.
#[derive(Clone, Debug)]
pub struct Diagram<A>
where A: DiagramData
{
    pub(crate) nodes: Vec<Option<A::Node>>,
    pub(crate) node_count: usize,
    pub(crate) wires: Vec<Option<Vec<A::Wire>>>,
    pub(crate) wire_count: usize,
    pub(crate) inputs: Vec<IONodeId>,
    pub(crate) outputs: Vec<IONodeId>,
    pub(crate) free: Vec<NodeId>,
    pub(crate) scalar: A::Scalar,
}

/// A diagram in the ZX-calculus.
///
/// See [`zx`][crate::graph::zx] for more information.
pub type DiagramZX = Diagram<ZX>;

/// A diagram in the ZH-calculus.
///
/// See [`zh`][crate::graph::zh] for more information.
pub type DiagramZH = Diagram<ZH>;

/// A diagram in the Clifford+*T*-restricted ZX-calculus.
///
/// See [`clifft`][crate::graph::clifft] for more information.
pub type DiagramCT = Diagram<CT>;

impl<A> Default for Diagram<A>
where A: DiagramData
{
    fn default() -> Self { Self::new() }
}

impl<A> Diagram<A>
where A: DiagramData
{
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
            scalar: A::Scalar::one(),
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
    where I: IntoIterator<Item = A::Node>
    {
        let mut dg = Self::new();
        nodes.into_iter()
            .for_each(|def| { dg.add_node(def); });
        dg
    }

    /// Return the global scalar on the diagram.
    pub fn scalar(&self) -> A::Scalar { self.scalar }

    /// Apply a mapping function to the global scalar.
    pub fn map_scalar<F>(&mut self, map: F)
    where F: FnOnce(A::Scalar) -> A::Scalar
    {
        self.scalar = map(self.scalar);
    }

    /// Return the number of nodes.
    pub fn count_nodes(&self) -> usize { self.node_count }

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

    /// Return the number of wires.
    pub fn count_wires(&self) -> usize { self.wire_count }

    /// Get the node associated with a particular ID if it exists.
    pub fn get_node(&self, id: NodeId) -> Option<&A::Node> {
        self.nodes.get(id).and_then(|mb_n| mb_n.as_ref())
    }

    pub(crate) fn get_node_mut(&mut self, id: NodeId) -> Option<&mut A::Node> {
        self.nodes.get_mut(id).and_then(|mb_n| mb_n.as_mut())
    }

    /// Return `true` if a node exists with the given ID.
    pub fn has_node(&self, id: NodeId) -> bool {
        self.nodes.get(id).is_some_and(|mb_n| mb_n.is_some())
    }

    /// Return the qubit index corresponding to a node ID, if the node exists
    /// and is an input.
    pub fn get_input_index(&self, id: NodeId) -> Option<QubitId> {
        self.inputs.iter().enumerate()
            .find_map(|(qid, ioid)| ioid.has_id(id).then_some(qid))
    }

    /// Get the node ID associated with the `q`-th input qubit index.
    pub fn get_input_id(&self, q: QubitId) -> Option<IONodeId> {
        self.inputs.get(q).copied()
    }

    /// Return the qubit index corresponding to a node ID, if the node exists
    /// and is an output.
    pub fn get_output_index(&self, id: NodeId) -> Option<QubitId> {
        self.outputs.iter().enumerate()
            .find_map(|(qid, ioid)| ioid.has_id(id).then_some(qid))
    }

    /// Get the node ID associated with the `q`-th output qubit index.
    pub fn get_output_id(&self, q: QubitId) -> Option<IONodeId> {
        self.outputs.get(q).copied()
    }

    // get the raw neighbor IDs of a node, if it exists
    pub(crate) fn get_neighbors_mut(&mut self, id: NodeId)
        -> Option<&mut Vec<A::Wire>>
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
    pub fn add_node(&mut self, node: A::Node) -> NodeId {
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

    /// Remove the node associated with a particular ID and return its data if
    /// it exists.
    ///
    /// This method also removes all wires with an endpoint at the node.
    pub fn remove_node(&mut self, id: NodeId) -> GraphResult<A::Node> {
        self.remove_node_nb(id).map(|(n, _)| n)
    }

    /// Like [`remove_node`][Self::remove_node], but returning a list of the
    /// node's neighbors alongside its data.
    pub fn remove_node_nb(&mut self, id: NodeId)
        -> GraphResult<(A::Node, Vec<A::Wire>)>
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
            let k = self.get_input_index(id).unwrap();
            self.inputs.remove(k);
        }
        if node.is_output() || self.outputs.iter().any(|ioid| ioid.has_id(id)) {
            let k = self.get_output_index(id).unwrap();
            self.outputs.remove(k);
        }
        Ok((node, nnb_of))
    }

    // remove the node associated with a particular ID and return its data and
    // neighbor indices, but *don't* remove the wires connected to it
    //
    // input/output book-keeping is still performed, however
    pub(crate) fn delete_node(&mut self, id: NodeId)
        -> Option<(A::Node, Vec<A::Wire>)>
    {
        let node = self.nodes.get_mut(id)?.take()?;
        let nnb = self.wires.get_mut(id)?.take()?;
        if node.is_input() {
            let k = self.get_input_index(id).unwrap();
            self.inputs.remove(k);
        }
        if node.is_output() {
            let k = self.get_output_index(id).unwrap();
            self.outputs.remove(k);
        }
        self.free.push(id);
        self.node_count -= 1;
        Some((node, nnb))
    }

    /// Add an empty wire between two nodes.
    ///
    /// Fails if one or neither of the nodes exist.
    pub fn add_wire(&mut self, a: NodeId, b: NodeId) -> GraphResult<()> {
        self.get_node(a)
            .ok_or(AddWireMissingNode(a))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(a).unwrap())
                    .then_some(())
                    .ok_or(AddWireConnectedIO(a))
            })?;
        self.get_node(b)
            .ok_or(AddWireMissingNode(b))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(b).unwrap())
                    .then_some(())
                    .ok_or(AddWireConnectedIO(b))
            })?;
        self.wires[a].as_mut().unwrap().push(A::Wire::new_empty(b));
        self.wires[b].as_mut().unwrap().push(A::Wire::new_empty(a));
        self.wire_count += 1;
        Ok(())
    }

    /// Add an empty wire with an attached [input][NodeData::new_input] to a
    /// pre-existing node and return the new input node's ID.
    ///
    /// The pre-existing node cannot already by an input or output.
    pub fn add_input_wire(&mut self, id: NodeId) -> GraphResult<NodeId> {
        self.get_node(id)
            .ok_or(AddWireMissingNode(id))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(AddWireConnectedIO(id))
            })?;
        let input_id = self.add_node(A::Node::new_input());
        self.add_wire(id, input_id)?;
        Ok(input_id)
    }

    /// Add a wire with an attached [output][NodeData::new_output] to a
    /// pre-existing node and return the new output node's ID.
    ///
    /// The pre-existing node cannot already be an input or output.
    pub fn add_output_wire(&mut self, id: NodeId) -> GraphResult<NodeId>
    {
        self.get_node(id)
            .ok_or(AddWireMissingNode(id))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(AddWireConnectedIO(id))
            })?;
        let output_id = self.add_node(A::Node::new_output());
        self.add_wire(id, output_id)?;
        Ok(output_id)
    }

    /// Remove all or at most a fixed number of wires between two nodes.
    ///
    /// Returns the number of wires removed.
    ///
    /// Fails if either node does not exist.
    pub fn remove_wires(&mut self, a: NodeId, b: NodeId, nwires: Option<usize>)
        -> GraphResult<usize>
    {
        self.has_node(a).then_some(()).ok_or(RemoveWireMissingNode(a))?;
        self.has_node(b).then_some(()).ok_or(RemoveWireMissingNode(b))?;

        let mut to_remove: Vec<usize> = Vec::new();

        let nnb_a = self.get_neighbors_mut(a).unwrap();
        let len_nnb_a = nnb_a.len();
        nnb_a.iter().enumerate()
            .filter_map(|(k, nid)| nid.has_id(b).then_some(k))
            .take(nwires.unwrap_or(len_nnb_a))
            .for_each(|k| { to_remove.push(k); });
        let mut removed_a = to_remove.len();
        to_remove.drain(..).rev()
            .for_each(|k| { nnb_a.swap_remove(k); });

        let nnb_b = self.get_neighbors_mut(b).unwrap();
        let len_nnb_b = nnb_b.len();
        nnb_b.iter().enumerate()
            .filter_map(|(k, nid)| nid.has_id(a).then_some(k))
            .take(nwires.unwrap_or(len_nnb_b))
            .for_each(|k| { to_remove.push(k); });
        let removed_b = to_remove.len();
        to_remove.drain(..).rev()
            .for_each(|k| { nnb_b.swap_remove(k); });

        debug_assert!(removed_a == removed_b || a == b);
        if a == b { removed_a /= 2; }
        self.wire_count -= removed_a;
        Ok(removed_a)
    }

    /// Remove tracking for the `q`-th input qubit index as such, returning its
    /// associated node ID if it has a state applied.
    ///
    /// If the input is [free][IONodeId::Free], no operation is performed and
    /// `None` is returned.
    pub fn forget_input(&mut self, q: QubitId) -> Option<NodeId> {
        if self.inputs.get(q).is_some_and(|ioid| ioid.is_state()) {
            Some(self.inputs.remove(q).id())
        } else {
            None
        }
    }

    /// Remove tracking for the `q`-th output qubit index as such, returning its
    /// associated node ID if it has an effect applied.
    ///
    /// If the output is [free][IONodeId::Free], no operation is performed and
    /// `None` is returned.
    pub fn forget_output(&mut self, q: QubitId) -> Option<NodeId> {
        if self.outputs.get(q).is_some_and(|ioid| ioid.is_state()) {
            Some(self.outputs.remove(q).id())
        } else {
            None
        }
    }

    /// Return an iterator over all nodes, visited in index order.
    ///
    /// The iterator item type is `(`[`NodeId`]`, &A::Node)`.
    pub fn nodes(&self) -> Nodes<'_, A::Node> {
        Nodes { len: self.node_count, iter: self.nodes.iter().enumerate() }
    }

    /// Return an iterator over all interior (non-input/output) nodes in the
    /// diagram, visited in index order.
    ///
    /// The iterator item type is `(`[`NodeId`]`, &A::Node)`.
    pub fn nodes_inner(&self) -> NodesInner<'_, A::Node> {
        NodesInner { iter: self.nodes.iter().enumerate() }
    }

    /// Return an iterator over all diagram input node IDs, visited in qubit
    /// index order.
    ///
    /// The iterator item type is `(`[`QubitId`]`, &`[`IONodeId`]`)`.
    pub fn inputs(&self) -> Inputs<'_> {
        Inputs { iter: self.inputs.iter().enumerate() }
    }

    /// Return an iterator over all diagram output node IDs, visited in qubit
    /// index order.
    ///
    /// The iterator item type is `(`[`QubitId`]`, &`[`IONodeId`]`)`.
    pub fn outputs(&self) -> Outputs<'_> {
        Outputs { iter: self.outputs.iter().enumerate() }
    }

    /// Return an iterator over all wires in the diagram, visited in mostly
    /// arbitrary order.
    ///
    /// Specifically, the iterator steps through wires as pairs of node IDs such
    /// that those whose left IDs are smaller in value are visited before those
    /// larger in value. In other words, the left ID monotonically increases
    /// over the course of iteration while right IDs are free to vary in
    /// arbitrary order. Each wire is only visited once and any data attached to
    /// the wire itself is stored in the right ID.
    ///
    /// The iterator item type is `(`[`NodeId`]`, &A::Wire)`.
    pub fn wires(&self) -> Wires<'_, A::Wire> {
        Wires {
            nb_group: None,
            len: self.wire_count,
            seen: Vec::new(),
            iter: self.wires.iter().enumerate()
        }
    }

    /// Return an iterator over all wires in the diagram, visited in mostly
    /// arbitrary order (see [`wires`][Self::wires]), with node data attached.
    ///
    /// the iterator item type is `((`[`NodeId`], &A::Node), (&A::Wire,
    /// &A::Node))`.
    pub fn wires_data(&self) -> WiresData<'_, A, A::Node, A::Wire> {
        WiresData { dg: self, iter: self.wires() }
    }

    /// Return an iterator over all wires connected to only interior
    /// (non-input/output) nodes in the diagram, visited in mostly arbitrary
    /// order (see [`wires`][Self::wires]).
    ///
    /// The iterator item type is `(`[`NodeId`], &A::Wire)`.
    pub fn wires_inner(&self) -> WiresInner<'_, A, A::Node, A::Wire> {
        WiresInner { dg: self, iter: self.wires() }
    }

    /// Return an iterator over all wires connected to only interior
    /// (non-input/output) nodes in the diagram, visited in mostly arbitrary
    /// order (see [`wires`][Self::wires]), with node data attached.
    ///
    /// The iterator item type is `((`[`NodeId`], &A::Node), (&A::Wire,
    /// &A::Node))`.
    pub fn wires_data_inner(&self) -> WiresDataInner<'_, A, A::Node, A::Wire> {
        WiresDataInner { iter: self.wires_data() }
    }

    /// Return the ID of an arbitrary neighbor of the given node ID, if it
    /// exists.
    pub fn get_neighbor(&self, id: NodeId) -> Option<&A::Wire> {
        self.wires.get(id)
            .and_then(|mb_nnb| mb_nnb.as_ref())
            .and_then(|nnb| nnb.first())
    }

    /// Return an iterator over the IDs of neighbors of the given node ID,
    /// visited in arbitrary order, if they exist.
    ///
    /// Note that this iterator will contain duplicate elements if there are
    /// multiple wires going to the same neighbor; if the node is connected to
    /// itself, each such connection will be counted only once.
    ///
    /// The iterator item type is [`NodeId`].
    pub fn neighbor_ids(&self, id: NodeId)
        -> Option<NeighborIds<'_, A::Wire>>
    {
        self.wires.get(id)
            .and_then(|mb_nnb| {
                mb_nnb.as_ref()
                    .map(|nnb| {
                        NeighborIds {
                            id,
                            switch: true,
                            iter: nnb.iter()
                        }
                    })
            })
    }

    /// Find the ID of the first neighbor of the given node ID that satisfies
    /// some predicate, if it exists.
    pub fn find_neighbor_id<F>(&self, id: NodeId, mut pred: F)
        -> Option<&A::Wire>
    where F: FnMut(&A::Wire) -> bool
    {
        self.neighbor_ids(id)
            .and_then(|mut nnb| nnb.find(|id| pred(id)))
    }

    /// Find the ID of the first neighbor of a given node ID that satisfies some
    /// predicate, if it exists, and map to some output.
    pub fn find_map_neighbor_id<F, T>(&self, id: NodeId, pred: F)
        -> Option<T>
    where F: FnMut(&A::Wire) -> Option<T>
    {
        self.neighbor_ids(id)
            .and_then(|mut nnb| nnb.find_map(pred))
    }

    /// Return an iterator over neighbors of the given node ID, visited in
    /// arbitrary order, if it exists.
    ///
    /// Note that this iterator will contain duplicate elements if there are
    /// multiple wires going to the same neighbor; if the node is connected to
    /// itself, each such connection will be counted only once.
    ///
    /// The iterator item type is `(&A::Wire, &A::Node)`.
    pub fn neighbors(&self, id: NodeId)
        -> Option<Neighbors<'_, A, A::Node, A::Wire>>
    {
        self.wires.get(id)
            .and_then(|mb_nnb| {
                mb_nnb.as_ref()
                    .map(|nnb| {
                        Neighbors {
                            dg: self,
                            id,
                            switch: true,
                            iter: nnb.iter(),
                        }
                    })
            })
    }

    /// Find the first neighbor of the given node ID that satisfies some
    /// predicate, if it and the ID exist.
    pub fn find_neighbor<F>(&self, id: NodeId, mut pred: F)
        -> Option<(&A::Wire, &A::Node)>
    where F: FnMut(&A::Wire, &A::Node) -> bool
    {
        self.neighbors(id)
            .and_then(|mut nnb| nnb.find(|(id, n)| pred(id, n)))
    }

    /// Find the first neighbor of the given node ID that satisfies some
    /// predicate, if it and the ID exist, and map to some output.
    pub fn find_map_neighbor<F, T>(&self, id: NodeId, mut pred: F)
        -> Option<T>
    where F: FnMut(&A::Wire, &A::Node) -> Option<T>
    {
        self.neighbors(id)
            .and_then(|mut nnb| nnb.find_map(|(id, n)| pred(id, n)))
    }

    /// Return an iterator over all interior (non-input/output) neighbors of the
    /// given node ID, visited in arbitrary order, if it exists.
    ///
    /// Note that this iterator will contains duplicate elements if there are
    /// multiple wires going to the same neighbor; if the node is connected to
    /// itself, each such connection will be counted only once.
    ///
    /// The iterator type is `(&A::Wire, &A::Node)`.
    pub fn neighbors_inner(&self, id: NodeId)
        -> Option<NeighborsInner<'_, A, A::Node, A::Wire>>
    {
        self.neighbors(id)
            .map(|iter| NeighborsInner { iter })
    }

    /// Find the first interior (non-input/output) neighbor of the given node ID
    /// that satisfies some predicate, if it and the ID exist.
    pub fn find_neighbor_inner<F>(&self, id: NodeId, mut pred: F)
        -> Option<(&A::Wire, &A::Node)>
    where F: FnMut(&A::Wire, &A::Node) -> bool
    {
        self.neighbors_inner(id)
            .and_then(|mut nnb| nnb.find(|(id, n)| pred(id, n)))
    }

    /// Find the first interior (non-input/output) neighbor of the given node ID
    /// that satisfies some predicate, if it and the ID exist, and map to some
    /// output.
    pub fn find_map_neighbor_inner<F, T>(&self, id: NodeId, mut pred: F)
        -> Option<T>
    where F: FnMut(&A::Wire, &A::Node) -> Option<T>
    {
        self.neighbors_inner(id)
            .and_then(|mut nnb| nnb.find_map(|(id, n)| pred(id, n)))
    }

    /// Replace two existing diagram inputs or outputs with a Bell state/effect.
    ///
    /// The inputs or outputs will be replaced with phaseless spiders (for
    /// book-keeping purposes) with the same IDs.
    ///
    /// Fails if the given node IDs do not exist, are not both inputs or
    /// outputs, or refer to the same qubit.
    pub fn apply_bell(&mut self, a: NodeId, b: NodeId) -> GraphResult<()> {
        if a == b { return Err(ApplyBellSameQubit); }
        self.get_node(a).zip(self.get_node(b))
            .is_some_and(|(na, nb)| {
                (na.is_input() && nb.is_input())
                    || (na.is_output() && nb.is_output())
            })
            .then_some(())
            .ok_or(ApplyBellNotIO(a, b))?;
        if self.nodes[a].as_ref().unwrap().is_input() {
            self.inputs.iter_mut()
                .find(|ioid| ioid.has_id(a))
                .unwrap()
                .make_state();
        } else {
            self.outputs.iter_mut()
                .find(|ioid| ioid.has_id(a))
                .unwrap()
                .make_state();
        }
        self.nodes[a] = Some(A::Node::new_id());
        if self.nodes[b].as_ref().unwrap().is_input() {
            self.inputs.iter_mut()
                .find(|ioid| ioid.has_id(b))
                .unwrap()
                .make_state();
        } else {
            self.outputs.iter_mut()
                .find(|ioid| ioid.has_id(b))
                .unwrap()
                .make_state();
        }
        self.nodes[b] = Some(A::Node::new_id());
        self.add_wire(a, b)?;
        Ok(())
    }

    /// Like [`apply_bell`][Self::apply_bell], but using input qubit indices
    /// rather than bare node IDs.
    pub fn apply_bell_input(&mut self, qa: QubitId, qb: QubitId)
        -> GraphResult<()>
    {
        self.get_input_id(qa)
            .ok_or(ApplyBellMissingQubit(qa))
            .and_then(|nida| {
                self.get_input_id(qb)
                    .ok_or(ApplyBellMissingQubit(qb))
                    .map(|nidb| (nida, nidb))
            })
            .and_then(|(nida, nidb)| self.apply_bell(nida.id(), nidb.id()))
    }

    /// Like [`apply_bell`][Self::apply_bell], but using output qubit indices
    /// rather than bare node IDs.
    pub fn apply_bell_output(&mut self, qa: QubitId, qb: QubitId)
        -> GraphResult<()>
    {
        self.get_output_id(qa)
            .ok_or(ApplyBellMissingQubit(qa))
            .and_then(|nida| {
                self.get_output_id(qb)
                    .ok_or(ApplyBellMissingQubit(qb))
                    .map(|nidb| (nida, nidb))
            })
            .and_then(|(nida, nidb)| self.apply_bell(nida.id(), nidb.id()))
    }

    /// Swap inputs and outputs and conjugate all nodes as well as the global
    /// scalar, consuming `self`.
    pub fn adjoint(mut self) -> Self {
        self.adjoint_mut();
        self
    }

    /// Swap inputs and outputs and conjugate all nodes as well as the global
    /// scalar, modifying `self` in place.
    pub fn adjoint_mut(&mut self) {
        self.nodes.iter_mut()
            .flatten()
            .for_each(|node| { node.conjugate(); });
        std::mem::swap(&mut self.inputs, &mut self.outputs);
        self.scalar = self.scalar.conj();
    }

    pub(crate) fn append(&mut self, other: Self) -> usize {
        let Self {
            mut nodes,
            node_count,
            mut wires,
            wire_count,
            mut inputs,
            mut outputs,
            mut free,
            scalar,
        } = other;
        let shift = self.nodes.len();
        wires.iter_mut()
            .flatten()
            .flat_map(|nnb| nnb.iter_mut())
            .for_each(|nb| { nb.shift_id(shift); });
        inputs.iter_mut()
            .for_each(|ioid| { ioid.shift_id(shift); });
        outputs.iter_mut()
            .for_each(|ioid| { ioid.shift_id(shift); });
        free.iter_mut()
            .for_each(|nid| { *nid += shift; });
        self.nodes.append(&mut nodes);
        self.node_count += node_count;
        self.wires.append(&mut wires);
        self.wire_count += wire_count;
        self.inputs.append(&mut inputs);
        self.outputs.append(&mut outputs);
        self.free.append(&mut free);
        self.scalar *= scalar;
        shift
    }

    /// Return the tensor product `self ⊗ rhs`, consuming both.
    ///
    /// All node IDs from `rhs` will be adjusted to avoid collision. The exact
    /// size of the shift depends on how much space has been allocated for nodes
    /// in `self` (both current and previous), and is returned alongside the
    /// output diagram.
    pub fn tensor(mut self, rhs: Self) -> (Self, usize) {
        let shift = self.tensor_with(rhs);
        (self, shift)
    }

    /// Compute the tensor product of `self ⊗ rhs`, consuming `rhs` and
    /// modifying `self` in place.
    ///
    /// All node IDs from `rhs` will be adjusted to avoid collision. The exact
    /// size of the shift depends on how much space has been allocated for nodes
    /// in `self` (both current and previous), and will be returned.
    pub fn tensor_with(&mut self, rhs: Self) -> usize {
        self.append(rhs)
    }

    /// Return the composition `self ∘ rhs`, attempting to match the outputs of
    /// `rhs` to the inputs of `self` in qubit index order, consuming both
    /// `self` and `rhs`.
    ///
    /// The IDs of all nodes from `rhs` will be adjusted to avoid collision. The
    /// exact size of the shift depends on how much space has been allocated for
    /// nodes in `self` (both current and previous), and will be returned
    /// alongside the output diagram.
    ///
    /// This operation will fail if the numbers of inputs and outputs are not
    /// equal, or if they are equal but some inputs or outputs are not free
    /// wires. `self` will not be modified in this case.
    ///
    /// See also [`compose_rev`][Self::compose_rev].
    pub fn compose(mut self, rhs: Self) -> GraphResult<(Self, usize)> {
        let shift = self.compose_with(rhs)?;
        Ok((self, shift))
    }

    /// Compute the composition `self ∘ rhs`, attempting to match the outputs of
    /// `rhs` to the inputs of `self` in qubit index order, consuming `rhs` and
    /// modifying `self` in place.
    ///
    /// The IDs of all nodes from `rhs` will be adjusted to avoid collision. The
    /// exact size of the shift depends on how much space has been allocated for
    /// nodes in `self` (both current and previous), and will be returned.
    ///
    /// This operation will fail if the numbers of inputs and outputs are not
    /// equal, or if they are equal but some inputs or outputs are not free
    /// wires. `self` will not be modified in this case.
    ///
    /// See also [`compose_with_rev`][Self::compose_with_rev].
    pub fn compose_with(&mut self, rhs: Self) -> GraphResult<usize> {
        // find (input/output, interior neighbor)
        let self_inputs: Vec<(NodeId, Option<A::Wire>)> =
            self.inputs()
            .filter_map(|(_, ioid)| ioid.free_id())
            .map(|nid| {
                let mb_interior =
                    self.neighbor_ids(nid)
                    .and_then(|mut nnb| nnb.next().cloned());
                (nid, mb_interior)
            })
            .collect();
        let rhs_outputs: Vec<(NodeId, Option<A::Wire>)> =
            rhs.outputs()
            .filter_map(|(_, ioid)| ioid.free_id())
            .map(|nid| {
                let mb_interior =
                    rhs.neighbor_ids(nid)
                    .and_then(|mut nnb| nnb.next().cloned());
                (nid, mb_interior)
            })
            .collect();
        (self_inputs.len() == rhs_outputs.len()).then_some(())
            .ok_or(NonMatchingIO(rhs_outputs.len(), self_inputs.len()))?;
        let shift = self.append(rhs);
        let iter = rhs_outputs.into_iter().zip(self_inputs);
        for ((out_id, mb_out_int), (in_id, mb_in_int)) in iter {
            match (mb_out_int, mb_in_int) {
                (Some(out_int), Some(in_int)) => {
                    if out_int.is_empty() && in_int.is_empty() {
                        self.remove_node(out_id + shift).unwrap();
                        self.remove_node(in_id).unwrap();
                        self.add_wire(out_int.id() + shift, in_int.id())
                            .unwrap();
                    } else {
                        let n = self.get_node_mut(in_id).unwrap();
                        *n = A::Node::new_id();
                        self.delete_node(out_id + shift).unwrap();
                        self.get_neighbors_mut(out_int.id() + shift).unwrap()
                            .iter_mut()
                            .for_each(|out_int_nb| {
                                if out_int_nb.has_id(out_id) {
                                    out_int_nb.map_id(|_| in_id);
                                }
                            });
                    }
                },
                (Some(out_int), None) => {
                    self.remove_node(out_id + shift).unwrap();
                    self.remove_node(in_id).unwrap();
                    self.add_output_wire(out_int.id() + shift).unwrap();
                },
                (None, Some(in_int)) => {
                    self.remove_node(out_id + shift).unwrap();
                    self.remove_node(in_id).unwrap();
                    self.add_input_wire(in_int.id()).unwrap();
                },
                (None, None) => {
                    self.remove_node(out_id + shift).unwrap();
                    self.remove_node(in_id).unwrap();
                    let new_in = self.add_node(A::Node::new_input());
                    let new_out = self.add_node(A::Node::new_output());
                    self.add_wire(new_in, new_out).unwrap();
                },
            }
        }
        Ok(shift)
    }

    /// Return the composition `rhs ∘ self`, attempting to match the outputs of
    /// `self` to the inputs of `rhs` in qubit index order, consuming both
    /// `self` and `rhs`.
    ///
    /// The IDs of all nodes from `rhs` will be adjusted to avoid collision. The
    /// exact size of the shift depends on how much space has been allocated for
    /// nodes in `self` (both current and previous), and will be returned
    /// alongside the output diagram.
    ///
    /// This operation will fail if the numbers of inputs and outputs are not
    /// equal, or if they are equal but some inputs or outputs are not free
    /// wires. `self` will not be modified in this case.
    ///
    /// See also [`compose`][Self::compose].
    pub fn compose_rev(mut self, rhs: Self)
        -> GraphResult<(Self, usize)>
    {
        let shift = self.compose_with_rev(rhs)?;
        Ok((self, shift))
    }

    /// Compute the composition `rhs ∘ self`, attempting to match the outputs of
    /// `self` to the inputs of `rhs` in qubit index order, consuming `rhs` and
    /// modifying `self` in place.
    ///
    /// The IDs of all nodes from `rhs` will be adjusted to avoid collision. The
    /// exact size of the shift depends on how much space has been allocated for
    /// nodes in `self` (both current and previous), and will be returned.
    ///
    /// This operation will fail if the numbers of inputs and outputs are not
    /// equal, or if they are equal but some inputs or outputs are not free
    /// wires. `self` will not be modified in this case.
    ///
    /// See also [`compose_with`][Self::compose_with].
    pub fn compose_with_rev(&mut self, rhs: Self) -> GraphResult<usize> {
        // (input/output, interior neighbor)
        let self_outputs: Vec<(NodeId, Option<A::Wire>)> =
            self.outputs()
            .filter_map(|(_, ioid)| ioid.free_id())
            .map(|nid| {
                let mb_interior =
                    self.neighbor_ids(nid)
                    .and_then(|mut nnb| nnb.next().cloned());
                (nid, mb_interior)
            })
            .collect();
        let rhs_inputs: Vec<(NodeId, Option<A::Wire>)> =
            rhs.inputs()
            .filter_map(|(_, ioid)| ioid.free_id())
            .map(|nid| {
                let mb_interior =
                    rhs.neighbor_ids(nid)
                    .and_then(|mut nnb| nnb.next().cloned());
                (nid, mb_interior)
            })
            .collect();
        (self_outputs.len() == rhs_inputs.len()).then_some(())
            .ok_or(NonMatchingIO(self_outputs.len(), rhs_inputs.len()))?;
        let shift = self.append(rhs);
        let it = self_outputs.into_iter().zip(rhs_inputs);
        for ((out_id, mb_out_int), (in_id, mb_in_int)) in it {
            match (mb_out_int, mb_in_int) {
                (Some(out_int), Some(in_int)) => {
                    if out_int.is_empty() && in_int.is_empty() {
                        self.remove_node(out_id).unwrap();
                        self.remove_node(in_id + shift).unwrap();
                        self.add_wire(out_int.id(), in_int.id() + shift)
                            .unwrap();
                    } else {
                        let n = self.get_node_mut(out_id).unwrap();
                        *n = A::Node::new_id();
                        self.delete_node(in_id + shift).unwrap();
                        self.get_neighbors_mut(in_int.id() + shift).unwrap()
                            .iter_mut()
                            .for_each(|in_int_nb| {
                                if in_int_nb.has_id(in_id) {
                                    in_int_nb.map_id(|_| out_id);
                                }
                            });
                    }
                },
                (Some(out_int), None) => {
                    self.remove_node(out_id).unwrap();
                    self.remove_node(in_id + shift).unwrap();
                    self.add_output_wire(out_int.id()).unwrap();
                },
                (None, Some(in_int)) => {
                    self.remove_node(out_id).unwrap();
                    self.remove_node(in_id + shift).unwrap();
                    self.add_input_wire(in_int.id() + shift).unwrap();
                },
                (None, None) => {
                    self.remove_node(out_id).unwrap();
                    self.remove_node(in_id + shift).unwrap();
                    let new_in = self.add_node(A::Node::new_input());
                    let new_out = self.add_node(A::Node::new_output());
                    self.add_wire(new_in, new_out).unwrap();
                },
            }
        }
        Ok(shift)
    }

}

// helper for tensor conversions -- this acts as a vendor for wire IDs so that,
// when constructing a tensor diagram, each wire can be guaranteed to have a
// unique ID
#[derive(Clone, Debug)]
pub(crate) struct WireStore(FxHashMap<(NodeId, NodeId), Vec<usize>>);

impl FromIterator<(NodeId, NodeId)> for WireStore {
    fn from_iter<I>(iter: I) -> Self
    where I: IntoIterator<Item = (NodeId, NodeId)>
    {
        // assume each wire is visited only once in the input iterator
        //
        // self-loops should be caught by the relevant Diagram methods, so we'll
        // skip those
        //
        // reflexivity in keys should be okay, but the same ID cannot be stored
        // on multiple stacks
        let mut store: FxHashMap<(NodeId, NodeId), Vec<usize>> =
            FxHashMap::default();
        for (id, (a, b)) in iter.into_iter().enumerate() {
            if a == b { continue; }
            if let Some(ids) = store.get_mut(&(a, b)) {
                ids.push(id);
            } else if let Some(ids) = store.get_mut(&(b, a)) {
                ids.push(id);
            } else {
                store.insert((a, b), vec![id]);
            }
        }
        Self(store)
    }
}

impl WireStore {
    // get a list of wire IDs for a (left, right) node ID pair
    pub(crate) fn get(&self, a: NodeId, b: NodeId) -> Option<&Vec<usize>> {
        self.0.get(&(a, b)).or_else(|| self.0.get(&(b, a)))
    }
}

/// Iterator over all nodes in a diagram, visited in index order.
///
/// The iterator item type is `(`[`NodeId`]`, &N)`.
#[derive(Clone, Debug)]
pub struct Nodes<'a, N> {
    len: usize,
    iter: std::iter::Enumerate<std::slice::Iter<'a, Option<N>>>
}

impl<'a, N> Iterator for Nodes<'a, N> {
    type Item = (NodeId, &'a N);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find_map(|(id, mb_n)| {
            mb_n.as_ref()
                .map(|n| {
                    self.len = self.len.saturating_sub(1);
                    (id, n)
                })
        })
    }
}

impl<'a, N> DoubleEndedIterator for Nodes<'a, N> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|(_, mb_n)| mb_n.is_some())
            .map(|(id, mb_n)| {
                self.len = self.len.saturating_sub(1);
                (id, mb_n.as_ref().unwrap())
            })
    }
}

impl<'a, N> ExactSizeIterator for Nodes<'a, N> {
    fn len(&self) -> usize { self.len }
}

impl<'a, N> std::iter::FusedIterator for Nodes<'a, N> { }

/// Iterator over all interior (non-input/output) nodes in a diagram, visited in
/// index order.
///
/// The iterator item type is `(`[`NodeId`]`, &N)`.
pub struct NodesInner<'a, N> {
    iter: std::iter::Enumerate<std::slice::Iter<'a, Option<N>>>
}

impl<'a, N> Iterator for NodesInner<'a, N>
where N: NodeData
{
    type Item = (NodeId, &'a N);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find_map(|(id, mb_n)| {
            mb_n.as_ref()
                .and_then(|n| n.is_generator().then_some((id, n)))
        })
    }
}

impl<'a, N> DoubleEndedIterator for NodesInner<'a, N>
where N: NodeData
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|(_, mb_n)| {
            mb_n.as_ref().is_some_and(|n| n.is_generator())
        })
        .map(|(id, mb_n)| (id, mb_n.as_ref().unwrap()))
    }
}

impl<'a, N> std::iter::FusedIterator for NodesInner<'a, N>
where N: NodeData
{ }

/// Iterator over all diagram input node IDs, visited in qubit index order.
///
/// The iterator item type is `(`[`QubitId`]`, &`[`IONodeId`]`)`.
pub struct Inputs<'a> {
    iter: std::iter::Enumerate<std::slice::Iter<'a, IONodeId>>
}

impl<'a> Iterator for Inputs<'a> {
    type Item = (QubitId, &'a IONodeId);

    fn next(&mut self) -> Option<Self::Item> { self.iter.next() }
}

impl<'a> DoubleEndedIterator for Inputs<'a> {
    fn next_back(&mut self) -> Option<Self::Item> { self.iter.next_back() }
}

impl<'a> ExactSizeIterator for Inputs<'a> {
    fn len(&self) -> usize { self.iter.len() }
}

impl<'a> std::iter::FusedIterator for Inputs<'a> { }

/// Iterator over all diagram output node IDs, visited in qubit index order.
///
/// The iterator item type is `(`[`QubitId`]`, &`[`IONodeId`]`)`.
pub struct Outputs<'a> {
    iter: std::iter::Enumerate<std::slice::Iter<'a, IONodeId>>
}

impl<'a> Iterator for Outputs<'a> {
    type Item = (QubitId, &'a IONodeId);

    fn next(&mut self) -> Option<Self::Item> { self.iter.next() }
}

impl<'a> DoubleEndedIterator for Outputs<'a> {
    fn next_back(&mut self) -> Option<Self::Item> { self.iter.next_back() }
}

impl<'a> ExactSizeIterator for Outputs<'a> {
    fn len(&self) -> usize { self.iter.len() }
}

impl<'a> std::iter::FusedIterator for Outputs<'a> { }

/// Iterator over all wires in a diagram, visited in mostly arbitrary order.
///
/// Specifically, the iterator steps through wires as pairs of node IDs such
/// that those whose left IDs are smaller in value are visited before those
/// larger in value. In other words, the left ID monotonically increases over
/// the course of iteration while right IDs are free to vary in arbitrary order.
/// Each wire is only visited once.
///
/// The iterator item type is `(`[`NodeId`]`, &W)`.
pub struct Wires<'a, W> {
    nb_group: Option<(NodeId, std::slice::Iter<'a, W>)>,
    len: usize,
    seen: Vec<(NodeId, NodeId)>,
    iter: std::iter::Enumerate<std::slice::Iter<'a, Option<Vec<W>>>>,
}

impl<'a, W> Iterator for Wires<'a, W>
where W: WireData
{
    type Item = (NodeId, &'a W);

    fn next(&mut self) -> Option<Self::Item> {
        // this is basically `flatten`, but with some extra processing to
        // hold left node IDs and account for right IDs being held in an Option
        //
        // this would be perfectly doable with a `flat_map`, but I wanted a
        // custom iterator type and didn't want to have it be generic over a
        // closure
        loop {
            if let Some((left, group_iter)) = self.nb_group.as_mut() {
                if let Some(right) = group_iter.next() {
                    let pair = (*left, right);
                    let pair_rev = (pair.1.id(), pair.0);
                    if self.seen.contains(&(*left, right.id())) {
                        continue;
                    } else {
                        self.seen.push(pair_rev);
                        self.len = self.len.saturating_sub(1);
                        return Some(pair);
                    }
                } else if let Some((new_left, mb_new_nnb)) = self.iter.next() {
                    self.nb_group =
                        mb_new_nnb.as_ref()
                        .map(|new_nnb| (new_left, new_nnb.iter()));
                    continue;
                } else {
                    return None;
                }
            } else if let Some((new_left, mb_new_nnb)) = self.iter.next() {
                self.nb_group =
                    mb_new_nnb.as_ref()
                    .map(|new_nnb| (new_left, new_nnb.iter()));
                continue;
            } else {
                return None;
            }
        }
    }
}

impl<'a, W> DoubleEndedIterator for Wires<'a, W>
where W: WireData
{
    fn next_back(&mut self) -> Option<Self::Item> {
        loop {
            if let Some((left, group_iter)) = self.nb_group.as_mut() {
                if let Some(right) = group_iter.next_back() {
                    let pair = (*left, right);
                    let pair_rev = (pair.1.id(), pair.0);
                    if self.seen.contains(&(*left, right.id())) {
                        continue;
                    } else {
                        self.seen.push(pair_rev);
                        self.len = self.len.saturating_sub(1);
                        return Some(pair);
                    }
                } else if let Some((new_left, mb_new_nnb)) =
                    self.iter.next_back()
                {
                    self.nb_group =
                        mb_new_nnb.as_ref()
                        .map(|new_nnb| (new_left, new_nnb.iter()));
                    continue;
                } else {
                    return None;
                }
            } else if let Some((new_left, mb_new_nnb)) = self.iter.next_back() {
                self.nb_group =
                    mb_new_nnb.as_ref()
                    .map(|new_nnb| (new_left, new_nnb.iter()));
                continue;
            } else {
                return None;
            }
        }
    }
}

impl<'a, W> ExactSizeIterator for Wires<'a, W>
where W: WireData
{
    fn len(&self) -> usize { self.len }
}

impl<'a, W> std::iter::FusedIterator for Wires<'a, W>
where W: WireData
{ }

/// Iterator over all wires in a diagram, visited in mostly arbitrary order,
/// with node data attached.
///
/// Specifically, the iterator steps through wires as pairs of node IDs such
/// that those whose left IDs are smaller in value are visited before those
/// larger in value. In other words, the left ID monotonically increases over
/// the course of iteration while right IDs are free to vary in arbitrary order.
/// Each wire is only visited once.
///
/// The iterator item type is `((`[`NodeId`]`, &N), (&W, &N))`.
pub struct WiresData<'a, A, N, W>
where A: DiagramData<Node = N, Wire = W>
{
    dg: &'a Diagram<A>,
    iter: Wires<'a, W>,
}

impl<'a, A, N: 'a, W> Iterator for WiresData<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    W: WireData,
{
    type Item = ((NodeId, &'a N), (&'a W, &'a N));

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
            .map(|(l, r)| {
                (
                    (l, self.dg.get_node(l).unwrap()),
                    (r, self.dg.get_node(r.id()).unwrap()),
                )
            })
    }
}

impl<'a, A, N: 'a, W> DoubleEndedIterator for WiresData<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    W: WireData,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.next_back()
            .map(|(l, r)| {
                (
                    (l, self.dg.get_node(l).unwrap()),
                    (r, self.dg.get_node(r.id()).unwrap()),
                )
            })
    }
}

impl<'a, A, N: 'a, W> ExactSizeIterator for WiresData<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    W: WireData,
{
    fn len(&self) -> usize { self.iter.len() }
}

impl<'a, A, N: 'a, W> std::iter::FusedIterator for WiresData<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    W: WireData,
{ }

/// Iterator over all wires connected to only interior (non-input/output) nodes
/// in a diagram, visited in mostly arbitrary order (see [`Wires`]).
///
/// The iterator item type is `(`[`NodeId`]`, &W)`.
pub struct WiresInner<'a, A, N, W>
where A: DiagramData<Node = N, Wire = W>
{
    dg: &'a Diagram<A>,
    iter: Wires<'a, W>,
}

impl<'a, A, N, W> Iterator for WiresInner<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    N: NodeData,
    W: WireData,
{
    type Item = (NodeId, &'a W);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(l, r)| {
            self.dg.get_node(*l).unwrap().is_generator()
                && self.dg.get_node(r.id()).unwrap().is_generator()
        })
    }
}

impl<'a, A, N, W> DoubleEndedIterator for WiresInner<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    N: NodeData,
    W: WireData,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|(l, r)| {
            self.dg.get_node(*l).unwrap().is_generator()
                && self.dg.get_node(r.id()).unwrap().is_generator()
        })
    }
}

impl<'a, A, N, W> std::iter::FusedIterator for WiresInner<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    N: NodeData,
    W: WireData,
{ }

/// Iterator over all wires connected to only interior (non-input/output) nodes
/// in a diagram, visited in mostly arbitrary order (see [`WiresData`]), with
/// node data attached.
///
/// The iterator item type is `((`[`NodeId`]`, &N), (&W, &N))`.
pub struct WiresDataInner<'a, A, N, W>
where A: DiagramData<Node = N, Wire = W>
{
    iter: WiresData<'a, A, N, W>
}

impl<'a, A, N: 'a, W> Iterator for WiresDataInner<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    N: NodeData,
    W: WireData,
{
    type Item = ((NodeId, &'a N), (&'a W, &'a N));

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|((_, l), (_, r))| {
            l.is_generator() && r.is_generator()
        })
    }
}

impl<'a, A, N: 'a, W> DoubleEndedIterator for WiresDataInner<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    N: NodeData,
    W: WireData,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|((_, l), (_, r))| {
            l.is_generator() && r.is_generator()
        })
    }
}

impl<'a, A, N: 'a, W> std::iter::FusedIterator for WiresDataInner<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    N: NodeData,
    W: WireData,
{ }

/// Iterator over neighbor IDs of a node, visited in arbitrary order.
///
/// Node that this iterator will contain duplicate elements if there are
/// multiple wires going to the same neighbor; if the node is connected to
/// itself, each such connection will be counted only once.
///
/// The iterator item type is `&W`.
pub struct NeighborIds<'a, W> {
    id: NodeId,
    switch: bool,
    iter: std::slice::Iter<'a, W>,
}

impl<'a, W> Iterator for NeighborIds<'a, W>
where W: WireData
{
    type Item = &'a W;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|id| {
            if id.id() == self.id {
                if self.switch {
                    self.switch = false;
                    true
                } else {
                    self.switch = true;
                    false
                }
            } else {
                true
            }
        })
    }
}

impl<'a, W> DoubleEndedIterator for NeighborIds<'a, W>
where W: WireData
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|id| {
            if id.id() == self.id {
                if self.switch {
                    self.switch = false;
                    true
                } else {
                    self.switch = true;
                    false
                }
            } else {
                true
            }
        })
    }
}

impl<'a, W> std::iter::FusedIterator for NeighborIds<'a, W>
where W: WireData
{ }

/// Iterator over neighbors of a node, visited in arbitrary order.
///
/// Note that this iterator will contain duplicate elements if there are
/// multiple wires going to the same neighbor; if the node is connected to
/// itself, each such connection will be counted only once.
///
/// The iterator item type is `(&W, &N)`.
pub struct Neighbors<'a, A, N, W>
where A: DiagramData<Node = N, Wire = W>
{
    dg: &'a Diagram<A>,
    id: NodeId,
    switch: bool,
    iter: std::slice::Iter<'a, W>,
}

impl<'a, A, N: 'a, W> Iterator for Neighbors<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    W: WireData,
{
    type Item = (&'a W, &'a N);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|id| {
            if id.id() == self.id {
                if self.switch {
                    self.switch = false;
                    true
                } else {
                    self.switch = true;
                    false
                }
            } else {
                true
            }
        })
        .map(|id| (id, self.dg.nodes[id.id()].as_ref().unwrap()))
    }
}

impl<'a, A, N: 'a, W> DoubleEndedIterator for Neighbors<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    W: WireData,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|id| {
            if id.id() == self.id {
                if self.switch {
                    self.switch = false;
                    true
                } else {
                    self.switch = true;
                    false
                }
            } else {
                true
            }
        })
        .map(|id| (id, self.dg.nodes[id.id()].as_ref().unwrap()))
    }
}

impl<'a, A, N: 'a, W> std::iter::FusedIterator for Neighbors<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    W: WireData,
{ }

/// Iterator over all interior (non-input/output) neighbors of a node, visited
/// in arbitrary order.
///
/// Note that this iterator will contain duplicate elements if there are
/// multiple wires going to the same neighbor; if the node is connected to
/// itself, each such connection will be counted only once.
///
/// The iterator item type is `(&W, &N)`.
pub struct NeighborsInner<'a, A, N, W>
where A: DiagramData<Node = N, Wire = W>
{
    iter: Neighbors<'a, A, N, W>
}

impl<'a, A, N: 'a, W> Iterator for NeighborsInner<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    N: NodeData,
    W: WireData,
{
    type Item = (&'a W, &'a N);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, n)| n.is_generator())
    }
}

impl<'a, A, N: 'a, W> DoubleEndedIterator for NeighborsInner<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    N: NodeData,
    W: WireData,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|(_, n)| n.is_generator())
    }
}

impl<'a, A, N: 'a, W> std::iter::FusedIterator for NeighborsInner<'a, A, N, W>
where
    A: DiagramData<Node = N, Wire = W>,
    N: NodeData,
    W: WireData,
{ }

/// Create a [graph-based diagram][Diagram] using an abbreviated syntax.
///
/// The first item (in angle brackets) sets the variant of the ZX-calculus in
/// use, and should be the name of a type implementing [`DiagramData`].
/// Following that, there are three distinct kinds of blocks, separated by a
/// single `+`. The first of these is headed with `nodes` and defines nodes in
/// the diagram with the syntax
/// ```text
/// <node_label> = <node_constructor> ( args... )
/// ```
/// where `<node_constructor>` is a variant or method of the node type
/// associated with the `DiagramData` implementation and `args...` is any
/// arguments to be passed to it. Note that `args...` may be left empty. Labels
/// given to nodes in this block are then used in the following blocks.
///
/// The second kind of block is headed by the name of a diagram method that
/// takes two node IDs as argument followed by zero additional arguments (e.g.
/// `add_wire` or `apply_bell`). The node IDs are identified using the syntax
/// ```text
/// <node1> -- <node2> [ -- <node3> ... ]
/// ```
/// Node labels can be chained to arbitrary depths, where the header method is
/// called for every (overlapping) consecutive pair.
///
/// The last kind of block is headed by the name of a diagram method that takes
/// a single node ID as an argument, followed by any number of additional
/// arguments (e.g. `apply_state`). The header method is called for every item
/// in the block following the syntax
/// ```text
/// <node_label> ( args... )
/// ```
///
/// All methods used as block headers must return `Result<_, `[`GraphError`]`>`.
///
/// All nodes' labels and corresponding IDs are returned in a
/// [`HashMap`][std::collections::HashMap]`<&'static `[`str`]`, `[`NodeId`]`>`.
///
/// The total return type is `Result<(`[`Diagram`]`,
/// `[`HashMap`][std::collections::HashMap]`<&'static `[`str`]`, `[`NodeId`]`>),
/// `[`GraphError`]`>`.
///
/// The normal usage
/// ```
/// # use zx_calc::graph::*;
/// use zx_calc::phase::Phase;
/// # fn main() -> Result<(), GraphError> {
/// // construct the single-qubit teleportation circuit
/// // post-selected measurement results
/// let m0 = Phase::pi();   // measure ∣1⟩ on qubit 0
/// let m1 = Phase::zero(); // measure ∣0⟩ on qubit 1
///
/// let mut diagram: Diagram<ZX> = Diagram::new();
/// let i0 = diagram.add_node(ZXNode::Input);
/// let i1 = diagram.add_node(ZXNode::Input);
/// let i2 = diagram.add_node(ZXNode::Input);
/// let cnot_z = diagram.add_node(ZXNode::z());
/// let cnot_x = diagram.add_node(ZXNode::x());
/// let xrot = diagram.add_node(ZXNode::X(m1));
/// let zrot = diagram.add_node(ZXNode::Z(m0));
/// let o0 = diagram.add_node(ZXNode::Output);
/// let o1 = diagram.add_node(ZXNode::Output);
/// let o2 = diagram.add_node(ZXNode::Output);
/// diagram.add_wire(i0, cnot_z)?;
/// diagram.add_wire_h(cnot_z, o0)?;
/// diagram.add_wire(cnot_z, cnot_x)?;
/// diagram.add_wire(i1, cnot_x)?;
/// diagram.add_wire(cnot_x, o1)?;
/// diagram.add_wire(i2, xrot)?;
/// diagram.add_wire(xrot, zrot)?;
/// diagram.add_wire(zrot, o2)?;
/// diagram.apply_state(o0, Spider::X(m0))?;
/// diagram.apply_state(o1, Spider::X(m1))?;
/// diagram.apply_bell(i1, i2)?;
/// # Ok(())
/// # }
/// ```
/// constructs the same diagram as
/// ```
/// # use zx_calc::graph::*;
/// use zx_calc::phase::Phase;
/// use zx_calc::graph_diagram;
/// // construct the single-qubit teleportation circuit
/// // post-selected measurement results
/// let m0 = Phase::pi();   // measure ∣1⟩ on qubit 0
/// let m1 = Phase::zero(); // measure ∣0⟩ on qubit 1
///
/// graph_diagram!(
///     <ZX>
///     nodes: {
///         i0 = input ( ),
///         i1 = input ( ),
///         i2 = input ( ),
///         cnot_z = z ( ),
///         cnot_x = x ( ),
///         xrot = X (m1),
///         zrot = Z (m0),
///         o0 = output ( ),
///         o1 = output ( ),
///         o2 = output ( ),
///     }
///     +
///     add_wire: {
///         i0 -- cnot_z,
///         i1 -- cnot_x -- o1,
///         i2 -- xrot -- zrot -- o2,
///     }
///     add_wire_h: { cnot_z -- o0 }
///     apply_bell: { i0 -- i1 }
///     +
///     apply_state: {
///         o0 (Spider::X(m0)),
///         o1 (Spider::X(m1)),
///     }
/// ).unwrap();
#[macro_export]
macro_rules! graph_diagram {
    (
        <$variant:ty>
        nodes : {
            $( $node_name:ident = $node:ident ( $( $arg:expr ),* $(,)? ) ),*
            $(,)?
        }
        +
        $(
            $bin_method:ident : {
                $( $node1_name:ident $( -- $nodek_name:ident )+ ),*
                $(,)?
            }
        )*
        +
        $(
            $uni_method:ident : {
                $( $node0_name:ident ( $( $uni_arg:expr ),* $(,)? ) ),*
                $(,)?
            }
        )*
    ) => {
        {
            let mut _diagram_ = $crate::graph::Diagram::<$variant>::new();
            $(
            let $node_name =
                _diagram_.add_node(
                    <$variant as $crate::graph::DiagramData>::Node::$node(
                        $( ($arg).into() ),*
                    )
                );
            )*
            Ok(())
            $($(.and_then(|_| {
                let mut _last_ = $node1_name;
                Ok(())
                $(.and_then(|_| {
                    let res = _diagram_.$bin_method(_last_, $nodek_name);
                    _last_ = $nodek_name;
                    res
                }))+
            }))*)*
            $(.and_then(|_| {
                Ok(())
                $(.and_then(|_| {
                    _diagram_.$uni_method($node0_name, $( ($uni_arg).into() ),*)
                }))*
            }))*
            .map(|_| {
                let _nodes_:
                    std::collections::HashMap<
                        &'static str,
                        $crate::graph::NodeId
                    > =
                    [$( (stringify!($node_name), $node_name) ),*]
                    .into_iter()
                    .collect();
                (_diagram_, _nodes_)
            })
        }
    }
}

pub use graph_diagram;

