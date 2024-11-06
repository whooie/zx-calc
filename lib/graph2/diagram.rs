use std::{ collections::VecDeque, fs, io::Write, path::Path };
use num_complex::Complex64 as C64;
use rustc_hash::FxHashMap;
use crate::{
    graph2::{
        GraphResult,
        Node,
        NodeData,
        NodeId,
        QubitId,
        Spider,
    },
    ketbra,
};

use crate::graph2::GraphError::*;

fn fst<T, U>(pair: (T, U)) -> T { pair.0 }

fn snd<T, U>(pair: (T, U)) -> U { pair.1 }

fn rev<T, U>(pair: (T, U)) -> (U, T) { (pair.1, pair.0) }

/// Represents a diagram in the ZX(H)-calculus.
///
/// Every node and edge is given a unique index for identification purposes.
/// Note that multiple edges may exist between two nodes.
///
/// Input and output qubits are numbered dynamically, but the order in which
/// they are added to the diagram is always preserved.
#[derive(Clone, Debug)]
pub struct Diagram {
    pub(crate) nodes: Vec<Option<Node>>,
    pub(crate) node_count: usize,
    pub(crate) wires: Vec<Option<Vec<NodeId>>>,
    pub(crate) wire_count: usize,
    pub(crate) inputs: Vec<NodeId>,
    pub(crate) outputs: Vec<NodeId>,
    pub(crate) free: Vec<NodeId>,
    pub(crate) scalar: C64,
}

impl Default for Diagram {
    fn default() -> Self { Self::new() }
}

impl Diagram {
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
            scalar: 1.0.into(),
        }
    }

    /// Create a new, disconnected diagram from a set of nodes.
    ///
    /// All nodes passed to this function are given unique IDs startinf from 0,
    /// i.e. if `n` nodes are passed to this function then the first node seen
    /// will have ID `0`, and the last will have ID `n - 1`. Input and output
    /// nodes will be given associated qubit numbers in the order they are seen,
    /// starting from zero.
    pub fn from_nodes<I>(nodes: I) -> Self
    where I: IntoIterator<Item = Node>
    {
        let mut dg = Self::new();
        nodes.into_iter()
            .for_each(|def| { dg.add_node(def); });
        dg
    }

    /// Return the global scalar on the diagram.
    pub fn scalar(&self) -> C64 { self.scalar }

    /// Apply a mapping function to the global scalar.
    pub fn map_scalar<F>(&mut self, map: F)
    where F: FnOnce(C64) -> C64
    {
        self.scalar = map(self.scalar);
    }

    /// Return `true` if the gobal scalar is zero.
    pub fn scalar_is_zero(&self) -> bool {
        const EPSILON: f64 = 1e-12;
        self.scalar.norm() < EPSILON
    }

    /// Return the number of nodes.
    pub fn count_nodes(&self) -> usize { self.node_count }

    /// Return the number of Z-spiders.
    pub fn count_z(&self) -> usize {
        self.nodes.iter().filter(|mb_n| mb_n.is_some_and(|n| n.is_z())).count()
    }

    /// Return the number of X-spiders.
    pub fn count_x(&self) -> usize {
        self.nodes.iter().filter(|mb_n| mb_n.is_some_and(|n| n.is_x())).count()
    }

    /// Return the number of H-boxes.
    pub fn count_h(&self) -> usize {
        self.nodes.iter().filter(|mb_n| mb_n.is_some_and(|n| n.is_h())).count()
    }

    /// Return the number of inputs.
    pub fn count_inputs(&self) -> usize { self.inputs.len() }

    /// Return the number of outputs.
    pub fn count_outputs(&self) -> usize { self.outputs.len() }

    /// Return the total number of spiders.
    pub fn count_spiders(&self) -> usize {
        self.nodes.iter()
            .filter(|mb_n| mb_n.is_some_and(|n| n.is_spider()))
            .count()
    }

    /// Return the number of wires.
    pub fn count_wires(&self) -> usize { self.wire_count }

    /// Get the node associated with a particular ID if it exists.
    pub fn get_node(&self, id: NodeId) -> Option<&Node> {
        self.nodes.get(id).and_then(|mb_n| mb_n.as_ref())
    }

    pub(crate) fn get_node_mut(&mut self, id: NodeId) -> Option<&mut Node> {
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
                        .find_map(|(qid, nid)| (*nid == id).then_some(qid))
                        .expect("bad book-keeping: unrecorded input")
                })
            })
    }

    // like `get_input_index`, but for use after the node has already been
    // removed -- this is purely a search through `self.inputs`
    pub(crate) fn find_input_index(&self, id: NodeId) -> Option<QubitId> {
        self.inputs.iter().enumerate()
            .find_map(|(qid, nid)| (*nid == id).then_some(qid))
    }

    /// Return the qubit index corresponding to a node ID, if the node exists
    /// and is an output.
    pub fn get_output_index(&self, id: NodeId) -> Option<QubitId> {
        self.get_node(id)
            .and_then(|n| {
                (!n.is_input()).then(|| {
                    self.outputs.iter().enumerate()
                        .find_map(|(qid, nid)| (*nid == id).then_some(qid))
                        .expect("bad book-keeping: unrecorded output")
                })
            })
    }

    // like `get_output_index`, but for use after the node has already been
    // removed -- this is purely a search through `self.outputs`
    pub(crate) fn find_output_index(&self, id: NodeId) -> Option<QubitId> {
        self.outputs.iter().enumerate()
            .find_map(|(qid, nid)| (*nid == id).then_some(qid))
    }

    // get the raw neighbor IDs of a node, if it exists
    pub(crate) fn get_neighbors_mut(&mut self, id: NodeId)
        -> Option<&mut Vec<NodeId>>
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
                            nnb.iter().filter(|id| **id == b).count();
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
    pub fn add_node(&mut self, node: Node) -> NodeId {
        let id = self.fresh_node_id();
        if node.is_input() {
            self.inputs.push(id);
        } else if node.is_output() {
            self.outputs.push(id);
        }
        let _ = self.nodes[id].insert(node);
        let _ = self.wires[id].insert(Vec::new());
        id
    }

    /// Remove the node associated with a particular ID and return its data if
    /// it exists.
    ///
    /// This method also removes all wires with an endpoint at the node.
    pub fn remove_node(&mut self, id: NodeId) -> GraphResult<Node> {
        // for brevity:
        // "mb" == "maybe" (type is Option<_>)
        // "nnb" == "neighbors" (type is Vec<_>)
        // "nb" == "neighbor" (type is NodeId)
        // "_of" == of the node being removed
        let node =
            self.nodes.get_mut(id)
            .ok_or(RemoveNodeMissingNode(id))?
            .take()
            .ok_or(RemoveNodeMissingNode(id))?;
        self.node_count -= 1;
        self.free.push(id);
        let nnb_of = self.wires[id].take().unwrap();
        let self_loops = nnb_of.iter().filter(|nb| **nb == id).count() / 2;
        self.wire_count -= nnb_of.len() - self_loops;
        for nb_of in nnb_of.into_iter() {
            while let Some(k) =
                self.wires[nb_of].as_ref()
                    .and_then(|nnb| {
                        nnb.iter().enumerate()
                            .find_map(|(k, nb)| (*nb == id).then_some(k))
                    })
            {
                self.wires[nb_of].as_mut().unwrap()
                    .swap_remove(k);
            }
        }
        if node.is_input() || self.inputs.contains(&id) {
            let k = self.find_input_index(id).unwrap();
            self.inputs.remove(k);
        }
        if node.is_output() || self.outputs.contains(&id) {
            let k = self.find_output_index(id).unwrap();
            self.outputs.remove(k);
        }
        Ok(node)
    }

    /// Like [`remove_node`][Self::remove_node], but returning a list of the
    /// node's neighbors alongside its data.
    pub fn remove_node_nb(&mut self, id: NodeId)
        -> GraphResult<(Node, Vec<NodeId>)>
    {
        let node =
            self.nodes.get_mut(id)
            .ok_or(RemoveNodeMissingNode(id))?
            .take()
            .ok_or(RemoveNodeMissingNode(id))?;
        self.node_count -= 1;
        self.free.push(id);
        let nnb_of = self.wires[id].take().unwrap();
        let self_loops = nnb_of.iter().filter(|nb| **nb == id).count() / 2;
        self.wire_count -= nnb_of.len() - self_loops;
        for &nb_of in nnb_of.iter() {
            while let Some(k) =
                self.wires[nb_of].as_ref()
                    .and_then(|nnb| {
                        nnb.iter().enumerate()
                            .find_map(|(k, nb)| (*nb == id).then_some(k))
                    })
            {
                self.wires[nb_of].as_mut().unwrap()
                    .swap_remove(k);
            }
        }
        if node.is_input() || self.inputs.contains(&id) {
            let k = self.find_input_index(id).unwrap();
            self.inputs.remove(k);
        }
        if node.is_output() || self.outputs.contains(&id) {
            let k = self.find_output_index(id).unwrap();
            self.outputs.remove(k);
        }
        Ok((node, nnb_of))
    }

    // remove the node associated with a particular ID and return its data and
    // neighbor indices, but *don't* remove the wires connected to it
    //
    // input/output book-keeping is still performed, however
    pub(crate) fn delete_node(&mut self, id: NodeId)
        -> Option<(Node, Vec<NodeId>)>
    {
        let node = self.nodes.get_mut(id)?.take()?;
        let nnb = self.wires.get_mut(id)?.take()?;
        if node.is_input() {
            let k = self.find_input_index(id).unwrap();
            self.inputs.remove(k);
        }
        if node.is_output() {
            let k = self.find_output_index(id).unwrap();
            self.outputs.remove(k);
        }
        self.free.push(id);
        self.node_count -= 1;
        Some((node, nnb))
    }

    // /// Replace the data at a given node ID, preserving all connected wires,
    // /// returning the previous data under the ID if it existed.
    // pub fn replace_node(&mut self, id: NodeId, mut data: Node) -> Option<Node> {
    //     if let Some(prev) = self.get_node_mut(id) {
    //         std::mem::swap(&mut data, prev);
    //         if data.is_input() && prev.is_output() {
    //             let k =
    //                 self.inputs.iter().enumerate()
    //                 .find_map(|(k, nid)| (*nid == id).then_some(k))
    //                 .expect("bad book-keeping: unrecorded input node");
    //             self.inputs.remove(k);
    //             self.outputs.push(id);
    //         } else if data.is_output() && prev.is_input() {
    //             let k =
    //                 self.outputs.iter().enumerate()
    //                 .find_map(|(k, nid)| (*nid == id).then_some(k))
    //                 .expect("bad book-keeping: unrecorded input node");
    //             self.outputs.remove(k);
    //             self.inputs.push(id);
    //         } else if data.is_input() && !prev.is_input() {
    //             let k =
    //                 self.inputs.iter().enumerate()
    //                 .find_map(|(k, nid)| (*nid == id).then_some(k))
    //                 .expect("bad book-keeping: unrecorded input node");
    //             self.inputs.remove(k);
    //         } else if data.is_output() && !prev.is_output() {
    //             let k =
    //                 self.inputs.iter().enumerate()
    //                 .find_map(|(k, nid)| (*nid == id).then_some(k))
    //                 .expect("bad book-keeping: unrecorded input node");
    //             self.outputs.remove(k);
    //         }
    //         Some(data)
    //     } else {
    //         None
    //     }
    // }

    /// Remove all nodes that have no wires going to any other node in the
    /// diagram, returning their IDs and data.
    pub fn remove_disconnected(&mut self) -> Vec<(NodeId, Node)> {
        // mutability rules mean we need to allocate
        let to_remove: Vec<NodeId> =
            self.wires.iter().enumerate()
            .filter_map(|(id, mb_nnb)| {
                mb_nnb.as_ref()
                    .is_some_and(|nnb| nnb.is_empty())
                    .then_some(id)
            })
            .collect();
        to_remove.into_iter()
            .flat_map(|id| self.delete_node(id).map(|(n, _)| (id, n)))
            .collect()
    }

    /// Add a wire between two nodes.
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
        self.wires[a].as_mut().unwrap().push(b);
        self.wires[b].as_mut().unwrap().push(a);
        self.wire_count += 1;
        Ok(())
    }

    /// Add a wire with an attached [`Input`][Node::Input] to a pre-existing
    /// node and return the new input node's ID.
    ///
    /// The pre-existing node cannot already be an input or output.
    pub fn add_input_wire(&mut self, id: NodeId) -> GraphResult<NodeId> {
        self.get_node(id)
            .ok_or(AddWireMissingNode(id))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(AddWireConnectedIO(id))
            })?;
        let input_id = self.add_node(Node::Input);
        self.add_wire(id, input_id)?;
        Ok(input_id)
    }

    /// Add a wire with an attached [`Output`][Node::Output] to a pre-existing
    /// node and return the new output node's ID.
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
        let output_id = self.add_node(Node::Output);
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
        self.has_node(a).then_some(())
            .ok_or(RemoveWireMissingNode(a))?;
        self.has_node(b).then_some(())
            .ok_or(RemoveWireMissingNode(b))?;

        let mut to_remove: Vec<usize> = Vec::new();

        let nnb_a = self.get_neighbors_mut(a).unwrap();
        let len_nnb_a = nnb_a.len();
        nnb_a.iter().enumerate()
            .filter_map(|(k, nid)| (*nid == b).then_some(k))
            .take(nwires.unwrap_or(len_nnb_a))
            .for_each(|k| { to_remove.push(k); });
        let mut removed_a = to_remove.len();
        to_remove.drain(..).rev()
            .for_each(|k| { nnb_a.swap_remove(k); });

        let nnb_b = self.get_neighbors_mut(b).unwrap();
        let len_nnb_b = nnb_b.len();
        nnb_b.iter().enumerate()
            .filter_map(|(k, nid)| (*nid == a).then_some(k))
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

    /// Get the node ID associated with the `q`-th input qubit index.
    pub fn input_id(&self, q: QubitId) -> Option<NodeId> {
        self.inputs.get(q).copied()
    }

    /// Get the node ID associated with the `q`-th output qubit index.
    pub fn output_id(&self, q: QubitId) -> Option<NodeId> {
        self.outputs.get(q).copied()
    }

    /// Replace an existing diagram input or output with a spider of arity 1.
    ///
    /// The new spider will have the same node ID as the input or output.
    ///
    /// Fails if the given node ID does not exist or is not an input or output.
    pub fn apply_state(&mut self, id: NodeId, spider: Spider)
        -> GraphResult<()>
    {
        if self.get_node(id).is_some_and(|n| !n.is_generator()) {
            let prev = self.get_node_mut(id).unwrap();
            *prev = spider.into();
            Ok(())
        } else {
            Err(ApplyStateNotIO(id))
        }
    }

    /// Like [`apply_state`][Self::apply_state], but using input qubit indices
    /// rather than bare node IDs.
    pub fn apply_state_input(&mut self, q: QubitId, spider: Spider)
        -> GraphResult<()>
    {
        self.input_id(q)
            .ok_or(ApplyStateMissingQubit(q))
            .and_then(|nid| self.apply_state(nid, spider))
    }

    /// Like [`apply_state`][Self::apply_state], but using output qubit indices
    /// rather than bare node IDs.
    pub fn apply_effect_output(&mut self, q: QubitId, spider: Spider)
        -> GraphResult<()>
    {
        self.output_id(q)
            .ok_or(ApplyStateMissingQubit(q))
            .and_then(|nid| self.apply_state(nid, spider))
    }

    /// Replace two existing diagram inputs or outputs with a Bell state/effect
    /// of optional phase. Returns the ID of the new spider node if a phase is
    /// passed.
    ///
    /// The inputs or outputs will be replaced with phaseless spiders (for
    /// book-keeping purposes) with the same IDs.
    ///
    /// Fails if the given node IDs do not exist, are not both inputs or
    /// outputs, or refer to the same qubit.
    pub fn apply_bell(&mut self, a: NodeId, b: NodeId, phase: Option<Spider>)
        -> GraphResult<Option<NodeId>>
    {
        if a == b { return Err(ApplyBellSameQubit); }
        self.get_node(a).zip(self.get_node(b))
            .is_some_and(|(na, nb)| {
                (na.is_input() && nb.is_input())
                    || (na.is_output() && nb.is_output())
            })
            .then_some(())
            .ok_or(ApplyBellNotIO(a, b))?;
        let _ = self.nodes[a].insert(Node::z());
        let _ = self.nodes[b].insert(Node::z());
        if let Some(spider) = phase {
            let spider_id = self.add_node(spider.into());
            self.add_wire(a, spider_id)?;
            self.add_wire(b, spider_id)?;
            Ok(Some(spider_id))
        } else {
            self.add_wire(a, b)?;
            Ok(None)
        }
    }

    /// Like [`apply_bell`][Self::apply_bell], but using input qubit indices
    /// rather than bare node IDs.
    pub fn apply_bell_input(
        &mut self,
        qa: QubitId,
        qb: QubitId,
        phase: Option<Spider>,
    ) -> GraphResult<Option<NodeId>>
    {
        self.input_id(qa)
            .ok_or(ApplyBellMissingQubit(qa))
            .and_then(|nida| {
                self.input_id(qb)
                    .ok_or(ApplyBellMissingQubit(qb))
                    .map(|nidb| (nida, nidb))
            })
            .and_then(|(nida, nidb)| self.apply_bell(nida, nidb, phase))
    }

    /// Like [`apply_bell`][Self::apply_bell], but using output qubit indices
    /// rather than bare node IDs.
    pub fn apply_bell_output(
        &mut self,
        qa: QubitId,
        qb: QubitId,
        phase: Option<Spider>,
    ) -> GraphResult<Option<NodeId>>
    {
        self.output_id(qa)
            .ok_or(ApplyBellMissingQubit(qa))
            .and_then(|nida| {
                self.output_id(qb)
                    .ok_or(ApplyBellMissingQubit(qb))
                    .map(|nidb| (nida, nidb))
            })
            .and_then(|(nida, nidb)| self.apply_bell(nida, nidb, phase))
    }

    /// Return an iterator over all nodes, visited in index order.
    ///
    /// The iterator item type is `(`[`NodeId`]`, &`[`Node`]`)`.
    pub fn nodes(&self) -> Nodes<'_> {
        Nodes { len: self.node_count, iter: self.nodes.iter().enumerate() }
    }

    /// Return an iterator over all interior (non-input/output) nodes in the
    /// diagram, visited in index order.
    ///
    /// The iterator item type is `(`[`NodeId`]`, &`[`Node`]`)`.
    pub fn nodes_inner(&self) -> NodesInner<'_> {
        NodesInner { iter: self.nodes.iter().enumerate() }
    }

    /// Return an iterator over all diagram input node IDs, visited in qubit
    /// index order.
    ///
    /// The iterator item type is `(`[`QubitId`]`, `[`NodeId`]`)`.
    pub fn inputs(&self) -> Inputs<'_> {
        Inputs { iter: self.inputs.iter().enumerate() }
    }

    /// Return an iterator over all diagram output node IDs, visited in qubit
    /// index order.
    ///
    /// The iterator item type is `(`[`QubitId`]`, `[`NodeId`]`)`.
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
    /// arbitrary order. Each wire is only visited once.
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`NodeId`]`)`.
    pub fn wires(&self) -> Wires<'_> {
        Wires {
            nb_group: None,
            len: self.wire_count,
            seen: Vec::with_capacity(self.wire_count),
            iter: self.wires.iter().enumerate() }
    }

    /// Return an iterator over all wires in the diagram, visited in mostly
    /// arbitrary order, with node data attached.
    ///
    /// Specifically, the iterator steps through wires as pairs of node IDs such
    /// that those whose left IDs are smaller in value are visited before those
    /// larger in value. In other words, the left ID monotonically increases
    /// over the course of iteration while right IDs are free to vary in
    /// arbitrary order. Each wire is only visited once.
    ///
    /// The iterator item type is `(`[`NodeData`]`, `[`NodeData`]`)`.
    pub fn wires_data(&self) -> WiresData<'_> {
        WiresData { dg: self, iter: self.wires() }
    }

    /// Return an iterator over all wires connected to only interior
    /// (non-input/output) nodes in the diagram, visited in mostly arbitrary
    /// order (see [`wires`][Self::wires]).
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`NodeId`]`)`.
    pub fn wires_inner(&self) -> WiresInner<'_> {
        WiresInner { dg: self, iter: self.wires() }
    }

    /// Return an iterator over all wires connected to only interior
    /// (non-input/output) nodes in the diagram, visited in mostly arbitrary
    /// order (see [`wires_data`][Self::wires_data]).
    ///
    /// The iterator item type is `(`[`NodeData`]`, `[`NodeData`]`)`.
    pub fn wires_data_inner(&self) -> WiresDataInner<'_> {
        WiresDataInner { iter: self.wires_data() }
    }

    /// Return the ID of an arbitrary neighbor of the given node ID, if it
    /// exists.
    pub fn get_neighbor_of(&self, id: NodeId) -> Option<NodeId> {
        self.wires.get(id)
            .and_then(|mb_nnb| mb_nnb.as_ref())
            .and_then(|nnb| nnb.first())
            .copied()
    }

    /// Return an iterator over the IDs of neighbors of the given node ID,
    /// visited in arbitrary order, if they exist.
    ///
    /// Note that this iterator will contain duplicate elements if there are
    /// multiple wires going to the same neighbor.
    ///
    /// The iterator item type is [`NodeId`].
    pub fn neighbor_ids_of(&self, id: NodeId) -> Option<NeighborIds<'_>> {
        self.wires.get(id)
            .and_then(|mb_nnb| {
                mb_nnb.as_ref()
                    .map(|nnb| {
                        NeighborIds {
                            id,
                            switch: true,
                            iter: nnb.iter(),
                        }
                    })
            })
    }

    /// Find the ID of the first neighbor of the given node ID that satisfies
    /// some predicate, if it exists.
    pub fn find_neighbor_id_of<F>(&self, id: NodeId, mut pred: F)
        -> Option<NodeId>
    where F: FnMut(NodeId) -> bool
    {
        self.neighbor_ids_of(id)
            .and_then(|mut nnb| nnb.find(|id| pred(*id)))
    }

    /// Find the ID of the first neighbor of the given node ID that satisfies
    /// some predicate, if it exists, and map to some output.
    pub fn find_map_neighbor_id_of<F, U>(&self, id: NodeId, pred: F)
        -> Option<U>
    where F: FnMut(NodeId) -> Option<U>
    {
        self.neighbor_ids_of(id)
            .and_then(|mut nnb| nnb.find_map(pred))
    }

    /// Return an iterator over neighbors of the given node ID, visited in
    /// arbitrary order, if it exists.
    ///
    /// Note that this iterator will contain duplicate elements if there are
    /// multiple wires going to the same neighbor.
    ///
    /// The iterator item type is `(`[`NodeId`]`, &`[`Node`]`)`.
    pub fn neighbors_of(&self, id: NodeId) -> Option<Neighbors<'_>> {
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
    pub fn find_neighbor_of<F>(&self, id: NodeId, mut pred: F)
        -> Option<(NodeId, &Node)>
    where F: FnMut(NodeId, &Node) -> bool
    {
        self.neighbors_of(id)
            .and_then(|mut nnb| nnb.find(|(id, n)| pred(*id, n)))
    }

    /// Find the first neighbor of the given node ID that satisfies some
    /// predicate, if it and the ID exist, and map to some output.
    pub fn find_map_neighbor_of<F, U>(&self, id: NodeId, mut pred: F)
        -> Option<U>
    where F: FnMut(NodeId, &Node) -> Option<U>
    {
        self.neighbors_of(id)
            .and_then(|mut nnb| nnb.find_map(|(id, n)| pred(id, n)))
    }

    /// Return an iterator over all interior (non-input/output) neighbors of the
    /// given node ID, visited in arbitrary order, if it exists.
    ///
    /// Note that this iterator will contain duplicate elements if there are
    /// multiple wires going to the same neighbor.
    ///
    /// The iterator item type is `(`[`NodeId`]`, &`[`Node`]`)`.
    pub fn neighbors_of_inner(&self, id: NodeId) -> Option<NeighborsInner<'_>> {
        self.neighbors_of(id)
            .map(|iter| NeighborsInner { iter })
    }

    /// Find the first interior (non-input/output) neighbor of the given node ID
    /// that satisfies some predicate, if it and the ID exist.
    pub fn find_neighbor_of_inner<F>(&self, id: NodeId, mut pred: F)
        -> Option<(NodeId, &Node)>
    where F: FnMut(NodeId, &Node) -> bool
    {
        self.neighbors_of_inner(id)
            .and_then(|mut nnb| nnb.find(|(id, n)| pred(*id, n)))
    }

    /// Find the first interior (non-input/output) neighbor of the given node ID
    /// that satisfies some predicate, if it and the ID exist, and map to some
    /// output.
    pub fn find_map_neighbor_of_inner<F, U>(&self, id: NodeId, mut pred: F)
        -> Option<U>
    where F: FnMut(NodeId, &Node) -> Option<U>
    {
        self.neighbors_of_inner(id)
            .and_then(|mut nnb| nnb.find_map(|(id, n)| pred(id, n)))
    }

    /// Swap inputs and outputs, flip the signs of all spiders' phases, and
    /// conjugate all H-boxes' arguments as well as the global scalar, consuming
    /// `self`.
    pub fn adjoint(mut self) -> Self {
        self.adjoint_mut();
        self
    }

    /// Swap inputs and outputs, flip the signs of all spiders' phases, and
    /// conjugate all H-box arguments as well as the global scalar, modifying
    /// `self` in place.
    pub fn adjoint_mut(&mut self) -> &mut Self {
        self.nodes.iter_mut()
            .flatten()
            .for_each(|node| { node.adjoint_mut(); });
        std::mem::swap(&mut self.inputs, &mut self.outputs);
        self.scalar = self.scalar.conj();
        self
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
            .for_each(|nb| { *nb += shift; });
        inputs.iter_mut()
            .for_each(|nid| { *nid += shift; });
        outputs.iter_mut()
            .for_each(|nid| { *nid += shift; });
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

    /// Return the tensor product of `self` and `other`, consuming both.
    ///
    /// All node IDs from `other` will be adjusted to avoid collision. The exact
    /// size of the shift depends on how much space has been allocated for nodes
    /// in `self` (both current and previous), and is returned alongside the
    /// output diagram.
    pub fn tensor(mut self, other: Self) -> (Self, usize) {
        let shift = self.tensor_with(other);
        (self, shift)
    }

    /// Compute the tensor product of `self` and `other`, consuming `other` and
    /// modifying `self` in place.
    ///
    /// All node IDs from `other` will be adjusted to avoid collision. The exact
    /// size of the shift depends on how much space has been allocated for nodes
    /// in `self` (both current and previous), and will be returned.
    pub fn tensor_with(&mut self, other: Self) -> usize {
        self.append(other)
    }

    /// Return the composition `self ∘ other`, attempting to match the outputs
    /// of `other` to the inputs of `self` in qubit index order, consuming both
    /// `self` and `other`.
    ///
    /// The IDs of all nodes from `other` will be adjusted to avoid collision.
    /// The exact size of the shift depends on how much space has been allocated
    /// for nodes in `self` (both current and previous), and will be returned
    /// alongside the output diagram.
    ///
    /// This operation will fail if the numbers of inputs and outputs are not
    /// equal, or if they are equal but some inputs or outputs are not free
    /// wires. `self` will not be modified in this case.
    ///
    /// See also [`compose_rev`][Self::compose_rev].
    pub fn compose(mut self, other: Self) -> GraphResult<(Self, usize)> {
        let shift = self.compose_with(other)?;
        Ok((self, shift))
    }

    /// Compute the composition `self ∘ other`, attempting to match the outputs
    /// of `other` to the inputs of `self` in qubit index order, consuming
    /// `other` and modifying `self` in place.
    ///
    /// The IDs of all nodes from `other` will be adjusted to avoid collision.
    /// The exact size of the shift depends on how much space has been allocated
    /// for nodes in `self` (both current and previous), and will be returned.
    ///
    /// This operation will fail if the numbers of inputs and outputs are not
    /// equal, or if they are equal but some inputs or outputs are not free
    /// wires. `self` will not be modified in this case.
    ///
    /// See also [`compose_with_rev`][Self::compose_with_rev].
    pub fn compose_with(&mut self, other: Self) -> GraphResult<usize> {
        // find (input/output, interior neighbor)
        let self_inputs: Vec<(NodeId, Option<NodeId>)> =
            self.inputs()
            .filter(|(_, nid)| self.nodes[*nid].as_ref().unwrap().is_input())
            .map(|(_, nid)| {
                let mb_interior =
                    self.neighbors_of(nid)
                    .and_then(|mut nnb| nnb.next().map(|(nb, _)| nb));
                (nid, mb_interior)
            })
            .collect();
        let other_outputs: Vec<(NodeId, Option<NodeId>)> =
            other.outputs()
            .filter(|(_, nid)| other.nodes[*nid].as_ref().unwrap().is_output())
            .map(|(_, nid)| {
                let mb_interior =
                    other.neighbors_of(nid)
                    .and_then(|mut nnb| nnb.next().map(|(nb, _)| nb));
                (nid, mb_interior)
            })
            .collect();
        (self_inputs.len() == other_outputs.len()).then_some(())
            .ok_or(NonMatchingIO(other_outputs.len(), self_inputs.len()))?;
        let shift = self.append(other);
        let it = other_outputs.into_iter().zip(self_inputs);
        for ((out_id, mb_out_int), (in_id, mb_in_int)) in it {
            self.remove_node(out_id + shift).unwrap();
            self.remove_node(in_id).unwrap();
            match (mb_out_int, mb_in_int) {
                (Some(out_int), Some(in_int)) => {
                    self.add_wire(out_int + shift, in_int).unwrap();
                },
                (Some(out_int), None) => {
                    self.add_output_wire(out_int + shift).unwrap();
                },
                (None, Some(in_int)) => {
                    self.add_input_wire(in_int).unwrap();
                },
                (None, None) => {
                    let new_in = self.add_node(Node::Input);
                    let new_out = self.add_node(Node::Output);
                    self.add_wire(new_in, new_out).unwrap();
                },
            }
        }
        Ok(shift)
    }

    /// Return the composition `other ∘ self`, attempting to match the outputs
    /// of `self` to the inputs of `other` in qubit index order, consuming both
    /// `self` and `other`.
    ///
    /// The IDs of all nodes from `other` will be adjusted to avoid collision.
    /// The exact size of the shift depends on how much space has been allocated
    /// for nodes in `self` (both current and previous), and will be returned
    /// alongside the output diagram.
    ///
    /// This operation will fail if the numbers of inputs and outputs are not
    /// equal, or if they are equal but some inputs or outputs are not free
    /// wires. `self` will not be modified in this case.
    ///
    /// See also [`compose`][Self::compose].
    pub fn compose_rev(mut self, other: Self)
        -> GraphResult<(Self, usize)>
    {
        let shift = self.compose_with_rev(other)?;
        Ok((self, shift))
    }

    /// Compute the composition `other ∘ self`, attempting to match the outputs
    /// of `self` to the inputs of `other` in qubit index order, consuming
    /// `other` and modifying `self` in place.
    ///
    /// The IDs of all nodes from `other` will be adjusted to avoid collision.
    /// The exact size of the shift depends on how much space has been allocated
    /// for nodes in `self` (both current and previous), and will be returned.
    ///
    /// This operation will fail if the numbers of inputs and outputs are not
    /// equal, or if they are equal but some inputs or outputs are not free
    /// wires. `self` will not be modified in this case.
    ///
    /// See also [`compose_with`][Self::compose_with].
    pub fn compose_with_rev(&mut self, other: Self) -> GraphResult<usize> {
        // (input/output, interior neighbor)
        let self_outputs: Vec<(NodeId, Option<NodeId>)> =
            self.outputs()
            .filter(|(_, nid)| self.nodes[*nid].as_ref().unwrap().is_output())
            .map(|(_, nid)| {
                let mb_interior =
                    self.neighbors_of(nid)
                    .and_then(|mut nnb| nnb.next().map(|(nb, _)| nb));
                (nid, mb_interior)
            })
            .collect();
        let other_inputs: Vec<(NodeId, Option<NodeId>)> =
            other.inputs()
            .filter(|(_, nid)| other.nodes[*nid].as_ref().unwrap().is_input())
            .map(|(_, nid)| {
                let mb_interior =
                    other.neighbors_of(nid)
                    .and_then(|mut nnb| nnb.next().map(|(nb, _)| nb));
                (nid, mb_interior)
            })
            .collect();
        (self_outputs.len() == other_inputs.len()).then_some(())
            .ok_or(NonMatchingIO(self_outputs.len(), other_inputs.len()))?;
        let shift = self.append(other);
        let it = self_outputs.into_iter().zip(other_inputs);
        for ((out_id, mb_out_int), (in_id, mb_in_int)) in it {
            self.remove_node(out_id).unwrap();
            self.remove_node(in_id + shift).unwrap();
            match (mb_out_int, mb_in_int) {
                (Some(out_int), Some(in_int)) => {
                    self.add_wire(out_int, in_int + shift).unwrap();
                },
                (Some(out_int), None) => {
                    self.add_output_wire(out_int).unwrap();
                },
                (None, Some(in_int)) => {
                    self.add_input_wire(in_int + shift).unwrap();
                },
                (None, None) => {
                    let new_in = self.add_node(Node::Input);
                    let new_out = self.add_node(Node::Output);
                    self.add_wire(new_in, new_out).unwrap();
                },
            }
        }
        Ok(shift)
    }

    // BFS explore from a given starting queue, accumulating a list of ketbra
    // elements
    fn ketbra_explore<'a>(
        &'a self,
        wire_nums: &WireStore,
        visited: &mut Vec<NodeId>,
        to_visit: &mut VecDeque<(NodeId, &'a Node)>,
        elements: &mut Vec<ketbra::Element>,
    ) {
        // upper bound on possible qubit indices
        let qcount = self.inputs.len().max(self.outputs.len());
        // want wire IDs line up with qubit indices when connected to
        // inputs/outputs, so shift any ID that could be a qubit ID by the total
        // wire count
        let nonq_wire = |id: usize| -> usize {
            if id < qcount { self.wire_count + id } else { id }
        };
        // self-loops are treated as inputs to a Bell effect, whose input wire
        // IDs are shifted to come after all existing (possibly shifted) wires
        let mut bell_wire = self.wire_count + qcount..;

        // loop variables
        let mut ins: Vec<usize>; // inputs to an element
        let mut outs: Vec<usize>; // outputs from an element
        let mut bell_wires: Vec<(usize, usize)> = Vec::new(); // for self-loop wires
        let mut empty_wires: Vec<(QubitId, QubitId)> = Vec::new(); // direct input-to-output
        while let Some((id, node)) = to_visit.pop_back() {
            visited.push(id);
            ins = Vec::new();
            outs = Vec::new();
            for (id2, node2) in self.neighbors_of(id).unwrap() {
                // self-loops
                if id2 == id {
                    let b1 = bell_wire.next().unwrap();
                    let b2 = bell_wire.next().unwrap();
                    outs.push(b1);
                    outs.push(b2);
                    bell_wires.push((b1, b2));
                    continue;
                }
                // empty wires
                if node.is_input() && node2.is_output() {
                    visited.push(id2);
                    let qin = self.get_input_index(id).unwrap();
                    let qout = self.get_output_index(id).unwrap();
                    empty_wires.push((qin, qout));
                    continue;
                }
                // everything else
                if !visited.contains(&id2)
                    && !to_visit.contains(&(id2, node2))
                {
                    to_visit.push_front((id2, node2));
                }
                if !node.is_generator() { break; }
                if node2.is_input() {
                    let qin = self.get_input_index(id2).unwrap();
                    ins.push(qin);
                } else if !node2.is_output() && visited.contains(&id2) {
                    wire_nums.get(id, id2).unwrap().iter()
                        .for_each(|wid| { ins.push(nonq_wire(*wid)); });
                } else if node2.is_output() {
                    let qout = self.get_output_index(id2).unwrap();
                    outs.push(qout);
                } else if !visited.contains(&id2) {
                    wire_nums.get(id, id2).unwrap().iter()
                        .for_each(|wid| { outs.push(nonq_wire(*wid)); });
                }
            }
            if node.is_generator() {
                elements.push(node.as_element(ins, outs));
            }
            if !bell_wires.is_empty() {
                for (b1, b2) in bell_wires.drain(..) {
                    elements.push(ketbra::Element::cap([b1, b2], None));
                }
            }
        }
        // empty wires may include swaps, which have to be dealt with
        //
        // individual swaps are be determined by sorting on input wire indices,
        // and then finding the swaps required to sort the output wire indices
        if !empty_wires.is_empty() {
            empty_wires.sort_by(|(qin_l, _), (qin_r, _)| qin_l.cmp(qin_r));
            let (idents_in, mut idents_out): (Vec<QubitId>, Vec<QubitId>) =
                empty_wires.into_iter().unzip();
            let mut swaps: Vec<ketbra::Element> = Vec::new();
            let mut mb_mismatch: Option<(usize, (&QubitId, &QubitId))>;
            loop {
                mb_mismatch =
                    idents_in.iter().zip(idents_out.iter()).enumerate()
                    .find(|(_, (qin, qout))| qin != qout);
                if let Some((k, (qin, qswap))) = mb_mismatch {
                    let (kswap, _) =
                        idents_out.iter().enumerate()
                        .find(|(_, qout)| qin == *qout)
                        .unwrap();
                    swaps.push(ketbra::Element::swap([*qin, *qswap]));
                    idents_out.swap(k, kswap);
                } else {
                    break;
                }
            }
            swaps.reverse();
            elements.append(&mut swaps);
        }
    }

    /// Convert `self` to a [`ketbra::Diagram`] representation.
    pub fn as_ketbra(&self) -> ketbra::Diagram {
        // assemble the ketbra diagram by iterating over nodes and analyzing
        // input/output wires relative to what's already been seen
        //
        // have to BFS explore starting from the input nodes in order to ensure
        // that input-adjencent nodes are placed first in the ketbra diagram
        //
        // we also want to have wire numbers in the ketbra diagram line up with
        // qubit indices for convenience, so if a given wire id (normally used
        // as-is for a ketbra wire index) coincides with a possible qubit index,
        // shift it by the maximum wire id in the (graph) diagram
        //
        // do this in two steps because nodes with paths to inputs/outputs need
        // to be visited in a special order, but everything else (i.e. part of a
        // scalar) doesn't
        //
        // self-wires are dealt with by adding two extra outgoing wires
        // immediately coupled to a spiderless Bell effect

        // ketbra accumulator
        let mut elements: Vec<ketbra::Element> =
            Vec::with_capacity(self.node_count);
        // have to index wires so that each one can have a unique ID in the
        // ketbra diagram
        let wire_nums: WireStore = self.wires().collect();

        // first step: all non-scalar nodes
        let mut visited: Vec<NodeId> = Vec::with_capacity(self.node_count);
        // init with input nodes to make sure they're seen first
        let mut to_visit: VecDeque<(NodeId, &Node)> =
            self.inputs()
            .map(|(_, nid)| (nid, self.get_node(nid).unwrap()))
            .collect();
        self.ketbra_explore(
            &wire_nums, &mut visited, &mut to_visit, &mut elements);

        // second step: all nodes that aren't part of a scalar
        // reset `to_visit` with everything not already visited
        to_visit
            = self.nodes_inner()
            .filter(|(id, _)| !visited.contains(id))
            .collect();
        self.ketbra_explore(
            &wire_nums, &mut visited, &mut to_visit, &mut elements);

        ketbra::Diagram::new(elements)
    }

    /// Find all nodes that are part of a scalar subgraph.
    ///
    /// A node is part of a scalar subgraph if there is no path from it to any
    /// `Input` or `Output` node. Note that the returned list of nodes may
    /// comprise more than one scalar.
    pub fn find_scalar_nodes(&self) -> Vec<(NodeId, &Node)> {
        // a node is part of a scalar subgraph if there is no path from it to
        // any Input or Output node
        //
        // all scalar subgraphs are found by starting with the set of all nodes
        // and removing those seen by BFS explorations starting at each of the
        // Input and Output nodes; everything that's left is part of a scalar
        // subgraph
        //
        // returned value may contain multiple disconnected scalar subgraphs

        let mut has_io: bool = false;
        let mut nodes: Vec<Option<&Node>> =
            self.nodes.iter()
            .map(|mb_node| {
                has_io = has_io || mb_node.is_some_and(|n| !n.is_generator());
                mb_node.as_ref()
            })
            .collect();
        if !has_io {
            return nodes.into_iter().enumerate()
                .filter_map(|(nid, mb_node)| mb_node.map(|n| (nid, n)))
                .collect();
        }
        let mut to_visit: VecDeque<NodeId>;
        let mut io_visited: Vec<NodeId> =
            Vec::with_capacity(self.inputs.len() + self.outputs.len());
        for io in self.inputs.iter().chain(&self.outputs).copied() {
            if io_visited.contains(&io) { continue; }
            to_visit = vec![io].into();
            while let Some(id) = to_visit.pop_back() {
                for (id2, node2) in self.neighbors_of(id).unwrap() {
                    if nodes[id2].is_some() {
                        to_visit.push_front(id2);
                    }
                    if node2.is_input() || node2.is_output() {
                        io_visited.push(id2);
                        nodes[id2] = None;
                    }
                }
                nodes[id] = None;
            }
        }
        nodes.into_iter().enumerate()
            .filter_map(|(nid, mb_node)| mb_node.map(|n| (nid, n)))
            .collect()
    }

    #[allow(unused_variables, unused_mut)]
    fn compute_scalar(&self, nodes: &[(NodeId, &Node)]) -> C64 {
        // compute the value of the scalar subgraph(s) by converting to a ketbra
        // diagram and contracting
        //
        // if `find_scalar_nodes` did its job right, the result of the
        // contraction is guaranteed to be an Element of a single term with no
        // inputs or outputs, the amplitude of which is the scalar
        //
        // conversion to a ketbra diagram is done by iterating over nodes and
        // analyzing input/output wires relative to what's already been seen
        //
        // this subgraph can be deformed arbitrarily, so there's no need to care
        // about the order of iteration
        //
        // self-wires are dealt with by adding two extra outgoing wires
        // immediately coupled to a spiderless Bell effect

        todo!()
    }

    /// Find all nodes that are part of a scalar subgraph and compute their
    /// total product.
    ///
    /// See also [`find_scalar_nodes`][Self::find_scalar_nodes].
    pub fn get_scalar(&self) -> C64 {
        let nodes = self.find_scalar_nodes();
        self.compute_scalar(&nodes)
    }

    /// Find, compute, and remove all scalars from `self`, returning their total
    /// product.
    ///
    /// See also [`find_scalar_nodes`][Self::find_scalar_nodes] and
    /// [`remove_scalar_nodes`][Self::remove_scalar_nodes].
    pub fn remove_scalars(&mut self) -> C64 {
        let nodes = self.find_scalar_nodes();
        let scalar = self.compute_scalar(&nodes);
        nodes.into_iter()
            .map(fst)
            .collect::<Vec<NodeId>>() // need to drop Node references to self
            .into_iter()
            .for_each(|id| { self.remove_node(id).unwrap(); });
        scalar
    }

    /// Like [`remove_scalars`][Self::remove_scalars], but do not actually
    /// compute any scalar values, only remove their nodes.
    ///
    /// See also [`find_scalar_nodes`][Self::find_scalar_nodes].
    pub fn remove_scalar_nodes(&mut self) {
        self.find_scalar_nodes()
            .into_iter()
            .map(fst)
            .collect::<Vec<NodeId>>() // need to drop Node references to self
            .into_iter()
            .for_each(|id| { self.remove_node(id).unwrap(); });
    }

    /// Return an object containing an encoding of `self` in the [DOT
    /// language][dot-lang].
    ///
    /// Rendering this object using the default formatter will result in a full
    /// DOT string representation of the diagram.
    ///
    /// [dot-lang]: https://en.wikipedia.org/wiki/DOT_(graph_description_language)
    pub fn to_graphviz(&self) -> GraphResult<tabbycat::Graph> {
        use tabbycat::*;
        use tabbycat::attributes::*;
        use crate::vizdefs::*;
        // initial declarations
        let mut statements =
            StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rankdir(RankDir::LR)),
            )
            .add_attr(
                AttrType::Node,
                AttrList::new()
                    .add_pair(fontname(FONT))
                    .add_pair(fontsize(FONTSIZE))
                    .add_pair(margin(NODE_MARGIN)),
            );
        // ensure all inputs are in a subgraph at the same rank
        let mut inputs_subgraph_stmt =
            StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rank(RankType::Source)),
            );
        let mut prev: Option<usize> = None;
        for (qid, &nid) in self.inputs.iter().enumerate() {
            let node = self.get_node(nid).unwrap();
            let attrs =
                if node.is_generator() {
                    node.graph_attrs()
                        .add_pair(xlabel(format!("In {}", qid)))
                } else {
                    AttrList::new()
                        .add_pair(label(format!("In {}", qid)))
                        .add_pair(shape(Shape::Plaintext))
                };
            inputs_subgraph_stmt =
                inputs_subgraph_stmt.add_node(nid.into(), None, Some(attrs));
            if let Some(pid) = prev {
                inputs_subgraph_stmt =
                    inputs_subgraph_stmt.add_edge(
                        Edge::head_node(pid.into(), Some(Port::compass(Compass::South)))
                        .line_to_node(nid.into(), Some(Port::compass(Compass::North)))
                        .add_attrpair(style(Style::Invisible))
                    );
            }
            prev = Some(nid);
        }
        statements =
            statements.add_subgraph(SubGraph::cluster(inputs_subgraph_stmt));
        // ensure all outputs are in a subgraph at the same rank
        let mut outputs_subgraph_stmt =
            StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rank(RankType::Sink)),
            );
        prev = None;
        for (qid, &nid) in self.outputs.iter().enumerate() {
            let node = self.get_node(nid).unwrap();
            let attrs =
                if node.is_generator() {
                    node.graph_attrs()
                        .add_pair(xlabel(format!("Out {}", qid)))
                } else {
                    AttrList::new()
                        .add_pair(label(format!("Out {}", qid)))
                        .add_pair(shape(Shape::Plaintext))
                };
            outputs_subgraph_stmt =
                outputs_subgraph_stmt.add_node(nid.into(), None, Some(attrs));
            if let Some(pid) = prev {
                outputs_subgraph_stmt =
                    outputs_subgraph_stmt.add_edge(
                        Edge::head_node(pid.into(), Some(Port::compass(Compass::South)))
                        .line_to_node(nid.into(), Some(Port::compass(Compass::North)))
                        .add_attrpair(style(Style::Invisible))
                    );
            }
            prev = Some(nid);
        }
        statements =
            statements.add_subgraph(SubGraph::cluster(outputs_subgraph_stmt));
        // add interior nodes
        for (id, node) in self.nodes_inner() {
            if self.inputs.contains(&id) || self.outputs.contains(&id) {
                continue;
            }
            let attrs = node.graph_attrs();
            statements = statements.add_node(id.into(), None, Some(attrs));
        }
        // add the overall scalar
        let mut a = self.scalar;
        a.re = (1e6 * a.re).round() / 1e6;
        a.im = (1e6 * a.im).round() / 1e6;
        let scalar_id = self.nodes.len();
        let attrs =
            AttrList::new()
            .add_pair(label(format!("{}", a)))
            .add_pair(shape(Shape::Rectangle))
            .add_pair(style(Style::Filled))
            .add_pair(fillcolor(H_COLOR));
        statements = statements.add_node(scalar_id.into(), None, Some(attrs));
        // add wires
        for (left, right) in self.wires() {
            statements =
                statements.add_edge(
                    Edge::head_node(left.into(), None)
                    .line_to_node(right.into(), None)
                );
        }
        let graphviz =
            GraphBuilder::default()
                .graph_type(GraphType::Graph)
                .strict(false)
                .id(Identity::quoted(""))
                .stmts(statements)
                .build()
                .unwrap();
        Ok(graphviz)
    }

    /// Like [`to_graphviz`][Self::to_graphviz], but render directly to a string
    /// and write it to `path`.
    pub fn save_graphviz<P>(&self, path: P) -> GraphResult<()>
    where P: AsRef<Path>
    {
        let graphviz = self.to_graphviz()?;
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

// helper for ketbra conversions -- this acts as a vendor for wire IDs so that,
// when constructing a ketbra diagram, each wire can be guaranteed to have a
// unique ID
#[derive(Clone, Debug)]
struct WireStore(FxHashMap<(NodeId, NodeId), Vec<usize>>);

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
    fn get(&self, a: NodeId, b: NodeId) -> Option<&Vec<usize>> {
        self.0.get(&(a, b)).or_else(|| self.0.get(&(b, a)))
    }
}

/// Iterator over all nodes in a diagram, visited in index order.
///
/// The iterator item type is `(`[`NodeId`]`, &`[`Node`]`)`.
#[derive(Clone, Debug)]
pub struct Nodes<'a> {
    len: usize,
    iter: std::iter::Enumerate<std::slice::Iter<'a, Option<Node>>>
}

impl<'a> Iterator for Nodes<'a> {
    type Item = (NodeId, &'a Node);

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

impl<'a> DoubleEndedIterator for Nodes<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|(_, mb_n)| mb_n.is_some())
            .map(|(id, mb_n)| {
                self.len = self.len.saturating_sub(1);
                (id, mb_n.as_ref().unwrap())
            })
    }
}

impl<'a> ExactSizeIterator for Nodes<'a> {
    fn len(&self) -> usize { self.len }
}

impl<'a> std::iter::FusedIterator for Nodes<'a> { }

/// Iterator over all interior (non-input/output) nodes in a diagram, visited in
/// index order.
///
/// The iterator item type is `(`[`NodeId`]`, &`[`Node`]`)`.
pub struct NodesInner<'a> {
    iter: std::iter::Enumerate<std::slice::Iter<'a, Option<Node>>>
}

impl<'a> Iterator for NodesInner<'a> {
    type Item = (NodeId, &'a Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find_map(|(id, mb_n)| {
            mb_n.as_ref()
                .and_then(|n| n.is_generator().then_some((id, n)))
        })
    }
}

impl<'a> DoubleEndedIterator for NodesInner<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|(_, mb_n)| mb_n.is_some_and(|n| n.is_generator()))
            .map(|(id, mb_n)| (id, mb_n.as_ref().unwrap()))
    }
}

impl<'a> std::iter::FusedIterator for NodesInner<'a> { }

/// Iterator over all diagram input node IDs, visited in qubit index order.
///
/// The iterator item type is `(`[`QubitId`]`, `[`NodeId`]`)`.
pub struct Inputs<'a> {
    iter: std::iter::Enumerate<std::slice::Iter<'a, NodeId>>
}

impl<'a> Iterator for Inputs<'a> {
    type Item = (NodeId, QubitId);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|(qid, nid)| (qid, *nid))
    }
}

impl<'a> DoubleEndedIterator for Inputs<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.next_back().map(|(qid, nid)| (qid, *nid))
    }
}

impl<'a> ExactSizeIterator for Inputs<'a> {
    fn len(&self) -> usize { self.iter.len() }
}

impl<'a> std::iter::FusedIterator for Inputs<'a> { }

/// Iterator over all diagram output node IDs, visited in qubit index order.
///
/// The iterator item type is `(`[`QubitId`]`, `[`NodeId`]`)`.
pub struct Outputs<'a> {
    iter: std::iter::Enumerate<std::slice::Iter<'a, NodeId>>
}

impl<'a> Iterator for Outputs<'a> {
    type Item = (NodeId, QubitId);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|(qid, nid)| (qid, *nid))
    }
}

impl<'a> DoubleEndedIterator for Outputs<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.next_back().map(|(qid, nid)| (qid, *nid))
    }
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
/// The iterator item type is `(`[`NodeId`]`, `[`NodeId`]`)`.
pub struct Wires<'a> {
    nb_group: Option<(NodeId, std::slice::Iter<'a, NodeId>)>,
    len: usize,
    seen: Vec<(NodeId, NodeId)>,
    iter: std::iter::Enumerate<std::slice::Iter<'a, Option<Vec<NodeId>>>>,
}

impl<'a> Iterator for Wires<'a> {
    type Item = (NodeId, NodeId);

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
                    let pair = (*left, *right);
                    let pair_rev = (*right, *left);
                    if self.seen.contains(&pair) {
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

impl<'a> DoubleEndedIterator for Wires<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        loop {
            if let Some((left, group_iter)) = self.nb_group.as_mut() {
                if let Some(right) = group_iter.next_back() {
                    let pair = (*left, *right);
                    let pair_rev = (*right, *left);
                    if self.seen.contains(&pair) {
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

impl<'a> ExactSizeIterator for Wires<'a> {
    fn len(&self) -> usize { self.len }
}

impl<'a> std::iter::FusedIterator for Wires<'a> { }

/// Iterator over all wires in a diagram, visited in mostly arbitrary order,
/// with node data attached.
///
/// Specifically, the iterator steps through wires as pairs of node IDs such
/// that those whose left IDs are smaller in value are visited before those
/// larger in value. In other wods, the left ID monotonically increases over the
/// course of iteration while right IDs are free to vary in arbitrary order.
/// Each wire is only visited once.
///
/// The iterator item type is `(`[`NodeData`]`, `[`NodeData`]`)`.
pub struct WiresData<'a> {
    dg: &'a Diagram,
    iter: Wires<'a>,
}

impl<'a> Iterator for WiresData<'a> {
    type Item = (NodeData<'a>, NodeData<'a>);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
            .map(|(l, r)| {
                (
                    (l, self.dg.get_node(l).unwrap()),
                    (r, self.dg.get_node(r).unwrap()),
                )
            })
    }
}

impl<'a> DoubleEndedIterator for WiresData<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.next_back()
            .map(|(l, r)| {
                (
                    (l, self.dg.get_node(l).unwrap()),
                    (r, self.dg.get_node(r).unwrap()),
                )
            })
    }
}

impl<'a> ExactSizeIterator for WiresData<'a> {
    fn len(&self) -> usize { self.iter.len() }
}

impl<'a> std::iter::FusedIterator for WiresData<'a> { }

/// Iterator over all wires connected to only interior (non-input/output) nodes
/// in a diagram, visited in mostly arbitrary order (see [`Wires`]).
///
/// The iterator item type is `(`[`NodeId`]`, `[`NodeId`]`)`.
pub struct WiresInner<'a> {
    dg: &'a Diagram,
    iter: Wires<'a>,
}

impl<'a> Iterator for WiresInner<'a> {
    type Item = (NodeId, NodeId);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(l, r)| {
            self.dg.get_node(*l).unwrap().is_generator()
                && self.dg.get_node(*r).unwrap().is_generator()
        })
    }
}

impl<'a> DoubleEndedIterator for WiresInner<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|(l, r)| {
            self.dg.get_node(*l).unwrap().is_generator()
                && self.dg.get_node(*r).unwrap().is_generator()
        })
    }
}

impl<'a> std::iter::FusedIterator for WiresInner<'a> { }

/// Iterator over all wires connected to only interior (non-input/output) nodes
/// in a diagram, visited in mostly arbitrary order (see [`WiresData`]), with
/// node data attached.
///
/// The iterator item type is `(`[`NodeData`]`, `[`NodeData`]`)`.
pub struct WiresDataInner<'a> {
    iter: WiresData<'a>
}

impl<'a> Iterator for WiresDataInner<'a> {
    type Item = (NodeData<'a>, NodeData<'a>);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|((_, l), (_, r))| {
            l.is_generator() && r.is_generator()
        })
    }
}

impl<'a> DoubleEndedIterator for WiresDataInner<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|((_, l), (_, r))| {
            l.is_generator() && r.is_generator()
        })
    }
}

impl<'a> std::iter::FusedIterator for WiresDataInner<'a> { }

/// Iterator over neighbor IDs of a node, visited in arbitrary order.
///
/// Node that this iterator will contain duplicate elements if there are
/// multiple wires going to the same neighbor; if the node is connected to
/// itself, each such connection will be counted only once.
///
/// The iterator item type is [`NodeId`].
pub struct NeighborIds<'a> {
    id: NodeId,
    switch: bool,
    iter: std::slice::Iter<'a, NodeId>,
}

impl<'a> Iterator for NeighborIds<'a> {
    type Item = NodeId;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|id| {
            if **id == self.id {
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
        .copied()
    }
}

impl<'a> DoubleEndedIterator for NeighborIds<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|id| {
            if **id == self.id {
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
        .copied()
    }
}

impl<'a> std::iter::FusedIterator for NeighborIds<'a> { }

/// Iterator over neighbors of a node, visited in arbitrary order.
///
/// Note that this iterator will contain duplicate elements if there are
/// multiple wires going to the same neighbor; if the node is connected to
/// itself, each such connection will be counted only once.
///
/// The iterator item type is `(`[`NodeId`]`, &`[`Node`]`)`.
pub struct Neighbors<'a> {
    dg: &'a Diagram,
    id: NodeId,
    switch: bool,
    iter: std::slice::Iter<'a, NodeId>,
}

impl<'a> Iterator for Neighbors<'a> {
    type Item = (NodeId, &'a Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|id| {
            if **id == self.id {
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
        .map(|id| (*id, self.dg.nodes[*id].as_ref().unwrap()))
    }
}

impl<'a> DoubleEndedIterator for Neighbors<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|id| {
            if **id == self.id {
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
        .map(|id| (*id, self.dg.nodes[*id].as_ref().unwrap()))
    }
}

impl<'a> std::iter::FusedIterator for Neighbors<'a> { }

/// Iterator over all interior (non-input/output) neighbors of a node, visited
/// in arbitrary order.
///
/// Note that this iterator will contain duplicate elements if there are
/// multiple wires going to the same neighbor; if the node is connected to
/// itself, each such connection will be counted only once.
///
/// The iterator item type is `(`[`NodeId`]`, &`[`Node`]`)`.
pub struct NeighborsInner<'a> {
    iter: Neighbors<'a>
}

impl<'a> Iterator for NeighborsInner<'a> {
    type Item = (NodeId, &'a Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, n)| n.is_generator())
    }
}

impl<'a> DoubleEndedIterator for NeighborsInner<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|(_, n)| n.is_generator())
    }
}

impl<'a> std::iter::FusedIterator for NeighborsInner<'a> { }

/// Create a [`Diagram`] using an abbreviated syntax.
///
/// The first block defines nodes with the syntax
/// ```text
/// <node_label> = <node_type> ( arg )
/// ```
/// where `<node_type>` is a variant or method of [`Node`], and `arg` is any
/// arguments to be passed to it. Note that `arg` may be left empty. Labels
/// given to nodes in this block are then used to define wire connections in the
/// second block with the syntax
/// ```text
/// <node1> -- <node2> [ -- <node3> ... ]
/// ```
/// All nodes' labels and corresponding IDs are returned in a
/// [`HashMap`][std::collections::HashMap]`<&'static `[`str`]`, `[`NodeId`]`>`.
///
/// The total return type is [`Result`]`<(`[`Diagram`]`,
/// `[`HashMap`][std::collections::HashMap]`<&'static `[`str`]`, `[`NodeId`]`>),
/// `[`crate::graph2::GraphError`]`>`.
///
/// The normal usage
/// ```
/// # use zx_calc::graph2::*;
///
/// # fn main() -> Result<(), GraphError> {
/// let mut diagram = Diagram::new();
///
/// let i = diagram.add_node(Node::Input);
/// let o = diagram.add_node(Node::Output);
/// let n1 = diagram.add_node(Node::z_pi());
/// let n2 = diagram.add_node(Node::x());
/// let n3 = diagram.add_node(Node::H((-1.0).into()));
/// diagram.add_wire(i, n1)?;
/// diagram.add_wire(o, n3)?;
/// diagram.add_wire(n1, n2)?;
/// diagram.add_wire(n1, n2)?;
/// diagram.add_wire(n2, n3)?;
/// # Ok(())
/// # }
/// ```
/// is equivalent to
/// ```
/// # use zx_calc::graph2::*;
/// use zx_calc::diagram;
///
/// diagram!(
///     nodes: {
///         i = input ( ),
///         o = output ( ),
///         n1 = z_pi ( ),
///         n2 = x ( ),
///         n3 = H (-1.0),
///     },
///     wires: {
///         i -- n1 -- n2 -- n3 -- o,
///         n1 -- n2,
///     },
/// );
/// ```
/// where in the macro usage, the diagram is returned alongside the hash map
/// ```text
/// {
///     "i" => 0,
///     "o" => 1,
///     "n1" => 2,
///     "n2" => 3,
///     "n3" => 4,
/// }
/// ```
#[macro_export]
macro_rules! diagram {
    (
        nodes: {
            $( $node_name:ident = $node_def:ident ( $( $arg:expr ),* $(,)? ) ),*
            $(,)?
        },
        wires: {
            $( $node1_name:ident $( -- $nodek_name:ident )+ ),* $(,)?
        } $(,)?
    ) => {
        {
            let mut diagram = $crate::graph2::Diagram::new();
            $(
                let $node_name =
                    diagram.add_node(
                        $crate::graph2::Node::$node_def($( ($arg).into() ),*)
                    );
            )*
            Ok(())
            $(
            .and_then(|_| {
                let mut last = $node1_name;
                Ok(())
                $(
                .and_then(|_| {
                    let res = diagram.add_wire(last, $nodek_name);
                    last = $nodek_name;
                    res
                })
                )+
            })
            )*
            .map(|_| {
                let nodes:
                    std::collections::HashMap<
                        &'static str,
                        $crate::graph2::NodeId
                    > =
                    [$( (stringify!($node_name), $node_name) ),*]
                    .into_iter()
                    .collect();
                (diagram, nodes)
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::phase::Phase;

    fn build_simple() -> Diagram {
        let mut diagram = Diagram::new();
        let z0 = diagram.add_node(Node::Z(Phase::new(1, 6))); // 0
        let z1 = diagram.add_node(Node::Z(Phase::new(1, 3))); // 1
        diagram.add_input_wire(z0).unwrap();                  // 2, (0, 2)
        diagram.add_input_wire(z0).unwrap();                  // 3, (0, 3)
        diagram.add_input_wire(z0).unwrap();                  // 4, (0, 4)
        let x0 = diagram.add_node(Node::X(Phase::new(7, 8))); // 5
        diagram.add_output_wire(z0).unwrap();                 // 6, (0, 6)
        diagram.add_output_wire(z1).unwrap();                 // 7, (1, 7)
        let h0 = diagram.add_node(Node::h());                 // 8
        diagram.add_wire(z0, z1).unwrap();                    // (0, 1)
        diagram.add_wire(z0, z1).unwrap();                    // (0, 1)
        diagram.add_wire(z0, h0).unwrap();                    // (0, 8)
        diagram.add_wire(h0, x0).unwrap();                    // (8, 5)
        diagram
    }

    #[test]
    fn counts() {
        let mut dg = build_simple();
        assert_eq!(dg.count_nodes(),   9);
        assert_eq!(dg.count_z(),       2);
        assert_eq!(dg.count_x(),       1);
        assert_eq!(dg.count_h(),       1);
        assert_eq!(dg.count_inputs(),  3);
        assert_eq!(dg.count_outputs(), 2);
        assert_eq!(dg.count_spiders(), 3);
        assert_eq!(dg.count_wires(),   9);
        dg.remove_node(3).unwrap();
        assert_eq!(dg.count_nodes(),   8);
        assert_eq!(dg.count_z(),       2);
        assert_eq!(dg.count_x(),       1);
        assert_eq!(dg.count_h(),       1);
        assert_eq!(dg.count_inputs(),  2);
        assert_eq!(dg.count_outputs(), 2);
        assert_eq!(dg.count_spiders(), 3);
        assert_eq!(dg.count_wires(),   8);
        dg.remove_node(8).unwrap();
        assert_eq!(dg.count_nodes(),   7);
        assert_eq!(dg.count_z(),       2);
        assert_eq!(dg.count_x(),       1);
        assert_eq!(dg.count_h(),       0);
        assert_eq!(dg.count_inputs(),  2);
        assert_eq!(dg.count_outputs(), 2);
        assert_eq!(dg.count_spiders(), 3);
        assert_eq!(dg.count_wires(),   6);
        dg.add_output_wire(1).unwrap();
        assert_eq!(dg.count_nodes(),   8);
        assert_eq!(dg.count_z(),       2);
        assert_eq!(dg.count_x(),       1);
        assert_eq!(dg.count_h(),       0);
        assert_eq!(dg.count_inputs(),  2);
        assert_eq!(dg.count_outputs(), 3);
        assert_eq!(dg.count_spiders(), 3);
        assert_eq!(dg.count_wires(),   7);

        let mut dg = Diagram::new();
        let z = dg.add_node(Node::z());
        dg.add_wire(z, z).unwrap();
        assert_eq!(dg.count_nodes(),   1);
        assert_eq!(dg.count_z(),       1);
        assert_eq!(dg.count_x(),       0);
        assert_eq!(dg.count_h(),       0);
        assert_eq!(dg.count_inputs(),  0);
        assert_eq!(dg.count_outputs(), 0);
        assert_eq!(dg.count_spiders(), 1);
        assert_eq!(dg.count_wires(),   1);
        dg.remove_node(z).unwrap();
        assert_eq!(dg.count_nodes(),   0);
        assert_eq!(dg.count_z(),       0);
        assert_eq!(dg.count_x(),       0);
        assert_eq!(dg.count_h(),       0);
        assert_eq!(dg.count_inputs(),  0);
        assert_eq!(dg.count_outputs(), 0);
        assert_eq!(dg.count_spiders(), 0);
        assert_eq!(dg.count_wires(),   0);
    }

    #[test]
    fn get_node() {
        let dg = build_simple();
        assert_eq!(dg.get_node(0), Some(&Node::Z(Phase::new(1, 6))));
        assert_eq!(dg.get_node(1), Some(&Node::Z(Phase::new(1, 3))));
        assert_eq!(dg.get_node(2), Some(&Node::Input));
        assert_eq!(dg.get_node(3), Some(&Node::Input));
        assert_eq!(dg.get_node(4), Some(&Node::Input));
        assert_eq!(dg.get_node(5), Some(&Node::X(Phase::new(7, 8))));
        assert_eq!(dg.get_node(6), Some(&Node::Output));
        assert_eq!(dg.get_node(7), Some(&Node::Output));
        assert_eq!(dg.get_node(8), Some(&Node::H((-1.0).into())));
        assert_eq!(dg.get_node(9), None);
    }

    #[test]
    fn arity() {
        let dg = build_simple();
        assert_eq!(dg.arity(0), Some(7));
        assert_eq!(dg.arity(1), Some(3));
        assert_eq!(dg.arity(2), Some(1));
        assert_eq!(dg.arity(3), Some(1));
        assert_eq!(dg.arity(4), Some(1));
        assert_eq!(dg.arity(5), Some(1));
        assert_eq!(dg.arity(6), Some(1));
        assert_eq!(dg.arity(7), Some(1));
        assert_eq!(dg.arity(8), Some(2));
        assert_eq!(dg.arity(9), None);
    }

    #[test]
    fn is_connected() {
        let mut dg = build_simple();
        dg.add_node(Node::x()); // 9
        assert_eq!(dg.is_connected(0), Some(true));
        assert_eq!(dg.is_connected(1), Some(true));
        assert_eq!(dg.is_connected(2), Some(true));
        assert_eq!(dg.is_connected(3), Some(true));
        assert_eq!(dg.is_connected(4), Some(true));
        assert_eq!(dg.is_connected(5), Some(true));
        assert_eq!(dg.is_connected(6), Some(true));
        assert_eq!(dg.is_connected(7), Some(true));
        assert_eq!(dg.is_connected(8), Some(true));
        assert_eq!(dg.is_connected(9), Some(false));
        assert_eq!(dg.is_connected(10), None);
    }

    #[test]
    fn mutual_arity() {
        let mut dg = build_simple();
        let z = dg.add_node(Node::z()); // 9
        dg.add_wire(z, z).unwrap();
        assert_eq!(dg.mutual_arity(0, 2), Some(1));
        assert_eq!(dg.mutual_arity(0, 3), Some(1));
        assert_eq!(dg.mutual_arity(0, 4), Some(1));
        assert_eq!(dg.mutual_arity(0, 6), Some(1));
        assert_eq!(dg.mutual_arity(1, 7), Some(1));
        assert_eq!(dg.mutual_arity(0, 1), Some(2));
        assert_eq!(dg.mutual_arity(0, 8), Some(1));
        assert_eq!(dg.mutual_arity(8, 5), Some(1));
        assert_eq!(dg.mutual_arity(0, 2), Some(1));
        assert_eq!(dg.mutual_arity(1, 2), Some(0));
        assert_eq!(dg.mutual_arity(1, 8), Some(0));
        assert_eq!(dg.mutual_arity(2, 3), Some(0));
        assert_eq!(dg.mutual_arity(z, z), Some(1));
        assert_eq!(dg.mutual_arity(0, 10), None);
        assert_eq!(dg.mutual_arity(10, 0), None);
        assert_eq!(dg.mutual_arity(10, 11), None);
    }

    #[test]
    fn remove_node() {
        let mut dg = build_simple();
        let z = dg.add_node(Node::z()); // 9
        dg.add_wire(z, z).unwrap();
        assert_eq!(dg.remove_node(3).unwrap(), Node::Input);
        assert_eq!(dg.remove_node(8).unwrap(), Node::H((-1.0).into()));
        assert_eq!(dg.remove_node(6).unwrap(), Node::Output);
        assert_eq!(dg.remove_node(z).unwrap(), Node::z());
        assert!(dg.remove_node(3).is_err());
        assert!(dg.remove_node(9).is_err());
        assert_eq!(dg.add_node(Node::h()), 9);
    }

    #[test]
    fn add_wire() {
        let mut dg = build_simple();
        assert!(dg.add_wire(0, 5).is_ok());
        assert!(dg.add_wire(1, 8).is_ok());
        assert!(dg.add_wire(0, 2).is_err());
        assert!(dg.add_wire(1, 6).is_err());
        assert!(dg.add_wire(3, 7).is_err());
        assert!(dg.add_wire(0, 9).is_err());
        assert!(dg.add_wire(9, 10).is_err());
    }

    #[test]
    fn add_input_wire() {
        let mut dg = build_simple();
        assert_eq!(dg.add_input_wire(5).unwrap(), 9);
        assert_eq!(dg.add_input_wire(0).unwrap(), 10);
        dg.remove_node(2).unwrap();
        assert_eq!(dg.add_input_wire(0).unwrap(), 2);
        assert_eq!(dg.count_inputs(), 5);
        assert!(dg.add_input_wire(11).is_err());
        assert!(dg.add_input_wire(2).is_err());
    }

    #[test]
    fn add_output_wire() {
        let mut dg = build_simple();
        assert_eq!(dg.add_output_wire(5).unwrap(), 9);
        assert_eq!(dg.add_output_wire(0).unwrap(), 10);
        dg.remove_node(6).unwrap();
        assert_eq!(dg.add_output_wire(0).unwrap(), 6);
        assert_eq!(dg.count_outputs(), 4);
        assert!(dg.add_output_wire(11).is_err());
        assert!(dg.add_output_wire(2).is_err());
    }

    #[test]
    fn remove_wires() {
        let mut dg = build_simple();
        assert_eq!(dg.count_wires(), 9);
        assert_eq!(dg.mutual_arity(0, 1), Some(2));
        (0..5).for_each(|_| { dg.add_wire(0, 1).unwrap(); });
        assert_eq!(dg.count_wires(), 14);
        assert_eq!(dg.mutual_arity(0, 1), Some(7));
        assert_eq!(dg.remove_wires(0, 1, Some(2)).unwrap(), 2);
        assert_eq!(dg.count_wires(), 12);
        assert_eq!(dg.mutual_arity(0, 1), Some(5));
        assert_eq!(dg.remove_wires(0, 1, None).unwrap(), 5);
        assert_eq!(dg.count_wires(), 7);
        assert_eq!(dg.mutual_arity(0, 1), Some(0));
        (0..3).for_each(|_| { dg.add_wire(0, 0).unwrap(); });
        assert_eq!(dg.count_wires(), 10);
        assert_eq!(dg.mutual_arity(0, 0), Some(3));
        assert_eq!(dg.remove_wires(0, 0, None).unwrap(), 3);
        assert_eq!(dg.count_wires(), 7);
        assert_eq!(dg.mutual_arity(0, 0), Some(0));
    }

    #[test]
    fn apply_state() {
        let mut dg = build_simple();

        assert_eq!(dg.count_nodes(),   9);
        assert_eq!(dg.count_inputs(),  3);
        assert_eq!(dg.count_outputs(), 2);
        assert_eq!(dg.count_z(),       2);
        assert_eq!(dg.count_x(),       1);

        dg.apply_state_input(1, Spider::Z(Phase::pi())).unwrap();
        assert!(dg.apply_state_input(3, Spider::Z(Phase::zero())).is_err());
        assert!(dg.apply_state(0, Spider::Z(Phase::zero())).is_err());
        assert!(dg.apply_state(10, Spider::Z(Phase::zero())).is_err());

        assert_eq!(dg.count_nodes(),   9);
        assert_eq!(dg.count_inputs(),  3);
        assert_eq!(dg.count_outputs(), 2);
        assert_eq!(dg.count_z(),       3);
        assert_eq!(dg.count_x(),       1);

        dg.apply_effect_output(1, Spider::X(Phase::new(1, 3))).unwrap();
        assert!(dg.apply_effect_output(2, Spider::X(Phase::zero())).is_err());

        assert_eq!(dg.count_nodes(),   9);
        assert_eq!(dg.count_inputs(),  3);
        assert_eq!(dg.count_outputs(), 2);
        assert_eq!(dg.count_z(),       3);
        assert_eq!(dg.count_x(),       2);
    }

    #[test]
    fn apply_bell() {
        let mut dg = build_simple();

        assert_eq!(dg.count_nodes(),   9);
        assert_eq!(dg.count_inputs(),  3);
        assert_eq!(dg.count_outputs(), 2);
        assert_eq!(dg.count_z(),       2);
        assert_eq!(dg.count_x(),       1);
        assert_eq!(dg.count_wires(),   9);

        assert!(dg.apply_bell_input(0, 2, None).unwrap().is_none());

        assert_eq!(dg.count_nodes(),   9);
        assert_eq!(dg.count_inputs(),  3);
        assert_eq!(dg.count_outputs(), 2);
        assert_eq!(dg.count_z(),       4);
        assert_eq!(dg.count_x(),       1);
        assert_eq!(dg.count_wires(),   10);

        assert_eq!(
            dg.apply_bell_output(0, 1, Some(Spider::X(Phase::pi()))).unwrap(),
            Some(9),
        );

        assert_eq!(dg.count_nodes(),   10);
        assert_eq!(dg.count_inputs(),  3);
        assert_eq!(dg.count_outputs(), 2);
        assert_eq!(dg.count_z(),       6);
        assert_eq!(dg.count_x(),       2);
        assert_eq!(dg.count_wires(),   12);

        assert!(dg.apply_bell(0, 1, None).is_err());
        assert!(dg.apply_bell(20, 21, None).is_err());
        assert!(dg.apply_bell(3, 21, None).is_err());
        assert!(dg.apply_bell(2, 4, None).is_err());
    }

    #[test]
    fn iter_nodes() {
        let mut dg = build_simple();

        let nodes: Vec<(NodeId, Node)> =
            dg.nodes().map(|(id, n)| (id, *n)).collect();
        let nodes_expected: Vec<(NodeId, Node)> =
            vec![
                (0, Node::Z(Phase::new(1, 6))),
                (1, Node::Z(Phase::new(1, 3))),
                (2, Node::Input),
                (3, Node::Input),
                (4, Node::Input),
                (5, Node::X(Phase::new(7, 8))),
                (6, Node::Output),
                (7, Node::Output),
                (8, Node::H((-1.0).into())),
            ];
        assert_eq!(nodes, nodes_expected);

        dg.remove_node(3).unwrap();
        dg.remove_node(8).unwrap();

        let nodes: Vec<(NodeId, Node)> =
            dg.nodes().map(|(id, n)| (id, *n)).collect();
        let nodes_expected: Vec<(NodeId, Node)> =
            vec![
                (0, Node::Z(Phase::new(1, 6))),
                (1, Node::Z(Phase::new(1, 3))),
                (2, Node::Input),
                (4, Node::Input),
                (5, Node::X(Phase::new(7, 8))),
                (6, Node::Output),
                (7, Node::Output),
            ];
        assert_eq!(nodes, nodes_expected);
    }

    #[test]
    fn iter_nodes_inner() {
        let mut dg = build_simple();

        let nodes_inner: Vec<(NodeId, Node)> =
            dg.nodes_inner().map(|(id, n)| (id, *n)).collect();
        let nodes_inner_expected: Vec<(NodeId, Node)> =
            vec![
                (0, Node::Z(Phase::new(1, 6))),
                (1, Node::Z(Phase::new(1, 3))),
                (5, Node::X(Phase::new(7, 8))),
                (8, Node::H((-1.0).into())),
            ];
        assert_eq!(nodes_inner, nodes_inner_expected);

        dg.remove_node(3).unwrap();
        dg.remove_node(8).unwrap();

        let nodes_inner: Vec<(NodeId, Node)> =
            dg.nodes_inner().map(|(id, n)| (id, *n)).collect();
        let nodes_inner_expected: Vec<(NodeId, Node)> =
            vec![
                (0, Node::Z(Phase::new(1, 6))),
                (1, Node::Z(Phase::new(1, 3))),
                (5, Node::X(Phase::new(7, 8))),
            ];
        assert_eq!(nodes_inner, nodes_inner_expected);
    }

    #[test]
    fn iter_inputs() {
        let mut dg = build_simple();

        let inputs: Vec<(QubitId, NodeId)> = dg.inputs().collect();
        let inputs_expected: Vec<(QubitId, NodeId)> =
            vec![ (0, 2), (1, 3), (2, 4) ];
        assert_eq!(inputs, inputs_expected);

        dg.remove_node(3).unwrap();
        dg.remove_node(8).unwrap();

        let inputs: Vec<(QubitId, NodeId)> = dg.inputs().collect();
        let inputs_expected: Vec<(QubitId, NodeId)> =
            vec![ (0, 2), (1, 4) ];
        assert_eq!(inputs, inputs_expected);
    }

    #[test]
    fn iter_outputs() {
        let mut dg = build_simple();

        let outputs: Vec<(QubitId, NodeId)> = dg.outputs().collect();
        let outputs_expected: Vec<(QubitId, NodeId)> =
            vec![ (0, 6), (1, 7) ];
        assert_eq!(outputs, outputs_expected);

        dg.remove_node(3).unwrap();
        dg.remove_node(8).unwrap();
        dg.remove_node(6).unwrap();

        let outputs: Vec<(QubitId, NodeId)> = dg.outputs().collect();
        let outputs_expected: Vec<(QubitId, NodeId)> =
            vec![ (0, 7) ];
        assert_eq!(outputs, outputs_expected);
    }

    #[test]
    fn iter_wires() {
        let mut dg = build_simple();

        let wires: Vec<(NodeId, NodeId)> = dg.wires().collect();
        let wires_expected: Vec<(NodeId, NodeId)> =
            vec![
                (0, 1),
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 4),
                (0, 6),
                (0, 8),
                (1, 7),
                (5, 8),
            ];
        assert_eq!(wires.len(), wires_expected.len());
        assert!(
            wires.into_iter()
                .all(|pair| {
                    wires_expected.contains(&pair)
                        || wires_expected.contains(&rev(pair))
                })
        );

        dg.remove_node(3).unwrap();
        dg.remove_node(8).unwrap();
        dg.remove_node(6).unwrap();

        let wires: Vec<(NodeId, NodeId)> = dg.wires().collect();
        let wires_expected: Vec<(NodeId, NodeId)> =
            vec![
                (0, 1),
                (0, 1),
                (0, 2),
                (0, 4),
                (1, 7),
            ];
        assert_eq!(wires.len(), wires_expected.len());
        assert!(
            wires.into_iter()
                .all(|pair| {
                    wires_expected.contains(&pair)
                        || wires_expected.contains(&rev(pair))
                })
        );
    }

    #[test]
    fn iter_wires_inner() {
        fn rev<T, U>(pair: (T, U)) -> (U, T) { (pair.1, pair.0) }

        let mut dg = build_simple();

        let wires_inner: Vec<(NodeId, NodeId)> = dg.wires_inner().collect();
        let wires_inner_expected: Vec<(NodeId, NodeId)> =
            vec![
                (0, 1),
                (0, 1),
                (0, 8),
                (5, 8),
            ];
        assert_eq!(wires_inner.len(), wires_inner_expected.len());
        assert!(
            wires_inner.into_iter()
                .all(|pair| {
                    wires_inner_expected.contains(&pair)
                        || wires_inner_expected.contains(&rev(pair))
                })
        );

        dg.remove_node(3).unwrap();
        dg.remove_node(8).unwrap();
        dg.remove_node(6).unwrap();

        let wires_inner: Vec<(NodeId, NodeId)> = dg.wires_inner().collect();
        let wires_inner_expected: Vec<(NodeId, NodeId)> =
            vec![
                (0, 1),
                (0, 1),
            ];
        assert_eq!(wires_inner.len(), wires_inner_expected.len());
        assert!(
            wires_inner.into_iter()
                .all(|pair| {
                    wires_inner_expected.contains(&pair)
                        || wires_inner_expected.contains(&rev(pair))
                })
        );
    }

    #[test]
    fn iter_neighbors_of() {
        let mut dg = build_simple();

        let neighbors_of: Vec<(NodeId, Node)> =
            dg.neighbors_of(0).unwrap().map(|(id, n)| (id, *n)).collect();
        let neighbors_of_expected: Vec<(NodeId, Node)> =
            vec![
                (2, Node::Input),
                (3, Node::Input),
                (4, Node::Input),
                (6, Node::Output),
                (1, Node::Z(Phase::new(1, 3))),
                (1, Node::Z(Phase::new(1, 3))),
                (8, Node::H((-1.0).into())),
            ];
        assert_eq!(neighbors_of.len(), neighbors_of_expected.len());
        assert!(
            neighbors_of.into_iter()
                .all(|x| neighbors_of_expected.contains(&x))
        );

        dg.remove_node(3).unwrap();
        dg.remove_node(8).unwrap();
        dg.remove_node(6).unwrap();

        let neighbors_of: Vec<(NodeId, Node)> =
            dg.neighbors_of(0).unwrap().map(|(id, n)| (id, *n)).collect();
        let neighbors_of_expected: Vec<(NodeId, Node)> =
            vec![
                (2, Node::Input),
                (4, Node::Input),
                (1, Node::Z(Phase::new(1, 3))),
                (1, Node::Z(Phase::new(1, 3))),
            ];
        assert_eq!(neighbors_of.len(), neighbors_of_expected.len());
        assert!(
            neighbors_of.into_iter()
                .all(|x| neighbors_of_expected.contains(&x))
        );

        let z = dg.add_node(Node::z_pi());
        dg.add_wire(z, z).unwrap();
        let neighbors_of: Vec<(NodeId, Node)> =
            dg.neighbors_of(z).unwrap().map(|(id, n)| (id, *n)).collect();
        let neighbors_of_expected: Vec<(NodeId, Node)> =
            vec![ (z, Node::z_pi()) ];
        assert_eq!(neighbors_of.len(), neighbors_of_expected.len());
        assert!(
            neighbors_of.into_iter()
                .all(|x| neighbors_of_expected.contains(&x))
        );
    }

    #[test]
    fn iter_neighbors_of_inner() {
        let mut dg = build_simple();

        let neighbors_of_inner: Vec<(NodeId, Node)> =
            dg.neighbors_of_inner(0).unwrap().map(|(id, n)| (id, *n)).collect();
        let neighbors_of_inner_expected: Vec<(NodeId, Node)> =
            vec![
                (1, Node::Z(Phase::new(1, 3))),
                (1, Node::Z(Phase::new(1, 3))),
                (8, Node::H((-1.0).into())),
            ];
        assert_eq!(neighbors_of_inner.len(), neighbors_of_inner_expected.len());
        assert!(
            neighbors_of_inner.into_iter()
                .all(|x| neighbors_of_inner_expected.contains(&x))
        );

        dg.remove_node(3).unwrap();
        dg.remove_node(8).unwrap();
        dg.remove_node(6).unwrap();

        let neighbors_of_inner: Vec<(NodeId, Node)> =
            dg.neighbors_of_inner(0).unwrap().map(|(id, n)| (id, *n)).collect();
        let neighbors_of_inner_expected: Vec<(NodeId, Node)> =
            vec![
                (1, Node::Z(Phase::new(1, 3))),
                (1, Node::Z(Phase::new(1, 3))),
            ];
        assert_eq!(neighbors_of_inner.len(), neighbors_of_inner_expected.len());
        assert!(
            neighbors_of_inner.into_iter()
                .all(|x| neighbors_of_inner_expected.contains(&x))
        );
    }

    #[test]
    fn adjoint() {
        let mut dg = build_simple();
        dg.adjoint_mut();

        let nodes: Vec<(NodeId, Node)> =
            dg.nodes().map(|(id, n)| (id, *n)).collect();
        let nodes_expected: Vec<(NodeId, Node)> =
            vec![
                (0, Node::Z(Phase::new(5, 6))),
                (1, Node::Z(Phase::new(2, 3))),
                (2, Node::Output),
                (3, Node::Output),
                (4, Node::Output),
                (5, Node::X(Phase::new(1, 8))),
                (6, Node::Input),
                (7, Node::Input),
                (8, Node::H((-1.0).into())),
            ];
        assert_eq!(nodes, nodes_expected);

        let wires: Vec<(NodeId, NodeId)> = dg.wires().collect();
        let wires_expected: Vec<(NodeId, NodeId)> =
            vec![
                (0, 1),
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 4),
                (0, 6),
                (0, 8),
                (1, 7),
                (5, 8),
            ];
        assert_eq!(wires.len(), wires_expected.len());
        assert!(
            wires.into_iter()
                .all(|pair| {
                    wires_expected.contains(&pair)
                        || wires_expected.contains(&rev(pair))
                })
        );
    }

    #[test]
    fn append() {
        let mut dg = build_simple();
        dg.append(build_simple());

        let nodes: Vec<(NodeId, Node)> =
            dg.nodes().map(|(id, n)| (id, *n)).collect();
        let nodes_expected: Vec<(NodeId, Node)> =
            vec![
                ( 0, Node::Z(Phase::new(1, 6))),
                ( 1, Node::Z(Phase::new(1, 3))),
                ( 2, Node::Input),
                ( 3, Node::Input),
                ( 4, Node::Input),
                ( 5, Node::X(Phase::new(7, 8))),
                ( 6, Node::Output),
                ( 7, Node::Output),
                ( 8, Node::H((-1.0).into())),
                ( 9, Node::Z(Phase::new(1, 6))),
                (10, Node::Z(Phase::new(1, 3))),
                (11, Node::Input),
                (12, Node::Input),
                (13, Node::Input),
                (14, Node::X(Phase::new(7, 8))),
                (15, Node::Output),
                (16, Node::Output),
                (17, Node::H((-1.0).into())),
            ];
        assert_eq!(nodes, nodes_expected);

        let wires: Vec<(NodeId, NodeId)> = dg.wires().collect();
        let wires_expected: Vec<(NodeId, NodeId)> =
            vec![
                ( 0,  1),
                ( 0,  1),
                ( 0,  2),
                ( 0,  3),
                ( 0,  4),
                ( 0,  6),
                ( 0,  8),
                ( 1,  7),
                ( 5,  8),
                ( 9, 10),
                ( 9, 10),
                ( 9, 11),
                ( 9, 12),
                ( 9, 13),
                ( 9, 15),
                ( 9, 17),
                (10, 16),
                (14, 17),
            ];
        assert_eq!(wires.len(), wires_expected.len());
        assert!(
            wires.into_iter()
                .all(|pair| {
                    wires_expected.contains(&pair)
                        || wires_expected.contains(&rev(pair))
                })
        );
    }

    #[test]
    fn tensor() {
        let mut dg = build_simple();
        dg.tensor_with(build_simple());

        let nodes: Vec<(NodeId, Node)> =
            dg.nodes().map(|(id, n)| (id, *n)).collect();
        let nodes_expected: Vec<(NodeId, Node)> =
            vec![
                ( 0, Node::Z(Phase::new(1, 6))),
                ( 1, Node::Z(Phase::new(1, 3))),
                ( 2, Node::Input),
                ( 3, Node::Input),
                ( 4, Node::Input),
                ( 5, Node::X(Phase::new(7, 8))),
                ( 6, Node::Output),
                ( 7, Node::Output),
                ( 8, Node::H((-1.0).into())),
                ( 9, Node::Z(Phase::new(1, 6))),
                (10, Node::Z(Phase::new(1, 3))),
                (11, Node::Input),
                (12, Node::Input),
                (13, Node::Input),
                (14, Node::X(Phase::new(7, 8))),
                (15, Node::Output),
                (16, Node::Output),
                (17, Node::H((-1.0).into())),
            ];
        assert_eq!(nodes, nodes_expected);

        let wires: Vec<(NodeId, NodeId)> = dg.wires().collect();
        let wires_expected: Vec<(NodeId, NodeId)> =
            vec![
                ( 0,  1),
                ( 0,  1),
                ( 0,  2),
                ( 0,  3),
                ( 0,  4),
                ( 0,  6),
                ( 0,  8),
                ( 1,  7),
                ( 5,  8),
                ( 9, 10),
                ( 9, 10),
                ( 9, 11),
                ( 9, 12),
                ( 9, 13),
                ( 9, 15),
                ( 9, 17),
                (10, 16),
                (14, 17),
            ];
        assert_eq!(wires.len(), wires_expected.len());
        assert!(
            wires.into_iter()
                .all(|pair| {
                    wires_expected.contains(&pair)
                        || wires_expected.contains(&rev(pair))
                })
        );
    }

    #[test]
    fn compose() {
        let mut d0 = build_simple();
        let mut dclone = d0.clone();
        let d1 = build_simple();
        assert!(dclone.compose_with(d1).is_err());
        assert_eq!(d0.nodes, dclone.nodes);
        assert_eq!(d0.node_count, dclone.node_count);
        assert_eq!(d0.wires, dclone.wires);
        assert_eq!(d0.wire_count, dclone.wire_count);
        assert_eq!(d0.inputs, dclone.inputs);
        assert_eq!(d0.outputs, dclone.outputs);
        assert_eq!(d0.free, dclone.free);

        let mut d1 = build_simple();
        d1.add_output_wire(1).unwrap();
        assert_eq!(d0.compose_with(d1).unwrap(), 9);

        let nodes: Vec<(NodeId, Node)> =
            d0.nodes().map(|(id, n)| (id, *n)).collect();
        let nodes_expected: Vec<(NodeId, Node)> =
            vec![
                ( 0, Node::Z(Phase::new(1, 6))),
                ( 1, Node::Z(Phase::new(1, 3))),
                ( 5, Node::X(Phase::new(7, 8))),
                ( 6, Node::Output),
                ( 7, Node::Output),
                ( 8, Node::H((-1.0).into())),
                ( 9, Node::Z(Phase::new(1, 6))),
                (10, Node::Z(Phase::new(1, 3))),
                (11, Node::Input),
                (12, Node::Input),
                (13, Node::Input),
                (14, Node::X(Phase::new(7, 8))),
                (17, Node::H((-1.0).into())),
            ];
        assert_eq!(nodes, nodes_expected);

        let wires: Vec<(NodeId, NodeId)> = d0.wires().collect();
        let wires_expected: Vec<(NodeId, NodeId)> =
            vec![
                ( 0,  1),
                ( 0,  1),
                ( 0,  6),
                ( 0,  8),
                ( 0,  9),
                ( 0, 10),
                ( 0, 10),
                ( 1,  7),
                ( 5,  8),
                ( 9, 10),
                ( 9, 10),
                ( 9, 11),
                ( 9, 12),
                ( 9, 13),
                ( 9, 17),
                (14, 17),
            ];
        assert_eq!(wires.len(), wires_expected.len());
        assert!(
            wires.into_iter()
                .all(|pair| {
                    wires_expected.contains(&pair)
                        || wires_expected.contains(&rev(pair))
                })
        );
    }

    #[test]
    fn compose_rev() {
        let mut d0 = build_simple();
        let mut dclone = d0.clone();
        let d1 = build_simple();
        assert!(dclone.compose_with_rev(d1).is_err());
        assert_eq!(d0.nodes, dclone.nodes);
        assert_eq!(d0.node_count, dclone.node_count);
        assert_eq!(d0.wires, dclone.wires);
        assert_eq!(d0.wire_count, dclone.wire_count);
        assert_eq!(d0.inputs, dclone.inputs);
        assert_eq!(d0.outputs, dclone.outputs);
        assert_eq!(d0.free, dclone.free);

        let d1 = build_simple();
        d0.add_output_wire(1).unwrap();
        assert_eq!(d0.compose_with_rev(d1).unwrap(), 10);

        let nodes: Vec<(NodeId, Node)> =
            d0.nodes().map(|(id, n)| (id, *n)).collect();
        let nodes_expected: Vec<(NodeId, Node)> =
            vec![
                ( 0, Node::Z(Phase::new(1, 6))),
                ( 1, Node::Z(Phase::new(1, 3))),
                ( 2, Node::Input),
                ( 3, Node::Input),
                ( 4, Node::Input),
                ( 5, Node::X(Phase::new(7, 8))),
                ( 8, Node::H((-1.0).into())),
                (10, Node::Z(Phase::new(1, 6))),
                (11, Node::Z(Phase::new(1, 3))),
                (15, Node::X(Phase::new(7, 8))),
                (16, Node::Output),
                (17, Node::Output),
                (18, Node::H((-1.0).into())),
            ];
        assert_eq!(nodes, nodes_expected);

        let wires: Vec<(NodeId, NodeId)> = d0.wires().collect();
        let wires_expected: Vec<(NodeId, NodeId)> =
            vec![
                ( 0,  1),
                ( 0,  1),
                ( 0,  2),
                ( 0,  3),
                ( 0,  4),
                ( 0,  8),
                ( 0, 10),
                ( 1, 10),
                ( 1, 10),
                ( 5,  8),
                (10, 11),
                (10, 11),
                (10, 16),
                (10, 18),
                (11, 17),
                (15, 18),
            ];
        println!("{:?}\n{:?}", wires, wires_expected);
        assert_eq!(wires.len(), wires_expected.len());
        assert!(
            wires.into_iter()
                .all(|pair| {
                    wires_expected.contains(&pair)
                        || wires_expected.contains(&rev(pair))
                })
        );
    }

    #[test]
    fn macro_build() -> GraphResult<()> {
        use std::collections::HashMap;
        let (diagram, ids): (Diagram, HashMap<&'static str, NodeId>) =
            diagram!(
                nodes: {
                    z0 = Z (Phase::new(1, 6)),
                    z1 = Z (Phase::new(1, 3)),
                    i0 = input (),
                    i1 = input (),
                    i2 = input (),
                    x0 = X (Phase::new(7, 8)),
                    o0 = output (),
                    o1 = output (),
                    h0 = h (),
                },
                wires: {
                    i0 -- z0,
                    i1 -- z0,
                    i2 -- z0 -- o0,
                    z1 -- o1,
                    z0 -- z1 -- z0,
                    z0 -- h0 -- x0,
                },
            )?;
        assert_eq!(ids.len(), 9);
        assert_eq!(ids["z0"], 0);
        assert_eq!(ids["z1"], 1);
        assert_eq!(ids["i0"], 2);
        assert_eq!(ids["i1"], 3);
        assert_eq!(ids["i2"], 4);
        assert_eq!(ids["x0"], 5);
        assert_eq!(ids["o0"], 6);
        assert_eq!(ids["o1"], 7);
        assert_eq!(ids["h0"], 8);

        let nodes: Vec<(NodeId, Node)> =
            diagram.nodes().map(|(id, n)| (id, *n)).collect();
        let wires: Vec<(NodeId, NodeId)> = diagram.wires().collect();

        let diagram_expected = build_simple();
        let nodes_expected: Vec<(NodeId, Node)> =
            diagram_expected.nodes().map(|(id, n)| (id, *n)).collect();
        let wires_expected: Vec<(NodeId, NodeId)> =
            diagram_expected.wires().collect();

        assert_eq!(nodes, nodes_expected);
        assert_eq!(wires.len(), wires_expected.len());
        assert!(
            wires.into_iter()
                .all(|(l, r)| {
                    wires_expected.contains(&(l, r))
                        || wires_expected.contains(&(r, l))
                })
        );

        Ok(())
    }

}

