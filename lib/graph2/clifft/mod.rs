//! Diagram variants specialized to Clifford+*T* quantum circuits.
//!
#![cfg_attr(
    feature = "doc-images",
    cfg_attr(
        all(),
        doc = ::embed_doc_image::embed_image!("hadamard_wire", "assets/hadamard_wire.svg"),
        doc = ::embed_doc_image::embed_image!("star_def", "assets/star_def.svg"),
        doc = ::embed_doc_image::embed_image!("star_wire", "assets/star_wire.svg"),
        doc = ::embed_doc_image::embed_image!("auto_wire_rules", "assets/auto_wire_rules.svg"),
    )
)]
#![cfg_attr(
    not(feature = "doc-images"),
    doc = "**Doc images are not enabled**. Compile with feature `doc-images` and Rust version >= 1.54 to enable."
)]
//!
//! Here, we restrict to only the "graph-like" ZX-calculus (i.e. only Z-spiders
//! and empty, binary H-boxes) with phases equal to an integer multiple of
//! *π*/4. Additionally, this module's tools are implemented with a greater
//! focus on performance than explanation, as with the rest of this crate.
//! Namely, the H-box is eliminated as an explicit generator, being replaced by
//! the "Hadamard wire", which is drawn as a dashed blue line.
//!
//! ![hadamard_wire][hadamard_wire]
//!
//! Further, the star generator is introduced to enable concise representations
//! of multi-control gates
//!
//! ![star_def][star_def]
//!
//! which is, like H-boxes, then replaced with a "star wire" as a dashed orange
//! line.
//!
//! ![star_wire][star_wire]
//!
//! One will likely want to deal with this diagram implementation mostly through
//! the high-level [`Circuit`][crate::circuit::CircuitDiagram] interface, since
//! many common circuit patterns become obfuscated when expressed under these
//! constraints.
//!
//! For more information on the above, see the paper whose methods are
//! implemented in this module, [arXiv:2307.01803][speedy-contraction].
//!
//! [speedy-contraction]: https://arxiv.org/abs/2307.01803

use std::{ /*collections::VecDeque,*/ fs, io::Write, path::Path };
use crate::{
    graph2::{
        GraphError,
        GraphResult,
        Diagram,
        DiagramData,
        IONodeId,
        NodeData,
        NodeId,
        QubitId,
        // WireData,
        // WireStore,
    },
};

use GraphError::*;

pub(crate) mod complex;
pub use complex::*;

pub(crate) mod tphase;
pub use tphase::*;

pub(crate) mod ctnode;
pub use ctnode::*;

pub(crate) mod ctwire;
pub use ctwire::*;

/// Marker type for diagrams in the Clifford+*T*-restricted ZX-calculus.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct CT;

impl DiagramData for CT {
    type Node = CTNode;
    type Wire = CTWire;
    type Scalar = Complex;
}

impl Diagram<CT> {
    /// Return the number of spiders.
    pub fn count_z(&self) -> usize {
        self.nodes.iter().filter(|mb_n| mb_n.is_some_and(|n| n.is_z())).count()
    }

    /// Return the number of Hadamard wires.
    pub fn count_h(&self) -> usize {
        self.wires().filter(|(_, w)| w.is_h()).count()
    }

    /// Return the number of star wires.
    pub fn count_s(&self) -> usize {
        self.wires().filter(|(_, w)| w.is_s()).count()
    }

    /// Add a Z-spider to the diagram and return its ID.
    ///
    /// `phase` defaults to 0.
    pub fn add_z(&mut self, phase: Option<TPhase>) -> NodeId {
        self.add_node(CTNode::Z(phase.unwrap_or(TPhase::T0)))
    }

    /// Get the number of empty wires attached to a node, if it exists.
    pub fn arity_e(&self, id: NodeId) -> Option<usize> {
        self.wires.get(id)
            .and_then(|mb_nnb| mb_nnb.as_ref())
            .map(|nnb| nnb.iter().filter(|w| w.is_e()).count())
    }

    /// Get the number of Hadamard wires attached to a node, if it exists.
    pub fn arity_h(&self, id: NodeId) -> Option<usize> {
        self.wires.get(id)
            .and_then(|mb_nnb| mb_nnb.as_ref())
            .map(|nnb| nnb.iter().filter(|w| w.is_h()).count())
    }

    /// Get the number of star wires attached to a node, if it exists.
    pub fn arity_s(&self, id: NodeId) -> Option<usize> {
        self.wires.get(id)
            .and_then(|mb_nnb| mb_nnb.as_ref())
            .map(|nnb| nnb.iter().filter(|w| w.is_s()).count())
    }

    /// Get the number of empty wires connecting two nodes, if they both exist.
    pub fn mutual_arity_e(&self, a: NodeId, b: NodeId) -> Option<usize> {
        self.has_node(b).then_some(())?;
        self.wires.get(a)
            .and_then(|mb_nnb| {
                mb_nnb.as_ref()
                    .map(|nnb| {
                        let ma =
                            nnb.iter().filter(|wire| wire.has_e_id(b)).count();
                        if a == b { ma / 2 } else { ma }
                    })
            })
    }

    /// Get the number of Hadamard wires connecting two nodes, if they both
    /// exist.
    pub fn mutual_arity_h(&self, a: NodeId, b: NodeId) -> Option<usize> {
        self.has_node(b).then_some(())?;
        self.wires.get(a)
            .and_then(|mb_nnb| {
                mb_nnb.as_ref()
                    .map(|nnb| {
                        let ma =
                            nnb.iter().filter(|wire| wire.has_h_id(b)).count();
                        if a == b { ma / 2 } else { ma }
                    })
            })
    }

    /// Get the number of star wires connecting two nodes, if they both exist.
    pub fn mutual_arity_s(&self, a: NodeId, b: NodeId) -> Option<usize> {
        self.has_node(b).then_some(())?;
        self.wires.get(a)
            .and_then(|mb_nnb| {
                mb_nnb.as_ref()
                    .map(|nnb| {
                        let ma =
                            nnb.iter().filter(|wire| wire.has_s_id(b)).count();
                        if a == b { ma / 2 } else { ma }
                    })
            })
    }

    /// Remove all or at most a fixed number of empty wires between two nodes.
    ///
    /// Returns the number of wires removed.
    ///
    /// Fails if either node does not exist.
    pub fn remove_wires_e(
        &mut self,
        a: NodeId,
        b: NodeId,
        nwires: Option<usize>,
    ) -> GraphResult<usize>
    {
        self.has_node(a).then_some(())
            .ok_or(RemoveWireMissingNode(a))?;
        self.has_node(b).then_some(())
            .ok_or(RemoveWireMissingNode(b))?;

        let mut to_remove: Vec<usize> = Vec::new();

        let nnb_a = self.get_neighbors_mut(a).unwrap();
        let len_nnb_a = nnb_a.len();
        nnb_a.iter().enumerate()
            .filter_map(|(k, nid)| nid.has_e_id(b).then_some(k))
            .take(nwires.unwrap_or(len_nnb_a))
            .for_each(|k| { to_remove.push(k); });
        let mut removed_a = to_remove.len();
        to_remove.drain(..).rev()
            .for_each(|k| { nnb_a.swap_remove(k); });

        let nnb_b = self.get_neighbors_mut(b).unwrap();
        let len_nnb_b = nnb_b.len();
        nnb_b.iter().enumerate()
            .filter_map(|(k, nid)| nid.has_e_id(a).then_some(k))
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

    /// Remove all or at most a fixed number of Hadamard wires between two
    /// nodes.
    ///
    /// Returns the number of wires removed.
    ///
    /// Fails if either node does not exist.
    pub fn remove_wires_h(
        &mut self,
        a: NodeId,
        b: NodeId,
        nwires: Option<usize>,
    ) -> GraphResult<usize>
    {
        self.has_node(a).then_some(())
            .ok_or(RemoveWireMissingNode(a))?;
        self.has_node(b).then_some(())
            .ok_or(RemoveWireMissingNode(b))?;

        let mut to_remove: Vec<usize> = Vec::new();

        let nnb_a = self.get_neighbors_mut(a).unwrap();
        let len_nnb_a = nnb_a.len();
        nnb_a.iter().enumerate()
            .filter_map(|(k, nid)| nid.has_h_id(b).then_some(k))
            .take(nwires.unwrap_or(len_nnb_a))
            .for_each(|k| { to_remove.push(k); });
        let mut removed_a = to_remove.len();
        to_remove.drain(..).rev()
            .for_each(|k| { nnb_a.swap_remove(k); });

        let nnb_b = self.get_neighbors_mut(b).unwrap();
        let len_nnb_b = nnb_b.len();
        nnb_b.iter().enumerate()
            .filter_map(|(k, nid)| nid.has_h_id(a).then_some(k))
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

    /// Remove all or at most a fixed number of star wires between two nodes.
    ///
    /// Returns the number of wires removed.
    ///
    /// Fails if either node does not exist.
    pub fn remove_wires_s(
        &mut self,
        a: NodeId,
        b: NodeId,
        nwires: Option<usize>,
    ) -> GraphResult<usize>
    {
        self.has_node(a).then_some(())
            .ok_or(RemoveWireMissingNode(a))?;
        self.has_node(b).then_some(())
            .ok_or(RemoveWireMissingNode(b))?;

        let mut to_remove: Vec<usize> = Vec::new();

        let nnb_a = self.get_neighbors_mut(a).unwrap();
        let len_nnb_a = nnb_a.len();
        nnb_a.iter().enumerate()
            .filter_map(|(k, nid)| nid.has_s_id(b).then_some(k))
            .take(nwires.unwrap_or(len_nnb_a))
            .for_each(|k| { to_remove.push(k); });
        let mut removed_a = to_remove.len();
        to_remove.drain(..).rev()
            .for_each(|k| { nnb_a.swap_remove(k); });

        let nnb_b = self.get_neighbors_mut(b).unwrap();
        let len_nnb_b = nnb_b.len();
        nnb_b.iter().enumerate()
            .filter_map(|(k, nid)| nid.has_s_id(a).then_some(k))
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

    /// Add a Hadamard wire between two nodes.
    ///
    /// Fails if one or neither of the nodes exist.
    pub fn add_wire_h(&mut self, a: NodeId, b: NodeId) -> GraphResult<()> {
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
        self.wires[a].as_mut().unwrap().push(CTWire::H(b));
        self.wires[b].as_mut().unwrap().push(CTWire::H(a));
        self.wire_count += 1;
        Ok(())
    }

    /// Add a star wire between two nodes.
    ///
    /// Fails if one or neither of the nodes exist.
    pub fn add_wire_s(&mut self, a: NodeId, b: NodeId) -> GraphResult<()> {
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
        self.wires[a].as_mut().unwrap().push(CTWire::S(b));
        self.wires[b].as_mut().unwrap().push(CTWire::S(a));
        self.wire_count += 1;
        Ok(())
    }

    /// Add a wire with an attached input spider state to a pre-existing node
    /// and return the new node's ID.
    ///
    /// The pre-existing node cannot already be an input or output.
    pub fn add_input_state_z(&mut self, id: NodeId, phase: TPhase)
        -> GraphResult<NodeId>
    {
        self.get_node(id)
            .ok_or(AddWireMissingNode(id))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(AddWireConnectedIO(id))
            })?;
        let input_id = self.add_node(CTNode::Z(phase));
        self.add_wire(id, input_id)?;
        self.inputs.push(IONodeId::State(input_id));
        Ok(input_id)
    }

    /// Add a Hadamard wire with an attached input spider state to a
    /// pre-existing node and return the new node's ID.
    ///
    /// The pre-existing node cannot already be an input or output.
    pub fn add_input_state_x(&mut self, id: NodeId, phase: TPhase)
        -> GraphResult<NodeId>
    {
        self.get_node(id)
            .ok_or(AddWireMissingNode(id))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(AddWireConnectedIO(id))
            })?;
        let input_id = self.add_node(CTNode::Z(phase));
        self.add_wire_h(id, input_id)?;
        self.inputs.push(IONodeId::State(input_id));
        Ok(input_id)
    }

    /// Add a wire with an attached input spider state to a pre-existing node
    /// and return the new node's ID.
    ///
    /// The pre-existing node cannot already be an input or output.
    pub fn add_output_z_effect(&mut self, id: NodeId, phase: TPhase)
        -> GraphResult<NodeId>
    {
        self.get_node(id)
            .ok_or(AddWireMissingNode(id))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(AddWireConnectedIO(id))
            })?;
        let output_id = self.add_node(CTNode::Z(phase));
        self.add_wire(id, output_id)?;
        self.outputs.push(IONodeId::State(output_id));
        Ok(output_id)
    }

    /// Add a Hadamard wire with an attached input spider state to a
    /// pre-existing node and return the new node's ID.
    ///
    /// The pre-existing node cannot already be an input or output.
    pub fn add_output_x_effect(&mut self, id: NodeId, phase: TPhase)
        -> GraphResult<NodeId>
    {
        self.get_node(id)
            .ok_or(AddWireMissingNode(id))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(AddWireConnectedIO(id))
            })?;
        let output_id = self.add_node(CTNode::Z(phase));
        self.add_wire_h(id, output_id)?;
        self.outputs.push(IONodeId::State(output_id));
        Ok(output_id)
    }

    /// Replace an existing diagram input or output with a spider of arity 1.
    ///
    /// The new spider will have the same node ID as the input or output.
    ///
    /// Fails if the given node ID does not exist or is not an input or output.
    pub fn apply_state_z(&mut self, id: NodeId, phase: TPhase)
        -> GraphResult<()>
    {
        if let Some(n) = self.get_node(id) {
            if n.is_input() {
                self.inputs.iter_mut()
                    .find(|ioid| ioid.has_id(id))
                    .unwrap()
                    .make_state();
                let prev = self.get_node_mut(id).unwrap();
                *prev = CTNode::Z(phase);
                Ok(())
            } else if n.is_output() {
                self.outputs.iter_mut()
                    .find(|ioid| ioid.has_id(id))
                    .unwrap()
                    .make_state();
                let prev = self.get_node_mut(id).unwrap();
                *prev = CTNode::Z(phase);
                Ok(())
            } else {
                Err(ApplyStateNotIO(id))
            }
        } else {
            Err(ApplyStateNotIO(id))
        }
    }

    /// Replace an existing diagram input or output with a spider of arity 1 on
    /// a Hadamard wire.
    ///
    /// The new spider will have the same node ID as the input or output.
    ///
    /// Fails if the given node ID does not exist or is not an input or output.
    pub fn apply_state_x(&mut self, id: NodeId, phase: TPhase)
        -> GraphResult<()>
    {
        if let Some(n) = self.get_node(id) {
            if n.is_input() {
                self.inputs.iter_mut()
                    .find(|ioid| ioid.has_id(id))
                    .unwrap()
                    .make_state();
                let prev = self.get_node_mut(id).unwrap();
                *prev = CTNode::Z(phase);
                let mb_int = self.wires[id].as_ref().unwrap().first().copied();
                if let Some(int_wire) = mb_int {
                    if int_wire.is_s() {
                        self.remove_wires_s(id, int_wire.id(), None).unwrap();
                        let new = self.add_node(CTNode::z0());
                        self.add_wire_s(new, int_wire.id()).unwrap();
                        self.add_wire_h(id, new).unwrap();
                    } else {
                        self.wires[id].as_mut().unwrap()
                            .first_mut().unwrap()
                            .toggle_h();
                    }
                }
                Ok(())
            } else if n.is_output() {
                self.outputs.iter_mut()
                    .find(|ioid| ioid.has_id(id))
                    .unwrap()
                    .make_state();
                let prev = self.get_node_mut(id).unwrap();
                *prev = CTNode::Z(phase);
                Ok(())
            } else {
                Err(ApplyStateNotIO(id))
            }
        } else {
            Err(ApplyStateNotIO(id))
        }
    }

    /// Like [`apply_state_z`][Self::apply_state], but using input qubit indices
    /// rather than bare node IDs.
    pub fn apply_state_z_input(&mut self, q: QubitId, phase: TPhase)
        -> GraphResult<()>
    {
        self.get_input_id(q)
            .ok_or(ApplyStateMissingQubit(q))
            .and_then(|nid| self.apply_state_z(nid.id(), phase))
    }

    /// Like [`apply_state_x`][Self::apply_state], but using input qubit indices
    /// rather than bare node IDs.
    pub fn apply_state_x_input(&mut self, q: QubitId, phase: TPhase)
        -> GraphResult<()>
    {
        self.get_input_id(q)
            .ok_or(ApplyStateMissingQubit(q))
            .and_then(|nid| self.apply_state_x(nid.id(), phase))
    }

    /// Like [`apply_state_z`][Self::apply_state], but using output qubit
    /// indices rather than bare node IDs.
    pub fn apply_z_effect_output(&mut self, q: QubitId, phase: TPhase)
        -> GraphResult<()>
    {
        self.get_output_id(q)
            .ok_or(ApplyStateMissingQubit(q))
            .and_then(|nid| self.apply_state_z(nid.id(), phase))
    }

    /// Like [`apply_state_x`][Self::apply_state], but using output qubit
    /// indices rather than bare node IDs.
    pub fn apply_x_effect_output(&mut self, q: QubitId, phase: TPhase)
        -> GraphResult<()>
    {
        self.get_output_id(q)
            .ok_or(ApplyStateMissingQubit(q))
            .and_then(|nid| self.apply_state_x(nid.id(), phase))
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

        // ensure all inputs are in a subgraph at the same rank, attaching the
        // overall scalar at the bottom
        let mut inputs_subgraph_stmt =
            StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rank(RankType::Source)),
            );
        let mut prev: Option<usize> = None;
        for (qid, ioid) in self.inputs.iter().enumerate() {
            let nid = ioid.id();
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
                        Edge::head_node(
                            pid.into(),
                            Some(Port::compass(Compass::South)),
                        )
                        .line_to_node(
                            nid.into(),
                            Some(Port::compass(Compass::North)),
                        )
                        .add_attrpair(style(Style::Invisible))
                    );
            }
            prev = Some(nid);
        }
        // add the overall scalar
        let scalar_id = self.nodes.len();
        let attrs =
            AttrList::new()
            .add_pair(label(format!("{}", self.scalar)))
            .add_pair(shape(Shape::Rectangle))
            .add_pair(style(Style::Filled))
            .add_pair(fillcolor(H_COLOR));
        inputs_subgraph_stmt =
            inputs_subgraph_stmt.add_node(scalar_id.into(), None, Some(attrs));
        if let Some(pid) = prev {
            inputs_subgraph_stmt =
                inputs_subgraph_stmt.add_edge(
                    Edge::head_node(
                        pid.into(),
                        Some(Port::compass(Compass::South)),
                    )
                    .line_to_node(
                        scalar_id.into(),
                        Some(Port::compass(Compass::North)),
                    )
                    .add_attrpair(style(Style::Invisible))
                );
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
        for (qid, ioid) in self.outputs.iter().enumerate() {
            let nid = ioid.id();
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
                        Edge::head_node(
                            pid.into(),
                            Some(Port::compass(Compass::South)),
                        )
                        .line_to_node(
                            nid.into(),
                            Some(Port::compass(Compass::North)),
                        )
                        .add_attrpair(style(Style::Invisible))
                    );
            }
            prev = Some(nid);
        }
        statements =
            statements.add_subgraph(SubGraph::cluster(outputs_subgraph_stmt));
        // add interior nodes
        for (id, node) in self.nodes_inner() {
            if self.inputs.iter().any(|ioid| ioid.has_id(id))
                || self.outputs.iter().any(|ioid| ioid.has_id(id))
            {
                continue;
            }
            let attrs = node.graph_attrs();
            statements = statements.add_node(id.into(), None, Some(attrs));
        }
        // add wires
        for (left, right) in self.wires() {
            match *right {
                CTWire::E(r) => {
                    statements =
                        statements.add_edge(
                            Edge::head_node(left.into(), None)
                            .line_to_node(r.into(), None)
                        );
                },
                CTWire::H(r) => {
                    statements =
                        statements.add_edge(
                            Edge::head_node(left.into(), None)
                            .line_to_node(r.into(), None)
                            .add_attrpair(style(Style::Dashed))
                            .add_attrpair(color(H_WIRE))
                        );
                },
                CTWire::S(r) => {
                    statements =
                        statements.add_edge(
                            Edge::head_node(left.into(), None)
                            .line_to_node(r.into(), None)
                            .add_attrpair(style(Style::Dashed))
                            .add_attrpair(color(STAR_WIRE))
                        );
                },
            }
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

