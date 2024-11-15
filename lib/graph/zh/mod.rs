//! Variant of the ZX-calculus based on generalization of the H-box.
//!
#![cfg_attr(
    feature = "doc-images",
    cfg_attr(
        all(),
        doc = ::embed_doc_image::embed_image!("h_box_def", "assets/h_box_def.svg"),
    )
)]
#![cfg_attr(
    not(feature = "doc-images"),
    doc = "**Doc images are not enabled**. Compile with feature `doc-images` and Rust version >= 1.54 to enable."
)]
//!
//! Here, the usual generator describing the Hadamard transformation is
//! generalized to arbitrary arity and given an argument, as follows:
//!
//! ![h_box_def][h_box_def]
//!
//! Note that in the case where *n* + *m* = 2 and *a* = –1, an extra scalar
//! factor of 1/√2 is added so that the H-box coincides with the usual
//! definition of the Hadamard gate.

use std::{ collections::VecDeque, fs, io::Write, path::Path };
use num_complex::Complex64 as C64;
use crate::{
    graph::{
        GraphError,
        GraphResult,
        Diagram,
        DiagramData,
        IONodeId,
        NodeId,
        QubitId,
        Spider,
        WireData,
        WireStore,
    },
    phase::Phase,
};

use GraphError::*;

pub(crate) mod zhnode;
pub use zhnode::*;

/// Marker type for diagrams in the ZH-calculus.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct ZH;

impl DiagramData for ZH {
    type Node = ZHNode;
    type Wire = NodeId;
    type Scalar = C64;
}

impl Diagram<ZH> {
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

    /// Return the total number of spiders.
    pub fn count_spiders(&self) -> usize {
        self.nodes.iter()
            .filter(|mb_n| mb_n.is_some_and(|n| n.is_spider()))
            .count()
    }

    /// Add a Z-spider to the diagram and return its ID.
    ///
    /// `phase` defaults to 0.
    pub fn add_z(&mut self, phase: Option<Phase>) -> NodeId {
        self.add_node(ZHNode::Z(phase.unwrap_or_else(Phase::zero)))
    }

    /// Add an X-spider to the diagram and return its ID.
    ///
    /// `phase` defaults to 0.
    pub fn add_x(&mut self, phase: Option<Phase>) -> NodeId {
        self.add_node(ZHNode::X(phase.unwrap_or_else(Phase::zero)))
    }

    /// Add an H-box to the diagram and return its ID.
    ///
    /// `arg` defaults to –1.
    pub fn add_h(&mut self, arg: Option<C64>) -> NodeId {
        self.add_node(ZHNode::H(arg.unwrap_or_else(|| -C64::from(1.0))))
    }

    /// Add a wire with an attached input spider state to a pre-existing node
    /// and return the new node's ID.
    ///
    /// The pre-existing node cannot already be an input or output.
    pub fn add_input_state(&mut self, id: NodeId, state: Spider)
        -> GraphResult<NodeId>
    {
        self.get_node(id)
            .ok_or(AddWireMissingNode(id))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(AddWireConnectedIO(id))
            })?;
        let input_id = self.add_node(state.into());
        self.add_wire(id, input_id)?;
        self.inputs.push(IONodeId::State(input_id));
        Ok(input_id)
    }

    /// Add a wire with an attached input spider state to a pre-existing node
    /// and return the new node's ID.
    ///
    /// The pre-existing node cannot already be an input or output.
    pub fn add_output_effect(&mut self, id: NodeId, effect: Spider)
        -> GraphResult<NodeId>
    {
        self.get_node(id)
            .ok_or(AddWireMissingNode(id))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(AddWireConnectedIO(id))
            })?;
        let output_id = self.add_node(effect.into());
        self.add_wire(id, output_id)?;
        self.outputs.push(IONodeId::State(output_id));
        Ok(output_id)
    }

    /// Replace an existing diagram input or output with a spider of arity 1.
    ///
    /// The new spider will have the same node ID as the input or output.
    ///
    /// Fails if the given node ID does not exist or is not an input or output.
    pub fn apply_state(&mut self, id: NodeId, spider: Spider)
        -> GraphResult<()>
    {
        if let Some(n) = self.get_node(id) {
            if n.is_input() {
                self.inputs.iter_mut()
                    .find(|ioid| ioid.has_id(id))
                    .unwrap()
                    .make_state();
                let prev = self.get_node_mut(id).unwrap();
                *prev = spider.into();
                Ok(())
            } else if n.is_output() {
                self.outputs.iter_mut()
                    .find(|ioid| ioid.has_id(id))
                    .unwrap()
                    .make_state();
                let prev = self.get_node_mut(id).unwrap();
                *prev = spider.into();
                Ok(())
            } else {
                Err(ApplyStateNotIO(id))
            }
        } else {
            Err(ApplyStateNotIO(id))
        }
    }

    /// Like [`apply_state`][Self::apply_state], but using input qubit indices
    /// rather than bare node IDs.
    pub fn apply_state_input(&mut self, q: QubitId, spider: Spider)
        -> GraphResult<()>
    {
        self.get_input_id(q)
            .ok_or(ApplyStateMissingQubit(q))
            .and_then(|nid| self.apply_state(nid.id(), spider))
    }

    /// Like [`apply_state`][Self::apply_state], but using output qubit indices
    /// rather than bare node IDs.
    pub fn apply_effect_output(&mut self, q: QubitId, spider: Spider)
        -> GraphResult<()>
    {
        self.get_output_id(q)
            .ok_or(ApplyStateMissingQubit(q))
            .and_then(|nid| self.apply_state(nid.id(), spider))
    }

    // BFS explore from a given starting queue, accumulating a list of tensor
    // elements
    fn element_explore<'a, T>(
        &'a self,
        wire_nums: &WireStore,
        visited: &mut Vec<NodeId>,
        to_visit: &mut VecDeque<(NodeId, &'a ZHNode)>,
        elements: &mut Vec<crate::tensor::Element<T>>,
    )
    where T: crate::tensor::ElementData
    {
        use crate::tensor;
        // upper bound on possible qubit indices
        let qcount = self.inputs.len().max(self.outputs.len());
        // want wire IDs to line up with qubit indices when connected to
        // inputs/outputs, so shift any ID that could be a qubit ID by the total
        // wire count
        let nonq_wire = |id: usize| -> usize {
            if id < qcount { self.wire_count + id } else { id }
        };
        // self-loops are treated as inputs to a Bell effect, whose input wire
        // IDs are shifted to come after all existing (possibly shifted) wires
        let mut bell_wire = self.wire_count + qcount..;

        // loop variables
        let mut ins: Vec<usize> = Vec::new(); // inputs to an element
        let mut outs: Vec<usize> = Vec::new(); // outputs from an element
        let mut bell_wires: Vec<(usize, usize)> = Vec::new(); // for self-loop wires
        let mut empty_wires: Vec<(QubitId, QubitId)> = Vec::new(); // direct input-to-output
        while let Some((id, node)) = to_visit.pop_back() {
            visited.push(id);
            for (id2, node2) in self.neighbors(id).unwrap() {
                let id2 = id2.id();
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
                    let qout = self.get_output_index(id2).unwrap();
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
                elements.push(node.to_element(&ins, &outs));
            }
            if !bell_wires.is_empty() {
                for (b1, b2) in bell_wires.drain(..) {
                    elements.push(tensor::Element::cap(b1, b2));
                }
            }
            ins.clear();
            outs.clear();
        }
        // empty wires may include swaps, which have to be dealt with
        //
        // individual swaps are be determined by sorting on input wire indices,
        // and then finding the swaps required to sort the output wire indices
        if !empty_wires.is_empty() {
            empty_wires.sort_by(|(qin_l, _), (qin_r, _)| qin_l.cmp(qin_r));
            let (idents_in, mut idents_out): (Vec<QubitId>, Vec<QubitId>) =
                empty_wires.into_iter().unzip();
            let mut swaps: Vec<tensor::Element<T>> = Vec::new();
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
                    swaps.push(tensor::Element::swap(*qin, *qswap));
                    idents_out.swap(k, kswap);
                } else {
                    break;
                }
            }
            swaps.reverse();
            elements.append(&mut swaps);
        }
    }

    /// Convert `self` to a [`tensor::Diagram`][crate::tensor::Diagram]
    /// representation.
    pub fn to_tensor<T>(&self) -> crate::tensor::Diagram<T>
    where T: crate::tensor::ElementData
    {
        use crate::tensor;
        // assemble the tensor diagram by iterating over nodes and analyzing
        // input/output wires relative to what's already been seen
        //
        // have to BFS explore starting from the input nodes in order to ensure
        // that input-adjencent nodes are placed first in the tensor diagram
        //
        // we also want to have wire numbers in the tensor diagram line up with
        // qubit indices for convenience, so if a given wire id (normally used
        // as-is for a tensor wire index) coincides with a possible qubit index,
        // shift it by the maximum wire id in the (graph) diagram
        //
        // do this in two steps because nodes with paths to inputs/outputs need
        // to be visited in a special order, but everything else (i.e. part of a
        // scalar) doesn't
        //
        // self-wires are dealt with by adding two extra outgoing wires
        // immediately coupled to a spiderless Bell effect

        // tensor accumulator
        let mut elements: Vec<tensor::Element<T>> = Vec::new();
        // have to index wires so that each one can have a unique ID in the
        // tensor diagram
        let wire_nums: WireStore =
            self.wires().map(|(l, r)| (l, r.id())).collect();

        // first step: all non-scalar nodes
        let mut visited: Vec<NodeId> = Vec::with_capacity(self.node_count);
        // init with input nodes to make sure they're seen first
        let mut to_visit: VecDeque<(NodeId, &ZHNode)> =
            self.inputs()
            .map(|(_, ioid)| ioid.id())
            .map(|nid| (nid, self.get_node(nid).unwrap()))
            .collect();
        self.element_explore(
            &wire_nums, &mut visited, &mut to_visit, &mut elements);

        // second step: all nodes that aren't part of a scalar
        // reset `to_visit` with everything not already visited
        to_visit
            = self.nodes_inner()
            .filter(|(id, _)| !visited.contains(id))
            .collect();
        self.element_explore(
            &wire_nums, &mut visited, &mut to_visit, &mut elements);

        tensor::Diagram::new(elements, Some(self.scalar))
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
            statements =
                statements.add_edge(
                    Edge::head_node(left.into(), None)
                    .line_to_node(right.id().into(), None)
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

