use std::io::Write;
use num_complex::Complex64 as C64;
use crate::{
    indexmap::IndexMap,
    ketbra2::{ Element2, Kind, KBError, KBResult },
};

use KBError::*;

/// Represents a series of [`Element2`]s as ket-bra-backed "slices" of a
/// ZX(H)-diagram.
///
/// Each `Element`'s inputs and outputs can span the entire lateral space of
/// wires.
#[derive(Clone, Debug)]
pub struct Diagram2 {
    pub(crate) elems: Vec<Element2>,
    pub(crate) scalar: C64,
}

impl FromIterator<Element2> for Diagram2 {
    fn from_iter<I>(iter: I) -> Self
    where I: IntoIterator<Item = Element2>
    {
        Self { elems: iter.into_iter().collect(), scalar: 1.0_f64.into() }
    }
}

impl From<Element2> for Diagram2 {
    fn from(elem: Element2) -> Self {
        Self { elems: vec![elem], scalar: 1.0_f64.into() }
    }
}

impl Diagram2 {
    /// Create a new diagram from a list of [`Element2`]s and an overall scalar,
    /// defaulting to 1.
    ///
    /// `Element`s should follow the left-to-right order in which they would be
    /// drawn in a real ZX(H)-diagram; i.e. the reverse of the order in which
    /// they appear on paper when dot-products are taken.
    pub fn new(elems: Vec<Element2>, scalar: Option<C64>) -> Self {
        Self { elems, scalar: scalar.unwrap_or(1.0.into()) }
    }

    /// Create a new diagram from a list of [`Element2`]s.
    ///
    /// `Element`s should follow the left-to-right order in which they would be
    /// drawn in a real ZX(H)-diagram; i.e. the reverse of the order in which
    /// they appear on paper when dot-products are taken.
    pub fn from_elems<I>(elems: I) -> Self
    where I: IntoIterator<Item = Element2>
    {
        elems.into_iter().collect()
    }

    /// Set an overall scalar factor.
    pub fn with_scalar(mut self, a: C64) -> Self {
        self.scalar = a;
        self
    }

    /// Set an overall scalar factor.
    pub fn set_scalar(&mut self, a: C64) -> &mut Self {
        self.scalar = a;
        self
    }

    /// Return the overall scalar factor.
    pub fn scalar(&self) -> C64 { self.scalar }

    /// Fold over all elements of the diagram with the [dot
    /// product][Element2::dot] and return the result as a single [`Element2`].
    pub fn contract(self) -> KBResult<Element2> {
        fn foldfn(mb_acc: Option<Element2>, elem: Element2)
            -> KBResult<Option<Element2>>
        {
            if let Some(acc) = mb_acc {
                acc.then(elem).map(Some)
            } else {
                Ok(Some(elem))
            }
        }

        self.elems.into_iter()
            .try_fold(None, foldfn)
            .map(|mb_res| {
                mb_res.unwrap_or_else(|| Element2::new_scalar(1.0))
                    .scalar_mul(self.scalar)
            })
    }

    /// Compose `self` with `rhs` by attaching the outputs of `rhs` to the
    /// inputs of `self`.
    ///
    /// See also [`compose_rev`][Self::compose_rev] for composition following a
    /// left-to-right diagram placement.
    pub fn compose(mut self, mut rhs: Self) -> Self {
        self.elems.append(&mut rhs.elems);
        self.scalar *= rhs.scalar;
        self
    }

    /// Compose `self` with `rhs` by attaching the outputs of `self` to the
    /// inputs of `rhs`.
    ///
    /// See also [`compose`][Self::compose] for composition following the usual
    /// `self ∘ rhs` operation.
    pub fn compose_rev(mut self, mut rhs: Self) -> Self {
        rhs.elems.append(&mut self.elems);
        rhs.scalar *= self.scalar;
        rhs
    }

    /// Return all input and output wire indices.
    pub fn ins_outs(&self) -> (Vec<usize>, Vec<usize>) {
        let mut ins: Vec<usize> = Vec::new();
        let mut outs: Vec<usize> = Vec::new();
        for element in self.elems.iter() {
            element.ins()
                .for_each(|k| {
                    if outs.contains(&k) {
                        vec_remove_elem(&mut outs, &k);
                    } else {
                        ins.push(k);
                    }
                });
            element.outs()
                .for_each(|k| { outs.push(k); });
        }
        ins.sort_unstable();
        outs.sort_unstable();
        (ins, outs)
    }

    // /// Convert `self` to a [`tensor::Diagram`][crate::tensor::Diagram]
    // /// representation.
    // pub fn as_tensor(&self) -> crate::tensor::Diagram {
    //     crate::tensor::Diagram::new(
    //         self.elems.iter().map(|elem| elem.as_tensor())
    //     ).with_scalar(self.scalar)
    // }

    /// Convert `self` to a [`graph::Diagram`][crate::graph2::Diagram]
    /// representation.
    ///
    /// Fails if any `Element` slices cannot be written as a pure generator or
    /// if any `Element` repeats a ket index without another `Element`'s
    /// matching bra index in between.
    pub fn as_graph(&self) -> KBResult<crate::graph2::Diagram> {
        use crate::graph2 as graph;
        let mut graph = graph::Diagram::new();
        let mut node_id: graph::NodeId;
        let (inputs, outputs) = self.ins_outs();
        let mut wires: IndexMap<graph::NodeId> =
            IndexMap::new(); // wire idx -> node idx
        let mut to_remove: Vec<graph::NodeId> = Vec::new(); // placeholder nodes

        // add inputs
        for idx in inputs.into_iter() {
            node_id = graph.add_node(graph::Node::Input);
            wires.insert(idx, node_id);
        }

        // add elements/wires
        for elem in self.elems.iter() {
            match elem.kind {
                Kind::Id(_) => { },
                Kind::Scalar(a) => { graph.map_scalar(|s| s * a); },
                Kind::Z(_) | Kind::X(_) | Kind::H(_) => {
                    let node =
                        match elem.kind {
                            Kind::Z(ph) => graph::Node::Z(ph),
                            Kind::X(ph) => graph::Node::X(ph),
                            Kind::H(a) => graph::Node::H(a),
                            _ => unreachable!(),
                        };
                    node_id = graph.add_node(node);
                    for in_wire in elem.ins() {
                        let prev_id = wires.remove(in_wire).unwrap();
                        graph.add_wire(prev_id, node_id).unwrap();
                    }
                    for out_wire in elem.outs() {
                        if wires.insert(out_wire, node_id).is_some() {
                            return Err(DuplicateKetKey(out_wire));
                        }
                    }
                },
                Kind::Swap(a, b) => {
                    wires.swap(a, b);
                },
                Kind::Cup(a, b) => {
                    let n = graph.add_node(graph::Node::z());
                    to_remove.push(n);
                    if wires.insert(a, n).is_some() {
                        return Err(DuplicateKetKey(a));
                    }
                    if wires.insert(b, n).is_some() {
                        return Err(DuplicateKetKey(b));
                    }
                },
                Kind::Cap(a, b) => {
                    let n1 = wires.remove(a).unwrap();
                    let n2 = wires.remove(b).unwrap();
                    graph.add_wire(n1, n2).unwrap();
                },
                Kind::Unknown => {
                    return Err(GraphConvNonGenerator);
                },
            }
        }

        // add outputs
        for idx in outputs.into_iter() {
            node_id = graph.add_node(graph::Node::Output);
            let prev_id = wires.remove(idx).unwrap();
            graph.add_wire(prev_id, node_id).unwrap();
        }

        // all nodes to remove here are placeholders for cups; they have exactly
        // two neighbors
        for remove in to_remove.into_iter() {
            let (_, nnb) = graph.remove_node_nb(remove).unwrap();
            graph.add_wire(nnb[0], nnb[1]).unwrap();
        }

        Ok(graph)
    }

    /// Return an object containing an encoding of `self` in the [dot
    /// language][dot-lang].
    ///
    /// Rendering this object using the default formatter will result in a full
    /// dot string representation of the diagram.
    ///
    /// If any `Element` slice cannot be written as a pure generator, it will be
    /// replaced with a block labeled `U_{k}` in the resulting graph, where `k`
    /// is an index corresponding to its position in `self`.
    ///
    /// [dot-lang]: https://en.wikipedia.org/wiki/DOT_(graph_description_language)
    pub fn to_graphviz(&self) -> KBResult<tabbycat::Graph> {
        use tabbycat::*;
        use tabbycat::attributes::*;
        use crate::vizdefs::*;
        let mut id_gen = 0..;
        let mut node_id: usize;

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
        let (inputs, outputs) = self.ins_outs();
        let mut wires: IndexMap<usize> = IndexMap::new(); // wire idx -> node id

        // ensure all inputs are in a subgraph at the same rank
        let mut inputs_subgraph_stmt =
            StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rank(RankType::Source)),
            );
        let mut prev: Option<usize> = None;
        for idx in inputs.into_iter() {
            node_id = id_gen.next().unwrap();
            wires.insert(idx, node_id);
            inputs_subgraph_stmt =
                inputs_subgraph_stmt.add_node(
                    node_id.into(),
                    None,
                    Some(
                        AttrList::new()
                            .add_pair(label(format!("In {}", idx)))
                            .add_pair(shape(Shape::Plaintext))
                    ),
                );
            if let Some(prev_idx) = prev {
                inputs_subgraph_stmt =
                    inputs_subgraph_stmt.add_edge(
                        Edge::head_node(
                            prev_idx.into(),
                            Some(Port::compass(Compass::South)),
                        )
                        .line_to_node(
                            node_id.into(),
                            Some(Port::compass(Compass::North)),
                        )
                        .add_attrpair(style(Style::Invisible))
                    );
            }
            prev = Some(node_id);
        }
        statements =
            statements.add_subgraph(SubGraph::cluster(inputs_subgraph_stmt));

        // add elements/wires
        let mut u_idx = 0..;
        for elem in self.elems.iter() {
            match elem.kind {
                Kind::Z(_)
                | Kind::X(_)
                | Kind::H(_)
                | Kind::Unknown
                => {
                    node_id = id_gen.next().unwrap();
                    let attrs =
                        match elem.kind {
                            Kind::Z(ph) => {
                                AttrList::new()
                                    .add_pair(label(ph.label()))
                                    .add_pair(shape(Shape::Circle))
                                    .add_pair(height(CIRCLE_HEIGHT))
                                    .add_pair(style(Style::Filled))
                                    .add_pair(fillcolor(Z_COLOR))
                            },
                            Kind::X(ph) => {
                                AttrList::new()
                                    .add_pair(label(ph.label()))
                                    .add_pair(shape(Shape::Circle))
                                    .add_pair(height(CIRCLE_HEIGHT))
                                    .add_pair(style(Style::Filled))
                                    .add_pair(fillcolor(X_COLOR))
                            },
                            Kind::H(a) => {
                                let a_label =
                                    if (a + 1.0).norm() < 1e-12 {
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
                            Kind::Unknown => {
                                let u_label =
                                    format!("U{}",
                                        subscript_str(u_idx.next().unwrap()));
                                AttrList::new()
                                    .add_pair(label(u_label))
                                    .add_pair(shape(Shape::Rectangle))
                                    .add_pair(style(Style::Filled))
                                    .add_pair(fillcolor(Color::White))
                            },
                            _ => unreachable!(),
                        };
                    statements =
                        statements.add_node(node_id.into(), None, Some(attrs));
                    for in_wire in elem.ins() {
                        let prev_id = wires.remove(in_wire).unwrap();
                        statements =
                            statements.add_edge(
                                Edge::head_node(prev_id.into(), None)
                                    .line_to_node(node_id.into(), None)
                            );
                    }
                    for out_wire in elem.outs() {
                        if wires.insert(out_wire, node_id).is_some() {
                            return Err(DuplicateKetKey(out_wire));
                        }
                    }
                },
                Kind::Scalar(mut a) => {
                    a.re = (1e6 * a.re).round() / 1e6;
                    a.im = (1e6 * a.im).round() / 1e6;
                    node_id = id_gen.next().unwrap();
                    let attrs =
                        AttrList::new()
                        .add_pair(label(format!("{}", a)))
                        .add_pair(shape(Shape::Rectangle))
                        .add_pair(style(Style::Filled))
                        .add_pair(fillcolor(H_COLOR));
                    statements =
                        statements.add_node(node_id.into(), None, Some(attrs));
                },
                Kind::Swap(a, b) => {
                    wires.swap(a, b);
                },
                Kind::Cup(a, b) => {
                    let n1 = id_gen.next().unwrap();
                    let attrs1 =
                        AttrList::new().add_pair(style(Style::Invisible));
                    statements =
                        statements.add_node(n1.into(), None, Some(attrs1));
                    wires.insert(a, n1);
                    let n2 = id_gen.next().unwrap();
                    let attrs2 =
                        AttrList::new().add_pair(style(Style::Invisible));
                    statements =
                        statements.add_node(n2.into(), None, Some(attrs2));
                    wires.insert(b, n2);
                    statements =
                        statements.add_edge(
                            Edge::head_node(n1.into(), None)
                                .line_to_node(n2.into(), None)
                        );
                },
                Kind::Cap(a, b) => {
                    let n1 = wires.remove(a).unwrap();
                    let n2 = wires.remove(b).unwrap();
                    statements =
                        statements.add_edge(
                            Edge::head_node(n1.into(), None)
                                .line_to_node(n2.into(), None)
                        );
                },
                Kind::Id(_) => { },
            }
        }

        // add the overall scalar
        let mut a = self.scalar;
        a.re = (1e6 * a.re).round() / 1e6;
        a.im = (1e6 * a.im).round() / 1e6;
        node_id = id_gen.next().unwrap();
        let attrs =
            AttrList::new()
            .add_pair(label(format!("{}", a)))
            .add_pair(shape(Shape::Rectangle))
            .add_pair(style(Style::Filled))
            .add_pair(fillcolor(H_COLOR));
        statements =
            statements.add_node(node_id.into(), None, Some(attrs));

        // ensure all outputs are in a subgraph at the same rank
        let mut outputs_subgraph_stmt =
            StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rank(RankType::Sink)),
            );
        let mut prev: Option<usize> = None;
        for idx in outputs.into_iter() {
            node_id = id_gen.next().unwrap();
            outputs_subgraph_stmt =
                outputs_subgraph_stmt.add_node(
                    node_id.into(),
                    None,
                    Some(
                        AttrList::new()
                            .add_pair(label(format!("Out {}", idx)))
                            .add_pair(shape(Shape::Plaintext))
                    ),
                );
            if let Some(prev_idx) = prev {
                outputs_subgraph_stmt =
                    outputs_subgraph_stmt.add_edge(
                        Edge::head_node(
                            prev_idx.into(),
                            Some(Port::compass(Compass::South)),
                        )
                        .line_to_node(
                            node_id.into(),
                            Some(Port::compass(Compass::North)),
                        )
                        .add_attrpair(style(Style::Invisible))
                    );
            }
            prev = Some(node_id);
            let prev_id = wires.remove(idx).unwrap();
            statements =
                statements.add_edge(
                    Edge::head_node(prev_id.into(), None)
                        .line_to_node(node_id.into(), None)
                );
        }
        statements =
            statements.add_subgraph(SubGraph::cluster(outputs_subgraph_stmt));

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
    pub fn save_graphviz<P>(&self, path: P) -> KBResult<()>
    where P: AsRef<std::path::Path>
    {
        let graphviz = self.to_graphviz()?;
        std::fs::OpenOptions::new()
            .write(true)
            .append(false)
            .create(true)
            .truncate(true)
            .open(path)?
            .write_all(format!("{}", graphviz).as_bytes())?;
        Ok(())
    }
}

fn vec_remove_elem<T>(v: &mut Vec<T>, target: &T) -> Option<T>
where T: PartialEq
{
    v.iter().enumerate()
        .find_map(|(k, t)| (t == target).then_some(k))
        .map(|k0| v.swap_remove(k0))
}

fn subscript_str(n: usize) -> String {
    format!("{}", n)
        .replace('0', "₀")
        .replace('1', "₁")
        .replace('2', "₂")
        .replace('3', "₃")
        .replace('4', "₄")
        .replace('5', "₅")
        .replace('6', "₆")
        .replace('7', "₇")
        .replace('8', "₈")
        .replace('9', "₉")
}

