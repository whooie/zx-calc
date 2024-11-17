use std::io::Write;
use num_complex::Complex64 as C64;
use crate::{
    graph,
    indexmap::IndexMap,
    tensor::{
        Element,
        ElementData,
        Kind,
        TensorError,
        TensorResult,
        dense::De,
        sparse::Sp,
    },
};
use TensorError::*;

/// Represents a series of [`Element`]s as "slices" of a ZX(H)-diagram.
///
/// Each `Element`'s inputs and outputs are allowed to span the entire lateral
/// space of wires, and the inputs and outputs to the diagram as a whole are
/// inferred from those of the elements.
///
/// As with `Element`, `Diagram`s are parameterized by the backing
/// representation of the underlying tensors in each slice.
#[derive(Clone, Debug)]
pub struct Diagram<A> {
    pub(crate) elems: Vec<Element<A>>,
    pub(crate) scalar: C64,
}

/// A [`Diagram`] with [dense][super::ElementDe] elements.
pub type DiagramDe = Diagram<De>;

/// A [`Diagram`] with [sparse][super::ElementSp] elements.
pub type DiagramSp = Diagram<Sp>;

impl<A> FromIterator<Element<A>> for Diagram<A> {
    fn from_iter<I>(iter: I) -> Self
    where I: IntoIterator<Item = Element<A>>
    {
        Self { elems: iter.into_iter().collect(), scalar: 1.0_f64.into() }
    }
}

impl<A> From<Element<A>> for Diagram<A> {
    fn from(elem: Element<A>) -> Self {
        Self { elems: vec![elem], scalar: 1.0_f64.into() }
    }
}

impl<A> Default for Diagram<A> {
    fn default() -> Self { Self { elems: Vec::new(), scalar: 1.0_f64.into() } }
}

impl<A> Diagram<A> {
    /// Create a new diagram from a list of [`Element`]s and an overall scalar,
    /// defaulting to 1.
    ///
    /// `Element`s should follow the left-to-right order in which they would be
    /// laid out diagrammatically; i.e. the reverse of the order in which dot
    /// products are taken.
    pub fn new(elems: Vec<Element<A>>, scalar: Option<C64>) -> Self {
        Self { elems, scalar: scalar.unwrap_or(1.0.into()) }
    }

    /// Create a new, empty diagram with overall scalar 1.
    pub fn empty() -> Self { Self::default() }

    /// Create a new diagram from a series of [`Element`]s.
    ///
    /// `Element`s should follow the left-to-right order in which they would be
    /// laid out diagrammatically; i.e. the reverse of the order in which dot
    /// products are taken.
    pub fn from_elems<I>(elems: I) -> Self
    where I: IntoIterator<Item = Element<A>>
    {
        elems.into_iter().collect()
    }

    /// Return an iterator over all elements.
    ///
    /// The iterator item type is `&`[`Element`]`<A>`.
    pub fn elems(&self) -> Elems<'_, A> { Elems { iter: self.elems.iter() } }

    /// Return the number of elements in `self`.
    pub fn len(&self) -> usize { self.elems.len() }

    /// Return `true` if `self` is the empty diagram; i.e. no elements *and*
    /// overall scalar equal to 1.
    pub fn is_empty(&self) -> bool {
        self.elems.is_empty() && (self.scalar - 1.0).norm() < 1e-12
    }

    /// Return `true` if `self` contains no elements.
    pub fn no_elems(&self) -> bool { self.elems.is_empty() }

    /// Set the overall scalar.
    pub fn with_scalar<C>(mut self, a: C) -> Self
    where C: Into<C64>
    {
        self.scalar = a.into();
        self
    }

    /// Set the overall scalar.
    pub fn set_scalar<C>(&mut self, a: C) -> &mut Self
    where C: Into<C64>
    {
        self.scalar = a.into();
        self
    }

    /// Apply a mapping function to the overall scalar.
    pub fn map_scalar<F>(&mut self, map: F) -> &mut Self
    where F: FnOnce(C64) -> C64
    {
        self.scalar = map(self.scalar);
        self
    }

    /// Return the overall scalar.
    pub fn scalar(&self) -> C64 { self.scalar }

    /// Return a reference to the `k`-th `Element`, if it exists.
    pub fn get(&self, k: usize) -> Option<&Element<A>> { self.elems.get(k) }

    /// Add a new `Element` onto the end of the diagram.
    pub fn push(&mut self, elem: Element<A>) { self.elems.push(elem); }

    /// Insert a new `Element` at the `k`-th position, shifting everything after
    /// it.
    pub fn insert(&mut self, k: usize, elem: Element<A>) {
        self.elems.insert(k, elem);
    }

    /// Remove the `Element` at the `k`-th position, returning it if it existed.
    pub fn remove(&mut self, k: usize) -> Option<Element<A>> {
        if k < self.elems.len() {
            Some(self.elems.remove(k))
        } else {
            None
        }
    }

    /// Compose `self` with `rhs` by attaching the outputs of `rhs` to the
    /// inputs of `self`.
    ///
    /// See also [`compose_rev`][Self::compose_rev] for composition following
    /// left-to-right diagram placement.
    pub fn compose(mut self, mut rhs: Self) -> Self {
        rhs.elems.append(&mut self.elems);
        rhs.scalar *= self.scalar;
        rhs
    }

    /// Compose `self` with `rhs` by attaching the outputs of `self` to the
    /// inputs of `rhs`.
    ///
    /// See also [`compose`][Self::compose] for composition following the usual
    /// `self ∘ rhs` operation.
    pub fn compose_rev(mut self, mut rhs: Self) -> Self {
        self.elems.append(&mut rhs.elems);
        self.scalar *= rhs.scalar;
        self
    }
}

impl<A> Diagram<A>
where A: ElementData
{
    /// If `self` contains only scalars, return their total product.
    pub fn as_scalar(&self) -> Option<C64> {
        self.elems.iter()
            .try_fold(
                C64::from(1.0),
                |acc, elem| elem.as_scalar().map(|a| acc * a),
            )
            .map(|a| self.scalar * a)
    }

    /// Remove scalar `Element`s, returning their total product *and* absorbing
    /// it into the global scalar.
    ///
    /// Returns `None` if no scalar `Element`s were found.
    pub fn remove_scalars(&mut self) -> Option<C64> {
        let mut to_remove: Vec<usize> = Vec::new();
        let res: Option<C64> =
            self.elems.iter().enumerate()
            .fold(
                None,
                |mb_acc, (k, elem)| {
                    match (mb_acc, elem.as_scalar()) {
                        (Some(acc), Some(a)) => {
                            to_remove.push(k);
                            Some(acc * a)
                        },
                        (None, Some(a)) => {
                            to_remove.push(k);
                            Some(a)
                        },
                        (Some(acc), None) => Some(acc),
                        (None, None) => None,
                    }
                },
            );
        if let Some(a) = res {
            to_remove.into_iter().rev()
                .for_each(|k| { self.elems.remove(k); });
            self.scalar *= a;
            Some(a)
        } else {
            None
        }
    }

    /// Return sorted lists of all input and output wire indices (i.e. ones that
    /// have no match within `self`).
    ///
    /// Note that the returned lists may contain duplicate elements if unmatched
    /// duplicate wire indices exist. In this case, subsequent calls to
    /// [`contract`][Self::contract] will return an error.
    pub fn ins_outs(&self) -> (Vec<usize>, Vec<usize>) {
        let mut ins: Vec<usize> = Vec::new();
        let mut outs: Vec<usize> = Vec::new();
        for elem in self.elems.iter() {
            elem.input_iter()
                .for_each(|k| {
                    if outs.contains(&k) {
                        vec_remove_elem(&mut outs, &k);
                    } else {
                        ins.push(k);
                    }
                });
            elem.output_iter()
                .for_each(|k| { outs.push(k); });
        }
        ins.sort_unstable();
        outs.sort_unstable();
        (ins, outs)
    }

    /// Search elements for any duplicated, unmatched wire indices, returning
    /// the wire index and the position of the element it belongs to.
    ///
    /// If this function returns `Some`, then subsequent calls to
    /// [`contract`][Self::contract] will return an error.
    pub fn find_dup(&self) -> Option<(usize, Dup)> {
        let mut ins: Vec<usize> = Vec::new();
        let mut outs: Vec<usize> = Vec::new();
        for (pos, elem) in self.elems.iter().enumerate() {
            for input in elem.input_iter() {
                if outs.contains(&input) {
                    vec_remove_elem(&mut outs, &input);
                } else if ins.contains(&input) {
                    return Some((pos, Dup::Bra(input)));
                } else {
                    ins.push(input);
                }
            }
            for output in elem.output_iter() {
                if outs.contains(&output) {
                    return Some((pos, Dup::Ket(output)));
                } else {
                    outs.push(output);
                }
            }
        }
        None
    }

    /// Reverse the order of elements, conjugating each one as well as the
    /// global scalar.
    pub fn adjoint_mut(&mut self) {
        self.scalar = self.scalar.conj();
        self.elems.reverse();
        self.elems.iter_mut()
            .for_each(|elem| { elem.adjoint_mut(); });
    }

    /// Return the adjoint of `self` with conjugated global scalar, and reversed
    /// and conjugated elements.
    pub fn adjoint(mut self) -> Self {
        self.adjoint_mut();
        self
    }

    /// Return the tensor product `self ⊗ other`, consuming both.
    ///
    /// All wire indices from `other` will be shifted to avoid collision. The
    /// exact size of the shift depends on the maximum wire index value
    /// currently in `self`. After shifting, all elements from `other` are
    /// appended to the end of `self` and global scalars are multiplied.
    pub fn tensor(mut self, mut other: Self) -> Self {
        let shift: usize =
            self.elems.iter()
            .map(|elem| {
                elem.input_iter().chain(elem.output_iter())
                    .max()
                    .map(|m| m + 1)
                    .unwrap_or(0)
            })
            .max()
            .unwrap_or(0);
        other.elems.iter_mut()
            .for_each(|elem| { elem.shift_indices(shift); });
        self.elems.append(&mut other.elems);
        self.scalar *= other.scalar;
        self
    }

    /// Fold over all elements of the diagram with the [dot
    /// product][Element::dot] and return the result as a single [`Element`].
    pub fn contract(self) -> Result<Element<A>, A::Error> {
        fn foldfn<B>(mb_acc: Option<Element<B>>, elem: Element<B>)
            -> Result<Option<Element<B>>, B::Error>
        where B: ElementData
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
                mb_res.unwrap_or_else(|| Element::scalar(1.0))
                    .scalar_mul(self.scalar)
            })
    }

    /// Convert `self` to a [`graph::Diagram`][crate::graph::Diagram]
    /// representation.
    ///
    /// Fails if any `Element` slices cannot be written as a pure generator or
    /// if any `Element` repeats a ket index without another `Element`'s
    /// matching bra index in between.
    pub fn to_graph(&self) -> TensorResult<graph::Diagram<graph::ZH>> {
        let mut graph = graph::Diagram::new();
        let mut node_id: graph::NodeId;
        let (inputs, outputs) = self.ins_outs();
        let mut wires: IndexMap<graph::NodeId> =
            IndexMap::new(); // wire idx -> node idx
        let mut to_remove: Vec<graph::NodeId> = Vec::new(); // placeholder nodes

        // add inputs
        for idx in inputs.into_iter() {
            node_id = graph.add_node(graph::ZHNode::Input);
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
                            Kind::Z(ph) => graph::ZHNode::Z(ph),
                            Kind::X(ph) => graph::ZHNode::X(ph),
                            Kind::H(a) => graph::ZHNode::H(a),
                            _ => unreachable!(),
                        };
                    node_id = graph.add_node(node);
                    for in_wire in elem.input_iter() {
                        let prev_id = wires.remove(in_wire).unwrap();
                        graph.add_wire(prev_id, node_id).unwrap();
                    }
                    for out_wire in elem.output_iter() {
                        if wires.insert(out_wire, node_id).is_some() {
                            return Err(DuplicateWire(out_wire));
                        }
                    }
                },
                Kind::Swap(a, b) => {
                    wires.swap(a, b);
                },
                Kind::Cup(a, b) => {
                    let n = graph.add_node(graph::ZHNode::z());
                    to_remove.push(n);
                    if wires.insert(a, n).is_some() {
                        return Err(DuplicateWire(a));
                    }
                    if wires.insert(b, n).is_some() {
                        return Err(DuplicateWire(b));
                    }
                },
                Kind::Cap(a, b) => {
                    let n1 = wires.remove(a).unwrap();
                    let n2 = wires.remove(b).unwrap();
                    graph.add_wire(n1, n2).unwrap();
                },
                Kind::Unknown => {
                    return Err(NonGenerator);
                },
            }
        }

        // add outputs
        for idx in outputs.into_iter() {
            node_id = graph.add_node(graph::ZHNode::Output);
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
    pub fn to_graphviz(&self) -> TensorResult<tabbycat::Graph> {
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

        // ensure all inputs are in a subgraph at the same rank, attaching the
        // overall scalar at the bottom
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
        inputs_subgraph_stmt =
            inputs_subgraph_stmt.add_node(node_id.into(), None, Some(attrs));
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
        statements =
            statements.add_subgraph(SubGraph::cluster(inputs_subgraph_stmt));

        // add elements/wires
        let mut u_idx = 0..;
        for elem in self.elems.iter() {
            match elem.kind {
                e if e.is_z() || e.is_x() || e.is_h() || e.is_unknown() => {
                    node_id = id_gen.next().unwrap();
                    let attrs: AttrList = match elem.kind {
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
                    for in_wire in elem.input_iter() {
                        let prev_id = wires.remove(in_wire).unwrap();
                        statements =
                            statements.add_edge(
                                Edge::head_node(prev_id.into(), None)
                                    .line_to_node(node_id.into(), None)
                            );
                    }
                    for out_wire in elem.output_iter() {
                        if wires.insert(out_wire, node_id).is_some() {
                            return Err(DuplicateWire(out_wire));
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
                _ => unreachable!(),
            }
        }

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
    pub fn save_graphviz<P>(&self, path: P) -> TensorResult<()>
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

/// A repeated, unmatched wire index in a [`Diagram`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Dup {
    /// The index is the output of an element.
    Ket(usize),
    /// The index is the input to an element.
    Bra(usize),
}

impl Dup {
    /// Return the bare wire index.
    pub fn wire_index(self) -> usize {
        match self {
            Self::Ket(k) => k,
            Self::Bra(k) => k,
        }
    }
}

impl Diagram<De> {
    fn to_sparse(&self) -> Diagram<Sp> {
        Diagram::<Sp>::from_elems(
            self.elems.iter()
                .map(|elem| elem.to_sparse())
        ).with_scalar(self.scalar)
    }
}

impl Diagram<Sp> {
    fn to_dense(&self) -> Diagram<De> {
        Diagram::<De>::from_elems(
            self.elems.iter()
                .map(|elem| elem.to_dense())
        ).with_scalar(self.scalar)
    }
}

/// Iterator over diagram elements.
///
/// The iterator item type is `&`[`Element`]`<A>`.
#[derive(Clone, Debug)]
pub struct Elems<'a, A> {
    iter: std::slice::Iter<'a, Element<A>>
}

impl<'a, A> Iterator for Elems<'a, A> {
    type Item = &'a Element<A>;

    fn next(&mut self) -> Option<Self::Item> { self.iter.next() }
}

impl<'a, A> DoubleEndedIterator for Elems<'a, A> {
    fn next_back(&mut self) -> Option<Self::Item> { self.iter.next_back() }
}

impl<'a, A> ExactSizeIterator for Elems<'a, A> {
    fn len(&self) -> usize { self.iter.len() }
}

impl<'a, A> std::iter::FusedIterator for Elems<'a, A> { }

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

/// Create a [tensor-based diagram][Diagram] using an abbreviated syntax.
///
/// The first item in angle brackets is the storage type for the diagram's
/// tensors, e.g. [`De`] or [`Sp`]. Everything that follows defines elements
/// in the diagram with the syntax
/// ```text
/// <element_type> ( args... )
/// ```
/// where `<element_type>` is a constructor method of [`Element`], and `args...`
/// is any arguments to be passed to it. Note that `args...` may be left empty.
///
/// This macro returns a single value of type `Diagram`.
///
/// The normal usage
/// ```
/// # use zx_calc::tensor::*;
/// let mut diagram = Diagram::<De>::empty();
///
/// diagram.push(Element::cup(1, 2));
/// diagram.push(Element::z_0([0], [0, 3]));
/// diagram.push(Element::x_0([1, 3], [1]));
/// diagram.push(Element::had(0));
/// diagram.push(Element::x_0([0], []));
/// diagram.push(Element::x_pi([1], []));
/// diagram.push(Element::x_pi([2], [2]));
/// diagram.push(Element::z_0([2], [2]));
/// ```
/// constructs the same diagram as
/// ```
/// # use zx_calc::tensor::*;
/// use zx_calc::tensor_diagram;
///
/// tensor_diagram!(
///     <De>
///     {
///         cup (1, 2),
///         z_0 ([0], [0, 3]),
///         x_0 ([1, 3], [1]),
///         had (0),
///         x_0 ([0], []),
///         x_pi ([1], []),
///         x_pi ([2], [2]),
///         z_0 ([2], [2]),
///     }
/// );
/// ```
#[macro_export]
macro_rules! tensor_diagram {
    (
        <$storage:ty>
        {
            $( $element_type:ident ( $( $arg:expr ),* $(,)? ) ),* $(,)?
        } $(,)?
    ) => {
        {
            let mut diagram = $crate::tensor::Diagram::<$storage>::empty();
            $(
            diagram.push(
                $crate::tensor::Element::<$storage>::$element_type($( $arg ),*)
            );
            )*
            diagram
        }
    }
}

pub use crate::tensor_diagram;

