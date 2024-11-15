use super::*;

/// Replace a phaseless binary spider or two adjacent binary H-boxes with an
/// empty wire.
///
/// ![identity][identity]
#[embed_doc_image::embed_doc_image("identity", "assets/rules/Identity.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Identity;

/// Nodes that can be simplified to an empty wire.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum IdentityKind {
    /// A phaseless binary spider surrounded by either empty or Hadamard wires.
    Spider(NodeId),
//  Spider(spider)
    /// Two adjacent binary H-boxes with default argument.
    HBox(NodeId, NodeId),
//  HBox(hbox1, hbox2)
}

/// Output of [`Identity::find`].
#[derive(Debug)]
pub struct IdentityData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) kind: IdentityKind,
}

impl RuleFinder<ZX> for Identity {
    type Output<'a> = IdentityData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if node.has_defarg()
                && dg.arity(id).unwrap() == 2
                && (
                    dg.arity_e(id).unwrap() == 2
                    || (
                        dg.arity_h(id).unwrap() == 2
                        // don't consider self-loops via Hadamard wires
                        && dg.mutual_arity(id, id).unwrap() == 0
                    )
                )
            {
                let kind = IdentityKind::Spider(id);
                return Some(IdentityData { dg, kind });
            }
        }
        None
    }
}

impl<'a> Rule<ZX> for IdentityData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, kind } = self;
        let IdentityKind::Spider(s) = kind else { unreachable!() };
        let mut nb_iter = dg.neighbor_ids(s).unwrap().copied();
        let nb1 = nb_iter.next().unwrap();
        if let Some(nb2) = nb_iter.next() {
            // not a self-neighbor
            dg.remove_node(s).unwrap();
            dg.add_wire(nb1.id(), nb2.id()).unwrap();
        } else {
            // self-neighbor
            assert_eq!(nb1.id(), s);
            dg.remove_node(s).unwrap();
            dg.scalar *= if nb1.is_h() { 0.0 } else { 2.0 };
        }
    }
}

impl RuleFinder<ZH> for Identity {
    type Output<'a> = IdentityData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if node.has_defarg() && dg.arity(id).unwrap() == 2 {
                if node.is_spider() {
                    let kind = IdentityKind::Spider(id);
                    return Some(IdentityData { dg, kind });
                } else if node.is_h() {
                    for (id2, node2) in dg.neighbors_inner(id).unwrap() {
                        if node2.is_h()
                            && node2.has_defarg()
                            && dg.arity(*id2).unwrap() == 2
                        {
                            let kind = IdentityKind::HBox(id, *id2);
                            return Some(IdentityData { dg, kind });
                        }
                    }
                }
            }
        }
        None
    }
}

impl<'a> Rule<ZH> for IdentityData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, kind } = self;
        match kind {
            IdentityKind::Spider(s) => {
                let mut nb_iter = dg.neighbor_ids(s).unwrap().copied();
                let nb1 = nb_iter.next().unwrap();
                if let Some(nb2) = nb_iter.next() {
                    // not a self-neighbor
                    dg.remove_node(s).unwrap();
                    dg.add_wire(nb1, nb2).unwrap();
                } else {
                    // self-neighbor
                    dg.remove_node(s).unwrap();
                    dg.scalar *= 2.0;
                }
            },
            IdentityKind::HBox(h1, h2) => {
                let mut nb_iter =
                    dg.neighbor_ids(h1).unwrap()
                    .chain(dg.neighbor_ids(h2).unwrap())
                    .filter(|id| **id != h1 && **id != h2)
                    .copied();
                if let Some(nb1) = nb_iter.next() {
                    let nb2 = nb_iter.next().unwrap();
                    dg.remove_node(h1).unwrap();
                    dg.remove_node(h2).unwrap();
                    if nb1 != h1 && nb1 != h2 && nb2 != h1 && nb2 != h2 {
                        dg.add_wire(nb1, nb2).unwrap();
                    }
                } else {
                    dg.remove_node(h1).unwrap();
                    dg.remove_node(h2).unwrap();
                    dg.scalar *= 2.0;
                }
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simplify_identity() {
        let mut dg: Diagram<ZH> = Diagram::new();
        let z0 = dg.add_node(ZHNode::z());
        let z1 = dg.add_node(ZHNode::z_pi());
        let z2 = dg.add_node(ZHNode::z());
        let z3 = dg.add_node(ZHNode::z());
        let z4 = dg.add_node(ZHNode::z());
        let z5 = dg.add_node(ZHNode::z());
        let x0 = dg.add_node(ZHNode::x());
        let x1 = dg.add_node(ZHNode::x_pi());
        let x2 = dg.add_node(ZHNode::x());
        let h0 = dg.add_node(ZHNode::h());
        let h1 = dg.add_node(ZHNode::h());
        let h2 = dg.add_node(ZHNode::H(1.0.into()));
        let h3 = dg.add_node(ZHNode::h());
        let h4 = dg.add_node(ZHNode::h());
        dg.add_input_wire(z0).unwrap();
        dg.add_output_wire(x0).unwrap();
        dg.add_wire(z0, z1).unwrap();
        dg.add_wire(z1, z2).unwrap();
        dg.add_wire(z2, z3).unwrap();
        dg.add_wire(z2, x0).unwrap();
        dg.add_wire(z4, z4).unwrap();
        dg.add_wire(z5, x1).unwrap();
        dg.add_wire(z5, x1).unwrap();
        dg.add_input_wire(h0).unwrap();
        dg.add_input_wire(h2).unwrap();
        dg.add_wire(h0, h1).unwrap();
        dg.add_wire(h1, h2).unwrap();
        dg.add_wire(h3, h4).unwrap();
        dg.add_wire(h3, x2).unwrap();
        dg.add_wire(x2, h4).unwrap();
        assert_eq!(dg.count_nodes(), 18);
        assert_eq!(dg.count_wires(), 16);
        assert_eq!(dg.count_z(),     6);
        assert_eq!(dg.count_x(),     3);
        assert_eq!(dg.count_h(),     5);
        assert_eq!(dg.simplify_rule_n(Identity, None), 7);
        assert_eq!(dg.count_nodes(), 9);
        assert_eq!(dg.count_wires(), 7);
        assert_eq!(dg.count_z(),     3);
        assert_eq!(dg.count_x(),     1);
        assert_eq!(dg.count_h(),     1);
    }
}

