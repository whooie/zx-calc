use super::*;

/// Either a phaseless spider with two wires or two adjacent H-boxes of arity 2
/// and default argument.
#[derive(Copy, Clone, Debug)]
pub struct Identity;

/// Output of [`Identity::find`].
#[derive(Debug)]
pub struct IdentityData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) kind: IdentityKind,
}

#[derive(Debug)]
pub(crate) enum IdentityKind {
    Spider(NodeId),
//  Spider(spider)
    HBox(NodeId, NodeId),
//  HBox(hbox1, hbox2)
}

impl RuleSeal for Identity { }
impl RuleFinder for Identity {
    type Output<'a> = IdentityData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if node.has_defarg() && dg.arity(id).unwrap() == 2 {
                if node.is_spider() {
                    let kind = IdentityKind::Spider(id);
                    return Some(IdentityData { dg, kind });
                } else if node.is_h() {
                    for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                        if node2.is_h()
                            && node2.has_defarg()
                            && dg.arity(id2).unwrap() == 2
                        {
                            let kind = IdentityKind::HBox(id, id2);
                            return Some(IdentityData { dg, kind });
                        }
                    }
                }
            }
        }
        None
    }
}

impl<'a> RuleSeal for IdentityData<'a> { }
impl<'a> Rule for IdentityData<'a> {
    fn simplify(self) {
        let Self { dg, kind } = self;
        match kind {
            IdentityKind::Spider(s) => {
                let mut nb_iter = dg.neighbors_of(s).unwrap().map(fst);
                let nb1 = nb_iter.next().unwrap();
                if let Some(nb2) = nb_iter.next() {
                    // not a dg-neighbor
                    dg.remove_node(s).unwrap();
                    dg.add_wire(nb1, nb2).unwrap();
                } else {
                    // dg-neighbor
                    dg.remove_node(s).unwrap();
                }
            },
            IdentityKind::HBox(h1, h2) => {
                let mut nb_iter =
                    dg.neighbors_of(h1).unwrap()
                    .chain(dg.neighbors_of(h2).unwrap())
                    .filter_map(|(id, _)| {
                        (id != h1 && id != h2).then_some(id)
                    });
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
        let mut dg = Diagram::new();
        let z0 = dg.add_node(Node::z());
        let z1 = dg.add_node(Node::z_pi());
        let z2 = dg.add_node(Node::z());
        let z3 = dg.add_node(Node::z());
        let z4 = dg.add_node(Node::z());
        let z5 = dg.add_node(Node::z());
        let x0 = dg.add_node(Node::x());
        let x1 = dg.add_node(Node::x_pi());
        let x2 = dg.add_node(Node::x());
        let h0 = dg.add_node(Node::h());
        let h1 = dg.add_node(Node::h());
        let h2 = dg.add_node(Node::H(1.0.into()));
        let h3 = dg.add_node(Node::h());
        let h4 = dg.add_node(Node::h());
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

