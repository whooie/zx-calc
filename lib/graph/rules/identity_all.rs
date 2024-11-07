use super::*;

/// Replace all phaseless binary spiders or adjacent binary H-boxes with empty
/// wires.
///
/// This is the comprehensive version of [`Identity`], searching more
/// efficiently for all such spiders or H-boxes.
///
/// ![identity][identity]
#[embed_doc_image::embed_doc_image("identity", "assets/rules/Identity.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct IdentityAll;

/// Output of [`IdentityAll::find`].
#[derive(Debug)]
pub struct IdentityAllData<'a> {
    pub(crate) dg: &'a mut Diagram,
    // all identities; guaranteed non-overlapping
    pub(crate) idents: Vec<IdentityKind>,
}

impl<'a> IdentityAllData<'a> {
    /// Return the number of identity groups found.
    ///
    /// All groups are guaranteed to be disconnected from each other.
    pub fn len(&self) -> usize { self.idents.len() }

    /// Return `true` if the number of identity groups is zero.
    pub fn is_empty(&self) -> bool { self.idents.is_empty() }

    /// Return a reference to all identity groups.
    pub fn groups(&self) -> &Vec<IdentityKind> { &self.idents }
}

impl RuleSeal for IdentityAll { }
impl RuleFinder for IdentityAll {
    type Output<'a> = IdentityAllData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let mut seen: Vec<NodeId> = Vec::with_capacity(dg.node_count);
        let mut idents: Vec<IdentityKind> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if seen.contains(&id) {
                continue;
            } else {
                seen.push(id);
            }
            if node.has_defarg() && dg.arity(id).unwrap() == 2 {
                if node.is_spider() {
                    idents.push(IdentityKind::Spider(id));
                } else if node.is_h() {
                    for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                        if seen.contains(&id2) { continue; }
                        if node2.is_h()
                            && node2.has_defarg()
                            && dg.arity(id2).unwrap() == 2
                        {
                            seen.push(id2);
                            idents.push(IdentityKind::HBox(id, id2));
                        }
                    }
                }
            }
        }
        if idents.is_empty() {
            None
        } else {
            Some(IdentityAllData { dg, idents })
        }
    }
}

impl<'a> RuleSeal for IdentityAllData<'a> { }
impl<'a> Rule for IdentityAllData<'a> {
    fn simplify(self) {
        let Self { dg, idents } = self;
        for kind in idents.into_iter() {
            match kind {
                IdentityKind::Spider(s) => {
                    let mut nb_iter = dg.neighbors_of(s).unwrap().map(fst);
                    let nb1 = nb_iter.next().unwrap();
                    if let Some(nb2) = nb_iter.next() {
                        // not a self-neighbor
                        dg.remove_node(s).unwrap();
                        dg.add_wire(nb1, nb2).unwrap();
                    } else {
                        // self-neighbor
                        dg.remove_node(s).unwrap();
                    }
                },
                IdentityKind::HBox(h1, h2) => {
                    let mut nb_iter =
                        dg.neighbor_ids_of(h1).unwrap()
                        .chain(dg.neighbor_ids_of(h2).unwrap())
                        .filter(|id| *id != h1 && *id != h2);
                    if let Some(nb1) = nb_iter.next() {
                        // not a self-loop
                        let nb2 = nb_iter.next().unwrap();
                        dg.remove_node(h1).unwrap();
                        dg.remove_node(h2).unwrap();
                        if nb1 != h1 && nb1 != h2 && nb2 != h1 && nb2 != h2 {
                            dg.add_wire(nb1, nb2).unwrap();
                        }
                    } else {
                        // self-loop
                        dg.remove_node(h1).unwrap();
                        dg.remove_node(h2).unwrap();
                    }
                }
            }
        }
    }
}

