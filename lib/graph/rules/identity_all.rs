use super::*;

/// Replace all phaseless binary spiders or adjacent binary H-boxes with empty
/// wires.
///
/// This is the comprehensive version of [`Identity`], searching more
/// efficiently for all such spiders or H-boxes.
///
/// ![identity][identity]
#[embed_doc_image::embed_doc_image("identity", "assets/rules/Identity.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct IdentityAll;

/// Output of [`IdentityAll::find`].
#[derive(Debug)]
pub struct IdentityAllData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    // all identities; guaranteed non-overlapping
    pub(crate) idents: Vec<IdentityKind>,
}

impl<'a, A> IdentityAllData<'a, A>
where A: DiagramData
{
    /// Return the number of identity groups found.
    ///
    /// All groups are guaranteed to be disconnected from each other.
    pub fn len(&self) -> usize { self.idents.len() }

    /// Return `true` if the number of identity groups is zero.
    pub fn is_empty(&self) -> bool { self.idents.is_empty() }

    /// Return a reference to all identity groups.
    pub fn groups(&self) -> &Vec<IdentityKind> { &self.idents }
}

impl RuleFinder<ZX> for IdentityAll {
    type Output<'a> = IdentityAllData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        let mut seen: Vec<NodeId> = Vec::new();
        let mut idents: Vec<IdentityKind> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if seen.contains(&id) {
                continue;
            } else {
                seen.push(id);
            }
            if node.has_defarg()
                && dg.arity(id).unwrap() == 2
                && dg.neighbor_ids(id).unwrap()
                    .all(|nb| !seen.contains(&nb.id()))
                && (
                    dg.arity_e(id).unwrap() == 2
                    || (
                        dg.arity_h(id).unwrap() == 2
                        && dg.mutual_arity(id, id).unwrap() == 0
                    )
                )
            {
                idents.push(IdentityKind::Spider(id));
            }
        }
        if idents.is_empty() {
            None
        } else {
            Some(IdentityAllData { dg, idents })
        }
    }
}

impl<'a> Rule<ZX> for IdentityAllData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, idents } = self;
        for kind in idents.into_iter() {
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
}

impl RuleFinder<ZH> for IdentityAll {
    type Output<'a> = IdentityAllData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let mut seen: Vec<NodeId> = Vec::new();
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
                    for (id2, node2) in dg.neighbors_inner(id).unwrap() {
                        if seen.contains(id2) { continue; }
                        if node2.is_h()
                            && node2.has_defarg()
                            && dg.arity(*id2).unwrap() == 2
                        {
                            seen.push(*id2);
                            idents.push(IdentityKind::HBox(id, *id2));
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

impl<'a> Rule<ZH> for IdentityAllData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, idents } = self;
        for kind in idents.into_iter() {
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
                        dg.scalar *= 2.0;
                    }
                }
            }
        }
    }
}

