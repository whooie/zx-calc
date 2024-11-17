use crate::phase::Phase;
use super::*;

/// Change the colors of all spiders surrounded by Hadamard wires.
///
/// This rule is a comprehensive version of [`ColorFlip`], searching more
/// efficiently for all spider neighborhoods that satisfy the conditions for the
/// rule.
///
/// ![color_flip][color_flip]
#[embed_doc_image::embed_doc_image("color_flip", "assets/rules/ColorFlip.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct ColorFlipAll;

/// A spider surrounded by Hadamard wires.
///
/// The first item gives the node ID of the central spider whose color is to be
/// flipped; the second is its phase; and the third is a list of all surrounding
/// H-boxes' IDs paired with their outer neighbors' IDs.
pub type ColorFlipGroup = (NodeId, Phase, Vec<(NodeId, NodeId)>);

/// Output of [`ColorFlipAll::find`].
#[derive(Debug)]
pub struct ColorFlipAllData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) groups: Vec<ColorFlipGroup>,
}

impl<'a, A> ColorFlipAllData<'a, A>
where A: DiagramData
{
    /// Return the number of flippable spiders found.
    pub fn len(&self) -> usize { self.groups.len() }

    /// Return `true` if the number of groups is zero.
    pub fn is_empty(&self) -> bool { self.groups.is_empty() }

    /// Return a reference to all groups.
    ///
    /// All central spiders and H-boxes are guaranteed to be disconnected from
    /// other groups' spiders and H-boxes. Outer neighbors may overlap, however.
    pub fn groups(&self) -> &Vec<ColorFlipGroup> { &self.groups }
}

impl RuleFinder<ZX> for ColorFlipAll {
    type Output<'a> = ColorFlipAllData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        let mut seen: Vec<NodeId> = Vec::new();
        let mut groups: Vec<ColorFlipGroup> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if seen.contains(&id) {
                continue;
            } else {
                seen.push(id);
            }
            if dg.neighbor_ids(id).unwrap().all(|w| w.is_h()) {
                dg.neighbor_ids(id).unwrap()
                    .for_each(|w| { seen.push(w.id()); });
                groups.push((id, node.phase().unwrap(), Vec::with_capacity(0)));
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(ColorFlipAllData { dg, groups })
        }
    }
}

impl<'a> Rule<ZX> for ColorFlipAllData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, groups } = self;
        for (s, ph, _) in groups.into_iter() {
            dg.get_neighbors_mut(s).unwrap()
                .iter_mut()
                .for_each(|wire| { wire.toggle(); });
            if dg.get_node(s).unwrap().is_z() {
                let n = dg.get_node_mut(s).unwrap();
                *n = ZXNode::X(ph);
            } else {
                let n = dg.get_node_mut(s).unwrap();
                *n = ZXNode::Z(ph);
            }
        }
    }
}

impl RuleFinder<ZH> for ColorFlipAll {
    type Output<'a> = ColorFlipAllData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let mut seen: Vec<NodeId> = Vec::new();
        let mut groups: Vec<ColorFlipGroup> = Vec::new();
        let mut hh: Vec<(NodeId, NodeId)> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if seen.contains(&id) {
                continue;
            } else {
                seen.push(id);
            }
            if node.is_spider() {
                let all_h =
                    dg.neighbors(id).unwrap()
                    .all(|(wire, node2)| {
                        let id2 = wire.id();
                        !seen.contains(&id2)
                            && node2.is_h()
                            && node2.has_defarg()
                            && dg.arity(id2).unwrap() == 2
                            && dg.mutual_arity(id, id2).unwrap() == 1
                    });
                if all_h {
                    dg.neighbor_ids(id).unwrap()
                        .for_each(|wire| {
                            let id2 = wire.id();
                            seen.push(id2);
                            let nb =
                                dg.find_neighbor_id(id2, |w| !w.has_id(id))
                                .unwrap();
                            hh.push((id2, nb.id()));
                        });
                    let ph = node.phase().unwrap();
                    let h = std::mem::take(&mut hh);
                    groups.push((id, ph, h));
                }
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(ColorFlipAllData { dg, groups })
        }
    }
}

impl<'a> Rule<ZH> for ColorFlipAllData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, groups } = self;
        for (s, ph, hh) in groups.into_iter() {
            for (h, nb) in hh.into_iter() {
                dg.remove_node(h).unwrap();
                dg.add_wire(s, nb).unwrap();
            }
            if dg.get_node(s).unwrap().is_z() {
                let n = dg.get_node_mut(s).unwrap();
                *n = ZHNode::X(ph);
            } else {
                let n = dg.get_node_mut(s).unwrap();
                *n = ZHNode::Z(ph);
            }
        }
    }
}

