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
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct ColorFlipAll;

/// A spider surrounded by Hadamard wires.
///
/// The first item gives the node ID of the central spider whose color is to be
/// flipped; the second is its phase; and the third is a list of all surrounding
/// H-boxes' IDs paired with their outer neighbors' IDs.
pub type ColorFlipGroup = (NodeId, Phase, Vec<(NodeId, NodeId)>);

/// Output of [`ColorFlipAll::find`].
#[derive(Debug)]
pub struct ColorFlipAllData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) groups: Vec<ColorFlipGroup>,
}

impl<'a> ColorFlipAllData<'a> {
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

impl RuleSeal for ColorFlipAll { }
impl RuleFinder for ColorFlipAll {
    type Output<'a> = ColorFlipAllData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let mut seen: Vec<NodeId> = Vec::with_capacity(dg.node_count);
        let mut groups: Vec<ColorFlipGroup> = Vec::new();
        let mut hh: Vec<(NodeId, NodeId)> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if seen.contains(&id) {
                continue;
            } else {
                seen.push(id);
            }
            if node.is_spider() {
                for (id2, node2) in dg.neighbors_of(id).unwrap() {
                    if !seen.contains(&id2)
                        && node2.is_h()
                        && node2.has_defarg()
                        && dg.arity(id2).unwrap() == 2
                        && dg.mutual_arity(id, id2).unwrap() == 1
                    {
                        seen.push(id2);
                        let nb =
                            dg.find_neighbor_id_of(id2, |id3| id3 != id)
                            .unwrap();
                        hh.push((id2, nb));
                    } else {
                        hh.clear();
                        break;
                    }
                }
                if !hh.is_empty() {
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

impl<'a> RuleSeal for ColorFlipAllData<'a> { }
impl<'a> Rule for ColorFlipAllData<'a> {
    fn simplify(self) {
        let Self { dg, groups } = self;
        for (s, ph, hh) in groups.into_iter() {
            for (h, nb) in hh.into_iter() {
                dg.remove_node(h).unwrap();
                dg.add_wire(s, nb).unwrap();
            }
            if dg.get_node(s).unwrap().is_z() {
                let n = dg.get_node_mut(s).unwrap();
                *n = Node::X(ph);
            } else {
                let n = dg.get_node_mut(s).unwrap();
                *n = Node::Z(ph);
            }
        }
    }
}

