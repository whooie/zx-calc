use crate::phase::Phase;
use super::*;

/// Replace a Z-spider surrounded by Hadamard wires with an X-spider (or
/// vice-versa).
///
/// ![color flip][color_flip]
#[embed_doc_image::embed_doc_image("color_clip", "assets/rules/ColorFlip.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct ColorFlip;

/// Output of [`ColorFlip::find`].
#[derive(Debug)]
pub struct ColorFlipData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s: (NodeId, Phase), // spider
    pub(crate) hh: Vec<(NodeId, NodeId)>, // [(h-box, neighbor)]
}

impl RuleSeal for ColorFlip { }
impl RuleFinder for ColorFlip {
    type Output<'a> = ColorFlipData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let mut hh: Vec<(NodeId, NodeId)> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() {
                for (id2, node2) in dg.neighbors_of(id).unwrap() {
                    if node2.is_h()
                        && node2.has_defarg()
                        && dg.arity(id2).unwrap() == 2
                        && dg.mutual_arity(id, id2).unwrap() == 1
                    {
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
                    let s = (id, node.phase().unwrap());
                    return Some(ColorFlipData { dg, s, hh });
                }
            }
        }
        None
    }
}

impl<'a> RuleSeal for ColorFlipData<'a> { }
impl<'a> Rule for ColorFlipData<'a> {
    fn simplify(self) {
        let Self { dg, s, hh } = self;
        for (h, nb) in hh.into_iter() {
            dg.remove_node(h).unwrap();
            dg.add_wire(s.0, nb).unwrap();
        }
        if dg.get_node(s.0).unwrap().is_z() {
            let n = dg.get_node_mut(s.0).unwrap();
            *n = Node::X(s.1);
        } else {
            let n = dg.get_node_mut(s.0).unwrap();
            *n = Node::Z(s.1);
        }
    }
}

