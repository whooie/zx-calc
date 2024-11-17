use crate::phase::Phase;
use super::*;

/// Replace a Z-spider surrounded by Hadamard wires with an X-spider (or
/// vice-versa).
///
/// ![color flip][color_flip]
#[embed_doc_image::embed_doc_image("color_flip", "assets/rules/ColorFlip.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct ColorFlip;

/// Output of [`ColorFlip::find`].
#[derive(Debug)]
pub struct ColorFlipData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) s: (NodeId, Phase), // spider
    pub(crate) hh: Vec<(NodeId, NodeId)>, // ZH: [(h-box, neighbor)]
                                          // ZX: [empty]
}

impl RuleFinder<ZX> for ColorFlip {
    type Output<'a> = ColorFlipData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if dg.neighbor_ids(id).unwrap().all(|w| w.is_h()) {
                let s = (id, node.phase().unwrap());
                return Some(ColorFlipData { dg, s, hh: Vec::with_capacity(0) });
            }
        }
        None
    }
}

impl<'a> Rule<ZX> for ColorFlipData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, s, hh: _ } = self;
        dg.get_neighbors_mut(s.0).unwrap()
            .iter_mut()
            .for_each(|wire| { wire.toggle(); });
        if dg.get_node(s.0).unwrap().is_z() {
            let n = dg.get_node_mut(s.0).unwrap();
            *n = ZXNode::X(s.1);
        } else {
            let n = dg.get_node_mut(s.0).unwrap();
            *n = ZXNode::Z(s.1);
        }
    }
}

impl RuleFinder<ZH> for ColorFlip {
    type Output<'a> = ColorFlipData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let mut hh: Vec<(NodeId, NodeId)> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() {
                for (id2, node2) in dg.neighbors(id).unwrap() {
                    let id2 = id2.id();
                    if node2.is_h()
                        && node2.has_defarg()
                        && dg.arity(id2).unwrap() == 2
                        && dg.mutual_arity(id, id2).unwrap() == 1
                    {
                        let nb =
                            dg.find_neighbor_id(id2, |w| !w.has_id(id))
                            .unwrap();
                        hh.push((id2, nb.id()));
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

impl<'a> Rule<ZH> for ColorFlipData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, s, hh } = self;
        for (h, nb) in hh.into_iter() {
            dg.remove_node(h).unwrap();
            dg.add_wire(s.0, nb).unwrap();
        }
        if dg.get_node(s.0).unwrap().is_z() {
            let n = dg.get_node_mut(s.0).unwrap();
            *n = ZHNode::X(s.1);
        } else {
            let n = dg.get_node_mut(s.0).unwrap();
            *n = ZHNode::Z(s.1);
        }
    }
}

