use crate::phase::Phase;
use super::*;

/// Apply [`HEuler`] in reverse together with [`ColorFlip`] in tandem, flipping
/// the color of a central spider surrounded by oppositely colored spiders with
/// phase ±π/2.
///
/// ![h_euler_color_flip][h_euler_color_flip]
#[embed_doc_image::embed_doc_image("h_euler_color_flip", "assets/rules/HEulerColorFlip.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HEulerColorFlip;

/// Output of [`HEulerColorFlip::find`].
#[derive(Debug)]
pub struct HEulerColorFlipData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s: (NodeId, Phase), // central spider
    // [(pi/2 spider, phase, outer neighbor)]
    pub(crate) ss_pi2: Vec<(NodeId, Phase, NodeId)>,
}

impl RuleSeal for HEulerColorFlip { }
impl RuleFinder for HEulerColorFlip {
    type Output<'a> = HEulerColorFlipData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let mut ss_pi2: Vec<(NodeId, Phase, NodeId)> = Vec::new();
        let pos_pi2 = Phase::pi2();
        let neg_pi2 = -Phase::pi2();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() {
                for (id2, node2) in dg.neighbors_of(id).unwrap() {
                    if node2.is_diff_color_and(
                        node, |_, ph2| ph2 == pos_pi2 || ph2 == neg_pi2)
                        && dg.arity(id2).unwrap() == 2
                        && dg.mutual_arity(id, id2).unwrap() == 1
                    {
                        let nb =
                            dg.find_neighbor_of(id2, |id3, _| id3 != id)
                            .unwrap();
                        if nb.1.is_same_color(node) {
                            let ph2 = node2.phase().unwrap();
                            ss_pi2.push((id2, ph2, nb.0));
                        } else {
                            ss_pi2.clear();
                            break;
                        }
                    } else {
                        ss_pi2.clear();
                        break;
                    }
                }
                if !ss_pi2.is_empty() {
                    let s = (id, node.phase().unwrap());
                    return Some(HEulerColorFlipData { dg, s, ss_pi2 });
                }
            }
        }
        None
    }
}

impl<'a> RuleSeal for HEulerColorFlipData<'a> { }
impl<'a> Rule for HEulerColorFlipData<'a> {
    fn simplify(self) {
        let Self { dg, s, ss_pi2 } = self;
        let mut ph_sum = Phase::zero();
        for (s_pi2, ph, nb) in ss_pi2.into_iter() {
            dg.remove_node(s_pi2).unwrap();
            dg.add_wire(s.0, nb).unwrap();
            dg.get_node_mut(nb).unwrap()
                .map_phase(|ph_nb| ph_nb - ph);
            ph_sum += ph;
        }
        if dg.get_node(s.0).unwrap().is_z() {
            let n = dg.get_node_mut(s.0).unwrap();
            *n = Node::X(s.1 - ph_sum);
        } else {
            let n = dg.get_node_mut(s.0).unwrap();
            *n = Node::Z(s.1 - ph_sum);
        }
    }
}

