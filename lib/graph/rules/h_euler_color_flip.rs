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
pub struct HEulerColorFlipData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) s: (NodeId, Phase), // central spider
    // [(pi/2 spider, phase, outer neighbor)]
    pub(crate) ss_pi2: Vec<(NodeId, Phase, NodeId)>,
}

impl RuleFinder<ZX> for HEulerColorFlip {
    type Output<'a> = HEulerColorFlipData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        let mut ss_pi2: Vec<(NodeId, Phase, NodeId)> = Vec::new();
        let pos_pi2 = Phase::pi2();
        let neg_pi2 = -Phase::pi2();
        for (id, node) in dg.nodes_inner() {
            for (wire, node2) in dg.neighbors(id).unwrap() {
                let id2 = wire.id();
                if wire.is_e()
                    && node2.is_diff_color_and(
                        node, |_, ph2| ph2 == pos_pi2 || ph2 == neg_pi2)
                    && dg.arity(id2).unwrap() == 2
                    && dg.arity_e(id2).unwrap() == 2
                    && dg.mutual_arity(id, id2).unwrap() == 1
                {
                    let nb =
                        dg.find_neighbor(id2, |w, _| !w.has_id(id))
                        .unwrap();
                    if nb.1.is_same_color(node) {
                        let ph2 = node2.phase().unwrap();
                        ss_pi2.push((id2, ph2, nb.0.id()));
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
        None
    }
}

impl<'a> Rule<ZX> for HEulerColorFlipData<'a, ZX> {
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
            *n = ZXNode::X(s.1 - ph_sum);
        } else {
            let n = dg.get_node_mut(s.0).unwrap();
            *n = ZXNode::Z(s.1 - ph_sum);
        }
        dg.scalar *= (ph_sum / 2).cis();
    }
}

impl RuleFinder<ZH> for HEulerColorFlip {
    type Output<'a> = HEulerColorFlipData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let mut ss_pi2: Vec<(NodeId, Phase, NodeId)> = Vec::new();
        let pos_pi2 = Phase::pi2();
        let neg_pi2 = -Phase::pi2();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() {
                for (id2, node2) in dg.neighbors(id).unwrap() {
                    let id2 = id2.id();
                    if node2.is_diff_color_and(
                        node, |_, ph2| ph2 == pos_pi2 || ph2 == neg_pi2)
                        && dg.arity(id2).unwrap() == 2
                        && dg.mutual_arity(id, id2).unwrap() == 1
                    {
                        let nb =
                            dg.find_neighbor(id2, |id3, _| *id3 != id)
                            .unwrap();
                        if nb.1.is_same_color(node) {
                            let ph2 = node2.phase().unwrap();
                            ss_pi2.push((id2, ph2, *nb.0));
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

impl<'a> Rule<ZH> for HEulerColorFlipData<'a, ZH> {
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
            *n = ZHNode::X(s.1 - ph_sum);
        } else {
            let n = dg.get_node_mut(s.0).unwrap();
            *n = ZHNode::Z(s.1 - ph_sum);
        }
        dg.scalar *= (ph_sum / 2).cis();
    }
}

