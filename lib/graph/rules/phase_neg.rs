use crate::phase::Phase;
use super::*;

/// Flip the sign of a spider's phase when it's surrounded by oppositely colored
/// binary spiders with phase π.
///
/// ![phase_neg][phase_neg]
#[embed_doc_image::embed_doc_image("phase_neg", "assets/rules/PhaseNeg.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct PhaseNeg;

/// Output of [`PhaseNeg::find`].
#[derive(Debug)]
pub struct PhaseNegData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s: NodeId, // center spider
    pub(crate) ss_pi: Vec<(NodeId, NodeId)>, // [(pi spider, neighbor)]
}

impl RuleSeal for PhaseNeg { }
impl RuleFinder for PhaseNeg {
    type Output<'a> = PhaseNegData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let pi = Phase::pi();
        let mut ss_pi: Vec<(NodeId, NodeId)> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() {
                for (id2, node2) in dg.neighbors_of(id).unwrap() {
                    if node.is_diff_color_and(node2, |_, ph2| ph2 == pi)
                        && dg.arity(id2).unwrap() == 2
                        && dg.mutual_arity(id, id2).unwrap() == 1
                    {
                        let nb =
                            dg.find_neighbor_id_of(id2, |id3| id3 != id)
                            .unwrap();
                        ss_pi.push((id2, nb));
                    } else {
                        ss_pi.clear();
                        break;
                    }
                }
                if !ss_pi.is_empty() {
                    return Some(PhaseNegData { dg, s: id, ss_pi });
                }
            }
        }
        None
    }
}

impl<'a> RuleSeal for PhaseNegData<'a> { }
impl<'a> Rule for PhaseNegData<'a> {
    fn simplify(self) {
        let Self { dg, s, ss_pi } = self;
        for (s_pi, nb) in ss_pi.into_iter() {
            dg.remove_node(s_pi).unwrap();
            dg.add_wire(s, nb).unwrap();
        }
        dg.get_node_mut(s).unwrap()
            .map_phase(|ph| -ph);
        let ph = dg.get_node(s).unwrap().phase().unwrap();
        dg.scalar *= ph.cis();
    }
}

