use crate::phase::Phase;
use super::*;

/// A specific case of [`HEuler`], removing Hadamard wires connecting a spider
/// to itself and adding π to the spider's phase for each wire removed.
///
/// ![h_loop][h_loop]
#[embed_doc_image::embed_doc_image("h_loop", "assets/rules/HLoop.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HLoop;

/// Output of [`HLoop::find`].
#[derive(Debug)]
pub struct HLoopData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s: NodeId, // spider
    pub(crate) hh: Vec<NodeId>, // h-boxes
}

impl RuleSeal for HLoop { }
impl RuleFinder for HLoop {
    type Output<'a> = HLoopData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let mut hh: Vec<NodeId> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() {
                for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                    if node2.is_h()
                        && node2.has_defarg()
                        && dg.arity(id2).unwrap() == 2
                        && dg.mutual_arity(id, id2).unwrap() == 2
                    {
                        hh.push(id2);
                    }
                }
                if hh.is_empty() {
                    return Some(HLoopData { dg, s: id, hh });
                }
            }
        }
        None
    }
}

impl<'a> RuleSeal for HLoopData<'a> { }
impl<'a> Rule for HLoopData<'a> {
    fn simplify(self) {
        let Self { dg, s, hh } = self;
        let nloops = hh.len();
        hh.into_iter()
            .for_each(|h| { dg.remove_node(h).unwrap(); });
        dg.get_node_mut(s).unwrap()
            .map_phase(|ph| ph + (nloops as i64) * Phase::pi());
        dg.scalar *= std::f64::consts::FRAC_1_SQRT_2.powi(nloops as i32);
    }
}

