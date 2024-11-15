use crate::phase::Phase;
use super::*;

/// Remove all unary X-spiders with phase π connected to an H-box.
///
/// ![h_absorb][h_absorb]
#[embed_doc_image::embed_doc_image("h_absorb", "assets/rules/HAbsorb.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HAbsorb;

/// Output of [`HAbsorb::find`].
#[derive(Debug)]
pub struct HAbsorbData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) states: Vec<NodeId>, // pi x-spiders
}

impl RuleFinder<ZH> for HAbsorb {
    type Output<'a> = HAbsorbData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let mut states: Vec<NodeId> = Vec::new();
        let pi = Phase::pi();
        for (id, node) in dg.nodes_inner() {
            if node.is_h() {
                dg.neighbors_inner(id).unwrap()
                    .filter_map(|(id2, node2)| {
                        (
                            node2.is_x_and(|ph| ph == pi)
                            && dg.arity(*id2).unwrap() == 1
                        ).then_some(*id2)
                    })
                    .for_each(|id2| { states.push(id2); });
                if !states.is_empty() {
                    return Some(HAbsorbData { dg, states });
                }
            }
        }
        None
    }
}

impl<'a> Rule<ZH> for HAbsorbData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, states } = self;
        let nstates = states.len();
        states.into_iter()
            .for_each(|id| { dg.remove_node(id).unwrap(); });
        dg.scalar *= std::f64::consts::SQRT_2.powi(nstates as i32);
    }
}

