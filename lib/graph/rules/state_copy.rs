use crate::phase::Phase;
use super::*;

/// Copy a unary spider representing a pure Z- or X-basis state through an
/// oppositely colored spider.
///
/// ![state_copy][state_copy]
#[embed_doc_image::embed_doc_image("state_copy", "assets/rules/StateCopy.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct StateCopy;

/// Output of [`StateCopy::find`].
#[derive(Debug)]
pub struct StateCopyData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) state: NodeId,
    pub(crate) spider: NodeId,
}

impl RuleFinder for StateCopy {
    type Output<'a> = StateCopyData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let zero = Phase::zero();
        let pi = Phase::pi();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() {
                for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                    if dg.arity(id2).unwrap() == 1
                        && node.is_diff_color_and(
                            node2, |_, ph2| ph2 == zero || ph2 == pi
                        )
                    {
                        let state = id2;
                        let spider = id;
                        return Some(StateCopyData { dg, state, spider});
                    }
                }
            }
        }
        None
    }
}

impl<'a> Rule for StateCopyData<'a> {
    fn simplify(self) {
        let Self { dg, state, spider } = self;
        let n = dg.remove_node(state).unwrap();
        let (state, nnb) = dg.remove_node_nb(spider).unwrap();
        let len_nnb = nnb.len();
        let state_phase = state.phase().unwrap();
        let spider_phase = dg.get_node(spider).unwrap().phase().unwrap();
        for nb in nnb.into_iter() {
            if nb == spider { continue; }
            let new = dg.add_node(n);
            dg.add_wire(new, nb).unwrap();
        }
        dg.scalar *=
            (
                (1.0 + spider_phase.cis())
                + state_phase.cis() * (1.0 - spider_phase.cis())
            ) * std::f64::consts::FRAC_1_SQRT_2.powi(len_nnb as i32 + 1);
    }
}

