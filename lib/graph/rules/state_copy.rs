use crate::phase::Phase;
use super::*;

/// Copy a unary spider representing a pure Z- or X-basis state through an
/// oppositely colored spider.
///
/// ![state_copy][state_copy]
#[embed_doc_image::embed_doc_image("state_copy", "assets/rules/StateCopy.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct StateCopy;

/// Output of [`StateCopy::find`].
#[derive(Debug)]
pub struct StateCopyData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) state: NodeId,
    pub(crate) spider: NodeId,
}

impl RuleFinder<ZX> for StateCopy {
    type Output<'a> = StateCopyData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            for (wire, node2) in dg.neighbors_inner(id).unwrap() {
                if dg.arity(wire.id()).unwrap() == 1
                    && (
                        (
                            wire.is_e()
                            && node.is_diff_color_and(
                                node2, |_, ph2| ph2.is_mult(2))
                        ) || (
                            wire.is_h()
                            && node.is_same_color_and(
                                node2, |_, ph2| ph2.is_mult(2))
                        )
                    )
                {
                    let state = wire.id();
                    let spider = id;
                    return Some(StateCopyData { dg, state, spider });
                }
            }
        }
        None
    }
}

impl<'a> Rule<ZX> for StateCopyData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, state, spider } = self;
        let (state_n, wire) = dg.remove_node_nb(state).unwrap();
        let (spider_n, nnb) = dg.remove_node_nb(spider).unwrap();
        let len_nnb = nnb.len();
        let state_phase = state_n.phase().unwrap();
        let spider_phase = spider_n.phase().unwrap();
        let new_state =
            match (state_n.is_z(), wire[0].is_h()) {
                (true, true) => ZXNode::X(state_phase),
                (true, false) => ZXNode::Z(state_phase),
                (false, true) => ZXNode::Z(state_phase),
                (false, false) => ZXNode::X(state_phase),
            };
        for nb in nnb.into_iter() {
            if nb.has_id(spider) { // self-wires are double-counted
                dg.scalar *= std::f64::consts::SQRT_2;
                continue;
            }
            let new = dg.add_node(new_state);
            if nb.is_e() {
                dg.add_wire(new, nb.id()).unwrap();
            } else {
                dg.add_wire_h(new, nb.id()).unwrap();
            }
        }
        dg.scalar *=
            (
                (1.0 + spider_phase.cis())
                + state_phase.cis() * (1.0 - spider_phase.cis())
            ) * std::f64::consts::FRAC_1_SQRT_2.powi(len_nnb as i32 + 1);
    }
}


impl RuleFinder<ZH> for StateCopy {
    type Output<'a> = StateCopyData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let zero = Phase::zero();
        let pi = Phase::pi();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() {
                for (id2, node2) in dg.neighbors_inner(id).unwrap() {
                    if dg.arity(*id2).unwrap() == 1
                        && node.is_diff_color_and(
                            node2, |_, ph2| ph2 == zero || ph2 == pi)
                    {
                        let state = *id2;
                        let spider = id;
                        return Some(StateCopyData { dg, state, spider });
                    }
                }
            }
        }
        None
    }
}

impl<'a> Rule<ZH> for StateCopyData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, state, spider } = self;
        let state_n = dg.remove_node(state).unwrap();
        let (spider_n, nnb) = dg.remove_node_nb(spider).unwrap();
        let len_nnb = nnb.len();
        let state_phase = state_n.phase().unwrap();
        let spider_phase = spider_n.phase().unwrap();
        for nb in nnb.into_iter() {
            if nb == spider {
                dg.scalar *= std::f64::consts::SQRT_2;
                continue;
            }
            let new = dg.add_node(state_n);
            dg.add_wire(new, nb).unwrap();
        }
        dg.scalar *=
            (
                (1.0 + spider_phase.cis())
                + state_phase.cis() * (1.0 - spider_phase.cis())
            ) * std::f64::consts::FRAC_1_SQRT_2.powi(len_nnb as i32 + 1);
    }
}

