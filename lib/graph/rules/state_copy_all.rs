use crate::phase::Phase;
use super::*;

/// Copy all unary spiders representing a pure Z- or X-basis state through any
/// oppositely colored spiders they're connected to.
///
/// This is the comprehensive version of [`StateCopy`], searching more
/// efficiently for all such pairs of spiders.
///
/// ![state_copy][state_copy]
#[embed_doc_image::embed_doc_image("state_copy", "assets/rules/StateCopy.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct StateCopyAll;

/// A unary spider with phase π or 0 connected to an oppositely colored spider.
///
/// The first item is the ID of the unary spider, and the second is the ID of
/// the other spider it's connected to.
pub type StateCopyPair = (NodeId, NodeId);

/// Output of [`StateCopyAll::find`].
#[derive(Debug)]
pub struct StateCopyAllData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) pairs: Vec<StateCopyPair>,
}

impl<'a, A> StateCopyAllData<'a, A>
where A: DiagramData
{
    /// Return the number of spider copy pairs found.
    pub fn len(&self) -> usize { self.pairs.len() }

    /// Return `true` if the number of spider copy pairs found is zero.
    pub fn is_empty(&self) -> bool { self.pairs.is_empty() }

    /// Return a reference to all found spider copy pairs.
    pub fn pairs(&self) -> &Vec<StateCopyPair> { &self.pairs }
}

impl RuleFinder<ZX> for StateCopyAll {
    type Output<'a> = StateCopyAllData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        let mut pairs: Vec<StateCopyPair> = Vec::new();
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
                    pairs.push((wire.id(), id));
                }
            }
        }
        if pairs.is_empty() {
            None
        } else {
            Some(StateCopyAllData { dg, pairs })
        }
    }
}

impl<'a> Rule<ZX> for StateCopyAllData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, pairs } = self;
        for (state, spider) in pairs.into_iter() {
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
}

impl RuleFinder<ZH> for StateCopyAll {
    type Output<'a> = StateCopyAllData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let zero = Phase::zero();
        let pi = Phase::pi();
        let mut pairs: Vec<StateCopyPair> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() {
                for (id2, node2) in dg.neighbors_inner(id).unwrap() {
                    if dg.arity(*id2).unwrap() == 1
                        && node.is_diff_color_and(
                            node2, |_, ph2| ph2 == zero || ph2 == pi)
                    {
                        pairs.push((*id2, id));
                    }
                }
            }
        }
        if pairs.is_empty() {
            None
        } else {
            Some(StateCopyAllData { dg, pairs })
        }
    }
}

impl<'a> Rule<ZH> for StateCopyAllData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, pairs } = self;
        for (state, spider) in pairs.into_iter() {
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
}

