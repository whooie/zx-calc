use crate::phase::Phase;
use super::*;

/// Convert between all Z- and X-spiders representing ±*y* states.
///
/// This is the comprehensive version of [`IState`], searching more efficiently
/// for all ±*y* states.
///
/// ![i_state][i_state]
#[embed_doc_image::embed_doc_image("i_state", "assets/rules/IState.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct IStateAll;

/// A unary spider representing a ±*y* state, optionally paired with an inner
/// spider.
///
/// The first two items are the unary spider's ID and phase, and the last is the
/// ID of the inner spider, if it exists. If `None`, the unary spider is
/// guaranteed to be an X-spider.
pub type IStatePair = (NodeId, Phase, Option<NodeId>);

/// Output of [`IStateAll::find`].
#[derive(Debug)]
pub struct IStateAllData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) pairs: Vec<IStatePair>,
}

impl<'a> IStateAllData<'a> {
    /// Return the number of ±*y* states found.
    pub fn len(&self) -> usize { self.pairs.len() }

    /// Return `true` if the number of ±*y* states found is zero.
    pub fn is_empty(&self) -> bool { self.pairs.is_empty() }

    /// Return a reference to all state-spider pairs.
    pub fn pairs(&self) -> &Vec<IStatePair> { &self.pairs }
}

impl RuleSeal for IStateAll { }
impl RuleFinder for IStateAll {
    type Output<'a> = IStateAllData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let pos_pi2 = Phase::pi2();
        let neg_pi2 = -Phase::pi2();
        let mut pairs: Vec<IStatePair> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider_and(|ph| ph == pos_pi2 || ph == neg_pi2)
                && dg.arity(id).unwrap() == 1
            {
                let mb_spider =
                    dg.find_map_neighbor_of(
                        id,
                        |id2, node2| node.is_diff_color(node2).then_some(id2)
                    );
                if mb_spider.is_some() || node.is_x() {
                    pairs.push((id, node.phase().unwrap(), mb_spider));
                }
            }
        }
        if pairs.is_empty() {
            None
        } else {
            Some(IStateAllData { dg, pairs })
        }
    }
}

impl<'a> RuleSeal for IStateAllData<'a> { }
impl<'a> Rule for IStateAllData<'a> {
    fn simplify(self) {
        let Self { dg, pairs } = self;
        for (state, ph0, mb_spider) in pairs.into_iter() {
            if let Some(inner) = mb_spider {
                dg.remove_node(state).unwrap();
                dg.get_node_mut(inner).unwrap()
                    .map_phase(|ph| ph - ph0);
            } else {
                let n = dg.get_node_mut(state).unwrap();
                *n = Node::Z(-ph0);
            }
            if ph0 == Phase::pi2() {
                dg.scalar *= Phase::pi4().cis();
            } else if ph0 == -Phase::pi2() {
                dg.scalar *= -Phase::pi4().cis();
            } else {
                unreachable!()
            }
        }
    }
}

