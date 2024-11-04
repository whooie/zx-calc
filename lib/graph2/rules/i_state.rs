use crate::phase::Phase;
use super::*;

/// Convert between Z- and X-spider representations of ±*y* states.
///
/// If a unary spider representing a ±*y* state is connected to another spider
/// of opposite color, convert it to the opposite-color representation and fuse
/// the spiders. Otherwise, convert an X-spider representation to a Z-spider
/// representation.
///
/// ![i_state][i_state]
#[embed_doc_image::embed_doc_image("i_state", "assets/rules/IState.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct IState;

/// Output of [`IState::find`].
#[derive(Debug)]
pub struct IStateData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) state: (NodeId, Phase),
    pub(crate) mb_spider: Option<NodeId>, // oppositely colored inner spider
}

impl RuleSeal for IState { }
impl RuleFinder for IState {
    type Output<'a> = IStateData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let pos_pi2 = Phase::pi2();
        let neg_pi2 = -Phase::pi2();
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
                    let state = (id, node.phase().unwrap());
                    return Some(IStateData { dg, state, mb_spider });
                }
            }
        }
        None
    }
}

impl<'a> RuleSeal for IStateData<'a> { }
impl<'a> Rule for IStateData<'a> {
    fn simplify(self) {
        let Self { dg, state, mb_spider } = self;
        if let Some(inner) = mb_spider {
            dg.remove_node(state.0).unwrap();
            dg.get_node_mut(inner).unwrap()
                .map_phase(|ph| ph - state.1);
        } else {
            let n = dg.get_node_mut(state.0).unwrap();
            *n = Node::Z(-state.1);
        }
        if state.1 == Phase::pi2() {
            dg.scalar *= Phase::pi4().cis();
        } else if state.1 == -Phase::pi2() {
            dg.scalar *= -Phase::pi4().cis();
        } else {
            unreachable!()
        }
    }
}

