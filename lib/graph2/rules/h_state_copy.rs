use crate::phase::Phase;
use super::*;

/// Z-spider and H-box version of [`HExplode`], converting a unary H-box or
/// Z-spider with phase π connected to a H-box into unary X-spiders with phase
/// π.
///
/// This is the H/Z-state counterpart to [`HExplode`].
///
/// ![h_state_copy][h_state_copy]
#[embed_doc_image::embed_doc_image("h_state_copy", "assets/rules/HStateCopy.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HStateCopy;

/// Output of [`HStateCopy::find`].
#[derive(Debug)]
pub struct HStateCopyData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) state: NodeId, // z-state or h-state
    pub(crate) h: NodeId, // h-box copier
}

impl RuleSeal for HStateCopy { }
impl RuleFinder for HStateCopy {
    type Output<'a> = HStateCopyData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        const EPSILON: f64 = 1e-12;
        let zero = Phase::zero();
        for (id, node) in dg.nodes_inner() {
            if node.is_h() {
                for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                    if dg.arity(id2).unwrap() == 1
                        && (
                            node2.is_z_and(|ph| ph == zero)
                            || node2.is_h_and(|a| a.norm() < EPSILON)
                        )
                    {
                        return Some(HStateCopyData { dg, state: id2, h: id });
                    }
                }
            }
        }
        None
    }
}

impl<'a> RuleSeal for HStateCopyData<'a> { }
impl<'a> Rule for HStateCopyData<'a> {
    fn simplify(self) {
        let Self { dg, state, h } = self;
        dg.remove_node(state).unwrap();
        let (_, nnb) = dg.remove_node_nb(h).unwrap();
        for nb in nnb.into_iter() {
            if nb == h { continue; }
            let x = dg.add_node(Node::X(Phase::pi()));
            dg.add_wire(x, nb).unwrap();
        }
    }
}

