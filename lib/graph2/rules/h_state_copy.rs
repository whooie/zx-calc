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
    pub(crate) h: (NodeId, bool), // h-box copier
}

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
                            || node2.is_h_and(|a| (a + 1.0).norm() < EPSILON)
                        )
                    {
                        let is_had =
                            dg.arity(id).unwrap() == 2 && node.has_defarg();
                        let h = (id, is_had);
                        return Some(HStateCopyData { dg, state: id2, h });
                    }
                }
            }
        }
        None
    }
}

impl<'a> Rule for HStateCopyData<'a> {
    fn simplify(self) {
        let Self { dg, state, h: (h, is_had) } = self;
        dg.remove_node(state).unwrap();
        if is_had {
            let n = dg.get_node_mut(h).unwrap();
            *n = Node::x_pi();
        } else {
            let (Node::H(a), nnb) = dg.remove_node_nb(h).unwrap()
                else { unreachable!() };
            let nnb_len = nnb.len();
            for nb in nnb.into_iter() {
                if nb == h {
                    dg.scalar *= 2.0;
                    continue;
                }
                let x = dg.add_node(Node::x_pi());
                dg.add_wire(x, nb).unwrap();
            }
            dg.scalar *=
                (1.0 - a) / std::f64::consts::SQRT_2.powi(nnb_len as i32);
        }
    }
}

