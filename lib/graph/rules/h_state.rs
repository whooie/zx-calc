use num_complex::Complex64 as C64;
use crate::phase::Phase;
use super::*;

/// Convert a unary H-box into a unary Z-spider.
///
/// ![h_state][h_state]
#[embed_doc_image::embed_doc_image("h_state", "assets/rules/HState.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HState;

/// Output of [`HState::find`].
#[derive(Debug)]
pub struct HStateData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) h: NodeId,
    pub(crate) arg: C64,
}

impl RuleFinder for HState {
    type Output<'a> = HStateData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        const EPSILON: f64 = 1e-12;
        const I: C64 = C64 { re: 0.0, im: 1.0 };
        for (id, node) in dg.nodes_inner() {
            if node.is_h() && dg.arity(id).unwrap() == 1 {
                let arg = node.arg().unwrap();
                if (arg + 1.0).norm() < EPSILON
                    || (arg - 1.0).norm() < EPSILON
                    || (arg + I).norm() < EPSILON
                    || (arg - I).norm() < EPSILON
                {
                    return Some(HStateData { dg, h: id, arg });
                }
            }
        }
        None
    }
}

impl<'a> Rule for HStateData<'a> {
    fn simplify(self) {
        const EPSILON: f64 = 1e-12;
        const I: C64 = C64 { re: 0.0, im: 1.0 };
        let Self { dg, h, arg } = self;
        let ph =
            if (arg + 1.0).norm() < EPSILON {
                Phase::pi()
            } else if (arg - 1.0).norm() < EPSILON {
                Phase::zero()
            } else if (arg + I).norm() < EPSILON {
                -Phase::pi2()
            } else {
                Phase::pi2()
            };
        let n = dg.get_node_mut(h).unwrap();
        *n = Node::Z(ph);
    }
}

