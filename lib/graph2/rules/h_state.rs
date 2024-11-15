use num_complex::Complex64 as C64;
use crate::{ c64_eq, phase::Phase };
use super::*;

/// Convert a unary H-box into a unary Z-spider.
///
/// ![h_state][h_state]
#[embed_doc_image::embed_doc_image("h_state", "assets/rules/HState.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HState;

/// Output of [`HState::find`].
#[derive(Debug)]
pub struct HStateData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) h: NodeId,
    pub(crate) arg: C64,
}

impl RuleFinder<ZH> for HState {
    type Output<'a> = HStateData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        const I: C64 = C64 { re: 0.0, im: 1.0 };
        for (id, node) in dg.nodes_inner() {
            if node.is_h() && dg.arity(id).unwrap() == 1
                && (
                    node.has_arg(1.0)
                    || node.has_arg(-1.0)
                    || node.has_arg(I)
                    || node.has_arg(-I)
                )
            {
                let arg = node.arg().unwrap();
                return Some(HStateData { dg, h: id, arg });
            }
        }
        None
    }
}

impl<'a> Rule<ZH> for HStateData<'a, ZH> {
    fn simplify(self) {
        const I: C64 = C64 { re: 0.0, im: 1.0 };
        let Self { dg, h, arg } = self;
        let ph =
            if c64_eq(arg, -1.0) {
                Phase::pi()
            } else if c64_eq(arg, 1.0) {
                Phase::zero()
            } else if c64_eq(arg, -I) {
                -Phase::pi2()
            } else {
                Phase::pi2()
            };
        let n = dg.get_node_mut(h).unwrap();
        *n = ZHNode::Z(ph);
    }
}

