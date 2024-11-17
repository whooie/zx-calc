use num_complex::Complex64 as C64;
use crate::{ c64_eq, phase::Phase };
use super::*;

/// Convert all unary H-boxes into unary Z-spiders.
///
/// This is the comprehensive version of [`HState`], searching more efficiently
/// for all unary H-boxes.
///
/// ![h_state][h_state]
#[embed_doc_image::embed_doc_image("h_state", "assets/rules/HState.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct HStateAll;

/// Output of [`HStateAll::find`].
#[derive(Debug)]
pub struct HStateAllData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) hh: Vec<(NodeId, C64)>,
}

impl<'a, A> HStateAllData<'a, A>
where A: DiagramData
{
    /// Return the number of unary H-boxes found with argument ±1.
    pub fn len(&self) -> usize { self.hh.len() }

    /// Return `true` if the number of H-boxes found is zero.
    pub fn is_empty(&self) -> bool { self.hh.is_empty() }

    /// Return a reference to all found H-boxes with their arguments.
    pub fn groups(&self) -> &Vec<(NodeId, C64)> { &self.hh }
}

impl RuleFinder<ZH> for HStateAll {
    type Output<'a> = HStateAllData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        const I: C64 = C64 { re: 0.0, im: 1.0 };
        let mut hh: Vec<(NodeId, C64)> = Vec::new();
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
                hh.push((id, arg));
            }
        }
        if hh.is_empty() {
            None
        } else {
            Some(HStateAllData { dg, hh })
        }
    }
}

impl<'a> Rule<ZH> for HStateAllData<'a, ZH> {
    fn simplify(self) {
        const I: C64 = C64 { re: 0.0, im: 1.0 };
        let Self { dg, hh } = self;
        for (h, arg) in hh.into_iter() {
            let ph =
                if c64_eq(arg, -1.0) {
                    Phase::pi()
                } else if c64_eq(arg, 1.0) {
                    Phase::zero()
                } else if c64_eq(arg, -I) {
                    -Phase::pi2()
                } else {
                    Phase::pi()
                };
            let n = dg.get_node_mut(h).unwrap();
            *n = ZHNode::Z(ph);
        }
    }
}

