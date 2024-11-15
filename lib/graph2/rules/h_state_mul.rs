use num_complex::Complex64 as C64;
use crate::phase::Phase;
use super::*;

/// Merge at unary H-boxes connected to the same Z-spider into one, multiplying
/// their arguments.
///
/// ![h_state_mul][h_state_mul]
#[embed_doc_image::embed_doc_image("h_state_mul", "assets/rules/HStateMul.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HStateMul;

/// Output of [`HStateMul::find`].
#[derive(Debug)]
pub struct HStateMulData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) s: NodeId, // spider
    pub(crate) hh: Vec<(NodeId, C64)>, // h-boxes
}

impl RuleFinder<ZH> for HStateMul {
    type Output<'a> = HStateMulData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let mut hh: Vec<(NodeId, C64)> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_z() {
                dg.neighbors_inner(id).unwrap()
                    .filter(|(id2, node2)| {
                        node2.is_h() && dg.arity(**id2).unwrap() == 1
                    })
                    .for_each(|(id2, node2)| {
                        hh.push((*id2, node2.arg().unwrap()));
                    });
            }
            if hh.len() > 1 {
                return Some(HStateMulData { dg, s: id, hh });
            } else if !hh.is_empty() {
                hh.clear();
            }
        }
        None
    }
}

impl<'a> Rule<ZH> for HStateMulData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, s, mut hh } = self;
        if dg.arity(s).unwrap() == hh.len() + 1
            && dg.get_node(s).unwrap().phase().unwrap() == Phase::zero()
        {
            let a_prod: C64 =
                hh.into_iter()
                .map(|(id, a_k)| {
                    dg.remove_node(id).unwrap();
                    a_k
                })
                .product();
            let n = dg.get_node_mut(s).unwrap();
            *n = ZHNode::H(a_prod);
        } else {
            let (h0, _) = hh.pop().unwrap();
            let a_prod: C64 =
                hh.into_iter()
                .map(|(id, a_k)| {
                    dg.remove_node(id).unwrap();
                    a_k
                })
                .product();
            dg.get_node_mut(h0).unwrap()
                .map_arg(|a0| a0 * a_prod);
        }
    }
}

