use crate::phase::Phase;
use super::*;

/// Propagate a unary X-spider with phase 0 through an H-box, rewriting as a
/// number of unary Z-spiders with phase 0.
///
/// ![h_explode][h_explode]
#[embed_doc_image::embed_doc_image("h_explode", "assets/rules/HExplode.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HExplode;

/// Output of [`HExplode::find`].
#[derive(Debug)]
pub struct HExplodeData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) x: NodeId,
    pub(crate) h: (NodeId, bool),
}

impl RuleSeal for HExplode { }
impl RuleFinder for HExplode {
    type Output<'a> = HExplodeData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let zero = Phase::zero();
        for (id, node) in dg.nodes_inner() {
            if node.is_x_and(|ph| ph == zero) && dg.arity(id).unwrap() == 1 {
                let mb_h =
                    dg.neighbors_of(id).unwrap()
                    .find(|(_, node2)| node2.is_h());
                if let Some((id2, node2)) = mb_h {
                    let is_had =
                        dg.arity(id2).unwrap() == 2 && node2.has_defarg();
                    let h = (id2, is_had);
                    return Some(HExplodeData { dg, x: id, h });
                }
            }
        }
        None
    }
}

impl<'a> RuleSeal for HExplodeData<'a> { }
impl<'a> Rule for HExplodeData<'a> {
    fn simplify(self) {
        let Self { dg, x, h: (h, is_had) } = self;
        dg.remove_node(x).unwrap();
        if is_had {
            let n = dg.get_node_mut(h).unwrap();
            *n = Node::z();
        } else {
            let (_, nnb) = dg.remove_node_nb(h).unwrap();
            for nb in nnb.into_iter() {
                if nb == h {
                    dg.scalar *= 2.0;
                    continue;
                }
                let z = dg.add_node(Node::z());
                dg.add_wire(z, nb).unwrap();
            }
            dg.scalar *= std::f64::consts::SQRT_2;
        }
    }
}

