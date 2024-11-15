use super::*;

/// Expand an H-box with argument 1 into unary Z-spiders with phase 0.
///
/// ![h_multi_state][h_multi_state]
#[embed_doc_image::embed_doc_image("h_multi_state", "assets/rules/HMultiState.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HMultiState;

/// Output of [`HMultiState::find`].
#[derive(Debug)]
pub struct HMultiStateData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) h: NodeId,
}

impl RuleFinder<ZH> for HMultiState {
    type Output<'a> = HMultiStateData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if node.has_arg(1.0) && dg.arity(id).unwrap() > 0 {
                return Some(HMultiStateData { dg, h: id });
            }
        }
        None
    }
}

impl<'a> Rule<ZH> for HMultiStateData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, h } = self;
        let (_, nnb) = dg.remove_node_nb(h).unwrap();
        let mut self_loops_x2: i32 = 0;
        nnb.into_iter()
            .for_each(|nb| {
                if nb == h {
                    self_loops_x2 += 1;
                } else {
                    let z = dg.add_node(ZHNode::z());
                    dg.add_wire(z, nb).unwrap();
                }
            });
        dg.scalar *= 2.0_f64.powi(self_loops_x2 / 2);
    }
}

