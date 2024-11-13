use super::*;

/// Expand an H-box with argument 1 into unary Z-spiders with phase 0.
///
/// ![h_multi_state][h_multi_state]
#[embed_doc_image::embed_doc_image("h_multi_state", "assets/rules/HMultiState.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HMultiState;

/// Output of [`HMultiState::find`].
#[derive(Debug)]
pub struct HMultiStateData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) h: NodeId,
}

// Note that self-wires are not important here, since we demand arg == 1 which
// is a fixed point of the usual self-wire argument transformation
//
//      a -> 1 + (a - 1) / 2^n
//
// (See `HSelfLoop`.)

impl RuleFinder for HMultiState {
    type Output<'a> = HMultiStateData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        const EPSILON: f64 = 1e-12;
        for (id, node) in dg.nodes_inner() {
            let non_self_arity =
                dg.arity(id).unwrap()
                - dg.mutual_arity(id, id).unwrap();
            if node.is_h_and(|arg| (arg + 1.0).norm() < EPSILON)
                && non_self_arity > 1
            {
                return Some(HMultiStateData { dg, h: id });
            }
        }
        None
    }
}

impl<'a> Rule for HMultiStateData<'a> {
    fn simplify(self) {
        let Self { dg, h } = self;
        let (_, nnb) = dg.remove_node_nb(h).unwrap();
        nnb.into_iter()
            .filter(|nb| *nb != h)
            .for_each(|nb| {
                let z = dg.add_node(Node::z());
                dg.add_wire(z, nb).unwrap();
            });
    }
}

