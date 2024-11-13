use super::*;

/// Expand all H-boxes with argument 1 into unary Z-spiders with phase 0.
///
/// This is the comprehensive version of [`HMultiState`], searching more
/// efficiently for all H-boxes satisfying the conditions for this rule.
///
/// ![h_multi_state][h_multi_state]
#[embed_doc_image::embed_doc_image("h_multi_state", "assets/rules/HMultiState.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HMultiStateAll;

/// Output of [`HMultiStateAll::find`].
#[derive(Debug)]
pub struct HMultiStateAllData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) hh: Vec<NodeId>,
}

impl<'a> HMultiStateAllData<'a> {
    /// Return the number of H-boxes found with argument 1 and more than one
    /// outgoing wire.
    pub fn len(&self) -> usize { self.hh.len() }

    /// Return `true` if the number of H-boxes found is zero.
    pub fn is_empty(&self) -> bool { self.hh.is_empty() }

    /// Return a reference to all found H-boxes.
    pub fn groups(&self) -> &Vec<NodeId> { &self.hh }
}

impl RuleFinder for HMultiStateAll {
    type Output<'a> = HMultiStateAllData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        const EPSILON: f64 = 1e-12;
        let mut hh: Vec<NodeId> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            let non_self_arity =
                dg.arity(id).unwrap()
                - dg.mutual_arity(id, id).unwrap();
            if node.is_h_and(|arg| (arg + 1.0).norm() < EPSILON)
                && non_self_arity > 1
            {
                hh.push(id);
            }
        }
        if hh.is_empty() {
            None
        } else {
            Some(HMultiStateAllData { dg, hh })
        }
    }
}

impl<'a> Rule for HMultiStateAllData<'a> {
    fn simplify(self) {
        let Self { dg, hh } = self;
        for h in hh.into_iter() {
            let (_, nnb) = dg.remove_node_nb(h).unwrap();
            nnb.into_iter()
                .filter(|nb| *nb != h)
                .for_each(|nb| {
                    let z = dg.add_node(Node::z());
                    dg.add_wire(z, nb).unwrap();
                });
        }
    }
}

