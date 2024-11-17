use super::*;

/// Expand all H-boxes with argument 1 into unary Z-spiders with phase 0.
///
/// This is the comprehensive version of [`HMultiState`], searching more
/// efficiently for all H-boxes satisfying the conditions for this rule.
///
/// ![h_multi_state][h_multi_state]
#[embed_doc_image::embed_doc_image("h_multi_state", "assets/rules/HMultiState.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct HMultiStateAll;

/// Output of [`HMultiStateAll::find`].
#[derive(Debug)]
pub struct HMultiStateAllData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) hh: Vec<NodeId>,
}

impl<'a, A> HMultiStateAllData<'a, A>
where A: DiagramData
{
    /// Return the number of H-boxes found with argument 1 and more than one
    /// outgoing wire.
    pub fn len(&self) -> usize { self.hh.len() }

    /// Return `true` if the number of H-boxes found is zero.
    pub fn is_empty(&self) -> bool { self.hh.is_empty() }

    /// Return a reference to all found H-boxes.
    pub fn groups(&self) -> &Vec<NodeId> { &self.hh }
}

impl RuleFinder<ZH> for HMultiStateAll {
    type Output<'a> = HMultiStateAllData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let mut hh: Vec<NodeId> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.has_arg(1.0) && dg.arity(id).unwrap() > 0 {
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

impl<'a> Rule<ZH> for HMultiStateAllData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, hh } = self;
        let mut self_loops_x2: i32 = 0;
        for h in hh.into_iter() {
            let (_, nnb) = dg.remove_node_nb(h).unwrap();
            nnb.into_iter()
                .for_each(|nb| {
                    if nb == h {
                        self_loops_x2 += 1;
                    } else {
                        let z = dg.add_node(ZHNode::z());
                        dg.add_wire(z, nb).unwrap();
                    }
                });
        }
        dg.scalar *= 2.0_f64.powi(self_loops_x2 / 2);
    }
}

