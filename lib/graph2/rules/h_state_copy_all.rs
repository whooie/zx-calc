use crate::phase::Phase;
use super::*;

/// Convert all unary H-boxes or Z-spiders with phase π connected to any H-box
/// into unary X-spiders with phase π.
///
/// This is the comprehensive version of [`HStateCopy`], searching more
/// efficiently for all H-boxes connected to any unary H-box or Z-spider with
/// phase π.
///
/// ![h_state_copy][h_state_copy]
#[embed_doc_image::embed_doc_image("h_state_copy", "assets/rules/HStateCopy.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HStateCopyAll;

/// A Z- or H-state/H-box pair used by [`HStateCopyAll`].
///
/// The first item is the ID of a unary Z-spider with phase π or unary H-box
/// with argument –1, and the second is the ID of the H-box it's connected to.
pub type StateHPair = (NodeId, NodeId);

/// Output of [`HStateCopyAll::find`].
#[derive(Debug)]
pub struct HStateCopyAllData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) groups: Vec<StateHPair>, // [(z/h state, h-box)]
}

impl<'a> HStateCopyAllData<'a> {
    /// Return the number of H-boxes found with connected Z- or H-states.
    pub fn len(&self) -> usize { self.groups.len() }

    /// Return `true` if the number of groups found is zero.
    pub fn is_empty(&self) -> bool { self.groups.is_empty() }

    /// Return a reference to all groups.
    pub fn groups(&self) -> &Vec<StateHPair> { &self.groups }
}

impl RuleSeal for HStateCopyAll { }
impl RuleFinder for HStateCopyAll {
    type Output<'a> = HStateCopyAllData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        const EPSILON: f64 = 1e-12;
        let zero = Phase::zero();
        let mut groups: Vec<StateHPair> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_h() {
                for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                    if dg.arity(id2).unwrap() == 1
                        && (
                            node2.is_z_and(|ph| ph == zero)
                            || node2.is_h_and(|a| a.norm() < EPSILON)
                        )
                    {
                        groups.push((id2, id));
                    }
                }
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(HStateCopyAllData { dg, groups })
        }
    }
}

impl<'a> RuleSeal for HStateCopyAllData<'a> { }
impl<'a> Rule for HStateCopyAllData<'a> {
    fn simplify(self) {
        let Self { dg, groups } = self;
        for (state, h) in groups.into_iter() {
            dg.remove_node(state).unwrap();
            let (_, nnb) = dg.remove_node_nb(h).unwrap();
            for nb in nnb.into_iter() {
                if nb == h { continue; }
                let x = dg.add_node(Node::X(Phase::pi()));
                dg.add_wire(x, nb).unwrap();
            }
        }
    }
}

