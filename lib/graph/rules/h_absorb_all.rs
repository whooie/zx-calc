use crate::phase::Phase;
use super::*;

/// Remove all unary X-spider with phase π connected to any H-box.
///
/// This rule is a comprehensive version of [`HAbsorb`], searching more
/// efficiently for all H-boxes and spiders that satisfy the conditions for the
/// rule.
///
/// ![h_absorb][h_absorb]
#[embed_doc_image::embed_doc_image("h_absorb", "assets/rules/HAbsorb.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HAbsorbAll;

/// Output of [`HAbsorbAll::find`].
#[derive(Debug)]
pub struct HAbsorbAllData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) groups: Vec<Vec<NodeId>>, // all pi x-spiders, grouped by h-box
}

impl<'a> HAbsorbAllData<'a> {
    /// Return the total number of unary X(π)-spiders to remove.
    pub fn len_all(&self) -> usize { self.groups.iter().map(|g| g.len()).sum() }

    /// Return the number of H-boxes with attached unary X(π)-spiders.
    pub fn len(&self) -> usize { self.groups.len() }

    /// Return `true` if the number of unary X(π)-spiders to remove is zero.
    pub fn is_empty(&self) -> bool { self.groups.is_empty() }

    /// Return a reference to all unary X(π)-spiders, grouped by H-box.
    pub fn groups(&self) -> &Vec<Vec<NodeId>> { &self.groups }
}

impl RuleFinder for HAbsorbAll {
    type Output<'a> = HAbsorbAllData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let mut groups: Vec<Vec<NodeId>> = Vec::new();
        let mut group: Vec<NodeId> = Vec::new();
        let pi = Phase::pi();
        // no possibility of group overlap here
        for (id, node) in dg.nodes_inner() {
            if node.is_h() {
                dg.neighbors_of_inner(id).unwrap()
                    .filter_map(|(id2, node2)| {
                        (
                            node2.is_x_and(|ph| ph == pi)
                            && dg.arity(id2).unwrap() == 1
                        ).then_some(id2)
                    })
                    .for_each(|id2| { group.push(id2); });
                if !group.is_empty() {
                    let g = std::mem::take(&mut group);
                    groups.push(g);
                }
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(HAbsorbAllData { dg, groups })
        }
    }
}

impl<'a> Rule for HAbsorbAllData<'a> {
    fn simplify(self) {
        let Self { dg, groups } = self;
        for states in groups.into_iter() {
            let nstates = states.len();
            states.into_iter()
                .for_each(|id| { dg.remove_node(id).unwrap(); });
            dg.scalar *= std::f64::consts::SQRT_2.powi(nstates as i32);
        }
    }
}

