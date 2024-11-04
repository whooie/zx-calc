use crate::phase::Phase;
use super::*;

/// Expand all unary X-spiders with phase 0 connected to all H-boxes into unary
/// Z-spiders with phase 0.
///
/// This rule is the comprehensive version of [`HExplode`], searching more
/// efficiently for all spiders and H-boxes satisfying the conditions for this
/// rule.
///
/// ![h_explode][h_explode]
#[embed_doc_image::embed_doc_image("h_explode", "assets/rules/HExplode.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HExplodeAll;

/// An X-state/H-box pair used by [`HExplodeAll`].
///
/// The first item is the ID of a unary X-spider with phase 0; the second is the
/// ID of the H-box it's connected to; and the last is `true` when the H-box has
/// arity 2 and argument –1.
pub type XHPair = (NodeId, NodeId, bool);

/// Output of [`HExplodeAll::find`].
#[derive(Debug)]
pub struct HExplodeAllData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) groups: Vec<XHPair>, // [(x spider, h-box)]
}

impl<'a> HExplodeAllData<'a> {
    /// Return the number of H-boxes found with connected X-spider states.
    pub fn len(&self) -> usize { self.groups.len() }

    /// Return `true` if the number of groups found is zero.
    pub fn is_empty(&self) -> bool { self.groups.is_empty() }

    /// Return a reference to all groups.
    pub fn groups(&self) -> &Vec<XHPair> { &self.groups }
}

impl RuleSeal for HExplodeAll { }
impl RuleFinder for HExplodeAll {
    type Output<'a> = HExplodeAllData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let mut groups: Vec<XHPair> = Vec::new();
        let zero = Phase::zero();
        for (id, node) in dg.nodes_inner() {
            if node.is_h() {
                let mb_x =
                    dg.neighbors_of(id).unwrap()
                    .find_map(|(id2, node2)| {
                        (
                            node2.is_x_and(|ph| ph == zero)
                            && dg.arity(id2).unwrap() == 1
                        ).then_some(id2)
                    });
                if let Some(x) = mb_x {
                    let is_had =
                        dg.arity(id).unwrap() == 2 && node.has_defarg();
                    groups.push((x, id, is_had));
                }
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(HExplodeAllData { dg, groups })
        }
    }
}

impl<'a> RuleSeal for HExplodeAllData<'a> { }
impl<'a> Rule for HExplodeAllData<'a> {
    fn simplify(self) {
        let Self { dg, groups } = self;
        for (x, h, is_had) in groups.into_iter() {
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
}

