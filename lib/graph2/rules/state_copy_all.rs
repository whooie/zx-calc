use crate::phase::Phase;
use super::*;

/// Copy all unary spiders representing a pure Z- or X-basis state through any
/// oppositely colored spiders they're connected to.
///
/// This is the comprehensive version of [`StateCopy`], searching more
/// efficiently for all such pairs of spiders.
///
/// ![state_copy][state_copy]
#[embed_doc_image::embed_doc_image("state_copy", "assets/rules/StateCopy.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct StateCopyAll;

/// A unary spider with phase π or 0 connected to an oppositely colored spider.
///
/// The first item is the ID of the unary spider, and the second is the ID of
/// the other spider it's connected to.
pub type StateCopyPair = (NodeId, NodeId);

/// Output of [`StateCopyAll::find`].
#[derive(Debug)]
pub struct StateCopyAllData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) pairs: Vec<StateCopyPair>,
}

impl<'a> StateCopyAllData<'a> {
    /// Return the number of spider copy pairs found.
    pub fn len(&self) -> usize { self.pairs.len() }

    /// Return `true` if the number of spider copy pairs found is zero.
    pub fn is_empty(&self) -> bool { self.pairs.is_empty() }

    /// Return a reference to all found spider copy pairs.
    pub fn pairs(&self) -> &Vec<StateCopyPair> { &self.pairs }
}

impl RuleSeal for StateCopyAll { }
impl RuleFinder for StateCopyAll {
    type Output<'a> = StateCopyAllData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let zero = Phase::zero();
        let pi = Phase::pi();
        let mut pairs: Vec<StateCopyPair> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() {
                for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                    if dg.arity(id2).unwrap() == 1
                        && node.is_diff_color_and(
                            node2, |_, ph2| ph2 == zero || ph2 == pi
                        )
                    {
                        pairs.push((id2, id));
                    }
                }
            }
        }
        if pairs.is_empty() {
            None
        } else {
            Some(StateCopyAllData { dg, pairs })
        }
    }
}

impl<'a> RuleSeal for StateCopyAllData<'a> { }
impl<'a> Rule for StateCopyAllData<'a> {
    fn simplify(self) {
        let Self { dg, pairs } = self;
        for (state, spider) in pairs.into_iter() {
            let n = dg.remove_node(state).unwrap();
            let (_, nnb) = dg.remove_node_nb(spider).unwrap();
            for nb in nnb.into_iter() {
                if nb == spider { continue; }
                let new = dg.add_node(n);
                dg.add_wire(new, nb).unwrap();
            }
        }
    }
}

