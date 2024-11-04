use crate::phase::Phase;
use super::*;

/// Remove all Hadamard wires connected any spider to itself, adding π to the
/// spiders' phases for each wire removed.
///
/// This is the comprehensive version of [`HLoop`], searching more effiicently
/// for all instances of the rule.
///
/// ![h_loop][h_loop]
#[embed_doc_image::embed_doc_image("h_loop", "assets/rules/HLoop.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HLoopAll;

/// A spider and all Hadamard wires connecting it to itself.
type HLoopGroup = (NodeId, Vec<NodeId>);

/// Output of [`HLoopAll::find`].
#[derive(Debug)]
pub struct HLoopAllData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) groups: Vec<HLoopGroup>,
}

impl<'a> HLoopAllData<'a> {
    /// Return the number of spiders found with at least one Hadamard wire
    /// self-loop.
    pub fn len(&self) -> usize { self.groups.len() }

    /// Return `true` if the number of spiders found is zero.
    pub fn is_empty(&self) -> bool { self.groups.is_empty() }

    /// Return a reference to all spiders and Hadamard self-loops.
    pub fn groups(&self) -> &Vec<HLoopGroup> { &self.groups }
}

impl RuleSeal for HLoopAll { }
impl RuleFinder for HLoopAll {
    type Output<'a> = HLoopAllData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let mut groups: Vec<HLoopGroup> = Vec::new();
        let mut hh: Vec<NodeId> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() {
                for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                    if node2.is_h()
                        && node2.has_defarg()
                        && dg.arity(id2).unwrap() == 2
                        && dg.mutual_arity(id, id2).unwrap() == 2
                    {
                        hh.push(id2);
                    }
                }
                if !hh.is_empty() {
                    let g = std::mem::take(&mut hh);
                    groups.push((id, g));
                }
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(HLoopAllData { dg, groups })
        }
    }
}

impl<'a> RuleSeal for HLoopAllData<'a> { }
impl<'a> Rule for HLoopAllData<'a> {
    fn simplify(self) {
        let Self { dg, groups } = self;
        for (s, hh) in groups.into_iter() {
            let nloops = hh.len();
            hh.into_iter()
                .for_each(|h| { dg.remove_node(h).unwrap(); });
            dg.get_node_mut(s).unwrap()
                .map_phase(|ph| ph + (nloops as i64) * Phase::pi());
            dg.scalar *= std::f64::consts::FRAC_1_SQRT_2.powi(nloops as i32);
        }
    }
}
