use num_complex::Complex64 as C64;
use crate::phase::Phase;
use super::*;

/// Merge all groups of unary H-boxes connected to Z-spiders, replacing each
/// with a single unary H-box of appropriate argument.
///
/// This is the comprehensive version of [`HStateMul`], searching more
/// efficiently for all Z-spiders connected to at least two unary H-boxes.
///
/// ![h_state_mul][h_state_mul]
#[embed_doc_image::embed_doc_image("h_state_mul", "assets/rules/HStateMul.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HStateMulAll;

/// A Z-spider and a group of attached unary H-boxes.
///
/// The first item is the ID of the spider, and the second is a list of IDs for
/// the H-boxes and their arguments. Each list is at least two elements long.
pub type HStateMulGroup = (NodeId, Vec<(NodeId, C64)>);

/// Output of [`HStateMulAll::find`].
#[derive(Debug)]
pub struct HStateMulAllData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) groups: Vec<HStateMulGroup>,
}

impl<'a, A> HStateMulAllData<'a, A>
where A: DiagramData
{
    /// Return the number of Z-spiders found with at least two unary H-boxes
    /// attached.
    pub fn len(&self) -> usize { self.groups.len() }

    /// Return `true` if the number of groups is zero.
    pub fn is_empty(&self) -> bool { self.groups.is_empty() }

    /// Return a reference to all groups.
    pub fn groups(&self) -> &Vec<HStateMulGroup> { &self.groups }
}

impl RuleFinder<ZH> for HStateMulAll {
    type Output<'a> = HStateMulAllData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let mut groups: Vec<HStateMulGroup> = Vec::new();
        let mut hh: Vec<(NodeId, C64)> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_z() {
                dg.neighbors_inner(id).unwrap()
                    .filter(|(id2, node2)| {
                        node2.is_h() && dg.arity(**id2).unwrap() == 1
                    })
                    .for_each(|(id2, node2)| {
                        hh.push((*id2, node2.arg().unwrap()));
                    });
            }
            if hh.len() > 1 {
                let g = std::mem::take(&mut hh);
                groups.push((id, g));
            } else if !hh.is_empty() {
                hh.clear();
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(HStateMulAllData { dg, groups })
        }
    }
}

impl<'a> Rule<ZH> for HStateMulAllData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, groups } = self;
        let zero = Phase::zero();
        for (s, mut hh) in groups.into_iter() {
            if dg.arity(s).unwrap() == hh.len() + 1
                && dg.get_node(s).unwrap().phase().unwrap() == zero
            {
                let a_prod: C64 =
                    hh.into_iter()
                    .map(|(id, a_k)| {
                        dg.remove_node(id).unwrap();
                        a_k
                    })
                    .product();
                let n = dg.get_node_mut(s).unwrap();
                *n = ZHNode::H(a_prod);
            } else {
                let (h0, _) = hh.pop().unwrap();
                let a_prod: C64 =
                    hh.into_iter()
                    .map(|(id, a_k)| {
                        dg.remove_node(id).unwrap();
                        a_k
                    })
                    .product();
                dg.get_node_mut(h0).unwrap()
                    .map_arg(|a0| a0 * a_prod);
            }
        }
    }
}

