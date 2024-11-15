use super::*;

/// Remove all empty wires connecting any spider to itself.
///
/// This is the comprehensive version of [`SpiderSelfLoop`], searching more
/// efficiently for all spiders with at least one self-loop.
///
/// ![spider_self_loop][spider_self_loop]
#[embed_doc_image::embed_doc_image("spider_self_loop", "assets/rules/SpiderSelfLoop.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct SpiderSelfLoopAll;

/// Output of [`SpiderSelfLoopAll::find`].
#[derive(Debug)]
pub struct SpiderSelfLoopAllData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) spiders: Vec<NodeId>,
}

impl<'a, A> SpiderSelfLoopAllData<'a, A>
where A: DiagramData
{
    /// Return the number of spiders found with self-loops.
    pub fn len(&self) -> usize { self.spiders.len() }

    /// Return `true` if the number of spiders found is zero.
    pub fn is_empty(&self) -> bool { self.spiders.is_empty() }

    /// Return a reference to all spiders found with self-loops.
    pub fn groups(&self) -> &Vec<NodeId> { &self.spiders }
}

impl RuleFinder<ZX> for SpiderSelfLoopAll {
    type Output<'a> = SpiderSelfLoopAllData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        let mut spiders: Vec<NodeId> = Vec::new();
        for (id, _node) in dg.nodes_inner() {
            if dg.mutual_arity_e(id, id).unwrap() > 0 {
                spiders.push(id);
            }
        }
        if spiders.is_empty() {
            None
        } else {
            Some(SpiderSelfLoopAllData { dg, spiders })
        }
    }
}

impl<'a> Rule<ZX> for SpiderSelfLoopAllData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, spiders } = self;
        spiders.into_iter()
            .for_each(|s| { dg.remove_wires_e(s, s, None).unwrap(); });
    }
}

impl RuleFinder<ZH> for SpiderSelfLoopAll {
    type Output<'a> = SpiderSelfLoopAllData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let mut spiders: Vec<NodeId> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() && dg.mutual_arity(id, id).unwrap() > 0 {
                spiders.push(id);
            }
        }
        if spiders.is_empty() {
            None
        } else {
            Some(SpiderSelfLoopAllData { dg, spiders })
        }
    }
}

impl<'a> Rule<ZH> for SpiderSelfLoopAllData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, spiders } = self;
        spiders.into_iter()
            .for_each(|s| { dg.remove_wires(s, s, None).unwrap(); });
    }
}

