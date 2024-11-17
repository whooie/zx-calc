use super::*;

/// Remove all empty wires connecting a spider to itself.
///
/// ![spider_self_loop][spider_self_loop]
#[embed_doc_image::embed_doc_image("spider_self_loop", "assets/rules/SpiderSelfLoop.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct SpiderSelfLoop;

/// Output of [`SpiderSelfLoop::find`].
#[derive(Debug)]
pub struct SpiderSelfLoopData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) s: NodeId, // spider
}

impl RuleFinder<ZX> for SpiderSelfLoop {
    type Output<'a> = SpiderSelfLoopData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        for (id, _node) in dg.nodes_inner() {
            if dg.mutual_arity_e(id, id).unwrap() > 0 {
                return Some(SpiderSelfLoopData { dg, s: id });
            }
        }
        None
    }
}

impl<'a> Rule<ZX> for SpiderSelfLoopData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, s } = self;
        dg.remove_wires_e(s, s, None).unwrap();
    }
}

impl RuleFinder<ZH> for SpiderSelfLoop {
    type Output<'a> = SpiderSelfLoopData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() && dg.mutual_arity(id, id).unwrap() > 0 {
                return Some(SpiderSelfLoopData { dg, s: id });
            }
        }
        None
    }
}

impl<'a> Rule<ZH> for SpiderSelfLoopData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, s } = self;
        dg.remove_wires(s, s, None).unwrap();
    }
}

