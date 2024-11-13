use super::*;

/// Remove all empty wires connecting a spider to itself.
///
/// ![spider_self_loop][spider_self_loop]
#[embed_doc_image::embed_doc_image("spider_self_loop", "assets/rules/SpiderSelfLoop.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct SpiderSelfLoop;

/// Output of [`SpiderSelfLoop::find`].
#[derive(Debug)]
pub struct SpiderSelfLoopData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s: NodeId, // spider
}

impl RuleFinder for SpiderSelfLoop {
    type Output<'a> = SpiderSelfLoopData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if (node.is_z() || node.is_x())
                && dg.mutual_arity(id, id).unwrap() > 0
            {
                return Some(SpiderSelfLoopData { dg, s: id });
            }
        }
        None
    }
}

impl<'a> Rule for SpiderSelfLoopData<'a> {
    fn simplify(self) {
        let SpiderSelfLoopData { dg, s } = self;
        dg.remove_wires(s, s, None).unwrap();
    }
}

