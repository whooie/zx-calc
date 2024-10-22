use super::*;

/// Any spider of either color with one or more self-loops.
#[derive(Copy, Clone, Debug)]
pub struct SpiderSelfLoop;

/// Output of [`SpiderSelfLoop::find`].
#[derive(Debug)]
pub struct SpiderSelfLoopData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s: NodeId, // spider
}

impl RuleSeal for SpiderSelfLoop { }
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

impl<'a> RuleSeal for SpiderSelfLoopData<'a> { }
impl<'a> Rule for SpiderSelfLoopData<'a> {
    fn simplify(self) {
        let SpiderSelfLoopData { dg, s } = self;
        dg.remove_wires(s, s, None).unwrap();
    }
}

