use super::*;

/// Remove all wires connecting an H-box to itself, adjusting its argument for
/// each wire removed.
///
/// The argument adjustment per wire takes *a* → (1 + *a*) / 2.
///
/// ![h_self_loop][h_self_loop]
#[embed_doc_image::embed_doc_image("h_self_loop", "assets/rules/HSelfLoop.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HSelfLoop;

/// Output of [`HSelfLoop::find`].
#[derive(Debug)]
pub struct HSelfLoopData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) h: NodeId,
}

impl RuleSeal for HSelfLoop { }
impl RuleFinder for HSelfLoop {
    type Output<'a> = HSelfLoopData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if node.is_h() && dg.mutual_arity(id, id).unwrap() > 0 {
                return Some(HSelfLoopData { dg, h: id });
            }
        }
        None
    }
}

impl<'a> RuleSeal for HSelfLoopData<'a> { }
impl<'a> Rule for HSelfLoopData<'a> {
    fn simplify(self) {
        let Self { dg, h } = self;
        let n = dg.remove_wires(h, h, None).unwrap() as i32;
        dg.get_node_mut(h).unwrap()
            .map_arg(|a| 1.0 + (a - 1.0) / 2.0_f64.powi(n));
    }
}

