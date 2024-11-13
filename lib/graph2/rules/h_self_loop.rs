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

impl<'a> Rule for HSelfLoopData<'a> {
    fn simplify(self) {
        let Self { dg, h } = self;
        if dg.arity(h).unwrap() == 2
            && dg.get_node(h).unwrap().has_defarg()
        {
            dg.remove_node(h).unwrap();
            dg.scalar = 0.0.into();
        } else {
            let n = dg.remove_wires(h, h, None).unwrap() as i32;
            dg.get_node_mut(h).unwrap()
                .map_arg(|a| 1.0 + (a - 1.0) / 2.0_f64.powi(n));
            dg.scalar *= 2.0_f64.powi(n);
        }
    }
}

