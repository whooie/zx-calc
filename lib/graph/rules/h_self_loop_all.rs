use super::*;

/// Remove all wires connecting all H-boxes to themselves, adjusting their
/// arguments accordingly.
///
/// This is the comprehensive version of [`HSelfLoop`], searching more
/// efficiently for all H-boxes with at least one wire connecting it to itself.
///
/// ![h_self_loop][h_self_loop]
#[embed_doc_image::embed_doc_image("h_self_loop", "assets/rules/HSelfLoop.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HSelfLoopAll;

/// Output of [`HSelfLoopAll::find`].
#[derive(Debug)]
pub struct HSelfLoopAllData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) hh: Vec<NodeId>,
}

impl<'a> HSelfLoopAllData<'a> {
    /// Return the number of H-boxes with self-loops found.
    pub fn len(&self) -> usize { self.hh.len() }

    /// Return `true` if the number of H-boxes is zero.
    pub fn is_empty(&self) -> bool { self.hh.is_empty() }

    /// Return a reference to all H-boxes.
    pub fn groups(&self) -> &Vec<NodeId> { &self.hh }
}

impl RuleFinder for HSelfLoopAll {
    type Output<'a> = HSelfLoopAllData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let mut hh: Vec<NodeId> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_h() && dg.mutual_arity(id, id).unwrap() > 0 {
                hh.push(id);
            }
        }
        if hh.is_empty() {
            None
        } else {
            Some(HSelfLoopAllData { dg, hh })
        }
    }
}

impl<'a> Rule for HSelfLoopAllData<'a> {
    fn simplify(self) {
        let Self { dg, hh } = self;
        for h in hh.into_iter() {
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
}

