use crate::phase::Phase;
use super::*;

/// A specific case of [`HEuler`], removing Hadamard wires connecting a spider
/// to itself and adding π to the spider's phase for each wire removed.
///
/// ![h_loop][h_loop]
#[embed_doc_image::embed_doc_image("h_loop", "assets/rules/HLoop.svg")]
#[derive(Copy, Clone, Debug)]
pub struct HLoop;

/// Output of [`HLoop::find`].
#[derive(Debug)]
pub struct HLoopData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s: NodeId, // spider
    pub(crate) h: NodeId, // h-box
}

impl RuleSeal for HLoop { }
impl RuleFinder for HLoop {
    type Output<'a> = HLoopData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if node.is_z() || node.is_x() {
                for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                    if node2.is_h()
                        && node2.has_defarg()
                        && dg.arity(id2).unwrap() == 2
                        && dg.mutual_arity(id, id2).unwrap() == 2
                    {
                        return Some(HLoopData { dg, s: id, h: id2 });
                    }
                }
            }
        }
        None
    }
}

impl<'a> RuleSeal for HLoopData<'a> { }
impl<'a> Rule for HLoopData<'a> {
    fn simplify(self) {
        let Self { dg, s, h } = self;
        dg.nodes[s].as_mut().unwrap()
            .map_phase(|ph| ph + Phase::pi());
        dg.remove_node(h).unwrap();
    }
}

