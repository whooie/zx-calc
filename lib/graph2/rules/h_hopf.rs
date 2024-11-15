use crate::phase::Phase;
use super::*;

/// H-box version of [`Hopf`], removing all but one wire between a multiply
/// connected Z-spider and H-box.
///
/// ![h_hopf][h_hopf]
#[embed_doc_image::embed_doc_image("h_hopf", "assets/rules/HHopf.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HHopf;

/// Output of [`HHopf::find`].
#[derive(Debug)]
pub struct HHopfData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) z: NodeId, // z-spider
    pub(crate) h: NodeId, // h-box
}

impl RuleFinder<ZH> for HHopf {
    type Output<'a> = HHopfData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if node.is_z() {
                for (id2, node2) in dg.neighbors_inner(id).unwrap() {
                    if node2.is_h() && dg.mutual_arity(id, *id2).unwrap() > 1 {
                        let id2 = *id2;
                        return Some(HHopfData { dg, z: id, h: id2 });
                    }
                }
            }
        }
        None
    }
}

impl<'a> Rule<ZH> for HHopfData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, z, h } = self;
        let w = dg.mutual_arity(z, h).unwrap();
        dg.remove_wires(z, h, Some(w - 1)).unwrap();
        dg.scalar *= std::f64::consts::SQRT_2.powi((w - 1) as i32);
        if dg.arity(z).unwrap() == 2
            && dg.get_node(z).unwrap().has_phase(Phase::zero())
        {
            let nb = *dg.find_neighbor_id(z, |id| *id != h).unwrap();
            dg.remove_node(z).unwrap();
            dg.add_wire(nb, h).unwrap();
        }
    }
}

