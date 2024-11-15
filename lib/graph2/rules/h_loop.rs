use crate::phase::Phase;
use super::*;

/// A specific case of [`HEuler`], removing Hadamard wires connecting a spider
/// to itself and adding π to the spider's phase for each wire removed.
///
/// ![h_loop][h_loop]
#[embed_doc_image::embed_doc_image("h_loop", "assets/rules/HLoop.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HLoop;

/// Output of [`HLoop::find`].
#[derive(Debug)]
pub struct HLoopData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) s: NodeId, // spider
    pub(crate) hh: Vec<NodeId>, // h-boxes (empty if ZX)
}

impl RuleFinder<ZX> for HLoop {
    type Output<'a> = HLoopData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        for (id, _node) in dg.nodes_inner() {
            if dg.mutual_arity_h(id, id).unwrap() > 0 {
                return Some(HLoopData { dg, s: id, hh: Vec::with_capacity(0) });
            }
        }
        None
    }
}

impl<'a> Rule<ZX> for HLoopData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, s, hh: _ } = self;
        let nloops = dg.remove_wires_h(s, s, None).unwrap();
        if nloops % 2 == 1 {
            dg.get_node_mut(s).unwrap()
                .map_phase(|ph| ph + Phase::pi());
        }
        dg.scalar *= std::f64::consts::FRAC_1_SQRT_2.powi(nloops as i32);
    }
}

impl RuleFinder<ZH> for HLoop {
    type Output<'a> = HLoopData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let mut hh: Vec<NodeId> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_spider() {
                for (id2, node2) in dg.neighbors_inner(id).unwrap() {
                    if node2.is_h()
                        && node2.has_defarg()
                        && dg.arity(*id2).unwrap() == 2
                        && dg.mutual_arity(id, *id2).unwrap() == 2
                    {
                        hh.push(*id2);
                    }
                }
                if hh.is_empty() {
                    return Some(HLoopData { dg, s: id, hh });
                }
            }
        }
        None
    }
}

impl<'a> Rule<ZH> for HLoopData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, s, hh } = self;
        let nloops = hh.len();
        hh.into_iter()
            .for_each(|h| { dg.remove_node(h).unwrap(); });
        if nloops % 2 == 1 {
            dg.get_node_mut(s).unwrap()
                .map_phase(|ph| ph + Phase::pi());
        }
        dg.scalar *= std::f64::consts::FRAC_1_SQRT_2.powi(nloops as i32);
    }
}

