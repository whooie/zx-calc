use crate::phase::Phase;
use super::*;

/// Run the H-box introduction rule in reverse, eliminating two identical
/// H-boxes in a loop between spiders.
///
/// ![h_intro][h_intro]
#[embed_doc_image::embed_doc_image("h_intro", "assets/rules/HIntro.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct HIntro;

/// Output of [`HIntro::find`].
#[derive(Debug)]
pub struct HIntroData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) x: NodeId,
    pub(crate) h1: (NodeId, bool), // (id, is_had)
    pub(crate) h2: (NodeId, bool), // (id, is_had)
    pub(crate) z: (NodeId, bool), // (id, phase == 0)
}

impl RuleFinder<ZH> for HIntro {
    type Output<'a> = HIntroData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let n0: ZHNodeSpec0 = |_, _, n| { n.is_z() };
        let n1: ZHNodeSpec = |dg, _, id, n| {
            n.is_x_and(|ph| ph == Phase::pi()) && dg.arity(*id).unwrap() == 2
        };
        let n2: ZHNodeSpec = |dg, _, id, n| {
            n.is_h() && dg.arity(*id).unwrap() == 2
        };
        let n3: ZHNodeSpec = |dg, p, id, n| {
            n.is_z() && dg.mutual_arity(*id, p[2].0).unwrap() == 1
        };
        let n4: ZHNodeSpec = |dg, p, id, n| {
            n.is_h()
                && dg.arity(*id).unwrap() == 2
                && dg.mutual_arity(*id, p[0].0).unwrap() == 1
                && n.has_arg(p[2].1.arg().unwrap())
        };
        let p = dg.find_path(dg.nodes_inner(), n0, &[n1, n2, n3, n4])?;
        let x = p[1].0;
        let h1 = (p[2].0, p[2].1.has_defarg());
        let h2 = (p[4].0, p[4].1.has_defarg());
        let z = (p[3].0, p[3].1.has_defarg());
        Some(HIntroData { dg, x, h1, h2, z })
    }
}

impl<'a> Rule<ZH> for HIntroData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, x, h1, h2, z } = self;
        if h1.1 && h2.1 {
            dg.remove_node(h1.0).unwrap();
            dg.remove_node(h2.0).unwrap();
            dg.get_node_mut(z.0).unwrap()
                .map_phase(|ph| ph + Phase::pi());
            dg.scalar /= 2.0;
        } else {
            dg.remove_node(x).unwrap();
            dg.remove_node(h2.0).unwrap();
            if dg.arity(z.0).unwrap() == 2 && z.1 {
                let h = dg.remove_node(h1.0).unwrap();
                let n = dg.get_node_mut(z.0).unwrap();
                *n = h;
            }
        }
    }
}

