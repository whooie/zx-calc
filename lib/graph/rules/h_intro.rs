use crate::phase::Phase;
use super::*;

/// Run the H-box introduction rule in reverse, eliminating two identical
/// H-boxes in a loop between spiders.
///
/// ![h_intro][h_intro]
#[embed_doc_image::embed_doc_image("h_intro", "assets/rules/HIntro.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HIntro;

/// Output of [`HIntro::find`].
#[derive(Debug)]
pub struct HIntroData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) x: NodeId,
    pub(crate) h1: (NodeId, bool), // (id, is_had)
    pub(crate) h2: (NodeId, bool), // (id, is_had)
    pub(crate) z: (NodeId, bool), // (id, phase == 0)
}

impl RuleFinder for HIntro {
    type Output<'a> = HIntroData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        const EPSILON: f64 = 1e-12;
        let n0: PathNodeSpec = |_, _, _, n| {
            n.is_z()
        };
        let n1: PathNodeSpec = |dg, _, id, n| {
            n.is_x_and(|ph| ph == Phase::pi()) && dg.arity(id).unwrap() == 2
        };
        let n2: PathNodeSpec = |dg, _, id, n| {
            n.is_h() && dg.arity(id).unwrap() == 2
        };
        let n3: PathNodeSpec = |dg, p, id, n| {
            n.is_z() && dg.mutual_arity(id, p[2].0).unwrap() == 1
        };
        let n4: PathNodeSpec = |dg, p, id, n| {
            n.is_h()
                && dg.arity(id).unwrap() == 2
                && dg.mutual_arity(id, p[0].0).unwrap() == 1
                && (n.arg().unwrap() - p[2].1.arg().unwrap()).norm() < EPSILON
        };
        let p = dg.find_path(dg.nodes_inner(), &[n0, n1, n2, n3, n4])?;
        let x = p[1].0;
        let h1 = (p[2].0, p[2].1.has_defarg());
        let h2 = (p[4].0, p[4].1.has_defarg());
        let z = (p[3].0, p[3].1.has_defarg());
        Some(HIntroData { dg, x, h1, h2, z })
    }
}

impl<'a> Rule for HIntroData<'a> {
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

