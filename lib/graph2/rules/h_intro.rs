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
    pub(crate) z1: NodeId,
    pub(crate) x: NodeId,
    pub(crate) h1: NodeId,
    pub(crate) h2: NodeId,
    pub(crate) z2: NodeId,
}

impl RuleSeal for HIntro { }
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
        let z1 = p[0].0;
        let x = p[1].0;
        let h1 = p[2].0;
        let h2 = p[4].0;
        let z2 = p[3].0;
        Some(HIntroData { dg, z1, x, h1, h2, z2 })
    }
}

impl<'a> RuleSeal for HIntroData<'a> { }
impl<'a> Rule for HIntroData<'a> {
    fn simplify(self) {
        let Self { dg, z1: _, x, h1, h2, z2 } = self;
        dg.remove_node(x).unwrap();
        dg.remove_node(h2).unwrap();
        if dg.arity(z2).unwrap() == 2 {
            let nb = dg.find_neighbor_id_of(z2, |id| id != h1).unwrap();
            dg.remove_node(z2).unwrap();
            dg.add_wire(h1, nb).unwrap();
        }
    }
}

