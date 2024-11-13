use crate::phase::Phase;
use super::*;

/// Replace an H-box on a single wire joining two spiders of the same color with
/// its Euler angle decomposition, fusing spiders afterward.
///
/// ![h_euler][h_euler]
#[embed_doc_image::embed_doc_image("h_euler", "assets/rules/HEuler.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HEuler;

/// Output of [`HEuler::find`].
#[derive(Debug)]
pub struct HEulerData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s1: NodeId, // spider 1
    pub(crate) h: NodeId, // h-box
    pub(crate) s2: NodeId, // spider 2
}

impl RuleFinder for HEuler {
    type Output<'a> = HEulerData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let n0: PathNodeSpec = |_, _, _, n| {
            n.is_spider()
        };
        let n1: PathNodeSpec = |dg, _, id, n| {
            n.is_h() && n.has_defarg() && dg.arity(id).unwrap() == 2
        };
        let n2: PathNodeSpec = |_, p, _, n| {
            n.is_same_color(p[0].1)
        };
        let path = dg.find_path(dg.nodes_inner(), &[n0, n1, n2])?;
        let (s1, h, s2) = (path[0].0, path[1].0, path[2].0);
        Some(HEulerData { dg, s1, h, s2 })
    }
}

impl<'a> Rule for HEulerData<'a> {
    fn simplify(self) {
        let Self { dg, s1, h, s2 } = self;
        dg.nodes[s1].as_mut().unwrap()
            .map_phase(|ph| ph + Phase::pi2());
        dg.nodes[s2].as_mut().unwrap()
            .map_phase(|ph| ph + Phase::pi2());
        if dg.nodes[s1].as_ref().unwrap().is_z() {
            let _ = dg.nodes[h].insert(Node::X(Phase::pi2()));
        } else {
            let _ = dg.nodes[h].insert(Node::Z(Phase::pi2()));
        }
        if dg.nodes[s1].as_ref().unwrap()
            .is_spider_and(|ph| ph == Phase::zero())
            && dg.arity(s1).unwrap() == 2
        {
            let mut nb_iter = dg.neighbors_of(s1).unwrap().map(fst);
            let nb1 = nb_iter.next().unwrap();
            let nb2 = nb_iter.next().unwrap();
            dg.remove_node(s1).unwrap();
            dg.add_wire(nb1, nb2).unwrap();
        }
        if dg.nodes[s2].as_ref().unwrap()
            .is_spider_and(|ph| ph == Phase::zero())
            && dg.arity(s2).unwrap() == 2
        {
            let mut nb_iter = dg.neighbors_of(s2).unwrap().map(fst);
            let nb1 = nb_iter.next().unwrap();
            let nb2 = nb_iter.next().unwrap();
            dg.remove_node(s2).unwrap();
            dg.add_wire(nb1, nb2).unwrap();
        }
        dg.scalar *= (-Phase::pi4()).cis();
    }
}

