use crate::phase::Phase;
use super::*;

/// Replace an H-box on a single wire joining two spiders of the same color with
/// its Euler angle decomposition, fusing spiders afterward.
///
/// ![h_euler][h_euler]
#[embed_doc_image::embed_doc_image("h_euler", "assets/rules/HEuler.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct HEuler;

/// Output of [`HEuler::find`].
#[derive(Debug)]
pub struct HEulerData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) s1: NodeId, // spider 1
    pub(crate) h: NodeId, // h-box (zero if ZX)
    pub(crate) s2: NodeId, // spider 2
}

impl RuleFinder<ZX> for HEuler {
    type Output<'a> = HEulerData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            for (wire, node2) in dg.neighbors_inner(id).unwrap() {
                if wire.is_h() && node2.is_same_color(node) {
                    let id2 = wire.id();
                    return Some(HEulerData { dg, s1: id, h: 0, s2: id2 });
                }
            }
        }
        None
    }
}

impl<'a> Rule<ZX> for HEulerData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, s1, h: _, s2 } = self;
        dg.nodes[s1].as_mut().unwrap()
            .map_phase(|ph| ph + Phase::pi2());
        dg.nodes[s2].as_mut().unwrap()
            .map_phase(|ph| ph + Phase::pi2());
        let new =
            if dg.nodes[s1].as_ref().unwrap().is_z() {
                dg.add_node(ZXNode::X(Phase::pi2()))
            } else {
                dg.add_node(ZXNode::Z(Phase::pi2()))
            };
        dg.remove_wires_h(s1, s2, Some(1)).unwrap();
        dg.add_wire(s1, new).unwrap();
        dg.add_wire(new, s2).unwrap();
        if dg.nodes[s1].as_ref().unwrap()
            .is_spider_and(|ph| ph == Phase::zero())
            && dg.arity(s1).unwrap() == 2
            && dg.arity_e(s1).unwrap() == 2
        {
            let mut nb_iter = dg.neighbors(s1).unwrap().map(|(nb, _)| nb.id());
            let nb1 = nb_iter.next().unwrap();
            let nb2 = nb_iter.next().unwrap();
            dg.remove_node(s1).unwrap();
            dg.add_wire(nb1, nb2).unwrap();
        }
        if dg.nodes[s2].as_ref().unwrap()
            .is_spider_and(|ph| ph == Phase::zero())
            && dg.arity(s2).unwrap() == 2
            && dg.arity_e(s2).unwrap() == 2
        {
            let mut nb_iter = dg.neighbors(s2).unwrap().map(|(nb, _)| nb.id());
            let nb1 = nb_iter.next().unwrap();
            let nb2 = nb_iter.next().unwrap();
            dg.remove_node(s2).unwrap();
            dg.add_wire(nb1, nb2).unwrap();
        }
        dg.scalar *= (-Phase::pi4()).cis();
    }
}

impl RuleFinder<ZH> for HEuler {
    type Output<'a> = HEulerData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let n0: ZHNodeSpec0 = |_, _, n| { n.is_spider() };
        let n1: ZHNodeSpec = |dg, _, id, n| {
            n.is_h() && n.has_defarg() && dg.arity(*id).unwrap() == 2
        };
        let n2: ZHNodeSpec = |_, p, _, n| {
            n.is_same_color(p[0].1)
        };
        let path = dg.find_path(dg.nodes_inner(), n0, &[n1, n2])?;
        let (s1, h, s2) = (path[0].0, path[1].0, path[2].0);
        Some(HEulerData { dg, s1, h, s2 })
    }
}

impl<'a> Rule<ZH> for HEulerData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, s1, h, s2 } = self;
        dg.nodes[s1].as_mut().unwrap()
            .map_phase(|ph| ph + Phase::pi2());
        dg.nodes[s2].as_mut().unwrap()
            .map_phase(|ph| ph + Phase::pi2());
        if dg.nodes[s1].as_ref().unwrap().is_z() {
            dg.nodes[h] = Some(ZHNode::X(Phase::pi2()));
        } else {
            dg.nodes[h] = Some(ZHNode::Z(Phase::pi2()));
        }
        if dg.nodes[s1].as_ref().unwrap()
            .is_spider_and(|ph| ph == Phase::zero())
            && dg.arity(s1).unwrap() == 2
        {
            let mut nb_iter = dg.neighbors(s1).unwrap().map(|(nb, _)| *nb);
            let nb1 = nb_iter.next().unwrap();
            let nb2 = nb_iter.next().unwrap();
            dg.remove_node(s1).unwrap();
            dg.add_wire(nb1, nb2).unwrap();
        }
        if dg.nodes[s2].as_ref().unwrap()
            .is_spider_and(|ph| ph == Phase::zero())
            && dg.arity(s2).unwrap() == 2
        {
            let mut nb_iter = dg.neighbors(s2).unwrap().map(|(nb, _)| *nb);
            let nb1 = nb_iter.next().unwrap();
            let nb2 = nb_iter.next().unwrap();
            dg.remove_node(s2).unwrap();
            dg.add_wire(nb1, nb2).unwrap();
        }
        dg.scalar *= (-Phase::pi4()).cis();
    }
}

