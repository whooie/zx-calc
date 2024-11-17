use crate::phase::Phase;
use super::*;

/// When two spiders of the same color are connected by a single binary spider
/// with phase π, commute the inner spider through the one of lesser arity and
/// fuse the outer spiders.
///
/// Empty self-wires are also removed.
///
/// ![pi_commute][pi_commute]
#[embed_doc_image::embed_doc_image("pi_commute", "assets/rules/PiCommute.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct PiCommute;

/// Output of [`PiCommute::find`].
#[derive(Debug)]
pub struct PiCommuteData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) s1: NodeId,
    pub(crate) s_pi: NodeId, // inner π-phase spider
    pub(crate) s2: NodeId,
}

impl RuleFinder<ZX> for PiCommute {
    type Output<'a> = PiCommuteData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        let n0: ZXNodeSpec0 = |_, _, _| { true };
        let n1: ZXNodeSpec = |dg, p, w, n| {
            (
                n.is_diff_color_and(p[0].1, |ph, _| ph == Phase::pi())
                && dg.arity_e(w.id()).unwrap() == 2
            ) || (
                n.is_same_color_and(p[0].1, |ph, _| ph == Phase::pi())
                && dg.arity_h(w.id()).unwrap() == 2
            )
        };
        let n2: ZXNodeSpec = |_, p, _, n| { n.is_same_color(p[0].1) };
        let p = dg.find_path(dg.nodes_inner(), n0, &[n1, n2])?;
        let s1 = p[0].0;
        let s_pi = p[1].0;
        let s2 = p[2].0;
        Some(PiCommuteData { dg, s1, s_pi, s2 })
    }
}

impl<'a> Rule<ZX> for PiCommuteData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, s1, s_pi, s2 } = self;
        dg.remove_wires_e(s1, s2, None).unwrap();
        dg.remove_node(s_pi).unwrap();
        let (smin, smax) =
            if dg.arity(s1).unwrap() < dg.arity(s2).unwrap() {
                (s1, s2)
            } else {
                (s2, s1)
            };
        let (nodemin, nnbmin) = dg.delete_node(smin).unwrap();
        let new_pi =
            if nodemin.is_z() {
                ZXNode::X(Phase::pi())
            } else {
                ZXNode::Z(Phase::pi())
            };
        let phmin = nodemin.phase().unwrap();
        dg.nodes[smax].as_mut().unwrap()
            .map_phase(|phmax| phmax - phmin);
        dg.scalar *= (-phmin).cis();
        let pi_spiders: Vec<NodeId> =
            nnbmin.into_iter()
            .filter_map(|nbmin| {
                if nbmin.has_h_id(smin) {
                    let snew1 = dg.add_node(new_pi);
                    let snew2 = dg.add_node(new_pi);
                    dg.add_wire_h(snew1, snew2).unwrap();
                    dg.add_wire(smax, snew1).unwrap();
                    dg.add_wire(smax, snew2).unwrap();
                    None
                } else if !nbmin.has_id(smin) {
                    let snew = dg.add_node(new_pi);
                    dg.add_wire(smax, snew).unwrap();
                    Some(snew)
                } else {
                    None
                }
            })
            .collect();
        dg.wires.iter_mut()
            .flatten()
            .flat_map(|nnb| nnb.iter_mut())
            .filter(|nb| nb.has_id(smin))
            .zip(pi_spiders)
            .for_each(|(nb, snew)| { nb.map_id(|_| snew); });
    }
}

impl RuleFinder<ZH> for PiCommute {
    type Output<'a> = PiCommuteData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let n0: ZHNodeSpec0 = |_, _, n| { n.is_spider() };
        let n1: ZHNodeSpec = |dg, p, id, n| {
            n.is_diff_color_and(p[0].1, |ph, _| ph == Phase::pi())
                && dg.arity(*id).unwrap() == 2
        };
        let n2: ZHNodeSpec = |_, p, _, n| { n.is_same_color(p[0].1) };
        let p = dg.find_path(dg.nodes_inner(), n0, &[n1, n2])?;
        let s1 = p[0].0;
        let s_pi = p[1].0;
        let s2 = p[2].0;
        Some(PiCommuteData { dg, s1, s_pi, s2 })
    }
}

impl<'a> Rule<ZH> for PiCommuteData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, s1, s_pi, s2 } = self;
        let n_pi = dg.remove_node(s_pi).unwrap();
        let (smin, smax) =
            if dg.arity(s1).unwrap() < dg.arity(s2).unwrap() {
                (s1, s2)
            } else {
                (s2, s1)
            };
        let (nodemin, nnbmin) = dg.delete_node(smin).unwrap();
        let phmin = nodemin.phase().unwrap();
        dg.nodes[smax].as_mut().unwrap()
            .map_phase(|phmax| phmax - phmin);
        dg.scalar *= (-phmin).cis();
        let pi_spiders: Vec<NodeId> =
            nnbmin.into_iter()
            .filter_map(|nbmin| {
                if nbmin != smin {
                    let snew = dg.add_node(n_pi);
                    dg.add_wire(smax, snew).unwrap();
                    Some(snew)
                } else {
                    None
                }
            })
            .collect();
        dg.wires.iter_mut()
            .flatten()
            .flat_map(|nnb| nnb.iter_mut())
            .filter(|nb| **nb == smin)
            .zip(pi_spiders)
            .for_each(|(nb, snew)| { *nb = snew; });
    }
}

