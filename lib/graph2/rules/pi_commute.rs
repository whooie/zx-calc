use crate::phase::Phase;
use super::*;

/// When two spiders of the same color are connected by a single binary spider
/// with phase π, commute the inner spider through the one of lesser arity and
/// fuse the outer spiders.
///
/// ![pi_commute][pi_commute]
#[embed_doc_image::embed_doc_image("pi_commute", "assets/rules/PiCommute.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct PiCommute;

/// Output of [`PiCommute::find`].
#[derive(Debug)]
pub struct PiCommuteData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s1: NodeId,
    pub(crate) s_pi: NodeId, // inner π-phase spider
    pub(crate) s2: NodeId,
}

impl RuleFinder for PiCommute {
    type Output<'a> = PiCommuteData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let n0: PathNodeSpec = |_, _, _, n| {
            n.is_spider()
        };
        let n1: PathNodeSpec = |dg, p, id, n| {
            n.is_diff_color_and(p[0].1, |ph, _| ph == Phase::pi())
                && dg.arity(id).unwrap() == 2
        };
        let n2: PathNodeSpec = |_, p, _, n| {
            n.is_same_color(p[0].1)
        };
        let p = dg.find_path(dg.nodes_inner(), &[n0, n1, n2])?;
        let s1 = p[0].0;
        let s_pi = p[1].0;
        let s2 = p[2].0;
        Some(PiCommuteData { dg, s1, s_pi, s2 })
    }
}

impl<'a> Rule for PiCommuteData<'a> {
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
        let mut self_loops: usize = 0;
        let pi_spiders: Vec<NodeId> =
            nnbmin.into_iter()
            .filter_map(|nbmin| {
                if nbmin != smin {
                    let snew = dg.add_node(n_pi);
                    dg.add_wire(smax, snew).unwrap();
                    dg.add_wire(snew, nbmin).unwrap();
                    Some(snew)
                } else {
                    self_loops += 1;
                    None
                }
            })
            .collect();
        self_loops /= 2; // self-loops are double-counted
        dg.wires.iter_mut()
            .flatten()
            .flat_map(|nnb| nnb.iter_mut())
            .filter(|nb| **nb == smax)
            .zip(pi_spiders)
            .for_each(|(nb, snew)| { *nb = snew; });
        (0..self_loops)
            .for_each(|_| { dg.add_wire(smax, smax).unwrap(); });
    }
}

