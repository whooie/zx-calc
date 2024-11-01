use crate::phase::Phase;
use super::*;

/// Combine the arguments of two H-boxes connecting an X- and a Z-spider into a
/// single average.
///
/// ![h_avg][h_avg]
#[embed_doc_image::embed_doc_image("h_avg", "assets/rules/HAvg.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HAvg;

/// Output of [`HAvg::find`].
#[derive(Debug)]
pub struct HAvgData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) h1: NodeId,
    pub(crate) h2: NodeId,
    pub(crate) x: NodeId,
    pub(crate) z: (NodeId, Phase, usize),
}

impl RuleSeal for HAvg { }
impl RuleFinder for HAvg {
    type Output<'a> = HAvgData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let n0: PathNodeSpec = |dg, _, id, n| {
            n.is_x_and(|ph| ph == Phase::pi()) && dg.arity(id).unwrap() == 2
        };
        let n1: PathNodeSpec = |dg, _, id, n| {
            n.is_h() && dg.arity(id).unwrap() == 2
        };
        let n2: PathNodeSpec = |_, _, _, n| {
            n.is_z_and(|ph| ph == Phase::pi())
        };
        let n3: PathNodeSpec = |dg, p, id, n| {
            n.is_h() && dg.arity(id).unwrap() == 2
                && dg.mutual_arity(id, p[0].0).unwrap() == 1
        };
        let p = dg.find_path(dg.nodes_inner(), &[n0, n1, n2, n3])?;
        let h1 = p[1].0;
        let h2 = p[3].0;
        let x = p[0].0;
        let z = (p[2].0, p[2].1.phase().unwrap(), dg.arity(p[2].0).unwrap());
        Some(HAvgData { dg, h1, h2, x, z })
    }
}

impl<'a> RuleSeal for HAvgData<'a> { }
impl<'a> Rule for HAvgData<'a> {
    fn simplify(self) {
        let Self { dg, h1, h2, x, z } = self;
        dg.remove_node(x).unwrap();
        let Node::H(a1) = dg.remove_node(h1).unwrap() else { unreachable!() };
        let Node::H(a2) = dg.remove_node(h2).unwrap() else { unreachable!() };
        if z.1 == Phase::zero() && z.2 == 3 {
            let n = dg.get_node_mut(z.0).unwrap();
            *n = Node::H((a1 + a2) / 2.0);
        } else {
            let h = dg.add_node(Node::H((a1 + a2) / 2.0));
            dg.add_wire(h, z.0).unwrap();
        }
    }
}

