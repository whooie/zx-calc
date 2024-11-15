use num_complex::Complex64 as C64;
use crate::{ c64_eq, phase::Phase };
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
pub struct HAvgData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) h1: (NodeId, C64),
    pub(crate) h2: (NodeId, C64),
    pub(crate) x: NodeId,
    pub(crate) z: (NodeId, Phase, usize),
}

impl RuleFinder<ZH> for HAvg {
    type Output<'a> = HAvgData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let n0: ZHNodeSpec0 = |dg, id, n| {
            n.is_x_and(|ph| ph == Phase::pi()) && dg.arity(id).unwrap() == 2
        };
        let n1: ZHNodeSpec = |dg, _, id, n| {
            n.is_h() && dg.arity(*id).unwrap() == 2
        };
        let n2: ZHNodeSpec = |_, _, _, n| {
            n.is_z_and(|ph| ph == Phase::pi())
        };
        let n3: ZHNodeSpec = |dg, p, id, n| {
            n.is_h() && dg.arity(*id).unwrap() == 2
                && dg.mutual_arity(*id, p[0].0).unwrap() == 1
        };
        let p = dg.find_path(dg.nodes_inner(), n0, &[n1, n2, n3])?;
        let h1 = (p[1].0, p[1].1.arg().unwrap());
        let h2 = (p[3].0, p[3].1.arg().unwrap());
        let x = p[0].0;
        let z = (p[2].0, p[2].1.phase().unwrap(), dg.arity(p[2].0).unwrap());
        Some(HAvgData { dg, h1, h2, x, z })
    }
}

impl<'a> Rule<ZH> for HAvgData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, h1, h2, x, z } = self;
        dg.remove_node(x).unwrap();
        let is_had1 = c64_eq(h1.1, -1.0);
        let is_had2 = c64_eq(h2.1, -1.0);
        if is_had1 && is_had2 {
            dg.remove_node(h1.0).unwrap();
            dg.remove_node(h2.0).unwrap();
            dg.get_node_mut(z.0).unwrap()
                .map_phase(|ph| ph + Phase::pi());
            return;
        }
        let rem_z = z.1 == Phase::zero() && z.2 == 3;
        let scalar =
            if is_had1 || is_had2 { std::f64::consts::SQRT_2 } else { 2.0 };
        let node =
            if rem_z {
                dg.remove_node(h1.0).unwrap();
                dg.remove_node(h2.0).unwrap();
                dg.get_node_mut(z.0).unwrap()
            } else {
                dg.remove_node(h2.0).unwrap();
                dg.get_node_mut(h1.0).unwrap()
            };
        *node = ZHNode::H((h1.1 + h2.1) / 2.0);
        dg.scalar *= scalar;
    }
}

