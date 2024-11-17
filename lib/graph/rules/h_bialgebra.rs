use num_complex::Complex64 as C64;
use crate::{ c64_eq, phase::Phase };
use super::*;

/// H-box version of [`Bialgebra`], rewriting two Z-spiders and two H-boxes in a
/// checkerboard square according to the (co)monoids they comprise.
///
/// Any identities left by the transformations below are automatically removed.
/// If both H-boxes are trinary with argument –1, the rule results in an
/// X-spider, and if both H-boxes are binary, [`HStateMul`] is automatically
/// applied.
///
/// ![h_bialgebra][h_bialgebra]
#[embed_doc_image::embed_doc_image("h_bialgebra", "assets/rules/HBialgebra.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct HBialgebra;

/// Output of [`HBialgebra::find`].
#[derive(Debug)]
pub struct HBialgebraData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) z1: (NodeId, Phase, usize), // (id, phase, arity)
    pub(crate) z2: (NodeId, Phase, usize),
    pub(crate) h1: (NodeId, C64,   usize),
    pub(crate) h2: (NodeId, C64,   usize),
}

impl RuleFinder<ZH> for HBialgebra {
    type Output<'a> = HBialgebraData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let n0: ZHNodeSpec0 = |_, _, n| { n.is_z() };
        let n1: ZHNodeSpec = |_, _, _, n| { n.is_h() };
        let n2: ZHNodeSpec = |_, _, _, n| { n.is_z() };
        let n3: ZHNodeSpec = |dg, p, id, n| {
            n.is_h() && dg.mutual_arity(p[0].0, *id).unwrap() > 0
        };
        let p = dg.find_path(dg.nodes_inner(), n0, &[n1, n2, n3])?;
        let z1 = (p[0].0, p[0].1.phase().unwrap(), dg.arity(p[0].0).unwrap());
        let z2 = (p[2].0, p[2].1.phase().unwrap(), dg.arity(p[2].0).unwrap());
        let h1 = (p[1].0, p[1].1.arg().unwrap(),   dg.arity(p[1].0).unwrap());
        let h2 = (p[3].0, p[3].1.arg().unwrap(),   dg.arity(p[3].0).unwrap());
        Some(HBialgebraData { dg, z1, z2, h1, h2 })
    }
}

impl<'a> Rule<ZH> for HBialgebraData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, z1, z2, h1, h2 } = self;
        // remove wires internal to the rule to make handling the Z side easier
        dg.remove_wires(z1.0, h1.0, Some(1)).unwrap();
        dg.remove_wires(z1.0, h2.0, Some(1)).unwrap();
        dg.remove_wires(z2.0, h1.0, Some(1)).unwrap();
        dg.remove_wires(z2.0, h2.0, Some(1)).unwrap();

        // handle the H side first by top-level cases; output a connection point
        // for use when handling the Z side
        let hz: NodeId;
        if h1.2 == 2 && h2.2 == 2 {
            // hstatemul case
            let new_h = dg.add_node(ZHNode::H(h1.1 * h2.1));
            // both original H-boxes should be disconnected by this point
            dg.delete_node(h1.0);
            dg.delete_node(h2.0);
            hz = new_h;
        } else if h1.2 == 3 && c64_eq(h1.1, -1.0)
            && h2.2 == 3 && c64_eq(h2.1, -1.0)
        {
            // colorflip case
            let x = dg.add_node(ZHNode::x());
            let h = dg.add_node(ZHNode::h());
            dg.add_wire(x, h).unwrap();
            // both original H-boxes should have only one neighbor by this point
            let nb1 = *dg.get_neighbor(h1.0).unwrap();
            let nb2 = *dg.get_neighbor(h2.0).unwrap();
            dg.add_wire(x, nb1).unwrap();
            dg.add_wire(x, nb2).unwrap();
            hz = h;
        } else {
            // general case
            let h = dg.add_node(ZHNode::h());
            let mh = dg.add_node(ZHNode::h());
            let z = dg.add_node(ZHNode::z());
            dg.add_wire(h, mh).unwrap();
            dg.add_wire(mh, z).unwrap();
            // both original H-boxes will remain; we just need to connect
            // everything up
            dg.add_wire(h1.0, z).unwrap();
            dg.add_wire(h2.0, z).unwrap();
            hz = h;
        }

        // now handle the Z side, removing the original spiders if they're
        // identities
        if z1.2 == 3 && z1.1 == Phase::zero() {
            // should only have one neighbor by this point
            let nb = *dg.get_neighbor(z1.0).unwrap();
            dg.add_wire(hz, nb).unwrap();
        } else {
            dg.add_wire(hz, z1.0).unwrap();
        }
        if z2.2 == 3 && z2.1 == Phase::zero() {
            // should only have one neighbor by this point
            let nb = *dg.get_neighbor(z2.0).unwrap();
            dg.add_wire(hz, nb).unwrap();
        } else {
            dg.add_wire(hz, z2.0).unwrap();
        }
    }
}

