use crate::phase::Phase;
use super::*;

/// Rewrite two Z-spiders and two X-spiders in a checkerboard square according
/// to the (co)monoids they comprise.
///
/// Any identities left by the transformations below are automatically removed.
///
/// ![bialgebra][bialgebra]
#[embed_doc_image::embed_doc_image("bialgebra", "assets/rules/Bialgebra.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Bialgebra;

/// Output of [`Bialgebra::find`].
#[derive(Debug)]
pub struct BialgebraData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) z1: (NodeId, Phase, usize), // (id, phase, arity)
    pub(crate) z2: (NodeId, Phase, usize),
    pub(crate) x1: (NodeId, Phase, usize),
    pub(crate) x2: (NodeId, Phase, usize),
}

impl RuleFinder for Bialgebra {
    type Output<'a> = BialgebraData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let n0: PathNodeSpec = |_, _, _, n| { n.is_z() };
        let n1: PathNodeSpec = |_, _, _, n| { n.is_x() };
        let n2: PathNodeSpec = |_, _, _, n| { n.is_z() };
        let n3: PathNodeSpec = |dg, p, id, n| {
            n.is_x() && dg.mutual_arity(p[0].0, id).unwrap() > 0
        };
        let p = dg.find_path(dg.nodes_inner(), &[n0, n1, n2, n3])?;
        let z1 = (p[0].0, p[0].1.phase().unwrap(), dg.arity(p[0].0).unwrap());
        let z2 = (p[2].0, p[2].1.phase().unwrap(), dg.arity(p[2].0).unwrap());
        let x1 = (p[1].0, p[1].1.phase().unwrap(), dg.arity(p[1].0).unwrap());
        let x2 = (p[3].0, p[3].1.phase().unwrap(), dg.arity(p[3].0).unwrap());
        Some(BialgebraData { dg, z1, z2, x1, x2 })
    }
}

impl<'a> Rule for BialgebraData<'a> {
    fn simplify(self) {
        let Self { dg, z1, z2, x1, x2 } = self;
        let pi = Phase::pi();
        dg.scalar *= std::f64::consts::FRAC_1_SQRT_2;
        let (new_z, move_ph_z) =
            if z1 == z2 && z1.1.is_mult(2) && z2.1.is_mult(2) {
                if z1.1 == pi { dg.scalar *= -1.0; }
                (dg.add_node(Node::Z(z1.1)), true)
            } else {
                (dg.add_node(Node::z()), false)
            };
        let (new_x, move_ph_x) =
            if x1 == x2 && x1.1.is_mult(2) && x2.1.is_mult(2) {
                if x1.1 == pi { dg.scalar *= -1.0; }
                (dg.add_node(Node::X(x1.1)), true)
            } else {
                (dg.add_node(Node::x()), false)
            };
        dg.add_wire(new_z, new_x).unwrap();

        // when handling nodes, pre-emptively remove identity spiders that would
        // have been left by this rule
        if z1.2 == 3 && (z1.1 == Phase::zero() || move_ph_z) {
            let nb = dg.find_neighbor_id_of(z1.0, |id| id != x1.0 && id != x2.0)
                .unwrap();
            dg.remove_node(z1.0).unwrap();
            dg.add_wire(nb, new_x).unwrap();
        } else {
            if move_ph_z {
                dg.get_node_mut(z1.0).unwrap()
                    .map_phase(|_| Phase::zero());
            }
            dg.remove_wires(z1.0, x1.0, Some(1)).unwrap();
            dg.remove_wires(z1.0, x2.0, Some(1)).unwrap();
            dg.add_wire(z1.0, new_x).unwrap();
        }
        if z2.2 == 3 && (z2.1 == Phase::zero() || move_ph_z) {
            let nb = dg.find_neighbor_id_of(z2.0, |id| id != x1.0 && id != x2.0)
                .unwrap();
            dg.remove_node(z2.0).unwrap();
            dg.add_wire(nb, new_x).unwrap();
        } else {
            if move_ph_z {
                dg.get_node_mut(z2.0).unwrap()
                    .map_phase(|_| Phase::zero());
            }
            dg.remove_wires(z2.0, x1.0, Some(1)).unwrap();
            dg.remove_wires(z2.0, x2.0, Some(1)).unwrap();
            dg.add_wire(z2.0, new_x).unwrap();
        }
        // z spiders may have already been removed
        if x1.2 == 3 && (x1.1 == Phase::zero() || move_ph_x) {
            let nb = dg.find_neighbor_id_of(x1.0, |id| id != z1.0 && id != z2.0)
                .unwrap();
            dg.remove_node(x1.0).unwrap();
            dg.add_wire(nb, new_z).unwrap();
        } else {
            if move_ph_x {
                dg.get_node_mut(x1.0).unwrap()
                    .map_phase(|_| Phase::zero());
            }
            dg.remove_wires(x1.0, z1.0, Some(1)).ok();
            dg.remove_wires(x1.0, z2.0, Some(1)).ok();
            dg.add_wire(x1.0, new_z).unwrap();
        }
        if x2.2 == 3 && (x2.1 == Phase::zero() || move_ph_x) {
            let nb = dg.find_neighbor_id_of(x2.0, |id| id != z1.0 && id != z2.0)
                .unwrap();
            dg.remove_node(x2.0).unwrap();
            dg.add_wire(nb, new_z).unwrap();
        } else {
            if move_ph_x {
                dg.get_node_mut(x2.0).unwrap()
                    .map_phase(|_| Phase::zero());
            }
            dg.remove_wires(x2.0, z1.0, Some(1)).ok();
            dg.remove_wires(x2.0, z2.0, Some(1)).ok();
            dg.add_wire(x2.0, new_z).unwrap();
        }
    }
}

