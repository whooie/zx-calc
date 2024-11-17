use num_complex::Complex64 as C64;
use crate::c64_eq;
use super::*;

/// Evaluate and remove a scalar subdiagram comprising one or two nodes.
///
/// ![scalar_pair][scalar_pair]
#[embed_doc_image::embed_doc_image("scalar_pair", "assets/rules/ScalarPair.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct ScalarPair;

/// Output of [`ScalarPair::find`].
#[derive(Debug)]
pub struct ScalarPairData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) n1: NodeId,
    pub(crate) n2: Option<(NodeId, usize)>, // (second node, mutual arity)
}

impl RuleFinder<ZX> for ScalarPair {
    type Output<'a> = ScalarPairData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        for (id, _node) in dg.nodes_inner() {
            let n = dg.arity(id).unwrap();
            if n == 0 || dg.mutual_arity(id, id).unwrap() == n {
                return Some(ScalarPairData { dg, n1: id, n2: None });
            } else if let Some(nb) =
                dg.find_neighbor_id(
                    id, |nb| dg.mutual_arity(id, nb.id()).unwrap() == n)
            {
                let n2 = Some((nb.id(), n));
                return Some(ScalarPairData { dg, n1: id, n2 });
            }
        }
        None
    }
}

impl<'a> Rule<ZX> for ScalarPairData<'a, ZX> {
    fn simplify(self) {
        use ZXNode::*;
        use std::f64::consts::FRAC_1_SQRT_2;
        let Self { dg, n1, n2 } = self;
        if let Some((n2, _)) = n2 {
            let (n1, wires) = dg.delete_node(n1).unwrap();
            let (n2, _    ) = dg.delete_node(n2).unwrap();
            let ne = wires.iter().filter(|w| w.is_e()).count();
            let nh = wires.iter().filter(|w| w.is_h()).count();
            let z: C64 =
                match (n1, n2) {
                    (Z(ph1), Z(ph2)) | (X(ph1), X(ph2)) if ne > 0 =>
                        if nh % 2 == 0 {
                            (1.0 + (ph1 + ph2).cis())
                                * FRAC_1_SQRT_2.powi(nh as i32)
                        } else {
                            (1.0 - (ph1 + ph2).cis())
                                * FRAC_1_SQRT_2.powi(nh as i32)
                        },
                    (Z(ph1), Z(ph2)) | (X(ph1), X(ph2)) =>
                        if nh % 2 == 0 {
                            (1.0 + ph1.cis()) * (1.0 + ph2.cis())
                                * FRAC_1_SQRT_2.powi(nh as i32)
                        } else {
                            ((1.0 + ph2.cis()) + ph1.cis() * (1.0 - ph2.cis()))
                                * FRAC_1_SQRT_2.powi(nh as i32)
                        },
                    (Z(ph1), X(ph2)) | (X(ph2), Z(ph1)) if nh > 0 =>
                        if ne % 2 == 0 {
                            (1.0 + (ph1 + ph2).cis())
                                * FRAC_1_SQRT_2.powi(ne as i32)
                        } else {
                            (1.0 - (ph1 + ph2).cis())
                                * FRAC_1_SQRT_2.powi(ne as i32)
                        },
                    (Z(ph1), X(ph2)) | (X(ph2), Z(ph1)) =>
                        if ne % 2 == 0 {
                            (1.0 + ph1.cis()) * (1.0 + ph2.cis())
                                * FRAC_1_SQRT_2.powi(ne as i32)
                        } else {
                            ((1.0 + ph2.cis()) + ph1.cis() * (1.0 - ph2.cis()))
                                * FRAC_1_SQRT_2.powi(ne as i32)
                        },
                    _ => unreachable!(),
                };
            dg.scalar *= z;
        } else {
            let (n1, wires) = dg.delete_node(n1).unwrap();
            let nh = wires.iter().filter(|w| w.is_h()).count() / 2;
            let z: C64 =
                match n1 {
                    Z(ph) | X(ph) =>
                        if nh % 2 == 0 {
                            (1.0 + ph.cis()) * FRAC_1_SQRT_2.powi(nh as i32)
                        } else {
                            (1.0 - ph.cis()) * FRAC_1_SQRT_2.powi(nh as i32)
                        },
                    _ => unreachable!(),
                };
            dg.scalar *= z;
        }
    }
}

impl RuleFinder<ZH> for ScalarPair {
    type Output<'a> = ScalarPairData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if node.is_generator() {
                let n = dg.arity(id).unwrap();
                if n == 0 {
                    return Some(ScalarPairData { dg, n1: id, n2: None });
                } else if let Some(id2) =
                    dg.find_neighbor_id(
                        id, |id2| dg.mutual_arity(id, *id2).unwrap() == n)
                {
                    let n2 = Some((*id2, n));
                    return Some(ScalarPairData { dg, n1: id, n2 });
                }
            }
        }
        None
    }
}

impl<'a> Rule<ZH> for ScalarPairData<'a, ZH> {
    fn simplify(self) {
        use ZHNode::*;
        use std::f64::consts::{ SQRT_2, FRAC_1_SQRT_2 };
        let Self { dg, n1, n2 } = self;
        if let Some((n2, nwires)) = n2 {
            let (n1, _) = dg.delete_node(n1).unwrap();
            let (n2, _) = dg.delete_node(n2).unwrap();
            let z: C64 =
                match (n1, n2) {
                    (Z(ph1), Z(ph2)) | (X(ph1), X(ph2)) =>
                        1.0 + (ph1 + ph2).cis(),

                    (Z(ph1), X(ph2)) | (X(ph2), Z(ph1)) =>
                        if nwires % 2 == 0 {
                            (1.0 + ph1.cis()) * (1.0 + ph2.cis())
                                * FRAC_1_SQRT_2.powi(nwires as i32)
                        } else {
                            ((1.0 + ph2.cis()) + ph1.cis() * (1.0 - ph2.cis()))
                                * FRAC_1_SQRT_2.powi(nwires as i32)
                        },

                    (Z(ph), H(a)) | (H(a), Z(ph)) =>
                        if c64_eq(a, -1.0) && nwires == 2 {
                            (1.0 - ph.cis()) * FRAC_1_SQRT_2
                        } else {
                            1.0 + a * ph.cis()
                        },

                    (X(ph), H(a)) | (H(a), X(ph)) =>
                        if c64_eq(a, -1.0) && nwires == 2 {
                            (1.0 - ph.cis()) * FRAC_1_SQRT_2
                        } else if nwires % 2 == 0 {
                            SQRT_2.powi(nwires as i32)
                            + (a - 1.0) * (1.0 + ph.cis())
                                * FRAC_1_SQRT_2.powi(nwires as i32)
                        } else {
                            SQRT_2.powi(nwires as i32)
                            + (a - 1.0) * (1.0 - ph.cis())
                                * FRAC_1_SQRT_2.powi(nwires as i32)
                        },

                    (H(a), H(b)) => {
                        let a_check = c64_eq(a, -1.0);
                        let b_check = c64_eq(b, -1.0);
                        match (a_check, b_check) {
                            (true, true) if nwires == 2 =>
                                2.0.into(),
                            (true, false) if nwires == 2 =>
                                2.0 * SQRT_2 - (b + 1.0) * FRAC_1_SQRT_2,
                            (false, true) if nwires == 2 =>
                                2.0 * SQRT_2 - (a + 1.0) * FRAC_1_SQRT_2,
                            _ =>
                                2.0_f64.powi(nwires as i32) - 1.0 + a * b,
                        }
                    },
                    _ => unreachable!(),
                };
            dg.scalar *= z;
        } else {
            let (n1, _) = dg.delete_node(n1).unwrap();
            let z: C64 =
                match n1 {
                    Z(ph) => 1.0 + ph.cis(),
                    X(ph) => 1.0 + ph.cis(),
                    H(a) => a,
                    _ => unreachable!(),
                };
            dg.scalar *= z;
        }
    }
}

