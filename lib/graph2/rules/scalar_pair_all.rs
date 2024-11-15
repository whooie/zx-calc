use num_complex::Complex64 as C64;
use crate::c64_eq;
use super::*;

/// Evaluate and remove all scalar subdiagrams comprising one or two nodes.
///
/// This is the comprehensive version of [`ScalarPair`], searching more
/// efficiently for all small scalar subdiagrams.
///
/// ![scalar_pair][scalar_pair]
#[embed_doc_image::embed_doc_image("scalar_pair", "assets/rules/ScalarPair.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct ScalarPairAll;

/// A scalar subdiagram comprising one or two nodes.
///
/// The first item is the ID of any generator (non-input/output) node, and the
/// second is the ID of a second node connected only to the first (which itself
/// is only connected to the second) and the number of connections. If the
/// second item is `None`, then the first node has arity 0.
pub type BasicScalar = (NodeId, Option<(NodeId, usize)>);

/// Output of [`ScalarPairAll::find`].
#[derive(Debug)]
pub struct ScalarPairAllData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) groups: Vec<BasicScalar>,
}

impl<'a, A> ScalarPairAllData<'a, A>
where A: DiagramData
{
    /// Return the number of basic scalars comprising either one or two nodes
    /// found.
    pub fn len(&self) -> usize { self.groups.len() }

    /// Return `true` if the number of groups found is zero.
    pub fn is_empty(&self) -> bool { self.groups.is_empty() }

    /// Return a reference to all groups found.
    pub fn groups(&self) -> &Vec<BasicScalar> { &self.groups }
}

impl RuleFinder<ZX> for ScalarPairAll {
    type Output<'a> = ScalarPairAllData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        let mut seen: Vec<NodeId> = Vec::new();
        let mut groups: Vec<BasicScalar> = Vec::new();
        for (id, _node) in dg.nodes_inner() {
            if seen.contains(&id) {
                continue;
            } else {
                seen.push(id);
            }
            let n = dg.arity(id).unwrap();
            if n == 0 || dg.mutual_arity(id, id).unwrap() == n {
                groups.push((id, None));
            } else if let Some(nb) =
                dg.find_neighbor_id(
                    id, |nb| dg.mutual_arity(id, nb.id()).unwrap() == n)
            {
                seen.push(nb.id());
                let n2 = Some((nb.id(), n));
                groups.push((id, n2));
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(ScalarPairAllData { dg, groups })
        }
    }
}

impl<'a> Rule<ZX> for ScalarPairAllData<'a, ZX> {
    fn simplify(self) {
        use ZXNode::*;
        use std::f64::consts::FRAC_1_SQRT_2;
        let Self { dg, groups } = self;
        for (n1, n2) in groups.into_iter() {
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
}

impl RuleFinder<ZH> for ScalarPairAll {
    type Output<'a> = ScalarPairAllData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let mut seen: Vec<NodeId> = Vec::new();
        let mut groups: Vec<BasicScalar> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if seen.contains(&id) {
                continue;
            } else {
                seen.push(id);
            }
            if node.is_generator() {
                let n = dg.arity(id).unwrap();
                if n == 0 {
                    groups.push((id, None));
                } else if let Some(id2) =
                    dg.find_neighbor_id(
                        id, |id2| dg.mutual_arity(id, *id2).unwrap() == n)
                {
                    seen.push(*id2);
                    let n2 = Some((*id2, n));
                    groups.push((id, n2));
                }
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(ScalarPairAllData { dg, groups })
        }
    }
}

impl<'a> Rule<ZH> for ScalarPairAllData<'a, ZH> {
    fn simplify(self) {
        use ZHNode::*;
        use std::f64::consts::{ SQRT_2, FRAC_1_SQRT_2 };
        let Self { dg, groups } = self;
        for (n1, n2) in groups.into_iter() {
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
}

