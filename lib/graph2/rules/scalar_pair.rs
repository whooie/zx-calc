use num_complex::Complex64 as C64;
use super::*;

/// Evaluate and remove a scalar subdiagram comprising one or two nodes.
///
/// ![scalar_pair][scalar_pair]
#[embed_doc_image::embed_doc_image("scalar_pair", "assets/rules/ScalarPair.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct ScalarPair;

/// Output of [`ScalarPair::find`].
#[derive(Debug)]
pub struct ScalarPairData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) n1: NodeId,
    pub(crate) n2: Option<(NodeId, usize)>, // (second node, mutual arity)
}

impl RuleSeal for ScalarPair { }
impl RuleFinder for ScalarPair {
    type Output<'a> = ScalarPairData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if node.is_generator() {
                let n = dg.arity(id).unwrap();
                if n == 0 {
                    return Some(ScalarPairData { dg, n1: id, n2: None });
                } else if let Some(id2) =
                    dg.find_neighbor_id_of(
                        id, |id2| dg.mutual_arity(id, id2).unwrap() == n
                    )
                {
                    let n2 = Some((id2, n));
                    return Some(ScalarPairData { dg, n1: id, n2 });
                }
            }
        }
        None
    }
}

impl<'a> RuleSeal for ScalarPairData<'a> { }
impl<'a> Rule for ScalarPairData<'a> {
    fn simplify(self) {
        use Node::*;
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
                            (1.0 + ph1.cis()) * (1.0 * ph2.cis())
                                / 2.0_f64.powi(nwires as i32 / 2)
                        } else {
                            ((1.0 + ph1.cis()) + ph2.cis() * (1.0 - ph1.cis()))
                                / 2.0_f64.powf(nwires as f64 / 2.0)
                        },

                    (Z(ph), H(a)) | (H(a), Z(ph)) =>
                        1.0 + a * ph.cis(),

                    (X(ph), H(a)) | (H(a), X(ph)) =>
                        if nwires % 2 == 0 {
                            2.0_f64.powi(nwires as i32) - 1.0 + a
                                + ph.cis() * (a - 1.0)
                        } else {
                            2.0_f64.powi(nwires as i32) - 1.0 + a
                                - ph.cis() * (a - 1.0)
                        },

                    (H(a), H(b)) =>
                        2.0_f64.powi(nwires as i32) - 1.0 + a * b,

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

