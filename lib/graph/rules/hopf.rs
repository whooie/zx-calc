use super::*;

/// Remove wires in pairs between two spiders of opposite colors.
///
/// ![hopf][hopf]
#[embed_doc_image::embed_doc_image("hopf", "assets/rules/Hopf.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Hopf;

/// Output of [`Hopf::find`].
#[derive(Debug)]
pub struct HopfData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s1: NodeId, // spider 1
    pub(crate) s2: NodeId, // spider 2
    pub(crate) n: usize, // mutual arity
}

impl RuleFinder for Hopf {
    type Output<'a> = HopfData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if node.is_z() {
                for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                    if node2.is_x() {
                        let n = dg.mutual_arity(id, id2).unwrap();
                        return Some(HopfData { dg, s1: id, s2: id2, n });
                    }
                }
            } else if node.is_x() {
                for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                    if node2.is_z() {
                        let n = dg.mutual_arity(id, id2).unwrap();
                        return Some(HopfData { dg, s1: id, s2: id2, n });
                    }
                }
            }
        }
        None
    }
}

impl<'a> Rule for HopfData<'a> {
    fn simplify(self) {
        let HopfData { dg, s1, s2, n } = self;
        dg.remove_wires(s1, s2, Some(2 * (n / 2))).unwrap();
        dg.scalar *= std::f64::consts::FRAC_1_SQRT_2.powi(2 * (n / 2) as i32)
    }
}

