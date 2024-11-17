use super::*;

/// Remove wires in pairs between two spiders of opposite colors.
///
/// ![hopf][hopf]
#[embed_doc_image::embed_doc_image("hopf", "assets/rules/Hopf.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct Hopf;

/// Output of [`Hopf::find`].
#[derive(Debug)]
pub struct HopfData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) s1: NodeId, // spider 1
    pub(crate) s2: NodeId, // spider 2
    pub(crate) n: usize, // mutual (empty) arity
}

impl RuleFinder<ZX> for Hopf {
    type Output<'a> = HopfData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            for (wire, node2) in dg.neighbors_inner(id).unwrap() {
                if node2.is_diff_color(node) {
                    let id2 = wire.id();
                    let n = dg.mutual_arity_e(id, id2).unwrap();
                    if n > 1 {
                        return Some(HopfData { dg, s1: id, s2: id2, n });
                    }
                }
            }
        }
        None
    }
}

impl<'a> Rule<ZX> for HopfData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, s1, s2, n } = self;
        dg.remove_wires_e(s1, s2, Some(2 * (n / 2))).unwrap();
        dg.scalar *= std::f64::consts::FRAC_1_SQRT_2.powi(2 * (n / 2) as i32);
    }
}

impl RuleFinder<ZH> for Hopf {
    type Output<'a> = HopfData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            for (id2, node2) in dg.neighbors_inner(id).unwrap() {
                if node2.is_diff_color(node) {
                    let id2 = *id2;
                    let n = dg.mutual_arity(id, id2).unwrap();
                    if n > 1 {
                        return Some(HopfData { dg, s1: id, s2: id2, n });
                    }
                }
            }
        }
        None
    }
}

impl<'a> Rule<ZH> for HopfData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, s1, s2, n } = self;
        dg.remove_wires(s1, s2, Some(2 * (n / 2))).unwrap();
        dg.scalar *= std::f64::consts::FRAC_1_SQRT_2.powi(2 * (n / 2) as i32);
    }
}

