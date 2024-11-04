use super::*;

/// Remove Hadamard wires in pairs between two spiders of the same color.
///
/// ![h2hopf][h2hopf]
#[embed_doc_image::embed_doc_image("h2hopf", "assets/rules/H2Hopf.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct H2Hopf;

/// Output of [`H2Hopf::find`].
#[derive(Debug)]
pub struct H2HopfData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s1: NodeId, // spider 1
    pub(crate) s2: NodeId, // spider 2
}

impl RuleSeal for H2Hopf { }
impl RuleFinder for H2Hopf {
    type Output<'a> = H2HopfData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<H2HopfData<'_>> {
        let n0: PathNodeSpec = |_, _, _, n| {
            n.is_z() || n.is_x()
        };
        let n1: PathNodeSpec = |dg, _, id, n| {
            n.is_h() && n.has_defarg() && dg.arity(id).unwrap() == 2
        };
        let n2: PathNodeSpec = |dg, p, id, n| {
            n.is_same_color(p[0].1)
                && dg.mutual_arity(p[0].0, id).unwrap() == 0
        };
        let n3: PathNodeSpec = |dg, p, id, n| {
            n.is_h() && n.has_defarg() && dg.arity(id).unwrap() == 2
                && dg.mutual_arity(id, p[0].0).unwrap() == 1
        };
        let path = dg.find_path(dg.nodes_inner(), &[n0, n1, n2, n3])?;
        let s1 = path[0].0;
        let s2 = path[2].0;
        Some(H2HopfData { dg, s1, s2 })
    }
}

impl<'a> RuleSeal for H2HopfData<'a> { }
impl<'a> Rule for H2HopfData<'a> {
    fn simplify(self) {
        let Self { dg, s1, s2 } = self;
        let hboxes: Vec<NodeId> =
            dg.neighbors_of(s1).unwrap()
            .filter_map(|(id, node)| {
                (
                    node.is_h()
                    && node.has_defarg()
                    && dg.arity(id).unwrap() == 2
                    && dg.mutual_arity(id, s2).unwrap() == 1
                ).then_some(id)
            })
            .collect();
        let n = hboxes.len();
        hboxes.into_iter().take(2 * (n / 2))
            .for_each(|h| { dg.remove_node(h).unwrap(); });
        dg.scalar *= std::f64::consts::FRAC_1_SQRT_2.powi(2 * (n / 2) as i32);
    }
}

