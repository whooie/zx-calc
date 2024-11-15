use super::*;

/// H-box version of [`Fuse`], combining two H-boxes connected by a single
/// Hadamard wire into one.
///
/// ![h_fuse][h_fuse]
#[embed_doc_image::embed_doc_image("h_fuse", "assets/rules/HFuse.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HFuse;

/// Output of [`HFuse::find`].
#[derive(Debug)]
pub struct HFuseData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) ha: NodeId, // h-box with argument
    pub(crate) h1: NodeId, // middle h-box
    pub(crate) h2: NodeId, // empty end h-box
}

impl RuleFinder<ZH> for HFuse {
    type Output<'a> = HFuseData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let n0: ZHNodeSpec0 = |_, _, n| { n.is_h() };
        let n1: ZHNodeSpec = |dg, _, id, n| {
            n.is_h() && n.has_defarg() && dg.arity(*id).unwrap() == 2
        };
        let n2: ZHNodeSpec = |_, _, _, n| { n.is_h() && n.has_defarg() };
        let p = dg.find_path(dg.nodes_inner(), n0, &[n1, n2])?;
        let ha = p[0].0;
        let h1 = p[1].0;
        let h2 = p[2].0;
        Some(HFuseData { dg, ha, h1, h2 })
    }
}

impl<'a> Rule<ZH> for HFuseData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, ha, h1, h2 } = self;
        // the scalar for the general case is √2; this gets reduced to 1 if ha
        // is binary with argument -1 or h2 is binary (with argument -1).
        let ha_is_had =
            dg.get_node(ha).unwrap().has_arg(1.0)
            && dg.arity(ha).unwrap() == 2;
        if !ha_is_had && dg.arity(h2).unwrap() != 2 {
            dg.scalar *= std::f64::consts::SQRT_2;
        }
        dg.remove_node(h1).unwrap();
        let (_, nnb) = dg.remove_node_nb(h2).unwrap();
        nnb.into_iter()
            .for_each(|nb| {
                if nb == h2 {
                    dg.add_wire(ha, ha).unwrap();
                } else {
                    dg.add_wire(ha, nb).unwrap();
                }
            });
    }
}

