use super::*;

/// When two spiders of opposite color are connected by a Hadamard wire, commute
/// the H-box through the spider of lesser arity and fuse the spiders.
///
/// Empty self-loops are also removed.
///
/// ![h_move][h_move]
#[embed_doc_image::embed_doc_image("h_move", "assets/rules/HMove.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HMove;

/// Output of [`HMove::find`].
#[derive(Debug)]
pub struct HMoveData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) s1: NodeId, // spider 1
    pub(crate) h: NodeId, // middle h-box (zero for ZX diagrams)
    pub(crate) s2: NodeId, // spider 2
}

impl RuleFinder<ZX> for HMove {
    type Output<'a> = HMoveData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            for (wire, node2) in dg.neighbors_inner(id).unwrap() {
                if wire.is_h() && node2.is_diff_color(node) {
                    let id2 = wire.id();
                    return Some(HMoveData { dg, s1: id, h: 0, s2: id2 });
                }
            }
        }
        None
    }
}

impl<'a> Rule<ZX> for HMoveData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, s1, h: _, s2 } = self;
        dg.remove_wires_e(s1, s2, None).unwrap();
        let (smin, smax) =
            if dg.arity(s1).unwrap() < dg.arity(s2).unwrap() {
                (s1, s2)
            } else {
                (s2, s1)
            };
        dg.get_neighbors_mut(smin).unwrap()
            .iter_mut()
            .for_each(|wire| { wire.toggle(); });
        let (nodemin, nnbmin) = dg.delete_node(smin).unwrap();
        let phmin = nodemin.phase().unwrap();
        dg.nodes[smax].as_mut().unwrap()
            .map_phase(|phmax| phmax + phmin);
        dg.wires.iter_mut()
            .flatten()
            .flat_map(|nnb| nnb.iter_mut())
            .for_each(|nb| { if nb.has_id(smin) { nb.map_id(|_| smax); } });
        let nnbmax = dg.wires[smax].as_mut().unwrap();
        nnbmin.into_iter()
            .for_each(|nb| { if !nb.has_e_id(smin) { nnbmax.push(nb); } });
    }
}

impl RuleFinder<ZH> for HMove {
    type Output<'a> = HMoveData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let n0: ZHNodeSpec0 = |_, _, n| { n.is_spider() };
        let n1: ZHNodeSpec = |dg, _, id, n| {
            n.is_h() && n.has_defarg() && dg.arity(*id).unwrap() == 2
        };
        let n2: ZHNodeSpec = |_, p, _, n| {
            (p[0].1.is_z() && n.is_x()) || (p[0].1.is_x() && n.is_z())
        };
        let path = dg.find_path(dg.nodes_inner(), n0, &[n1, n2])?;
        let s1 = path[0].0;
        let h = path[1].0;
        let s2 = path[2].0;
        Some(HMoveData { dg, s1, h, s2 })
    }
}

impl<'a> Rule<ZH> for HMoveData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, s1, h, s2 } = self;
        dg.remove_node(h).unwrap();
        let (smin, smax) =
            if dg.arity(s1).unwrap() < dg.arity(s2).unwrap() {
                (s1, s2)
            } else {
                (s2, s1)
            };
        let (nodemin, nnbmin) = dg.delete_node(smin).unwrap();
        let phmin = nodemin.phase().unwrap();
        dg.nodes[smax].as_mut().unwrap()
            .map_phase(|phmax| phmax + phmin);
        let hboxes: Vec<NodeId> =
            nnbmin.into_iter()
            .filter_map(|nbmin| {
                if nbmin != smin {
                    let hnew = dg.add_node(ZHNode::h());
                    dg.add_wire(smax, hnew).unwrap();
                    Some(hnew)
                } else {
                    None
                }
            })
            .collect();
        dg.wires.iter_mut()
            .flatten()
            .flat_map(|nnb| nnb.iter_mut())
            .filter(|nb| **nb == smin)
            .zip(hboxes)
            .for_each(|(nb, hnew)| { *nb = hnew; });
    }
}

