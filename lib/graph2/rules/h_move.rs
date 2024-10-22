use super::*;

/// Two spiders of arbitrary arity and phase with opposite colors sandwiching an
/// H-box with default argument on a single wire.
#[derive(Copy, Clone, Debug)]
pub struct HMove;

/// Output of [`HMove::find`].
#[derive(Debug)]
pub struct HMoveData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s1: NodeId, // spider 1
    pub(crate) h: NodeId, // middle h-box
    pub(crate) s2: NodeId, // spider 2
}

impl RuleSeal for HMove { }
impl RuleFinder for HMove {
    type Output<'a> = HMoveData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        let n0: PathNodeSpec = |_, _, _, n| {
            n.is_spider()
        };
        let n1: PathNodeSpec = |dg, _, id, n| {
            n.is_h() && n.has_defarg() && dg.arity(id).unwrap() == 2
        };
        let n2: PathNodeSpec = |_, p, _, n| {
            (p[0].1.is_z() && n.is_x()) || (p[0].1.is_x() && n.is_z())
        };
        let path = dg.find_path(dg.nodes_inner(), &[n0, n1, n2])?;
        let s1 = path[0].0;
        let h = path[1].0;
        let s2 = path[2].0;
        Some(HMoveData { dg, s1, h, s2 })
    }
}

impl<'a> RuleSeal for HMoveData<'a> { }
impl<'a> Rule for HMoveData<'a> {
    fn simplify(self) {
        let Self { dg, s1, h, s2 } = self;
        let _ = dg.remove_node(h).unwrap();
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
        let mut self_loops: usize = 0;
        let hboxes: Vec<NodeId> =
            nnbmin.into_iter()
            .filter_map(|nbmin| {
                if nbmin != smin {
                    let hnew = dg.add_node(Node::h());
                    dg.add_wire(smax, hnew).unwrap();
                    Some(hnew)
                } else {
                    self_loops += 1;
                    None
                }
            })
            .collect();
        self_loops /= 2; // self-loops are double-counted
        dg.wires.iter_mut()
            .flatten()
            .flat_map(|nnb| nnb.iter_mut())
            .filter(|nb| **nb == smax)
            .zip(hboxes)
            .for_each(|(nb, hnew)| { *nb = hnew; });
        (0..self_loops)
            .for_each(|_| { dg.add_wire(smax, smax).unwrap(); });
    }
}

