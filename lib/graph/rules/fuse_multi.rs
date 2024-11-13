use crate::phase::Phase;
use super::*;

/// Fuse groups of multiple connected spiders with the same color.
///
/// This rule lies between [`Fuse`] and [`FuseAll`] in terms of efficiency: Like
/// `FuseAll`, it finds more than two connected spiders to fuse together, but
/// like `Fuse` it won't always perform all available fuse operations. All
/// self-loops within the subgraphs induced by each group are also removed.
///
/// ![fuse][fuse]
#[embed_doc_image::embed_doc_image("fuse", "assets/rules/Fuse.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct FuseMulti;

/// Output of [`FuseMulti::find`].
#[derive(Debug)]
pub struct FuseMultiData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) init: NodeId,
    pub(crate) others: Vec<NodeId>,
    pub(crate) ph_sum: Phase,
}

impl<'a> FuseMultiData<'a> {
    /// Return the number of spiders found.
    ///
    /// Guaranteed to be at least 2.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize { 1 + self.others.len() }
}

impl RuleFinder for FuseMulti {
    type Output<'a> = FuseMultiData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        // DFS to find a group of connected spiders of the same color
        let mut others: Vec<NodeId> = Vec::new();
        let mut ph_sum: Phase;
        let mut visited: Vec<NodeId> = Vec::new();
        let mut to_visit: Vec<(NodeId, &Node)> = Vec::new();
        for (id0, n0) in dg.nodes_inner() {
            if !n0.is_spider() || visited.contains(&id0) { continue; }
            ph_sum = Phase::zero();
            visited.push(id0);
            dg.neighbors_of_inner(id0).unwrap()
                .filter(|(id, n)| {
                    n.is_same_color(n0) && !visited.contains(id)
                })
                .for_each(|(id, n)| { to_visit.push((id, n)); });
            while let Some((id_next, n_next)) = to_visit.pop() {
                others.push(id_next);
                ph_sum += n_next.phase().unwrap();
                visited.push(id_next);
                dg.neighbors_of_inner(id_next).unwrap()
                    .filter(|(id, n)| {
                        n.is_same_color(n0) && !visited.contains(id)
                    })
                    .for_each(|(id, n)| { to_visit.push((id, n)); });
            }
            if !others.is_empty() {
                return Some(FuseMultiData { dg, init: id0, others, ph_sum });
            }
        }
        None
    }
}

impl<'a> Rule for FuseMultiData<'a> {
    fn simplify(self) {
        let Self { dg, init: id0, others, ph_sum } = self;
        let mut all_nb: Vec<NodeId> = Vec::new();
        // deal with wires first -- we want to use `delete_node` here and manage
        // connections manually, so we need to keep an accurate count of how
        // many wires are disappearing
        //
        // namely, we want to remove all wires within the group, as well as any
        // self-loops on any group nodes, including `init`
        //
        // every such wire will be double-counted
        let mut remove_wire_x2: usize = 0;

        // deal with the base node first, converting all wires within the group
        // to self-wires so they can all be removed in one pass
        let nnb0 = dg.wires[id0].as_mut().unwrap();
        for nb in nnb0.iter_mut() {
            if *nb == id0 {
                remove_wire_x2 += 1;
            } else if others.contains(nb) {
                remove_wire_x2 += 1;
                *nb = id0;
            }
        }
        // remove wires
        while let Some(k) =
            nnb0.iter().enumerate()
            .find_map(|(j, id)| (*id == id0).then_some(j))
        {
            nnb0.swap_remove(k);
        }

        // now deal with other nodes -- removal of their wires is taken care of
        // by `delete_node`, so we just have to count them and only keep track
        // of neighbors outside the group
        for id1 in others.iter() {
            for nb in dg.delete_node(*id1).unwrap().1.into_iter() {
                if nb == id0 || others.contains(&nb) {
                    remove_wire_x2 += 1;
                } else {
                    all_nb.push(nb);
                }
            }
        }
        // nnb0 was dropped because we needed to call `delete_node`
        dg.wires[id0].as_mut().unwrap().append(&mut all_nb);

        // take care of external connections and total wire count
        dg.wires.iter_mut().flatten()
            .flat_map(|nnb| nnb.iter_mut())
            .for_each(|nb| { if others.contains(nb) { *nb = id0; } });
        assert!(remove_wire_x2 % 2 == 0);
        dg.wire_count -= remove_wire_x2 / 2;

        dg.nodes[id0].as_mut().unwrap()
            .map_phase(|ph0| ph0 + ph_sum);
    }
}

