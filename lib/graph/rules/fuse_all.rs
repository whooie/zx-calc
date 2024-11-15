use crate::phase::Phase;
use super::*;

/// Fuse all groups of connected spiders with the same color.
///
/// This rule is a comprehensive version of [`Fuse`], searching more efficiently
/// for all groups of spiders that satsfy the conditions for the rule. All
/// self-loops within the subgraphs induced by each group are also removed.
///
/// ![fuse][fuse]
#[embed_doc_image::embed_doc_image("fuse", "assets/rules/Fuse.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct FuseAll;

/// A fusable group of spiders.
///
/// The first item gives the node ID of a "base" spider that will still exist
/// after fusing; the second item is a list of all other spiders in the group;
/// and the third is the sum of all non-base spiders' phases.
pub type SpiderGroup = (NodeId, Vec<NodeId>, Phase);
//                     (initial spider, other spiders, sum of others' phases)

/// Output of [`FuseAll::find`].
#[derive(Debug)]
pub struct FuseAllData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    // all groups of spiders; groups are guaranteed disconnected from each other
    pub(crate) groups: Vec<SpiderGroup>,
}

impl<'a, A> FuseAllData<'a, A>
where A: DiagramData
{
    /// Return the number of connected groups found.
    ///
    /// All groups are guaranteed to be disconnected from each other and contain
    /// at least two spiders.
    pub fn len(&self) -> usize { self.groups.len() }

    /// Return `true` if the number of groups is zero.
    pub fn is_empty(&self) -> bool { self.groups.is_empty() }

    /// Return a reference to all groups.
    pub fn groups(&self) -> &Vec<SpiderGroup> { &self.groups }
}

impl RuleFinder<ZX> for FuseAll {
    type Output<'a> = FuseAllData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        // DFS to find all groups of connected spiders of the same color,
        // starting at every non-io node and ignoring ones already seen
        let mut groups: Vec<SpiderGroup> = Vec::new();
        let mut group: Vec<NodeId> = Vec::new();
        let mut ph_sum: Phase;
        let mut visited: Vec<NodeId> = Vec::new();
        let mut to_visit: Vec<(NodeId, &ZXNode)> = Vec::new();
        for (id0, n0) in dg.nodes_inner() {
            if visited.contains(&id0) { continue; }
            ph_sum = Phase::zero();
            visited.push(id0);
            dg.neighbors_inner(id0).unwrap()
                .filter(|(w, n)| {
                    w.is_e()
                        && n.is_same_color(n0)
                        && !visited.contains(&w.id())
                })
                .for_each(|(w, n)| { to_visit.push((w.id(), n)); });
            while let Some((id_next, n_next)) = to_visit.pop() {
                group.push(id_next);
                ph_sum += n_next.phase().unwrap();
                visited.push(id_next);
                dg.neighbors_inner(id_next).unwrap()
                    .filter(|(w, n)| {
                        w.is_e()
                            && n.is_same_color(n0)
                            && !visited.contains(&w.id())
                    })
                    .for_each(|(w, n)| { to_visit.push((w.id(), n)); });
            }
            if !group.is_empty() {
                let g = std::mem::take(&mut group);
                groups.push((id0, g, ph_sum));
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(FuseAllData { dg, groups })
        }
    }
}

impl<'a> Rule<ZX> for FuseAllData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, groups } = self;
        let mut all_nb: Vec<ZXWire> = Vec::new();
        for (id0, others, ph_sum) in groups.into_iter() {
            // deal with wires first -- we want to use `delete_node` here and
            // manage connections manually, so we need to keep an accurate count
            // of how many wires are disappearing
            //
            // namely, we want to remove all wires within the group, as well as
            // any self-loops on any group nodes including the "base" node, id0
            //
            // *but only for empty wires*
            //
            // every wire removed will be double-counted
            let mut remove_wire_x2: usize = 0;

            // deal with the base node first, converting all wires within the
            // group to self-wires so they can all be removed in one pass
            let nnb0 = dg.wires[id0].as_mut().unwrap();
            for nb in nnb0.iter_mut() {
                match nb {
                    ZXWire::E(id) => {
                        if *id == id0 {
                            remove_wire_x2 += 1;
                        } else if others.contains(id) {
                            remove_wire_x2 += 1;
                            *id = id0;
                        }
                    },
                    ZXWire::H(id) => {
                        if others.contains(id) { *id = id0; }
                    },
                }
            }
            // remove wires
            while let Some(k) =
                nnb0.iter().enumerate()
                .find_map(|(j, nb)| nb.has_id(id0).then_some(j))
            {
                nnb0.swap_remove(k);
            }

            // now deal with other nodes -- removal of their wires is taken care
            // of by `delete_node`, so we just have to count them and only keep
            // track of neighbors outside the group
            for id1 in others.iter() {
                for nb in dg.delete_node(*id1).unwrap().1.into_iter() {
                    if nb.is_e()
                        && (nb.has_id(id0) || others.contains(&nb.id()))
                    {
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
                .for_each(|nb| {
                    if others.iter().any(|id| nb.has_id(*id)) {
                        nb.map_id(|_| id0);
                    }
                });
            assert!(remove_wire_x2 % 2 == 0);
            dg.wire_count -= remove_wire_x2 / 2;

            dg.nodes[id0].as_mut().unwrap()
                .map_phase(|ph0| ph0 + ph_sum);
        }
    }
}

impl RuleFinder<ZH> for FuseAll {
    type Output<'a> = FuseAllData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        // DFS to find all groups of connected spiders of the same color,
        // starting at every non-io node and ignoring ones already seen
        let mut groups: Vec<SpiderGroup> = Vec::new();
        let mut group: Vec<NodeId> = Vec::new();
        let mut ph_sum: Phase;
        let mut visited: Vec<NodeId> = Vec::new();
        let mut to_visit: Vec<(NodeId, &ZHNode)> = Vec::new();
        for (id0, n0) in dg.nodes_inner() {
            if !n0.is_spider() || visited.contains(&id0) { continue; }
            ph_sum = Phase::zero();
            visited.push(id0);
            dg.neighbors_inner(id0).unwrap()
                .filter(|(id, n)| {
                    n.is_same_color(n0) && !visited.contains(id)
                })
                .for_each(|(id, n)| { to_visit.push((*id, n)); });
            while let Some((id_next, n_next)) = to_visit.pop() {
                group.push(id_next);
                ph_sum += n_next.phase().unwrap();
                visited.push(id_next);
                dg.neighbors_inner(id_next).unwrap()
                    .filter(|(id, n)| {
                        n.is_same_color(n0) && !visited.contains(id)
                    })
                    .for_each(|(id, n)| { to_visit.push((*id, n)); });
            }
            if !group.is_empty() {
                let g = std::mem::take(&mut group);
                groups.push((id0, g, ph_sum));
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(FuseAllData { dg, groups })
        }
    }
}

impl<'a> Rule<ZH> for FuseAllData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, groups } = self;
        let mut all_nb: Vec<NodeId> = Vec::new();
        for (id0, others, ph_sum) in groups.into_iter() {
            // deal with wires first -- we want to use `delete_node` here and
            // manage connections manually, so we need to keep an accurate count
            // of how many wires are disappearing
            //
            // namely, we want to remove all wires within the group, as well as
            // any self-loops on any group nodes including the "base" node, id0
            //
            // every such wire will be double-counted
            let mut remove_wire_x2: usize = 0;

            // deal with the base node first, converting all wires within the
            // group to self-wires so they can all be removed in one pass
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

            // now deal with other nodes -- removal of their wires is taken care
            // of by `delete_node`, so we just have to count them and only keep
            // track of neighbors outside the group
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
}

