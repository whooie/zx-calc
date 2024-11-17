use crate::phase::Phase;
use super::*;

/// Flip the signs of all spiders' phases when they're surrounded by oppositely
/// colored binary spiders with phase π.
///
/// This is the comprehensive version of [`PhaseNeg`], searching more
/// efficiently for all such spiders.
///
/// ![phase_neg][phase_neg]
#[embed_doc_image::embed_doc_image("phase_neg", "assets/rules/PhaseNeg.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct PhaseNegAll;

/// A spider surrounded by oppositely colored spiders with phase π.
///
/// The first item gives the node ID of the central spider, and the second is a
/// list of all surrounding spiders' IDs paired with their outer neighbors' IDs.
pub type PhaseNegGroup<W> = (NodeId, Vec<(NodeId, W)>);

/// Output of [`PhaseNegAll::find`].
#[derive(Debug)]
pub struct PhaseNegAllData<'a, A>
where A: DiagramData
{
    pub(crate) dg: &'a mut Diagram<A>,
    pub(crate) groups: Vec<PhaseNegGroup<A::Wire>>,
}

impl<'a, A> PhaseNegAllData<'a, A>
where A: DiagramData
{
    /// Return the number of phase-flippable spiders found.
    pub fn len(&self) -> usize { self.groups.len() }

    /// Return `true` if the number of groups found is zero.
    pub fn is_empty(&self) -> bool { self.groups.is_empty() }

    /// Return a reference to all groups.
    ///
    /// All central spiders and surrounding π-phase spiders are guaranteed to be
    /// disconnected from other groups' spiders. Outer neighbors may overlap,
    /// however.
    pub fn groups(&self) -> &Vec<PhaseNegGroup<A::Wire>> { &self.groups }
}

impl RuleFinder<ZX> for PhaseNegAll {
    type Output<'a> = PhaseNegAllData<'a, ZX>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        let pi = Phase::pi();
        let mut seen: Vec<NodeId> = Vec::new();
        let mut groups: Vec<PhaseNegGroup<ZXWire>> = Vec::new();
        let mut ss_pi: Vec<(NodeId, ZXWire)> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if seen.contains(&id) {
                continue;
            } else {
                seen.push(id);
            }
            for (wire, node2) in dg.neighbors(id).unwrap() {
                if wire.is_e()
                    && !seen.contains(&wire.id())
                    && node.is_diff_color_and(node2, |_, ph2| ph2 == pi)
                    && dg.arity(wire.id()).unwrap() == 2
                    && dg.mutual_arity_e(id, wire.id()).unwrap() == 1
                {
                    seen.push(wire.id());
                    let nb =
                        dg.find_neighbor_id(wire.id(), |nb| !nb.has_id(id))
                        .unwrap();
                    ss_pi.push((wire.id(), *nb));
                } else {
                    ss_pi.clear();
                    break;
                }
            }
            if !ss_pi.is_empty() {
                let g = std::mem::take(&mut ss_pi);
                groups.push((id, g));
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(PhaseNegAllData { dg, groups })
        }
    }
}

impl<'a> Rule<ZX> for PhaseNegAllData<'a, ZX> {
    fn simplify(self) {
        let Self { dg, groups } = self;
        for (s, ss_pi) in groups.into_iter() {
            for (s_pi, nb) in ss_pi.into_iter() {
                dg.remove_node(s_pi).unwrap();
                match nb {
                    ZXWire::E(id) => { dg.add_wire(s, id).unwrap(); },
                    ZXWire::H(id) => { dg.add_wire_h(s, id).unwrap(); },
                }
            }
            dg.get_node_mut(s).unwrap()
                .map_phase(|ph| -ph);
            let ph = dg.get_node(s).unwrap().phase().unwrap();
            dg.scalar *= ph.cis();
        }
    }
}


impl RuleFinder<ZH> for PhaseNegAll {
    type Output<'a> = PhaseNegAllData<'a, ZH>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        let pi = Phase::pi();
        let mut seen: Vec<NodeId> = Vec::new();
        let mut groups: Vec<PhaseNegGroup<NodeId>> = Vec::new();
        let mut ss_pi: Vec<(NodeId, NodeId)> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if seen.contains(&id) {
                continue;
            } else {
                seen.push(id);
            }
            if node.is_spider() {
                for (id2, node2) in dg.neighbors(id).unwrap() {
                    if !seen.contains(id2)
                        && node.is_diff_color_and(node2, |_, ph2| ph2 == pi)
                        && dg.arity(*id2).unwrap() == 2
                        && dg.mutual_arity(id, *id2).unwrap() == 1
                    {
                        seen.push(*id2);
                        let nb =
                            dg.find_neighbor_id(*id2, |id3| *id3 != id)
                            .unwrap();
                        ss_pi.push((*id2, *nb));
                    } else {
                        ss_pi.clear();
                        break;
                    }
                }
                if !ss_pi.is_empty() {
                    let g = std::mem::take(&mut ss_pi);
                    groups.push((id, g));
                }
            }
        }
        if groups.is_empty() {
            None
        } else {
            Some(PhaseNegAllData { dg, groups })
        }
    }
}

impl<'a> Rule<ZH> for PhaseNegAllData<'a, ZH> {
    fn simplify(self) {
        let Self { dg, groups } = self;
        for (s, ss_pi) in groups.into_iter() {
            for (s_pi, nb) in ss_pi.into_iter() {
                dg.remove_node(s_pi).unwrap();
                dg.add_wire(s, nb).unwrap();
            }
            dg.get_node_mut(s).unwrap()
                .map_phase(|ph| -ph);
            let ph = dg.get_node(s).unwrap().phase().unwrap();
            dg.scalar *= ph.cis();
        }
    }
}

