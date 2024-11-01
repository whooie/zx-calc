use super::*;

/// Replace two connected spiders of the same color with one, adding their
/// phases.
///
/// ![fuse][fuse]
#[embed_doc_image::embed_doc_image("fuse", "assets/rules/Fuse.svg")]
#[derive(Copy, Clone, Debug)]
pub struct Fuse;

/// Output of [`Fuse::find`].
#[derive(Debug)]
pub struct FuseData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) s1: NodeId, // spider 1
    pub(crate) s2: NodeId, // spider 2
}

impl<'a> FuseData<'a> {
    /// Return the node IDs of the two spiders.
    ///
    /// The left ID will still exist after fusing, while the right will be
    /// removed.
    pub fn ids(&self) -> (NodeId, NodeId) { (self.s1, self.s2) }
}

impl RuleSeal for Fuse { }
impl RuleFinder for Fuse {
    type Output<'a> = FuseData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        for (id, node) in dg.nodes_inner() {
            if node.is_z() {
                for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                    if id2 == id { continue; }
                    if node2.is_z() {
                        return Some(FuseData { dg, s1: id, s2: id2 });
                    }
                }
            } else if node.is_x() {
                for (id2, node2) in dg.neighbors_of_inner(id).unwrap() {
                    if id2 == id { continue; }
                    if node2.is_x() {
                        return Some(FuseData { dg, s1: id, s2: id2 });
                    }
                }
            }
        }
        None
    }
}

impl<'a> RuleSeal for FuseData<'a> { }
impl<'a> Rule for FuseData<'a> {
    fn simplify(self) {
        let Self { dg, s1, s2 } = self;
        dg.remove_wires(s1, s2, None).unwrap();
        let (node2, nnb2) = dg.delete_node(s2).unwrap();
        let ph2 = node2.phase().unwrap();
        dg.nodes[s1].as_mut().unwrap()
            .map_phase(|ph1| ph1 + ph2);
        dg.wires.iter_mut()
            .flatten()
            .flat_map(|nnb| nnb.iter_mut())
            .for_each(|nb| { if *nb == s2 { *nb = s1; } });
        let nnb1 = dg.wires[s1].as_mut().unwrap();
        nnb2.into_iter()
            .for_each(|nb| { if nb != s2 { nnb1.push(nb); } });
    }
}

