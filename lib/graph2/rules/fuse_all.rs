use crate::phase::Phase;
use super::*;

/// All cliques of spiders satisfying the conditions for a
/// [`Fuse`][super::Fuse].
#[derive(Clone, Debug)]
pub struct FuseAll(Vec<Vec<(NodeId, Phase)>>);
//                FuseAll([clique: [spider, phase]])


