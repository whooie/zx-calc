use super::*;

/// Four spiders in a checkerboard square, each with arity 3, such that the
/// phases of the spiders of each color all have idential phases that are an
/// integer multiple of π.
#[derive(Copy, Clone, Debug)]
pub struct BitBialgebra(NodeId, NodeId, NodeId, NodeId);
//         BitBialgebra(z-spider1, z-spider2, x-spider1, x-spider2)


