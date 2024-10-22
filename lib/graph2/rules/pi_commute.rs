use super::*;

/// Two spiders of arbitrary arity and phase with the same color sandwiching an
/// oppositely colored spider of arity 2 with phase π.
#[derive(Copy, Clone, Debug)]
pub struct PiCommute(NodeId, NodeId, NodeId);
//         PiCommute(z/x-spider1, x/z-spider, z/x-spider2)


