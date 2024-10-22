use super::*;

/// Two spiders of arbitrary arity and phase π with the same color sandwiching
/// an oppositely colored spider of arity 2 and arbitrary phase.
#[derive(Copy, Clone, Debug)]
pub struct PiBookend(NodeId, NodeId, NodeId);
//         PiBookend(z/x-spider1, x/z-spider, z/x-spider2)


