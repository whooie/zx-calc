use super::*;

/// An arbitrary number `n` of H-boxes, each with arity 1, connected to a
/// phaseless Z-spider of arity `n + 1`.
#[derive(Clone, Debug)]
pub struct HStateMul(Vec<NodeId>, NodeId);
//         HStateMul(hboxes, spider)


