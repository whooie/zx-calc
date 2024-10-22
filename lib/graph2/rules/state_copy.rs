use super::*;

/// A spider of arbitrary arity and phase connected to a single oppositely
/// colored state/effect with phase equal to an integer multiple of π.
#[derive(Clone, Debug)]
pub struct StateCopy(NodeId, NodeId);
//         StateCopy(state, spider)


