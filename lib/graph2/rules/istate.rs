use super::*;

/// A unary spider of arbitrary color with phase ±π/2 optionally connected to
/// another spider of arbitrary color.
#[derive(Copy, Clone, Debug)]
pub struct IState(NodeId, Option<NodeId>);
//         IState(state, maybe inner spider)


