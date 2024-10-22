use super::*;

/// An H-box of arbitrary arity and argument connected to any number of
/// phaseless X-states/effects.
#[derive(Clone, Debug)]
pub struct HExplode(Vec<NodeId>, NodeId);
//         HExplode(0pi x-states, hbox)


