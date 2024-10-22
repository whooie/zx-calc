use super::*;

/// An H-box of arbitrary arity connected to any number of X-states/effects,
/// each with phase π.
#[derive(Clone, Debug)]
pub struct HAbsorb(Vec<NodeId>);
//         HAbsorb(pi x-states)


