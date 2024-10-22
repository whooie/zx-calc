use super::*;

/// H-box of arbitrary arity and argument connected to an H-box of aribtrary
/// arity and default argument via a single wire carrying a single H-box of
/// default argument.
#[derive(Copy, Clone, Debug)]
pub struct HFuse(NodeId, NodeId, NodeId);
//         HFuse(arg hbox, mid hbox, hbox)


