use super::*;

/// An X-spider of arity 2 and phase π connected to a phaseless Z-spider of
/// arity 3 over two wires, each with an H-box of arbitrary argument.
#[derive(Copy, Clone, Debug)]
pub struct HAvg(NodeId, NodeId, NodeId, NodeId);
//         HAvg(pi x-spider, hbox1, hbox2, z-spider)


