use super::*;

/// Two phaseless Z-spiders, each of arity 3, connected by two wires, both with
/// an H-box of identical, arbitrary argument, one with an X-spider of phase π,
/// all of arity 2.
#[derive(Copy, Clone, Debug)]
pub struct HIntro(NodeId, NodeId, NodeId, NodeId, NodeId);
//         HIntro(z-spider1, pi x-spider, hbox1, hbox2, z-spider2)


