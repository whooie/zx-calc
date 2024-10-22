use crate::phase::Phase;
use super::*;

/// A central spider of arbitrary phase and arity at least 1, connected to only
/// oppositely colored spiders of arity 2 and phase ±π/2, which are themselves
/// connected to spiders the same color as the central spider with arbitrary
/// phase and arity.
#[derive(Clone, Debug)]
pub struct HEulerColorFlip(NodeId, Vec<(NodeId, Phase, NodeId)>);
//         HEulerColorFlip(center, [(pi/2 spider, phase, neighbor spider)])


