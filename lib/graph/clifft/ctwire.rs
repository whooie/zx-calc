use crate::graph::{ WireData, NodeId };

/// A wire in a Clifford+*T* diagram connecting to a certain node ID.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum CTWire {
    /// A normal, empty wire.
    E(NodeId),
    /// A Hadamard wire.
    H(NodeId),
    /// A star wire.
    S(NodeId),
}

impl CTWire {
    /// Return `true` if `self` is `E`.
    pub fn is_e(&self) -> bool { matches!(self, Self::E(_)) }

    /// Return `true` if `self` is `E` and the underlying node ID satisfies some
    /// predicate.
    pub fn is_e_and<F>(&self, pred: F) -> bool
    where F: FnOnce(NodeId) -> bool
    {
        match self {
            Self::E(id) => pred(*id),
            _ => false,
        }
    }

    /// Return `true` if `self` is `E` with an underlying node ID equal to `id`.
    pub fn has_e_id(&self, id: NodeId) -> bool {
        match self {
            Self::E(k) => *k == id,
            _ => false,
        }
    }

    pub(crate) fn make_e(&mut self) {
        match self {
            Self::E(_) => { },
            Self::H(id) | Self::S(id) => { *self = Self::E(*id); },
        }
    }

    /// Return `true` if `self` is `H`.
    pub fn is_h(&self) -> bool { matches!(self, Self::H(_)) }

    /// Return `true` if `self` is `H` and the underlying node ID satisfies some
    /// predicate.
    pub fn is_h_and<F>(&self, pred: F) -> bool
    where F: FnOnce(NodeId) -> bool
    {
        match self {
            Self::H(id) => pred(*id),
            _ => false,
        }
    }

    /// Return `true` if `self` is `H` with an underlying node ID equal to `id`.
    pub fn has_h_id(&self, id: NodeId) -> bool {
        match self {
            Self::H(k) => *k == id,
            _ => false,
        }
    }

    pub(crate) fn make_h(&mut self) {
        match self {
            Self::H(_) => { },
            Self::E(id) | Self::S(id) => { *self = Self::H(*id); },
        }
    }

    pub(crate) fn toggle_h(&mut self) {
        match self {
            Self::E(id) => { *self = Self::H(*id); },
            Self::H(id) => { *self = Self::E(*id); },
            Self::S(_) => { },
        }
    }

    /// Return `true` if `self` is `S`.
    pub fn is_s(&self) -> bool { matches!(self, Self::S(_)) }

    /// Return `true` if `self` is `S` and the underlying node ID satisfies some
    /// predicate.
    pub fn is_s_and<F>(&self, pred: F) -> bool
    where F: FnOnce(NodeId) -> bool
    {
        match self {
            Self::S(id) => pred(*id),
            _ => false,
        }
    }

    /// Return `true` if `self` is `S` with an underlying node ID equal to `id`.
    pub fn has_s_id(&self, id: NodeId) -> bool {
        match self {
            Self::S(k) => *k == id,
            _ => false,
        }
    }

    pub(crate) fn make_s(&mut self) {
        match self {
            Self::S(_) => { },
            Self::E(id) | Self::H(id) => { *self = Self::S(*id); },
        }
    }

    pub(crate) fn toggle_s(&mut self) {
        match self {
            Self::E(id) => { *self = Self::S(*id); },
            Self::S(id) => { *self = Self::E(*id); },
            Self::H(_) => { },
        }
    }

    /// Return `true` if `self` has an underlying node ID equal to `id`.
    pub fn has_id(&self, id: NodeId) -> bool {
        match self {
            Self::E(k) | Self::H(k) | Self::S(k) => *k == id
        }
    }

    /// Return the underlying node ID if `self` is `E`.
    pub fn e_id(&self) -> Option<NodeId> {
        match self {
            Self::E(id) => Some(*id),
            _ => None,
        }
    }

    /// Return the underlying node ID if `self` is `H`.
    pub fn h_id(&self) -> Option<NodeId> {
        match self {
            Self::H(id) => Some(*id),
            _ => None,
        }
    }

    /// Return the underlying node ID if `self` is `S`.
    pub fn s_id(&self) -> Option<NodeId> {
        match self {
            Self::S(id) => Some(*id),
            _ => None,
        }
    }

    /// Return the underlying node ID.
    pub fn id(&self) -> NodeId {
        match self {
            Self::E(id) | Self::H(id) | Self::S(id) => *id
        }
    }

    pub(crate) fn shift_id(&mut self, sh: usize) {
        match self {
            Self::E(id) | Self::H(id) | Self::S(id) => { *id += sh; }
        }
    }
}

impl WireData for CTWire {
    fn id(&self) -> NodeId {
        match self {
            Self::E(id) | Self::H(id) | Self::S(id) => *id,
        }
    }

    fn map_id<F>(&mut self, map: F)
    where F: FnOnce(NodeId) -> NodeId
    {
        match self {
            Self::E(id) | Self::H(id) | Self::S(id) => *id = map(*id),
        }
    }

    fn new_empty(id: NodeId) -> Self { Self::E(id) }

    fn is_empty(&self) -> bool { matches!(self, Self::E(_)) }
}

