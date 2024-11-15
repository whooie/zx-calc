use crate::{
    graph::{ NodeData, clifft::TPhase },
    tensor::{ Element, ElementData },
};

/// A single node in a Clifford+*T* diagram and its data.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum CTNode {
    /// A Z-spider, parameterized by a real phase.
    Z(TPhase),
    /// Termination of a wire as an input to the diagram.
    Input,
    /// Termination of a wire as an output of the diagram.
    Output,
}

impl CTNode {
    /// Create a new, phaseless Z-spider.
    pub fn z0() -> Self { Self::Z(TPhase::T0) }

    /// Create a new Z-spider with phase *π*/4.
    pub fn z1() -> Self { Self::Z(TPhase::T1) }

    /// Create a new Z-spider with phase *π*/2.
    pub fn z2() -> Self { Self::Z(TPhase::T2) }

    /// Create a new Z-spider with phase 3*π*/4.
    pub fn z3() -> Self { Self::Z(TPhase::T3) }

    /// Create a new Z-spider with phase *π*.
    pub fn z4() -> Self { Self::Z(TPhase::T4) }

    /// Create a new Z-spider with phase 5*π*/4.
    pub fn z5() -> Self { Self::Z(TPhase::T5) }

    /// Create a new Z-spider with phase 3*π*/2.
    pub fn z6() -> Self { Self::Z(TPhase::T6) }

    /// Create a new Z-spider with phase 7*π*/4.
    pub fn z7() -> Self { Self::Z(TPhase::T7) }

    /// Create a new input.
    pub fn input() -> Self { Self::Input }

    /// Create a new output.
    pub fn output() -> Self { Self::Output }

    pub(crate) fn map_phase<F>(&mut self, f: F)
    where F: FnOnce(TPhase) -> TPhase
    {
        if let Self::Z(ph) = self { *ph = f(*ph); }
    }

    pub(crate) fn adjoint_mut(&mut self) {
        match self {
            Self::Z(ph) => { *self = CTNode::Z(-*ph); },
            Self::Input => { *self = CTNode::Output; },
            Self::Output => { *self = CTNode::Input; },
        }
    }

    /// Return `true` if `self` is `Z`.
    pub fn is_z(&self) -> bool { matches!(self, Self::Z(_)) }

    /// Return `true` if `self` is `Z` and the phase satisfies some predicate.
    pub fn is_z_and<F>(&self, pred: F) -> bool
    where F: FnOnce(TPhase) -> bool
    {
        match self {
            Self::Z(ph) => pred(*ph),
            _ => false,
        }
    }

    /// Return `true` if `self` is `Input`.
    pub fn is_input(&self) -> bool { matches!(self, Self::Input) }

    /// Return `true` if `self` is `Output`.
    pub fn is_output(&self) -> bool { matches!(self, Self::Output) }

    /// Return `true` if `self` is `Z` and has zero phase. `Input` and `Output`
    /// return `false`.
    pub fn has_defarg(&self) -> bool {
        match self {
            Self::Z(ph) => *ph == TPhase::T0,
            Self::Input | Self::Output => false,
        }
    }

    /// If `self` is `Z`, return the associated phase.
    pub fn phase(&self) -> Option<TPhase> {
        match self {
            Self::Z(ph) => Some(*ph),
            _ => None,
        }
    }

    /// Return `true` if `self` is `Z` with the given phase, modulo 2π.
    pub fn has_phase(&self, phase: TPhase) -> bool {
        match self {
            Self::Z(ph) => *ph == phase,
            _ => false,
        }
    }

    /// Return `true` if `self` and `other` are both `Z`.
    pub fn is_zz(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (Self::Z(_), Self::Z(_))
        )
    }

    /// Return `true` if `self` and `other` are both `Z` and their phases
    /// satisfy a predicate.
    pub fn is_zz_and<F>(&self, other: &Self, pred: F) -> bool
    where F: FnOnce(TPhase, TPhase) -> bool
    {
        match (self, other) {
            (Self::Z(ph1), Self::Z(ph2))
                => pred(*ph1, *ph2),
            _ => false,
        }
    }

    pub(crate) fn to_element<A>(self, ins: &[usize], outs: &[usize])
        -> Element<A>
    where A: ElementData
    {
        match self {
            Self::Z(ph) => Element::z(
                ins.iter().copied(), outs.iter().copied(), Some(ph.into())),
            Self::Input | Self::Output => unreachable!()
        }
    }

    pub(crate) fn graph_attrs(&self) -> tabbycat::AttrList {
        use tabbycat::*;
        use tabbycat::attributes::*;
        use crate::vizdefs::*;
        match self {
            Self::Z(ph) => {
                let ph_label = ph.label();
                AttrList::new()
                    .add_pair(label(ph_label))
                    .add_pair(shape(Shape::Circle))
                    .add_pair(height(CIRCLE_HEIGHT))
                    .add_pair(style(Style::Filled))
                    .add_pair(fillcolor(Z_COLOR))
            },
            Self::Input => {
                AttrList::new()
                    .add_pair(label("Input"))
                    .add_pair(shape(Shape::Plaintext))
            },
            Self::Output => {
                AttrList::new()
                    .add_pair(label("Output"))
                    .add_pair(shape(Shape::Plaintext))
            },
        }
    }
}

impl NodeData for CTNode {
    fn new_id() -> Self { Self::z0() }

    fn new_input() -> Self { Self::Input }

    fn new_output() -> Self { Self::Output }

    fn is_input(&self) -> bool { matches!(self, Self::Input) }

    fn is_output(&self) -> bool { matches!(self, Self::Output) }

    fn conjugate(&mut self) {
        match self {
            Self::Z(ph) => { *ph = -*ph; },
            Self::Input => { *self = Self::Output; },
            Self::Output => { *self = Self::Input; },
        }
    }
}

