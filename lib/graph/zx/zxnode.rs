use crate::{
    graph::{ NodeData, Spider },
    phase::Phase,
    tensor::{ Element, ElementData },
};

/// A single node in a ZX-diagram and its data.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum ZXNode {
    /// A Z-spider, parameterized by a real phase.
    Z(Phase),
    /// An X-spider, parameterized by a real phase.
    X(Phase),
    /// Termination of a wire as an input to the diagram.
    Input,
    /// Termination of a wire as an output of the diagram.
    Output,
}

impl ZXNode {
    /// Create a new, phaseless Z-spider.
    pub fn z() -> Self { Self::Z(Phase::zero()) }

    /// Create a new Z-spider with π phase.
    pub fn z_pi() -> Self { Self::Z(Phase::pi()) }

    /// Create a new Z-spider with π/2 phase.
    pub fn z_pi2() -> Self { Self::Z(Phase::pi2()) }

    /// Create a new Z-spider with π/4 phase.
    pub fn z_pi4() -> Self { Self::Z(Phase::pi4()) }

    /// Create a new Z-spider with π/8 phase.
    pub fn z_pi8() -> Self { Self::Z(Phase::pi8()) }

    /// Create a new Z-spider with phase `(a / b) × 2π`.
    pub fn z_frac(a: i64, b: i64) -> Self { Self::Z(Phase::new(a, b)) }

    /// Create a new, phaseless X-spider.
    pub fn x() -> Self { Self::X(Phase::zero()) }

    /// Create a new X-spider with π phase.
    pub fn x_pi() -> Self { Self::X(Phase::pi()) }

    /// Create a new X-spider with π/2 phase.
    pub fn x_pi2() -> Self { Self::X(Phase::pi2()) }

    /// Create a new X-spider with π/4 phase.
    pub fn x_pi4() -> Self { Self::X(Phase::pi4()) }

    /// Create a new X-spider with π/8 phase.
    pub fn x_pi8() -> Self { Self::X(Phase::pi8()) }

    /// Cre ate a new X-spider with phase `(a / b) × 2π`.
    pub fn x_frac(a: i64, b: i64) -> Self { Self::X(Phase::new(a, b)) }

    /// Create a new input.
    pub fn input() -> Self { Self::Input }

    /// Create a new output.
    pub fn output() -> Self { Self::Output }

    pub(crate) fn map_phase<F>(&mut self, f: F)
    where F: FnOnce(Phase) -> Phase
    {
        match self {
            Self::Z(ph) => { *ph = f(*ph); },
            Self::X(ph) => { *ph = f(*ph); },
            _ => { },
        }
    }

    pub(crate) fn adjoint_mut(&mut self) {
        match self {
            Self::Z(ph) => { *self = ZXNode::Z(-*ph); },
            Self::X(ph) => { *self = ZXNode::X(-*ph); },
            Self::Input => { *self = ZXNode::Output; },
            Self::Output => { *self = ZXNode::Input; },
        }
    }

    /// Return `true` if `self` is `Z`.
    pub fn is_z(&self) -> bool { matches!(self, Self::Z(_)) }

    /// Return `true` if `self` is `Z` and the phase satisfies some predicate.
    pub fn is_z_and<F>(&self, pred: F) -> bool
    where F: FnOnce(Phase) -> bool
    {
        match self {
            Self::Z(ph) => pred(*ph),
            _ => false,
        }
    }

    /// Return `true` if `self` is `X`.
    pub fn is_x(&self) -> bool { matches!(self, Self::X(_)) }

    /// Return `true` if `self` is `X` and the phase satisfies some predicate.
    pub fn is_x_and<F>(&self, pred: F) -> bool
    where F: FnOnce(Phase) -> bool
    {
        match self {
            Self::X(ph) => pred(*ph),
            _ => false,
        }
    }

    /// Return `true` if `self` is `Input`.
    pub fn is_input(&self) -> bool { matches!(self, Self::Input) }

    /// Return `true` if `self` is `Output`.
    pub fn is_output(&self) -> bool { matches!(self, Self::Output) }

    /// Return `true` if `self` is `Z` or `X`.
    pub fn is_spider(&self) -> bool { matches!(self, Self::Z(_) | Self::X(_)) }

    /// Return `true` if `self` is `Z` or `X` and the phase satisfies some
    /// predicate.
    pub fn is_spider_and<F>(&self, pred: F) -> bool
    where F: FnOnce(Phase) -> bool
    {
        match self {
            Self::Z(ph) => pred(*ph),
            Self::X(ph) => pred(*ph),
            _ => false,
        }
    }

    /// Return `true` if `self` is a generator (`Z`, `X`, or `H`).
    pub fn is_generator(&self) -> bool {
        matches!(self, Self::Z(_) | Self::X(_))
    }

    /// Return `true` if `self` is `Z`, `X`, or `H` and has the default argument
    /// value (`0` for spider phases and `-1` for H-boxes). `Input` and `Output`
    /// return `false`.
    pub fn has_defarg(&self) -> bool {
        match self {
            Self::Z(ph) | Self::X(ph) => *ph == Phase::zero(),
            Self::Input | Self::Output => false,
        }
    }

    /// If `self` is `Z` or `X`, return the associated phase.
    pub fn phase(&self) -> Option<Phase> {
        match self {
            Self::Z(ph) | Self::X(ph) => Some(*ph),
            _ => None,
        }
    }

    /// Return `true` if `self` is a spider with the given phase, modulo 2π.
    pub fn has_phase(&self, phase: Phase) -> bool {
        match self {
            Self::Z(ph) | Self::X(ph) => *ph == phase,
            _ => false,
        }
    }

    /// Return `true` if `self` and `other` are both spiders with the same
    /// color.
    pub fn is_same_color(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (Self::Z(_), Self::Z(_)) | (Self::X(_), Self::X(_))
        )
    }

    /// Return `true` if `self` and `other` are both spiders with the same color
    /// and their phases satisfy a predicate.
    pub fn is_same_color_and<F>(&self, other: &Self, pred: F) -> bool
    where F: FnOnce(Phase, Phase) -> bool
    {
        match (self, other) {
            (Self::Z(ph1), Self::Z(ph2)) | (Self::X(ph1), Self::X(ph2))
                => pred(*ph1, *ph2),
            _ => false,
        }
    }

    /// Return `true` if `self` and `other` are both spiders with different
    /// color.
    pub fn is_diff_color(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (Self::Z(_), Self::X(_)) | (Self::X(_), Self::Z(_))
        )
    }

    /// Return `true` if `self` and `other` are both spiders with different
    /// color and their phases satisfy a predicate.
    pub fn is_diff_color_and<F>(&self, other: &Self, pred: F) -> bool
    where F: FnOnce(Phase, Phase) -> bool
    {
        match (self, other) {
            (Self::Z(ph1), Self::X(ph2)) | (Self::X(ph1), Self::Z(ph2))
                => pred(*ph1, *ph2),
            _ => false,
        }
    }

    pub(crate) fn to_element<A>(self, ins: &[usize], outs: &[usize])
        -> Element<A>
    where A: ElementData
    {
        match self {
            Self::Z(ph) =>
                Element::z(ins.iter().copied(), outs.iter().copied(), Some(ph)),
            Self::X(ph) =>
                Element::x(ins.iter().copied(), outs.iter().copied(), Some(ph)),
            Self::Input | Self::Output => unreachable!()
        }
    }

    pub(crate) fn graph_attrs(&self) -> tabbycat::AttrList {
        use tabbycat::*;
        use tabbycat::attributes::*;
        use crate::vizdefs::*;
        match self {
            ZXNode::Z(ph) => {
                let ph_label = ph.label();
                AttrList::new()
                    .add_pair(label(ph_label))
                    .add_pair(shape(Shape::Circle))
                    .add_pair(height(CIRCLE_HEIGHT))
                    .add_pair(style(Style::Filled))
                    .add_pair(fillcolor(Z_COLOR))
            },
            ZXNode::X(ph) => {
                let ph_label = ph.label();
                AttrList::new()
                    .add_pair(label(ph_label))
                    .add_pair(shape(Shape::Circle))
                    .add_pair(height(CIRCLE_HEIGHT))
                    .add_pair(style(Style::Filled))
                    .add_pair(fillcolor(X_COLOR))
            },
            ZXNode::Input => {
                AttrList::new()
                    .add_pair(label("Input"))
                    .add_pair(shape(Shape::Plaintext))
            },
            ZXNode::Output => {
                AttrList::new()
                    .add_pair(label("Output"))
                    .add_pair(shape(Shape::Plaintext))
            },
        }
    }
}

impl NodeData for ZXNode {
    fn new_id() -> Self { Self::z() }

    fn new_input() -> Self { Self::Input }

    fn new_output() -> Self { Self::Output }

    fn is_input(&self) -> bool { matches!(self, Self::Input) }

    fn is_output(&self) -> bool { matches!(self, Self::Output) }

    fn conjugate(&mut self) { self.adjoint_mut(); }
}

impl From<Spider> for ZXNode {
    fn from(spider: Spider) -> Self {
        match spider {
            Spider::Z(ph) => Self::Z(ph),
            Spider::X(ph) => Self::X(ph),
        }
    }
}

