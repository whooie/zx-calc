use num_complex::Complex64 as C64;
use num_traits::One;
use crate::{ ketbra, phase::Phase };

/// A single node in a diagram and its data.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Node {
    /// A Z-spider, parameterized by a real number `φ` representing the phase
    /// `π × φ`.
    Z(Phase),
    /// An X-spider, parameterized by a real number `φ` representing the phase
    /// `π × φ`.
    X(Phase),
    /// An H-box, parameterized by a general complex number.
    H(C64),
    /// Termination of a wire as an input to the diagram.
    Input,
    /// Termination of a wire as an output of the diagram.
    Output,
}

impl Node {
    /// Create a new, phaseless Z-spider.
    pub fn z() -> Self { Self::Z(Phase::zero()) }

    /// Create a new Z-spider with phase `(a / b) × 2π`.
    pub fn new_z(a: i64, b: i64) -> Self { Self::Z(Phase::new(a, b)) }

    /// Create a new Z-spider with π phase.
    pub fn z_pi() -> Self { Self::Z(Phase::pi()) }

    /// Create a new Z-spider with π/2 phase.
    pub fn z_pi2() -> Self { Self::Z(Phase::pi2()) }

    /// Create a new Z-spider with π/4 phase.
    pub fn z_pi4() -> Self { Self::Z(Phase::pi4()) }

    /// Create a new Z-spider with π/8 phase.
    pub fn z_pi8() -> Self { Self::Z(Phase::pi8()) }

    /// Create a new Z-spider with phase `2π / n`.
    pub fn z_frac(n: i64) -> Self { Self::Z(Phase::frac(n)) }

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

    /// Cre ate a new X-spider with phase `2π / n`.
    pub fn x_frac(n: i64) -> Self { Self::X(Phase::frac(n)) }

    /// Create a new H-box with default argument, `-1`.
    pub fn h() -> Self { Self::H(-C64::one()) }

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
            Self::Z(ph) => { *self = Node::Z(-*ph); },
            Self::X(ph) => { *self = Node::X(-*ph); },
            Self::H(a) => { *self = Node::H(a.conj()); },
            Self::Input => { *self = Node::Output; },
            Self::Output => { *self = Node::Input; },
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

    /// Return `true` if `self` is `H`.
    pub fn is_h(&self) -> bool { matches!(self, Self::H(_)) }

    /// Return `true` if `self` is `H` and the argument satisfies some
    /// predicate.
    pub fn is_h_and<F>(&self, pred: F) -> bool
    where F: FnOnce(C64) -> bool
    {
        match self {
            Self::H(a) => pred(*a),
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
        matches!(self, Self::Z(_) | Self::X(_) | Self::H(_))
    }

    /// Return `true` if `self` is `Z`, `X`, or `H` and has the default argument
    /// value (`0` for spider phases and `-1` for H-boxes). `Input` and `Output`
    /// return `false`.
    pub fn has_defarg(&self) -> bool {
        match self {
            Self::Z(ph) | Self::X(ph) => *ph == Phase::zero(),
            Self::H(a) => (*a + C64::one()).norm() < 1e-15,
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

    pub(crate) fn as_element<I, J>(&self, ins: I, outs: J) -> ketbra::Element
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        match self {
            Self::Z(ph) => ketbra::Element::Z(ins, outs, Some((*ph).into())),
            Self::X(ph) => ketbra::Element::X(ins, outs, Some((*ph).into())),
            Self::H(a) => ketbra::Element::H(ins, outs, Some(*a)),
            Self::Input => ketbra::Element::zero(),
            Self::Output => ketbra::Element::zero(),
        }
    }

    pub(crate) fn graph_attrs(&self) -> tabbycat::AttrList {
        use tabbycat::*;
        use tabbycat::attributes::*;
        use crate::vizdefs::*;
        match self {
            Node::Z(ph) => {
                let ph_label = ph.label();
                AttrList::new()
                    .add_pair(label(ph_label))
                    .add_pair(shape(Shape::Circle))
                    .add_pair(height(CIRCLE_HEIGHT))
                    .add_pair(style(Style::Filled))
                    .add_pair(fillcolor(Z_COLOR))
            },
            Node::X(ph) => {
                let ph_label = ph.label();
                AttrList::new()
                    .add_pair(label(ph_label))
                    .add_pair(shape(Shape::Circle))
                    .add_pair(height(CIRCLE_HEIGHT))
                    .add_pair(style(Style::Filled))
                    .add_pair(fillcolor(X_COLOR))
            },
            Node::H(a) => {
                let a_label =
                    if self.has_defarg() {
                        "".to_string()
                    } else {
                        format!("{}", a)
                    };
                AttrList::new()
                    .add_pair(label(a_label))
                    .add_pair(shape(Shape::Square))
                    .add_pair(height(SQUARE_HEIGHT))
                    .add_pair(style(Style::Filled))
                    .add_pair(fillcolor(H_COLOR))
            },
            Node::Input => {
                AttrList::new()
                    .add_pair(label("Input"))
                    .add_pair(shape(Shape::Plaintext))
            },
            Node::Output => {
                AttrList::new()
                    .add_pair(label("Output"))
                    .add_pair(shape(Shape::Plaintext))
            },
        }
    }
}

/// A single spider.
#[derive(Copy, Clone, Debug)]
pub enum Spider {
    Z(Phase),
    X(Phase),
}

impl Spider {
    /// A phaseless Z-spider.
    pub fn z() -> Self { Self::Z(Phase::zero()) }

    /// A phaseless X-spider.
    pub fn x() -> Self { Self::X(Phase::zero()) }
}

impl From<Spider> for Node {
    fn from(spider: Spider) -> Self {
        match spider {
            Spider::Z(ph) => Self::Z(ph),
            Spider::X(ph) => Self::X(ph),
        }
    }
}


