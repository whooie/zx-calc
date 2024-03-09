//! Graph-based tools to create and compute end products of diagrams in the
//! ZX(H)-calculus.
//!
//! Diagrams are represented as an undirected, unweighted graph with data
//! attached to its nodes.

use std::{
    collections::VecDeque,
    f64::consts::{ PI, FRAC_PI_2 as PI2, TAU },
    fs,
    io::Write,
    ops::{ Deref, DerefMut },
    path::Path,
};
use num_complex::Complex64 as C64;
use num_rational::Rational64 as R64;
use rustc_hash::{ FxHashMap as HashMap, FxHashSet as HashSet };
use thiserror::Error;
use crate::ketbra;

#[derive(Debug, Error)]
pub enum GraphError {
    #[error("error adding wire: missing node {0}")]
    AddWireMissingNode(usize),

    #[error("error adding wire: input/output node {0} is already connected")]
    AddWireConnectedInputOutput(usize),

    #[error("error applying state/effect: node {0} is not an input/output")]
    ApplyStateNotInputOutput(usize),

    #[error("error applying bell state/effect: nodes {0}, {1} are not both inputs/outputs")]
    ApplyBellStateNotInputOutput(usize, usize),

    #[error("error applying bell state/effect: input/output node {0} is disconnected")]
    ApplyBellStateDisconnected(usize),

    #[error("error constructing GraphViz representation: {0}")]
    GraphVizError(String),

    #[error("I/O error: {0}")]
    IOError(#[from] std::io::Error),
}
pub type GraphResult<T> = Result<T, GraphError>;

/// A definition for a single node in a diagram.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum NodeDef {
    /// A Z-spider, parameterized by a real phase.
    Z(f64),
    /// An X-spider, parameterized by a real phase.
    X(f64),
    /// An H-box, parameterized by a general complex number.
    H(C64),
    /// Termination of a wire as an input to the diagram.
    Input,
    /// Termination of a wire as an output of the diagram.
    Output,
}

impl NodeDef {
    /// Create a new, phaseless Z-spider.
    pub fn z() -> Self { Self::Z(0.0) }

    /// Create a new, phaseless X-spider.
    pub fn x() -> Self { Self::X(0.0) }

    /// Create a new H-box with default argument.
    pub fn h() -> Self { Self::H((-1.0).into()) }

    /// Create a new input.
    pub fn input() -> Self { Self::Input }

    /// Create a new output.
    pub fn output() -> Self { Self::Output }
}

/// A single node in a diagram.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Node {
    /// A Z-spider, parameterized by a real phase.
    Z(f64),
    /// An X-spider, parameterized by a real phase.
    X(f64),
    /// An H-box, parameterized by a general complex number.
    H(C64),
    /// Termination of a wire as an input to the diagram with qubit index.
    Input(QubitId),
    /// Termination of a wire as an output to the diagram with qubit index.
    Output(QubitId),
}

impl From<Node> for NodeDef {
    fn from(node: Node) -> Self {
        match node {
            Node::Z(ph) => Self::Z(ph),
            Node::X(ph) => Self::X(ph),
            Node::H(a) => Self::H(a),
            Node::Input(_) => Self::Input,
            Node::Output(_) => Self::Output,
        }
    }
}

impl Node {
    fn as_element<I, J>(&self, ins: I, outs: J) -> ketbra::Element
    where
        I: IntoIterator<Item = usize>,
        J: IntoIterator<Item = usize>,
    {
        match self {
            Self::Z(ph) => ketbra::Element::Z(ins, outs, Some(*ph)),
            Self::X(ph) => ketbra::Element::X(ins, outs, Some(*ph)),
            Self::H(a) => ketbra::Element::H(ins, outs, Some(*a)),
            Self::Input(_) => ketbra::Element::zero(),
            Self::Output(_) => ketbra::Element::zero(),
        }
    }

    /// Convert `self` to a [`NodeDef`], discarding qubit IDs from `Input` and
    /// `Output`.
    pub fn as_def(&self) -> NodeDef {
        match *self {
            Self::Z(ph) => NodeDef::Z(ph),
            Self::X(ph) => NodeDef::X(ph),
            Self::H(a) => NodeDef::H(a),
            Self::Input(_) => NodeDef::Input,
            Self::Output(_) => NodeDef::Output,
        }
    }

    /// Convert `self` to a [`NodeDef`], discarding qubit IDs from `Input` and
    /// `Output`, and applying a color swap to spiders.
    pub fn as_def_color_flip(&self) -> NodeDef {
        match *self {
            Self::Z(ph) => NodeDef::X(ph),
            Self::X(ph) => NodeDef::Z(ph),
            Self::H(a) => NodeDef::H(a),
            Self::Input(_) => NodeDef::Input,
            Self::Output(_) => NodeDef::Output,
        }
    }

    /// Return `true` if `self` is `Z`.
    pub fn is_z(&self) -> bool { matches!(self, Self::Z(_)) }

    /// Return `true` if `self` is `Z` and the phase satisfies some predicate.
    pub fn is_z_and<F>(&self, pred: F) -> bool
    where F: Fn(f64) -> bool
    {
        match *self {
            Self::Z(ph) => pred(ph),
            _ => false,
        }
    }

    /// Return `true` if `self` is `X`.
    pub fn is_x(&self) -> bool { matches!(self, Self::X(_)) }

    /// Return `true` if `self` is `X` and the phase satisfies some predicate.
    pub fn is_x_and<F>(&self, pred: F) -> bool
    where F: Fn(f64) -> bool
    {
        match *self {
            Self::X(ph) => pred(ph),
            _ => false,
        }
    }

    /// Return `true` if `self` is `H`.
    pub fn is_h(&self) -> bool { matches!(self, Self::H(_)) }

    /// Return `true` if `self` is `H`, and the argument satisfies some
    /// predicate.
    pub fn is_h_and<F>(&self, pred: F) -> bool
    where F: Fn(C64) -> bool
    {
        match *self {
            Self::H(a) => pred(a),
            _ => false,
        }
    }

    /// Return `true` if `self` is `Input`.
    pub fn is_input(&self) -> bool { matches!(self, Self::Input(_)) }

    /// Return `true` if `self` is `Input` and the qubit ID satisfies some
    /// predicate.
    pub fn is_input_and<F>(&self, pred: F) -> bool
    where F: Fn(QubitId) -> bool
    {
        match *self {
            Self::Input(id) => pred(id),
            _ => false,
        }
    }

    /// Return `true` if `self` is `Output`.
    pub fn is_output(&self) -> bool { matches!(self, Self::Output(_)) }

    /// Return `true` if `self` is `Output` and the qubit ID satisfies some
    /// predicate.
    pub fn is_output_and<F>(&self, pred: F) -> bool
    where F: Fn(QubitId) -> bool
    {
        match *self {
            Self::Output(id) => pred(id),
            _ => false,
        }
    }

    /// Return `true` if `self` is `Z` or `X`.
    pub fn is_spider(&self) -> bool { self.is_z() || self.is_x() }

    /// Return `true` if `self` is a spider whose phase satisfies a predicate.
    pub fn is_spider_and<F>(&self, pred: F) -> bool
    where F: Fn(f64) -> bool
    {
        match *self {
            Self::Z(ph) => pred(ph),
            Self::X(ph) => pred(ph),
            _ => false,
        }
    }

    /// Return `true` if `self` is a generator (`Z`, `X`, or `H`).
    pub fn is_generator(&self) -> bool {
        self.is_z() || self.is_x() || self.is_h()
    }

    /// Return `true` if `self` is `Z`, `X`, or `H` and has the default argument
    /// value. `Input` and `Output` return `false`.
    pub fn has_defarg(&self) -> bool {
        match *self {
            Self::Z(ph) => phase_eq(ph, 0.0),
            Self::X(ph) => phase_eq(ph, 0.0),
            Self::H(a) => (a + 1.0).norm() < EPSILON,
            Self::Input(_) => false,
            Self::Output(_) => false,
        }
    }

    /// If `self` is `Z` or `X`, return the associated phase.
    pub fn phase(&self) -> Option<f64> {
        match *self {
            Self::Z(ph) => Some(ph),
            Self::X(ph) => Some(ph),
            _ => None,
        }
    }

    /// Return `true` if `self` is a spider with the given phase, modulo 2π.
    pub fn has_phase(&self, phase: f64) -> bool {
        match *self {
            Self::Z(ph) => phase_eq(ph, phase),
            Self::X(ph) => phase_eq(ph, phase),
            _ => false,
        }
    }

    /// Return `true` if `self` and `other` are both spiders with the same
    /// color.
    pub fn is_same_color(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (Self::Z(_), Self::Z(_)) | (Self::X(_), Self::X(_)),
        )
    }

    /// Return `true` if `self` and `other` are both spiders with the same color
    /// and their phases satisfy a predicate.
    pub fn is_same_color_and<F>(&self, other: &Self, pred: F) -> bool
    where F: Fn(f64, f64) -> bool
    {
        match (self, other) {
            (Self::Z(ph1), Self::Z(ph2)) => pred(*ph1, *ph2),
            (Self::X(ph1), Self::X(ph2)) => pred(*ph1, *ph2),
            _ => false,
        }
    }

    /// Return `true` if `self` and `other` are both spiders with different
    /// color.
    pub fn is_diff_color(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (Self::Z(_), Self::X(_)) | (Self::X(_), Self::Z(_)),
        )
    }

    /// Return `true` if `self` and `other` are both spiders with different
    /// color and their phases satisfy a predicate.
    pub fn is_diff_color_and<F>(&self, other: &Self, pred: F) -> bool
    where F: Fn(f64, f64) -> bool
    {
        match (self, other) {
            (Self::Z(ph1), Self::X(ph2)) => pred(*ph1, *ph2),
            (Self::X(ph1), Self::Z(ph2)) => pred(*ph1, *ph2),
            _ => false,
        }
    }

    fn graph_attrs(&self) -> tabbycat::AttrList {
        use tabbycat::*;
        use tabbycat::attributes::*;
        use crate::vizdefs::*;
        match self {
            Node::Z(ph) => {
                let ph_label = format_phase(*ph);
                AttrList::new()
                    .add_pair(label(ph_label))
                    .add_pair(shape(Shape::Circle))
                    .add_pair(height(CIRCLE_HEIGHT))
                    .add_pair(style(Style::Filled))
                    .add_pair(fillcolor(Z_COLOR))
            },
            Node::X(ph) => {
                let ph_label = format_phase(*ph);
                AttrList::new()
                    .add_pair(label(ph_label))
                    .add_pair(shape(Shape::Circle))
                    .add_pair(height(CIRCLE_HEIGHT))
                    .add_pair(style(Style::Filled))
                    .add_pair(fillcolor(X_COLOR))
            },
            Node::H(a) => {
                let a_label
                    = if self.has_defarg() {
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
            Node::Input(QubitId(qid)) => {
                AttrList::new()
                    .add_pair(label(format!("In {}", qid)))
                    .add_pair(shape(Shape::Plaintext))
            },
            Node::Output(QubitId(qid)) => {
                AttrList::new()
                    .add_pair(label(format!("Out {}", qid)))
                    .add_pair(shape(Shape::Plaintext))
            },
        }
    }
}

/// A single spider.
#[derive(Copy, Clone, Debug)]
pub enum Spider {
    Z(f64),
    X(f64),
}

impl Spider {
    /// A phaseless Z-spider.
    pub fn z() -> Self { Self::Z(0.0) }

    /// A phaseless X-spider.
    pub fn x() -> Self { Self::X(0.0) }
}

impl From<Spider> for Node {
    fn from(spider: Spider) -> Self {
        match spider {
            Spider::Z(ph) => Self::Z(ph),
            Spider::X(ph) => Self::X(ph),
        }
    }
}

impl From<Spider> for NodeDef {
    fn from(spider: Spider) -> Self {
        match spider {
            Spider::Z(ph) => Self::Z(ph),
            Spider::X(ph) => Self::X(ph),
        }
    }
}

const EPSILON: f64 = 1e-15;

/// Returns `true` if `ph1` and `ph2` are equal up to [`EPSILON`], modulo 2π.
pub(crate) fn phase_eq(ph1: f64, ph2: f64) -> bool {
    (ph1 - ph2).rem_euclid(TAU).abs() < EPSILON
}

/// Return a string representation `ph`, reduced modulo π and with the π divided
/// out.
pub(crate) fn format_phase(ph: f64) -> String {
    let ph = ph % TAU;
    if phase_eq(ph, 0.0) {
        "".to_string()
    } else if ph.rem_euclid(PI) < EPSILON {
        if phase_eq(ph, PI) {
            "π".to_string()
        } else {
            format!("{:.0}π", ph / PI)
        }
    } else if let Some(r) = R64::approximate_float(ph / PI) {
        let numer = *r.numer();
        let denom = *r.denom();
        if denom <= 1000 {
            if numer == 1 {
                format!("π/{}", denom)
            } else if numer == -1 {
                format!("–π/{}", denom)
            } else {
                format!("({})π", r)
            }
        } else {
            format!("{}π", ph / PI)
        }
    } else {
        format!("{}π", ph / PI)
    }
}

macro_rules! isomorphism {
    (
        $docstring:literal,
        $name:ident ($iso_to:ident),
        derive: { $($derive:ident),* $(,)? },
        from: { $($from:ident),* $(,)? } $(,)?
    ) => {
        #[doc = $docstring]
        #[derive($($derive),*)]
        pub struct $name(pub $iso_to);

        impl From<$iso_to> for $name {
            fn from(x: $iso_to) -> Self { Self(x) }
        }

        impl From<$name> for $iso_to {
            fn from(x: $name) -> Self { x.0 }
        }

        impl AsRef<$iso_to> for $name {
            fn as_ref(&self) -> &$iso_to { &self.0 }
        }

        impl AsMut<$iso_to> for $name {
            fn as_mut(&mut self) -> &mut $iso_to { &mut self.0 }
        }

        impl Deref for $name {
            type Target = $iso_to;

            fn deref(&self) -> &Self::Target { &self.0 }
        }

        impl DerefMut for $name {
            fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
        }

        $(
            impl From<$from> for $name {
                fn from(x: $from) -> Self { Self(x.into()) }
            }
        )*
    }
}

macro_rules! copy_isomorphism {
    (
        $name:ident ($iso_to:ident),
        from: { $($from:ident),* $(,)? } $(,)?
    ) => {
        impl From<&$iso_to> for $name {
            fn from(x: &$iso_to) -> Self { Self(*x) }
        }

        impl From<&$name> for $name {
            fn from(x: &$name) -> Self { *x }
        }

        impl From<&$name> for $iso_to {
            fn from(x: &$name) -> Self { x.0 }
        }

        $(
            impl From<&$from> for $name {
                fn from(x: &$from) -> Self { Self((*x).into()) }
            }
        )*
    }
}

isomorphism!(
    "Sugared `usize` representing the ID of a single node.",
    NodeId (usize),
    derive: { Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash },
    from: { u8, u16 },
);
copy_isomorphism!(
    NodeId (usize),
    from: { u8, u16 },
);

isomorphism!(
    "Sugared `usize` representing the ID of a single wire.",
    WireId (usize),
    derive: { Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash },
    from: { u8, u16 },
);
copy_isomorphism!(
    WireId (usize),
    from: { u8, u16 },
);

isomorphism!(
    "Sugared `usize` representing the ID of a single qubit.",
    QubitId (usize),
    derive: { Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash },
    from: { u8, u16 },
);
copy_isomorphism!(
    QubitId (usize),
    from: { u8, u16 },
);

/// Represents the connection of one node to another.
///
/// Note that multiple wires are undirected and, since multiple edges are
/// allowed between nodes, a particular `Wire` value may not be unique in a
/// diagram even though node IDs are guaranteed to be unique.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct Wire(pub NodeId, pub NodeId);

impl<A, B> From<(A, B)> for Wire
where
    A: Into<NodeId>,
    B: Into<NodeId>,
{
    fn from(ab: (A, B)) -> Self {
        let (a, b) = ab;
        Self(a.into(), b.into())
    }
}

impl<A, B> From<Wire> for (A, B)
where
    NodeId: Into<A>,
    NodeId: Into<B>,
{
    fn from(wire: Wire) -> Self {
        let Wire(a, b) = wire;
        (a.into(), b.into())
    }
}

impl Wire {
    pub fn has_node<Id>(&self, id: Id) -> bool
    where Id: Into<NodeId>
    {
        let id = id.into();
        self.0 == id || self.1 == id
    }

    pub fn other_end_of<Id>(&self, start: Id) -> Option<NodeId>
    where Id: Into<NodeId>
    {
        let s = start.into();
        match *self {
            Wire(a, b) if a == s => Some(b),
            Wire(a, b) if b == s => Some(a),
            _ => None,
        }
    }
}

/// Like [`Wire`], but holding the data associated with the wire's endpoints
/// instead of their IDs.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct WireEnds(pub Node, pub Node);

/// Either a phaseless spider with two wires or two adjacent H-boxes of arity 2
/// and default argument.
#[derive(Copy, Clone, Debug)]
enum Identity {
    Spider(NodeId),
//  Spider(spider)
    HBox(NodeId, NodeId),
//  HBox(hbox1, hbox2)
}

/// Two spiders of the same color and arbitrary arity and phase with any number
/// of adjoining wires.
#[derive(Copy, Clone, Debug)]
struct Fuse(NodeId, NodeId);
//     Fuse(spider1, spider2)

/// Two spiders of opposite color and arbitrary arity and phase with exactly two
/// adjoining wires.
#[derive(Copy, Clone, Debug)]
struct Hopf(NodeId, NodeId);
//     Hopf(spider1, spider2)

/// Two spiders of the same color and arbitrary arity and phase with exactly two
/// adjoining wires, each with a single H-box with default argument.
#[derive(Copy, Clone, Debug)]
struct H2Hopf(NodeId, NodeId);
//     H2Hopf(spider1, spider2)

/// A spider of arbitrary arity and phase with a single self-loop containing a
/// single H-box with default argument.
#[derive(Copy, Clone, Debug)]
struct HLoop(NodeId, NodeId);
//     HLoop(spider, hbox)

/// Two spiders of arbitrarity arity and phase with the same color sandwiching
/// an H-box with default argument on a single wire.
#[derive(Copy, Clone, Debug)]
struct HEuler(NodeId, NodeId, NodeId);
//     HEuler(spider1, hbox, spider2)

/// Two spiders of arbitrary arity and phase with opposite colors sandwiching an
/// H-box with default argument on a single wire.
#[derive(Copy, Clone, Debug)]
struct HMove(NodeId, NodeId, NodeId);
//     HMove(spider1, hbox, spider2)

/// A spider of arbitrary arity and phase surrounded by a number of H-boxes of
/// default argument and arity 2.
#[derive(Copy, Clone, Debug)]
struct ColorFlip(NodeId);
//     ColorFlip(spider)

/// Two spiders of arbitrary arity and phase with the same color sandwiching an
/// oppositely colored spider of arity 2 with phase π on a single wire.
#[derive(Clone, Debug)]
struct PiCommute(NodeId, NodeId, NodeId);
//     PiCommute(z/x spider1, x/z spider, z/x spider2)

/// A spider of arbitrary arity and phase surrounded by a number of oppositely
/// colored spiders if arity 2 and phase π.
#[derive(Copy, Clone, Debug)]
struct PhaseNeg(NodeId);
//     PhaseNeg(spider)

/// `n` phaseless Z-spiders with all-to-all connectivity to `m` phaseless
/// X-spiders, where each Z-spider has arity `m + 1` and each X-spider has arity
/// `n + 1`.
#[derive(Clone, Debug)]
struct BiAlgebra(Vec<NodeId>, Vec<NodeId>);
//     BiAlgebra(z spiders,   x spiders)

/// Four spiders in a checkerboard square, each with arity 3, such that the
/// phases of the spiders of each color all have idential phases that are an
/// integer multiple of π.
#[derive(Copy, Clone, Debug)]
struct BitBiAlgebra(NodeId, NodeId, NodeId, NodeId);
//     BitBiAlgebra(z-spider1, z-spider2, x-spider1, x-spider2)

/// A unary spider of arbitrary color with phase ±π/2 optionally connected to
/// another spider of arbitrary color.
#[derive(Copy, Clone, Debug)]
struct IState(NodeId, Option<NodeId>);
//     IState(state, maybe inner spider)

/// A spider of arbitrary arity and phase connected to a single oppositely
/// colored state/effect with phase equal to an integer multiple of π.
#[derive(Clone, Debug)]
struct StateCopy(NodeId, NodeId);
//     StateCopy(state, spider)

/// H-box of arbitrary arity and argument connected to an H-box of aribtrary
/// arity and default argument via a single wire carrying a single H-box of
/// default argument.
#[derive(Copy, Clone, Debug)]
struct HFuse(NodeId, NodeId, NodeId);
//     HFuse(arg hbox, mid hbox, hbox)

/// An H-box of arbitrary arity connected to any number of X-states/effects,
/// each with phase π.
#[derive(Clone, Debug)]
struct HAbsorb(Vec<NodeId>);
//     HAbsorb(pi x-states)

/// An H-box of arbitrary arity and argument connected to a single phaseless
/// X-state/effect.
#[derive(Clone, Debug)]
struct HExplode(NodeId, NodeId);
//     HExplode(state,  hbox)

/// An H-box of arity 1 with argument ±1.
#[derive(Copy, Clone, Debug)]
struct HState(NodeId);
//     HState(hstate)

/// A unary H-box of argument -1 or Z-spider of phase π connected to an H-box of
/// arbitrary arity and argument -1.
#[derive(Copy, Clone, Debug)]
struct HStateCopy(NodeId, NodeId);
//     HStateCopy(state,  hcopier)

/// An H-box of arbitrary arity with argument 1.
#[derive(Copy, Clone, Debug)]
struct HMultiState(NodeId);
//     HMultiState(hbox)

/// A phaseless Z-spider and an H-box of default argument, each of arity 3,
/// connected by exactly two adjoining wires.
#[derive(Copy, Clone, Debug)]
struct HHopf(NodeId, NodeId);
//     HHopf(z-spider, hbox)

/// `n` phaseless Z-spiders with all-to-all connectivity to `m` default-argument
/// H-boxes, where each Z-spider has arity `m + 1` and each H-box has arity `n +
/// 1`.
#[derive(Clone, Debug)]
struct HBiAlgebra(Vec<NodeId>, Vec<NodeId>);

/// An X-spider of arity 2 and phase π connected to a phaseless Z-spider of
/// arity 3 over two wires, each with an H-box of arbitrary argument.
#[derive(Copy, Clone, Debug)]
struct HAvg(NodeId, NodeId, NodeId, NodeId);
//     HAvg(pi x-spider, hbox1, hbox2, z-spider)

/// Two phaseless Z-spiders of arbitrary arities connected to each other by an
/// arbitrary number of wires, each with an H-box of arbitrary argument.
#[derive(Clone, Debug)]
struct HMul(NodeId, Vec<NodeId>, NodeId);
//     HMul(spider1, hboxes, spider2)

/// An arbitrary number `n` of H-boxes, each with arity 1, connected to a
/// phaseless Z-spider of arity `n + 1`.
#[derive(Clone, Debug)]
struct HStateMul(Vec<NodeId>, NodeId);
//     HStateMul(hboxes, spider)

/// Two phaseless Z-spiders, each of arity 3, connected by two wires, both with
/// an H-box of identical, arbitrary argument, one with an X-spider of phase π,
/// all of arity 2.
#[derive(Clone, Debug)]
struct HIntro(NodeId, NodeId, NodeId, NodeId, NodeId);
//     HIntro(z-spider1, pi x-spider, hbox1, hbox2, z-spider2)

type NodePath = [(NodeId, Node)];
type PathNodeSpec = fn(&Diagram, &NodePath, &NodeId, &Node) -> bool;

/// Represents a diagram in the calculus.
///
/// Every node and edge is given a unique index necessary to identify it when
/// adding or removing components of the graph. Note that multiple edges may
/// exist between two nodes.
///
/// All accessors are O(1).
#[derive(Clone, Debug)]
pub struct Diagram {
    nodes: HashMap<NodeId, Node>,
    wires: HashMap<WireId, Wire>,
    wires_from: HashMap<NodeId, HashSet<WireId>>,
    wires_between: HashMap<Wire, HashSet<WireId>>,
    node_id: usize,
    input_counter: usize,
    input_gens: HashMap<NodeId, QubitId>,
    output_counter: usize,
    output_gens: HashMap<NodeId, QubitId>,
    edge_id: usize,
}

impl Default for Diagram {
    fn default() -> Self { Self::new() }
}

impl Diagram {
    /// Create a new, empty diagram.
    pub fn new() -> Self {
        Self {
            nodes: HashMap::default(),
            wires: HashMap::default(),
            wires_from: HashMap::default(),
            wires_between: HashMap::default(),
            node_id: 0,
            input_counter: 0,
            input_gens: HashMap::default(),
            output_counter: 0,
            output_gens: HashMap::default(),
            edge_id: 0,
        }
    }

    /// Create a new, disconnected diagram given a set of nodes.
    ///
    /// All nodes passed to this function are given unique IDs starting from 0,
    /// i.e. if `n` nodes are passed to this function then the first node seen
    /// will have ID `0`, and the last will have ID `n - 1`. Input and output
    /// nodes will be given associated qubit numbers in the order they are seen.
    pub fn from_nodes<I>(nodes: I) -> Self
    where I: IntoIterator<Item = NodeDef>
    {
        let mut dg = Self::new();
        nodes.into_iter()
            .for_each(|def| { dg.add_node(def); });
        dg
    }

    /// Return the number of nodes.
    pub fn count_nodes(&self) -> usize { self.nodes.len() }

    /// Return the number of Z-spiders.
    pub fn count_z(&self) -> usize {
        self.nodes.iter().filter(|(_, node)| node.is_z()).count()
    }

    /// Return the number of X-spiders.
    pub fn count_x(&self) -> usize {
        self.nodes.iter().filter(|(_, node)| node.is_x()).count()
    }

    /// Return the number of inputs.
    pub fn count_inputs(&self) -> usize { self.inputs().count() }

    /// Return the number of outputs.
    pub fn count_outputs(&self) -> usize { self.outputs().count() }

    /// Return the total number of spiders.
    pub fn count_spiders(&self) -> usize {
        self.nodes.iter().filter(|(_, node)| node.is_spider()).count()
    }

    /// Return the number of H-boxes.
    pub fn count_h(&self) -> usize {
        self.nodes.iter().filter(|(_, node)| node.is_h()).count()
    }

    /// Return the number of wires.
    pub fn count_wires(&self) -> usize { self.wires.len() }

    /// Get the node associated with a particular ID if it exists.
    pub fn get_node<Id>(&self, node_id: Id) -> Option<&Node>
    where Id: Into<NodeId>
    {
        self.nodes.get(&node_id.into())
    }

    /// Get the number of wires attached to a node, if it exists.
    pub fn arity<Id>(&self, node_id: Id) -> Option<usize>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        self.nodes.contains_key(&id)
            .then(|| {
                self.wires_from.get(&id)
                    .map(|wires| wires.len())
                    .unwrap_or(0)
            })
    }

    /// Return `true` if a node has neighbors, if it exists.
    pub fn is_connected<Id>(&self, node_id: Id) -> Option<bool>
    where Id: Into<NodeId>
    {
        self.wires_from.get(&node_id.into()).map(|wires| !wires.is_empty())
    }

    /// Get the number of wires connecting two nodes, if they both exist.
    pub fn mutual_arity<A, B>(&self, a: A, b: B) -> Option<usize>
    where
        A: Into<NodeId>,
        B: Into<NodeId>,
    {
        let a = a.into();
        let b = b.into();
        (self.nodes.contains_key(&a) && self.nodes.contains_key(&b))
            .then(|| {
                self.wires_between.get(&Wire(a, b))
                    .map(|wires| wires.len())
                    .unwrap_or(0)
            })
    }

    /// Add a node to the diagram and return its ID.
    pub fn add_node(&mut self, node_def: NodeDef) -> NodeId {
        let id = self.node_id;
        let node
            = match node_def {
                NodeDef::Z(ph) => Node::Z(ph),
                NodeDef::X(ph) => Node::X(ph),
                NodeDef::H(a) => Node::H(a),
                NodeDef::Input => {
                    let node = Node::Input(self.input_counter.into());
                    self.input_counter += 1;
                    node
                },
                NodeDef::Output => {
                    let node = Node::Output(self.output_counter.into());
                    self.output_counter += 1;
                    node
                },
            };
        self.nodes.insert(id.into(), node);
        self.wires_from.insert(id.into(), HashSet::default());
        self.node_id += 1;
        id.into()
    }

    fn remove_io_gen<Id>(&mut self, node_id: Id)
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        self.input_gens.remove(&id);
        self.output_gens.remove(&id);
    }

    /// Remove the node associated with a particular ID and return its data if
    /// it exists.
    ///
    /// This method also removes all wires with an endpoint at the node.
    pub fn remove_node<Id>(&mut self, node_id: Id) -> Option<Node>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        if let Some(wires_from) = self.wires_from(id) {
            let wire_ids: Vec<WireId> = wires_from.map(|(id, _)| id).collect();
            wire_ids.iter().for_each(|id| { self.wires.remove(id); });
            self.wires_from.values_mut()
                .for_each(|wires_from| {
                    wire_ids.iter().for_each(|id| { wires_from.remove(id); });
                });
            self.wires_between.values_mut()
                .for_each(|wires_btw| {
                    wire_ids.iter().for_each(|id| { wires_btw.remove(id); });
                });
            self.remove_io_gen(id);
            self.nodes.remove(&id)
        } else {
            None
        }
    }

    /// Remove all nodes that have no wires going to any other node in the
    /// diagram and return their IDs and data.
    pub fn remove_disconnected(&mut self) -> Vec<(NodeId, Node)> {
        let remove_ids: Vec<NodeId>
            = self.wires_from.iter()
            .filter_map(|(node_id, wires_from)| {
                wires_from.is_empty().then_some(*node_id)
            })
            .collect();
        remove_ids.into_iter()
            .filter_map(|id| self.remove_node(id).map(|node| (id, node)))
            .collect()
    }

    /// Get the wire associated with a particular ID if it exists.
    pub fn get_wire<Id>(&self, edge_id: Id) -> Option<&Wire>
    where Id: Into<WireId>
    {
        self.wires.get(&edge_id.into())
    }

    /// Add a wire between two nodes and return its ID.
    ///
    /// Fails if one or both of the nodes don't exist.
    pub fn add_wire<A, B>(&mut self, a: A, b: B) -> GraphResult<WireId>
    where
        A: Into<NodeId>,
        B: Into<NodeId>,
    {
        let a = a.into();
        let b = b.into();
        self.nodes.get(&a)
            .ok_or(GraphError::AddWireMissingNode(a.0))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(a).unwrap())
                    .then_some(())
                    .ok_or(GraphError::AddWireConnectedInputOutput(a.0))
            })?;
        self.nodes.get(&b)
            .ok_or(GraphError::AddWireMissingNode(b.0))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(b).unwrap())
                    .then_some(())
                    .ok_or(GraphError::AddWireConnectedInputOutput(b.0))
            })?;
        let id = self.edge_id;
        self.wires.insert(id.into(), Wire(a, b));
        if let Some(wires_from_a) = self.wires_from.get_mut(&a) {
            wires_from_a.insert(id.into());
        } else {
            self.wires_from.insert(a, [id.into()].into_iter().collect());
        }
        if let Some(wires_from_b) = self.wires_from.get_mut(&b) {
            wires_from_b.insert(id.into());
        } else {
            self.wires_from.insert(b, [id.into()].into_iter().collect());
        }
        if let Some(wires_btw) = self.wires_between.get_mut(&Wire(a, b)) {
            wires_btw.insert(id.into());
        } else {
            self.wires_between.insert(
                Wire(a, b), [id.into()].into_iter().collect());
        }
        if let Some(wires_btw) = self.wires_between.get_mut(&Wire(b, a)) {
            wires_btw.insert(id.into());
        } else {
            self.wires_between.insert(
                Wire(b, a), [id.into()].into_iter().collect());
        }
        self.edge_id += 1;
        self.remove_io_gen(a);
        self.remove_io_gen(b);
        Ok(id.into())
    }

    /// Add a wire with an attached [`Input`][Node::Input] to a pre-existing
    /// node and return the IDs of both the input and the wire as well as the
    /// index of the input qubit.
    ///
    /// The pre-existing node cannot already be an input or output.
    pub fn add_input_wire<Id>(&mut self, node_id: Id)
        -> GraphResult<(NodeId, WireId, QubitId)>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        self.nodes.get(&id)
            .ok_or(GraphError::AddWireMissingNode(id.0))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(GraphError::AddWireConnectedInputOutput(id.0))
            })?;
        let input_id = self.add_node(NodeDef::Input);
        let wire_id = self.add_wire(id, input_id)?;
        let qubit_id = self.input_counter - 1;
        Ok((input_id, wire_id, qubit_id.into()))
    }

    /// Add a wire with an attached [`Output`][Node::Output] to a pre-existing
    /// node and return the IDs of both the output and the wire as well as the
    /// index of the output qubit.
    ///
    /// The pre-existing node cannot already be an input or output.
    pub fn add_output_wire<Id>(&mut self, node_id: Id)
        -> GraphResult<(NodeId, WireId, QubitId)>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        self.nodes.get(&id)
            .ok_or(GraphError::AddWireMissingNode(id.0))
            .and_then(|node| {
                (node.is_generator() || !self.is_connected(id).unwrap())
                    .then_some(())
                    .ok_or(GraphError::AddWireConnectedInputOutput(id.0))
            })?;
        let output_id = self.add_node(NodeDef::Output);
        let wire_id = self.add_wire(id, output_id)?;
        let qubit_id = self.output_counter - 1;
        Ok((output_id, wire_id, qubit_id.into()))
    }

    /// Replace an existing diagram input or output with a spider of arity 1.
    /// Returns the ID of the new node.
    ///
    /// Fails if the given node ID does not exist or is not an input/output.
    ///
    /// Recall that, in the ZH-calculus, terminal H-boxes with arguments are
    /// equal to X-spiders with a corresponding phase, up to an overall scalar.
    pub fn apply_state<Id>(&mut self, node_id: Id, spider: Spider)
        -> GraphResult<NodeId>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        if self.nodes.get(&id).is_some_and(|n| !n.is_generator()) {
            if let Some((int_id, _)) = self.neighbors_of(id).unwrap().next() {
                let new_id = self.add_node(spider.into());
                self.add_wire(int_id, new_id)?;
                match self.remove_node(id) {
                    Some(Node::Input(qid))
                        => { self.input_gens.insert(new_id, qid); },
                    Some(Node::Output(qid))
                        => { self.output_gens.insert(new_id, qid); },
                    _ => unreachable!(),
                }
                Ok(new_id)
            } else {
                self.remove_node(id);
                let new_id = self.add_node(spider.into());
                Ok(new_id)
            }
        } else {
            Err(GraphError::ApplyStateNotInputOutput(id.0))
        }
    }

    /// Replace two existing diagram inputs or outputs with a Bell state/effect
    /// of optional phase. Returns the ID of the new spider node if a phase is
    /// passed.
    ///
    /// Fails if the given node IDs do not exist, are not input/outputs or have
    /// arity 0.
    ///
    /// *Panics if `a` and `b` are equal.*
    pub fn apply_bell<A, B>(&mut self, a: A, b: B, phase: Option<Spider>)
        -> GraphResult<Option<NodeId>>
    where
        A: Into<NodeId>,
        B: Into<NodeId>,
    {
        let a = a.into();
        let b = b.into();
        if a == b {
            panic!(
                "Diagram::apply_bell: \
                cannot apply both connections to the same node"
            );
        }
        self.nodes.get(&a)
            .and_then(|na| self.nodes.get(&b).map(|nb| (na, nb)))
            .is_some_and(|(na, nb)| {
                (na.is_input() && nb.is_input())
                    || (na.is_output() && nb.is_output())
            })
            .then_some(())
            .ok_or(GraphError::ApplyBellStateNotInputOutput(a.0, b.0))?;
        let int_a: NodeId
            = self.neighbors_of(a)
            .and_then(|mut neighbors| neighbors.next())
            .ok_or(GraphError::ApplyBellStateDisconnected(a.0))?
            .0;
        let id_a: NodeId = self.add_node(NodeDef::Z(0.0));
        self.add_wire(id_a, int_a).ok();
        let int_b: NodeId
            = self.neighbors_of(b)
            .and_then(|mut neighbors| neighbors.next())
            .ok_or(GraphError::ApplyBellStateDisconnected(b.0))?
            .0;
        let id_b: NodeId = self.add_node(NodeDef::Z(0.0));
        self.add_wire(id_b, int_b).ok();
        let maybe_spider_id
            = if let Some(spider) = phase {
                let new_id = self.add_node(spider.into());
                self.add_wire(id_a, new_id)?;
                self.add_wire(id_b, new_id)?;
                Some(new_id)
            } else {
                self.add_wire(id_a, id_b)?;
                None
            };
        match self.remove_node(a) {
            Some(Node::Input(qid)) => {
                self.input_gens.insert(id_a, qid);
            },
            Some(Node::Output(qid)) => {
                self.output_gens.insert(id_a, qid);
            },
            _ => unreachable!(),
        }
        match self.remove_node(b) {
            Some(Node::Input(qid)) => {
                self.input_gens.insert(id_b, qid);
            },
            Some(Node::Output(qid)) => {
                self.output_gens.insert(id_b, qid);
            },
            _ => unreachable!(),
        }
        Ok(maybe_spider_id)
    }

    /// Remove the wire associated with a particular ID and return its endpoints
    /// if it exists.
    ///
    /// Note that this will not remove any nodes that become disconnected as a
    /// result of this operation. To remove an input/output node (which have
    /// enforced arity 1), use [`Self::remove_node`] instead.
    pub fn remove_wire<Id>(&mut self, wire_id: Id) -> Option<Wire>
    where Id: Into<WireId>
    {
        let id = wire_id.into();
        if self.wires.contains_key(&id) {
            self.wires_from.values_mut()
                .for_each(|wires_from| { wires_from.remove(&id); });
            self.wires_between.values_mut()
                .for_each(|wires_between| { wires_between.remove(&id); });
            if let Some(wire) = self.wires.remove(&id) {
                self.remove_io_gen(wire.0);
                self.remove_io_gen(wire.1);
                Some(wire)
            } else {
                None
            }
        } else {
            None
        }
    }

    /// Return an iterator over all nodes in the diagram, visited in an
    /// arbitrary order.
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
    pub fn nodes(&self) -> Nodes<'_> {
        Nodes { iter: self.nodes.iter() }
    }

    /// Return an iterator over all interior (non-input/output) nodes in the
    /// diagram, visited in an arbitrary order.
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
    pub fn nodes_inner(&self) -> NodesInner<'_> {
        NodesInner { iter: self.nodes() }
    }

    /// Return an iterator over all diagram inputs, visited in an arbitrary
    /// order.
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
    pub fn inputs(&self) -> Inputs<'_> {
        Inputs { iter: self.nodes.iter() }
    }

    /// Return an iterator over all diagram outputs, visited in an arbitrary
    /// order.
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
    pub fn outputs(&self) -> Outputs<'_> {
        Outputs { iter: self.nodes.iter() }
    }

    /// Return an iterator over all wires in the diagram, visited in an
    /// arbitrary order.
    ///
    /// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
    pub fn wires(&self) -> Wires<'_> {
        Wires { iter: self.wires.iter() }
    }

    /// Return an iterator over all interior (not connected to an input/output
    /// node) wires in the diagram, visited in an arbitrary order.
    ///
    /// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
    pub fn wires_inner(&self) -> WiresInner<'_> {
        WiresInner { diagram: self, iter: self.wires() }
    }

    /// Return an iterator over all wires in the diagram with their endpoints'
    /// associated data, visited in an arbitrary order.
    ///
    /// The iterator item type is `(`[`WireId`]`, `[`WireEnds`]`)`.
    pub fn wires_ends(&self) -> WiresEnds<'_> {
        WiresEnds { diagram: self, iter: self.wires.iter() }
    }

    /// Return an iterator over all interior (not connected to an input/output
    /// node) wires in the diagram with their associated data, visited in an
    /// arbitrary order.
    ///
    /// The iterator item type is `(`[`WireId`]`, `[`WireEnds`]`)`.
    pub fn wires_ends_inner(&self) -> WiresEndsInner<'_> {
        WiresEndsInner { iter: self.wires_ends() }
    }

    /// Return an iterator over all wires with an endpoint at `node_id`, visited
    /// in an arbitrary order, if the node exists.
    ///
    /// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
    pub fn wires_from<Id>(&self, node_id: Id) -> Option<WiresFrom<'_>>
    where Id: Into<NodeId>
    {
        self.wires_from.get(&node_id.into())
            .map(|wire_ids| {
                WiresFrom { diagram: self, iter: wire_ids.iter() }
            })
    }

    /// Return an iterator over all interior (not connected to an input/output
    /// node) wires with an endpoint at `node_id`, visited in an arbitrary
    /// order, if the node exists and is not an input or output node.
    ///
    /// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
    pub fn wires_from_inner<Id>(&self, node_id: Id)
        -> Option<WiresFromInner<'_>>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        if let Some(node) = self.nodes.get(&id) {
            node.is_generator()
                .then(|| WiresFromInner { iter: self.wires_from(id).unwrap() })
        } else {
            None
        }
    }

    /// Return an iterator over all wires with endpoints at `a` and `b`, visited
    /// in an arbitrary order, if both nodes exist and there is at least one
    /// wire between them.
    ///
    /// The iterator item type is [`WireId`].
    pub fn wires_between<A, B>(&self, a: A, b: B) -> Option<WiresBetween<'_>>
    where
        A: Into<NodeId>,
        B: Into<NodeId>,
    {
        self.wires_between.get(&(a, b).into())
            .and_then(|wire_ids| {
                (!wire_ids.is_empty())
                    .then_some(WiresBetween { iter: wire_ids.iter() })
            })
    }

    /// Return an iterator over all interior (not connected to an input/output
    /// node) wires with endpoints at `a` and `b`, visited in an arbitrary order
    /// if both nodes exist, neither are input or output nodes, and there is at
    /// least one wire between them.
    ///
    /// The iterator item type is [`WireId`].
    pub fn wires_between_inner<A, B>(&self, a: A, b: B)
        -> Option<WiresBetweenInner<'_>>
    where
        A: Into<NodeId>,
        B: Into<NodeId>,
    {
        let a = a.into();
        let b = b.into();
        (
            self.nodes.get(&a).is_some_and(|node| node.is_generator())
            && self.nodes.get(&b).is_some_and(|node| node.is_generator())
        ).then(|| {
            self.wires_between(a, b)
                .map(|iter| WiresBetweenInner { diagram: self, iter })
        })
        .flatten()
    }

    /// Return an iterator over all nodes neighboring the node at `node_id`,
    /// visited in an arbitrary order, if it exists.
    ///
    /// Note that this iterator will contain duplicate elements if there are
    /// multiple wires going to the same node.
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
    pub fn neighbors_of<Id>(&self, node_id: Id) -> Option<Neighbors<'_>>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        self.wires_from.get(&id)
            .map(|wire_ids| {
                Neighbors { diagram: self, source: id, iter: wire_ids.iter() }
            })
    }

    /// Return an iterator over all interior (non-input/output) nodes
    /// neighboring the node at `node_id`, visited in an arbitrary order, if it
    /// exists and is not an input or output node.
    ///
    /// Note that this iterator will contain duplicate elements if there are
    /// multiple wires going to the same node.
    ///
    /// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
    pub fn neighbors_of_inner<Id>(&self, node_id: Id)
        -> Option<NeighborsInner<'_>>
    where Id: Into<NodeId>
    {
        let id = node_id.into();
        if self.nodes.get(&id).is_some_and(|node| node.is_generator()) {
            Some(NeighborsInner { iter: self.neighbors_of(id).unwrap() })
        } else {
            None
        }
    }

    /// Recursively find a finite sequence of nodes in the diagram satisfying a
    /// sequence of predicates, given a locally defined set of nodes to search
    /// over at each point in the path. The initial call to this function can
    /// pass any iterator as the set of nodes to search, but recursive calls
    /// will only use the set of neighbors of the last matching path node. Nodes
    /// with IDs matching any of those already in the path are automatically
    /// skipped. Any matching path can only have as many nodes as the number of
    /// predicates provided. Nodes matching a predicate are pushed onto the
    /// accumulator, which is returned as `Some` if all predicates are matched,
    /// otherwise `None` is returned instead.
    fn find_path<I>(
        &self,
        nodes: I,
        specs: &[PathNodeSpec],
        acc: Vec<(NodeId, Node)>,
    ) -> Option<Vec<(NodeId, Node)>>
    where I: Iterator<Item = (NodeId, Node)>
    {
        if let Some(spec) = specs.first() {
            for (id, node) in nodes {
                if acc.iter().any(|(path_id, _)| *path_id == id) {
                    continue;
                }
                if spec(self, &acc, &id, &node) {
                    let mut new_acc = acc.clone();
                    new_acc.push((id, node));
                    let rec: Option<Vec<(NodeId, Node)>>
                        = self.find_path(
                            self.neighbors_of(id).unwrap(),
                            &specs[1..],
                            new_acc.clone(),
                        );
                    if let Some(path) = rec {
                        return Some(path);
                    } else {
                        continue;
                    }
                }
            }
            None
        } else {
            Some(acc)
        }
    }

    fn find_identity(&self) -> Option<Identity> {
        for (id, node) in self.nodes_inner() {
            if node.has_defarg() && self.arity(id).unwrap() == 2 {
                if let Node::Z(_) = node {
                    return Some(Identity::Spider(id));
                } else if let Node::X(_) = node {
                    return Some(Identity::Spider(id));
                } else if let Node::H(_) = node {
                    for (id2, node2) in self.neighbors_of(id).unwrap() {
                        if
                            node2.is_h()
                            && node2.has_defarg()
                            && self.arity(id2).unwrap() == 2
                        {
                            return Some(Identity::HBox(id, id2));
                        }
                    }
                }
            }
        }
        None
    }

    /// Remove all identities from the diagram.
    ///
    /// The following structures are counted as identities:
    /// - a Z- or X-spider with two wires and zero phase
    /// - adjacent H-boxes, each with two wires and argument -1
    pub fn simplify_identity(&mut self) -> bool {
        if let Some(identity) = self.find_identity() {
            match identity {
                Identity::Spider(s) => {
                    let neighbors: Vec<NodeId>
                        = self.neighbors_of(s).unwrap()
                        .map(|(id, _)| id)
                        .collect();
                    self.remove_node(s);
                    self.add_wire(neighbors[0], neighbors[1]).ok();
                },
                Identity::HBox(h1, h2) => {
                    let neighbors: Vec<NodeId>
                        = self.neighbors_of(h1).unwrap()
                        .chain(self.neighbors_of(h2).unwrap())
                        .filter_map(|(id, _)| {
                            (id != h1 && id != h2).then_some(id)
                        })
                        .collect();
                    self.remove_node(h1);
                    self.remove_node(h2);
                    self.add_wire(neighbors[0], neighbors[1]).ok();
                },
            }
            true
        } else {
            false
        }
    }

    fn find_fuse(&self) -> Option<Fuse> {
        for (id, node) in self.nodes_inner() {
            if node.is_z() {
                for (id2, node2) in self.neighbors_of_inner(id).unwrap() {
                    if id2 == id {
                        continue;
                    }
                    if node2.is_z() {
                        return Some(Fuse(id, id2));
                    }
                }
            } else if node.is_x() {
                for (id2, node2) in self.neighbors_of_inner(id).unwrap() {
                    if id2 == id {
                        continue;
                    }
                    if node2.is_x() {
                        return Some(Fuse(id, id2));
                    }
                }
            }
        }
        None
    }

    /// Fuse adjacent spiders of the same color.
    pub fn simplify_fuse(&mut self) -> bool {
        if let Some(Fuse(s1, s2)) = self.find_fuse() {
            let neighbors2: Vec<NodeId>
                = self.neighbors_of(s2).unwrap()
                .filter_map(|(id, _)| (id != s1).then_some(id))
                .collect();
            let self_loops: usize
                = neighbors2.iter()
                .filter(|id| **id == s2)
                .count();
            let new
                = match (self.nodes[&s1], self.nodes[&s2]) {
                    (Node::Z(ph1), Node::Z(ph2))
                        => Node::Z((ph1 + ph2) % TAU),
                    (Node::X(ph1), Node::X(ph2))
                        => Node::X((ph1 + ph2) % TAU),
                    _ => unreachable!(),
                };
            self.remove_io_gen(s1);
            self.remove_node(s2);
            if let Some(node) = self.nodes.get_mut(&s1) { *node = new; }
            neighbors2.into_iter()
                .for_each(|id| { self.add_wire(s1, id).ok(); });
            (0..self_loops)
                .for_each(|_| { self.add_wire(s1, s1).ok(); });
            true
        } else {
            false
        }
    }

    fn find_hopf(&self) -> Option<Hopf> {
        for (id, node) in self.nodes_inner() {
            if node.is_z() {
                for (id2, node2) in self.neighbors_of_inner(id).unwrap() {
                    if
                        node2.is_x()
                        && self.mutual_arity(id, id2).unwrap() >= 2
                    {
                        return Some(Hopf(id, id2));
                    }
                }
            } else if node.is_x() {
                for (id2, node2) in self.neighbors_of_inner(id).unwrap() {
                    if
                        node2.is_z()
                        && self.mutual_arity(id, id2).unwrap() >= 2
                    {
                        return Some(Hopf(id, id2));
                    }
                }
            }
        }
        None
    }

    /// Apply the Hopf rule, disconnecting two spiders of opposite color
    /// joined by two wires.
    pub fn simplify_hopf(&mut self) -> bool {
        if let Some(Hopf(s1, s2)) = self.find_hopf() {
            self.remove_io_gen(s1);
            self.remove_io_gen(s2);
            let wires: Vec<WireId>
                = self.wires_between(s1, s2).unwrap()
                .collect();
            let n = wires.len();
            wires.into_iter().take(2 * (n / 2))
                .for_each(|id| { self.remove_wire(id); });
            true
        } else {
            false
        }
    }

    fn find_h2hopf(&self) -> Option<H2Hopf> {
        let n0 = |_: &Diagram, _: &NodePath, _: &NodeId, n: &Node| {
            n.is_z() || n.is_x()
        };
        let n1 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_h() && n.has_defarg() && dg.arity(id).unwrap() == 2
        };
        let n2 = |_: &Diagram, p: &NodePath, _: &NodeId, n: &Node| {
            n.is_same_color(&p[0].1)
        };
        let n3 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            n.is_h() && n.has_defarg() && dg.arity(id).unwrap() == 2
                && dg.mutual_arity(id, p[0].0).is_some_and(|deg| deg == 1)
        };
        self.find_path(self.nodes_inner(), &[n0, n1, n2, n3], Vec::new())
            .map(|path| H2Hopf(path[0].0, path[2].0))
    }

    /// Apply a variant of the Hopf rule, disconnecting two spiders of arbitrary
    /// arity and phase joined by two wires, each with a single H-box of default
    /// argument.
    pub fn simplify_h2hopf(&mut self) -> bool {
        if let Some(H2Hopf(s1, s2)) = self.find_h2hopf() {
            self.remove_io_gen(s1);
            self.remove_io_gen(s2);
            let hboxes: Vec<NodeId>
                = self.neighbors_of(s1).unwrap()
                .filter_map(|(id, node)| {
                    (
                        node.is_h()
                        && node.has_defarg()
                        && self.arity(id).unwrap() == 2
                        && self.mutual_arity(id, s2).unwrap() == 1
                    ).then_some(id)
                })
                .collect();
            let n = hboxes.len();
            hboxes.into_iter().take(2 * (n / 2))
                .for_each(|h| { self.remove_node(h); });
            true
        } else {
            false
        }
    }

    fn find_hloop(&self) -> Option<HLoop> {
        for (id, node) in self.nodes_inner() {
            if node.is_z() || node.is_x() {
                for (id2, node2) in self.neighbors_of_inner(id).unwrap() {
                    if
                        node2.is_h()
                        && node2.has_defarg()
                        && self.arity(id2).unwrap() == 2
                        && self.mutual_arity(id, id2).unwrap() == 2
                    {
                        return Some(HLoop(id, id2));
                    }
                }
            }
        }
        None
    }

    /// Simplify a spider connected to an H-box in a single loop by removing the
    /// loop and adding π to the spider's phase.
    pub fn simplify_hloop(&mut self) -> bool {
        if let Some(HLoop(s, h)) = self.find_hloop() {
            let new
                = match self.nodes[&s] {
                    Node::Z(ph) => Node::Z((ph + PI) % TAU),
                    Node::X(ph) => Node::X((ph + PI) % TAU),
                    _ => unreachable!(),
                };
            self.remove_io_gen(s);
            self.remove_node(h);
            if let Some(node) = self.nodes.get_mut(&s) { *node = new; }
            true
        } else {
            false
        }
    }

    fn find_h_euler(&self) -> Option<HEuler> {
        let n0 = |_: &Diagram, _: &NodePath, _: &NodeId, n: &Node| {
            n.is_spider()
        };
        let n1 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_h() && n.has_defarg() && dg.arity(id).unwrap() == 2
        };
        let n2 = |_: &Diagram, p: &NodePath, _: &NodeId, n: &Node| {
            n.is_same_color(&p[0].1)
        };
        self.find_path(self.nodes_inner(), &[n0, n1, n2], Vec::new())
            .map(|path| HEuler(path[0].0, path[1].0, path[2].0))
    }

    /// Simplify two spiders of the same color sandwiching an H-box on a single
    /// wire by rewriting the H-box as its Euler angles and fusing spiders
    /// appropriately.
    pub fn simplify_h_euler(&mut self) -> bool {
        fn unfuse(dg: &mut Diagram, s: NodeId, h: NodeId) -> NodeId {
            let (renew, new_def)
                = match dg.nodes[&s] {
                    Node::Z(ph) => (
                        Node::Z((ph - PI2).rem_euclid(TAU)),
                        NodeDef::Z(PI2),
                    ),
                    Node::X(ph) => (
                        Node::X((ph - PI2).rem_euclid(TAU)),
                        NodeDef::X(PI2),
                    ),
                    _ => unreachable!(),
                };
            if let Some(node) = dg.nodes.get_mut(&s) { *node = renew; }
            let new = dg.add_node(new_def);
            dg.add_wire(s, new).ok();
            dg.add_wire(new, h).ok();
            let prev_h_wire = dg.wires_between(s, h).unwrap().next().unwrap();
            dg.remove_wire(prev_h_wire);
            new
        }

        if let Some(HEuler(s1, h, s2)) = self.find_h_euler() {
            let unfuse1
                = self.arity(s1).unwrap() != 2
                || !phase_eq(self.nodes[&s1].phase().unwrap(), PI2);
            let s1 = if unfuse1 { unfuse(self, s1, h) } else { s1 };
            let unfuse2
                = self.arity(s2).unwrap() != 2
                || !phase_eq(self.nodes[&s2].phase().unwrap(), PI2);
            let s2 = if unfuse2 { unfuse(self, s2, h) } else { s2 };

            let neighbors: Vec<NodeId>
                = self.neighbors_of(s1).unwrap()
                .chain(self.neighbors_of(s2).unwrap())
                .filter_map(|(id, _)| (id != h).then_some(id))
                .collect();
            let new
                = if self.nodes[&s1].is_z() {
                    Node::X(-PI / 2.0)
                } else {
                    Node::Z(-PI / 2.0)
                };
            self.remove_node(s1);
            self.remove_node(s2);
            if let Some(node) = self.nodes.get_mut(&h) { *node = new; }
            self.remove_io_gen(h);
            self.add_wire(neighbors[0], h).ok();
            self.add_wire(neighbors[1], h).ok();
            true
        } else {
            false
        }
    }

    fn find_hmove(&self) -> Option<HMove> {
        let n0 = |_: &Diagram, _: &NodePath, _: &NodeId, n: &Node| {
            n.is_z() || n.is_x()
        };
        let n1 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            n.is_h() && n.has_defarg()
                && dg.mutual_arity(id, p[0].0).unwrap() == 1
        };
        let n2 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            ((p[0].1.is_z() && n.is_x()) || (p[0].1.is_x() && n.is_z()))
                && dg.mutual_arity(id, p[1].0).unwrap() == 1
        };
        self.find_path(self.nodes_inner(), &[n0, n1, n2], Vec::new())
            .map(|path| HMove(path[0].0, path[1].0, path[2].0))
    }

    /// Simplify two spiders of opposite color sandwiching an H-box on a single
    /// wire by moving the H-box through the spider with fewer wires and fusing.
    pub fn simplify_hmove(&mut self) -> bool {
        if let Some(HMove(s1, h, s2)) = self.find_hmove() {
            self.remove_io_gen(s1);
            self.remove_io_gen(s2);
            self.remove_node(h);
            let neighbors1: Vec<NodeId>
                = self.neighbors_of(s1).unwrap()
                .filter_map(|(id, _)| (id != s2).then_some(id))
                .collect();
            let neighbors2: Vec<NodeId>
                = self.neighbors_of(s2).unwrap()
                .filter_map(|(id, _)| (id != s1).then_some(id))
                .collect();
            let direct_wires = self.mutual_arity(s1, s2).unwrap();
            let (s, add_h_wire, add_wire, self_loops)
                : (NodeId, Vec<NodeId>, Vec<NodeId>, usize)
                = if neighbors1.len() <= neighbors2.len() {
                    let maybe_node2 = self.remove_node(s2);
                    match (self.nodes.get_mut(&s1), maybe_node2) {
                        (Some(node), Some(Node::X(ph2))) if node.is_z() => {
                            let Node::Z(ph1) = node else { unreachable!() };
                            *node = Node::Z((*ph1 + ph2) % TAU);
                        },
                        (Some(node), Some(Node::Z(ph2))) if node.is_x() => {
                            let Node::X(ph1) = node else { unreachable!() };
                            *node = Node::X((*ph1 + ph2) % TAU);
                        },
                        _ => unreachable!(),
                    }
                    let self_loops
                        = neighbors2.iter().filter(|id| **id == s2).count();
                    (s1, neighbors1, neighbors2, self_loops)
                } else {
                    let maybe_node1 = self.remove_node(s1);
                    match (self.nodes.get_mut(&s2), maybe_node1) {
                        (Some(node), Some(Node::X(ph1))) if node.is_z() => {
                            let Node::Z(ph2) = node else { unreachable!() };
                            *node = Node::Z((ph1 + *ph2) % TAU);
                        },
                        (Some(node), Some(Node::Z(ph1))) if node.is_x() => {
                            let Node::X(ph2) = node else { unreachable!() };
                            *node = Node::X((ph1 + *ph2) % TAU);
                        },
                        _ => unreachable!(),
                    }
                    let self_loops
                        = neighbors1.iter().filter(|id| **id == s1).count();
                    (s2, neighbors2, neighbors1, self_loops)
                };
            for id in add_h_wire.into_iter() {
                if id == s { continue; }
                let wires: Vec<WireId>
                    = self.wires_between(s, id).unwrap()
                    .collect();
                for wire_id in wires.into_iter() {
                    self.remove_wire(wire_id);
                    let h = self.add_node(NodeDef::h());
                    self.add_wire(s, h).ok();
                    self.add_wire(h, id).ok();
                }
            }
            for _ in 0..direct_wires {
                let h = self.add_node(NodeDef::h());
                self.add_wire(s, h).ok();
                self.add_wire(s, h).ok();
            }
            add_wire.into_iter().for_each(|id| { self.add_wire(s, id).ok(); });
            (0..self_loops).for_each(|_| { self.add_wire(s, s).ok(); });
            true
        } else {
            false
        }
    }

    fn find_color_flip(&self) -> Option<ColorFlip> {
        for (id, node) in self.nodes_inner() {
            if node.is_spider() {
                let is_color_flip
                    = self.neighbors_of(id).unwrap()
                    .all(|(id2, node2)| {
                        id == id2
                        || (
                            node2.is_h()
                                && node2.has_defarg()
                                && self.arity(id2).unwrap() == 2
                                && self.mutual_arity(id, id2).unwrap() == 1
                        )
                    });
                if is_color_flip {
                    return Some(ColorFlip(id));
                }
            }
        }
        None
    }

    fn do_color_flip(&mut self, color_flip: ColorFlip) {
        let ColorFlip(s) = color_flip;
        self.remove_io_gen(s);
        let neighbors: Vec<(NodeId, bool)>
            = self.neighbors_of(s).unwrap()
            .filter_map(|(id, node)| {
                (id != s)
                    .then_some(
                        (id, node.is_h() && self.arity(id).unwrap() == 2)
                    )
            })
            .collect();
        let new
            = match self.nodes[&s] {
                Node::Z(ph) => Node::X(ph),
                Node::X(ph) => Node::Z(ph),
                _ => unreachable!(),
            };
        for (id, is_h) in neighbors.into_iter() {
            if is_h {
                let maybe_h_neighbor
                    = self.neighbors_of(id).unwrap()
                    .find(|(id, _)| *id != s);
                if let Some((h_neighbor, _)) = maybe_h_neighbor {
                    self.remove_node(id);
                    self.add_wire(s, h_neighbor).ok();
                } else {
                    self.remove_node(id);
                    self.add_wire(s, s).ok();
                }
            } else {
                let wires: Vec<WireId>
                    = self.wires_between(s, id).unwrap()
                    .collect();
                for wire_id in wires.into_iter() {
                    self.remove_wire(wire_id);
                    let h = self.add_node(NodeDef::h());
                    self.add_wire(s, h).ok();
                    self.add_wire(h, id).ok();
                }
            }
        }
        if let Some(node) = self.nodes.get_mut(&s) { *node = new; }
    }

    /// Flip the color of the spider connected to the largest number of H-boxes
    /// with default argument and arity 2 if doing so would decrease the local
    /// number of such H-boxes.
    pub fn reduce_h(&mut self) -> bool {
        let spider: Option<NodeId>
            = self.nodes_inner()
            .filter_map(|(id, node)| {
                if node.is_spider() {
                    let score: isize
                        = self.neighbors_of(id).unwrap()
                        .map(|(id2, node2)| {
                            if
                                node2.is_h()
                                && node2.has_defarg()
                                && self.arity(id2).unwrap() == 2
                            {
                                -1
                            } else if id2 == id {
                                0
                            } else {
                                1
                            }
                        })
                        .sum();
                    (score < 0).then_some((id, score))
                } else {
                    None
                }
            })
            .min_by(|(_, score_l), (_, score_r)| score_l.cmp(score_r))
            .map(|(id, _)| id);
        if let Some(id) = spider {
            self.do_color_flip(ColorFlip(id));
            true
        } else {
            false
        }
    }

    /// Simplify a spider surrounded by H-boxes of default argument and arity 2
    /// by flipping the color of the spider and removing all H-boxes.
    pub fn simplify_color_flip(&mut self) -> bool {
        if let Some(color_flip) = self.find_color_flip() {
            self.do_color_flip(color_flip);
            true
        } else {
            false
        }
    }

    fn find_pi_commute(&self) -> Option<PiCommute> {
        let n0 = |_: &Diagram, _: &NodePath, _: &NodeId, n: &Node| {
            n.is_z() || n.is_x()
        };
        let n1 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            p[0].1.is_diff_color_and(n, |_, ph| phase_eq(ph, PI))
                && dg.mutual_arity(id, p[0].0).unwrap() == 1
        };
        let n2 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            p[0].1.is_same_color(n)
                && dg.mutual_arity(id, p[1].0).unwrap() == 1
                && dg.mutual_arity(id, p[0].0).unwrap() == 0
        };
        self.find_path(self.nodes_inner(), &[n0, n1, n2], Vec::new())
            .map(|path| PiCommute(path[0].0, path[1].0, path[2].0))
    }

    /// Simplify by applying the π-commute rule as an all-spider version of
    /// [`Self::simplify_hmove`].
    pub fn simplify_pi_commute(&mut self) -> bool {
        if let Some(PiCommute(s1, pi, s2)) = self.find_pi_commute() {
            self.remove_io_gen(s1);
            self.remove_io_gen(s2);
            let pi_def
                = match self.nodes[&pi] {
                    Node::Z(_) => NodeDef::Z(PI),
                    Node::X(_) => NodeDef::X(PI),
                    _ => unreachable!(),
                };
            self.remove_node(pi);
            let neighbors1: Vec<NodeId>
                = self.neighbors_of(s1).unwrap()
                .map(|(id, _)| id)
                .collect();
            let neighbors2: Vec<NodeId>
                = self.neighbors_of(s2).unwrap()
                .map(|(id, _)| id)
                .collect();
            let (s, add_pi_wire, add_wire, self_loops)
                : (NodeId, Vec<NodeId>, Vec<NodeId>, usize)
                = if neighbors1.len() <= neighbors2.len() {
                    let maybe_node2 = self.remove_node(s2);
                    match (self.nodes.get_mut(&s1), maybe_node2) {
                        (Some(node), Some(Node::Z(ph2))) if node.is_z() => {
                            let Node::Z(ph1) = node else { unreachable!() };
                            *node = Node::Z((ph2 - *ph1) % TAU);
                        },
                        (Some(node), Some(Node::X(ph2))) if node.is_x() => {
                            let Node::X(ph1) = node else { unreachable!() };
                            *node = Node::X((ph2 - *ph1) % TAU);
                        },
                        _ => unreachable!(),
                    }
                    let self_loops
                        = neighbors2.iter().filter(|id| **id == s2).count();
                    (s1, neighbors1, neighbors2, self_loops)
                } else {
                    let maybe_node1 = self.remove_node(s1);
                    match (self.nodes.get_mut(&s2), maybe_node1) {
                        (Some(node), Some(Node::Z(ph1))) if node.is_z() => {
                            let Node::Z(ph2) = node else { unreachable!() };
                            *node = Node::Z((ph1 - *ph2) % TAU);
                        },
                        (Some(node), Some(Node::X(ph1))) if node.is_x() => {
                            let Node::X(ph2) = node else { unreachable!() };
                            *node = Node::X((ph1 - *ph2) % TAU);
                        },
                        _ => unreachable!(),
                    }
                    let self_loops
                        = neighbors1.iter().filter(|id| **id == s1).count();
                    (s2, neighbors2, neighbors1, self_loops)
                };
            for id in add_pi_wire.into_iter() {
                if id == s { continue; }
                let wires: Vec<WireId>
                    = self.wires_between(s, id).unwrap()
                    .collect();
                for wire_id in wires.into_iter() {
                    self.remove_wire(wire_id);
                    let pi = self.add_node(pi_def);
                    self.add_wire(s, pi).ok();
                    self.add_wire(pi, id).ok();
                }
            }
            add_wire.into_iter().for_each(|id| { self.add_wire(s, id).ok(); });
            (0..self_loops).for_each(|_| { self.add_wire(s, s).ok(); });
            true
        } else {
            false
        }
    }

    fn find_phase_neg(&self) -> Option<PhaseNeg> {
        for (id, node) in self.nodes_inner() {
            if node.is_spider() {
                let is_phase_neg
                    = self.neighbors_of(id).unwrap()
                    .all(|(id2, node2)| {
                        id == id2
                        || (
                            node.is_diff_color_and(
                                &node2, |_, ph| phase_eq(ph, PI))
                            && self.arity(id2).unwrap() == 2
                            && self.mutual_arity(id, id2).unwrap() == 1
                        )
                    });
                if is_phase_neg {
                    return Some(PhaseNeg(id));
                }
            }
        }
        None
    }

    fn do_phase_neg(&mut self, phase_neg: PhaseNeg) {
        let PhaseNeg(s) = phase_neg;
        self.remove_io_gen(s);
        let pi_def
            = match self.nodes[&s] {
                Node::Z(_) => NodeDef::X(PI),
                Node::X(_) => NodeDef::Z(PI),
                _ => unreachable!(),
            };
        let neighbors: Vec<(NodeId, bool)>
            = self.neighbors_of(s).unwrap()
            .filter_map(|(id, node)| {
                if id != s {
                    let is_pi
                        = self.nodes[&s].is_diff_color_and(
                            &node, |_, ph| phase_eq(ph, PI))
                        && self.arity(id).unwrap() == 2;
                    Some((id, is_pi))
                } else {
                    None
                }
            })
            .collect();
        let new
            = match self.nodes[&s] {
                Node::Z(ph) => Node::Z((ph + PI) % TAU),
                Node::X(ph) => Node::X((ph + PI) % TAU),
                _ => unreachable!(),
            };
        for (id, is_pi) in neighbors.into_iter() {
            if is_pi {
                let maybe_pi_neighbor
                    = self.neighbors_of(id).unwrap()
                    .find(|(id, _)| *id != s);
                if let Some((pi_neighbor, _)) = maybe_pi_neighbor {
                    self.remove_node(id);
                    self.add_wire(s, pi_neighbor).ok();
                } else {
                    self.remove_node(id);
                    self.add_wire(s, s).ok();
                }
            } else {
                let wires: Vec<WireId>
                    = self.wires_between(s, id).unwrap()
                    .collect();
                for wire_id in wires.into_iter() {
                    self.remove_wire(wire_id);
                    let pi = self.add_node(pi_def);
                    self.add_wire(s, pi).ok();
                    self.add_wire(pi, id).ok();
                }
            }
        }
        if let Some(node) = self.nodes.get_mut(&s) { *node = new; }
    }

    /// Flip the sign of the phase of the spider connected to the largest number
    /// of oppositely colored spiders with phase π and arity 2 if doing so would
    /// decrease the local number of such spiders.
    pub fn reduce_pi(&mut self) -> bool {
        let center: Option<NodeId>
            = self.nodes_inner()
            .filter_map(|(id, node)| {
                if node.is_spider() {
                    let score: isize
                        = self.neighbors_of(id).unwrap()
                        .map(|(id2, node2)| {
                            let can_use
                                = node.is_diff_color_and(
                                    &node2, |_, ph| phase_eq(ph, PI))
                                && self.arity(id2).unwrap() == 2;
                            if can_use {
                                -1
                            } else if id2 == id {
                                0
                            } else {
                                1
                            }
                        })
                        .sum();
                    (score < 0).then_some((id, score))
                } else {
                    None
                }
            })
            .min_by(|(_, score_l), (_, score_r)| score_l.cmp(score_r))
            .map(|(id, _)| id);
        if let Some(id) = center {
            self.do_phase_neg(PhaseNeg(id));
            true
        } else {
            false
        }
    }

    /// Simplify a spider surrounded by oppositely colored spiders with phase π
    /// and arity 2 by flipping the sign of the phase of the central spider and
    /// removing all π-spiders.
    pub fn simplify_phase_neg(&mut self) -> bool {
        if let Some(phase_neg) = self.find_phase_neg() {
            self.do_phase_neg(phase_neg);
            true
        } else {
            false
        }
    }

    fn find_bialgebra(&self) -> Option<BiAlgebra> {
        // the problem is to find a maximal biclique in the diagram (comprising
        // Z- and X-spiders), subject to the constraint that each spider must
        // also have exactly one wire connecting it to any node outside the
        // biclique
        todo!()
    }

    /// Simplify by applying the bialgebra rule.
    #[allow(unused_variables)]
    pub fn simplify_bialgebra(&mut self) -> bool {
        if let Some(bialgebra) = self.find_bialgebra() {
            todo!()
        } else {
            false
        }
    }

    fn find_bit_bialgebra(&self) -> Option<BitBiAlgebra> {
        let n0 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_z_and(|ph| ph.rem_euclid(PI) < EPSILON)
                && dg.arity(id).unwrap() == 3
        };
        let n1 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            n.is_x_and(|ph| ph.rem_euclid(PI) < EPSILON)
                && dg.arity(id).unwrap() == 3
                && dg.mutual_arity(p[0].0, id).unwrap() == 1
        };
        let n2 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            n.is_z_and(|ph| ph.rem_euclid(PI) < EPSILON)
                && dg.arity(id).unwrap() == 3
                && dg.mutual_arity(p[1].0, id).unwrap() == 1
                && n.is_same_color_and(&p[0].1, phase_eq)
        };
        let n3 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            n.is_x_and(|ph| ph.rem_euclid(PI) < EPSILON)
                && dg.arity(id).unwrap() == 3
                && dg.mutual_arity(p[2].0, id).unwrap() == 1
                && dg.mutual_arity(p[0].0, id).unwrap() == 1
                && n.is_same_color_and(&p[1].1, phase_eq)
        };
        self.find_path(self.nodes_inner(), &[n0, n1, n2, n3], Vec::new())
            .map(|p| BitBiAlgebra(p[0].0, p[2].0, p[1].0, p[3].0))
    }

    /// Simplify by applying the 2-2 special case of the biagebra rule where
    /// the spiders have phases equal to integer multiples of π.
    pub fn simplify_bit_bialgebra(&mut self) -> bool {
        if let Some(BitBiAlgebra(z1, z2, x1, x2)) = self.find_bit_bialgebra() {
            let z_neighbors: (NodeId, NodeId)
                = (
                    self.neighbors_of(z1).unwrap()
                        .find(|(id, _)| *id != x1 && *id != x2).unwrap().0,
                    self.neighbors_of(z2).unwrap()
                        .find(|(id, _)| *id != x1 && *id != x2).unwrap().0,
                );
            let x_neighbors: (NodeId, NodeId)
                = (
                    self.neighbors_of(x1).unwrap()
                        .find(|(id, _)| *id != z1 && *id != z2).unwrap().0,
                    self.neighbors_of(x2).unwrap()
                        .find(|(id, _)| *id != z1 && *id != z2).unwrap().0,
                );
            let Some(Node::Z(z_ph))
                = self.remove_node(z1) else { unreachable!() };
            self.remove_node(z2);
            let Some(Node::X(x_ph))
                = self.remove_node(x1) else { unreachable!() };
            self.remove_node(x2);
            let z = self.add_node(NodeDef::Z(z_ph));
            self.add_wire(z, x_neighbors.0).ok();
            self.add_wire(z, x_neighbors.1).ok();
            let x = self.add_node(NodeDef::X(x_ph));
            self.add_wire(x, z_neighbors.0).ok();
            self.add_wire(x, z_neighbors.1).ok();
            self.add_wire(z, x).ok();
            true
        } else {
            false
        }
    }

    fn find_istate(&self) -> Option<IState> {
        for (id, node) in self.nodes_inner() {
            if
                node.is_spider_and(|ph| {
                    phase_eq(ph, PI2) || phase_eq(ph, -PI2) })
                && self.arity(id).unwrap() == 1
            {
                let maybe_spider
                    = self.neighbors_of(id).unwrap()
                    .next()
                    .and_then(|(id2, node2)| {
                        node2.is_spider().then_some(id2)
                    });
                if node.is_x() || maybe_spider.is_some() {
                    return Some(IState(id, maybe_spider));
                }
            }
        }
        None
    }

    /// Simplify by converting a state/effect spider with phase ±π/2 to the
    /// appropriate color and fusing if it is connected to another spider,
    /// otherwise converting it to Z.
    pub fn simplify_istate(&mut self) -> bool {
        if let Some(IState(state, maybe_spider)) = self.find_istate() {
            if let Some(spider_id) = maybe_spider {
                self.remove_io_gen(spider_id);
                let state_node = self.remove_node(state).unwrap();
                match (state_node, self.nodes.get_mut(&spider_id).unwrap()) {
                    (Node::Z(ph1), node) => match node {
                        Node::Z(ph2) => {
                            *node = Node::Z((*ph2 + ph1).rem_euclid(TAU));
                        },
                        Node::X(ph2) => {
                            *node = Node::X((*ph2 - ph1).rem_euclid(TAU));
                        },
                        _ => unreachable!(),
                    },
                    (Node::X(ph1), node) => match node {
                        Node::Z(ph2) => {
                            *node = Node::Z((*ph2 - ph1).rem_euclid(TAU));
                        },
                        Node::X(ph2) => {
                            *node = Node::X((*ph2 + ph1).rem_euclid(TAU));
                        },
                        _ => unreachable!(),
                    },
                    _ => unreachable!(),
                }
            } else if let Some(node) = self.nodes.get_mut(&state) {
                if let Node::X(ph) = node {
                    *node = Node::Z((-*ph).rem_euclid(TAU));
                }
            }
            true
        } else {
            false
        }
    }

    fn find_state_copy(&self) -> Option<StateCopy> {
        let n0 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_spider_and(|ph| ph.rem_euclid(PI) < EPSILON)
                && dg.arity(id).unwrap() == 1
        };
        let n1 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            n.is_diff_color(&p[0].1)
                && dg.arity(id).unwrap() >= 2
        };
        self.find_path(self.nodes_inner(), &[n0, n1], Vec::new())
            .map(|p| StateCopy(p[0].0, p[1].0))
    }

    /// Simplify by copying a state or effect through an oppositely colored
    /// spider.
    pub fn simplify_state_copy(&mut self) -> bool {
        if let Some(StateCopy(state, spider)) = self.find_state_copy() {
            let new_def = self.nodes[&state].as_def();
            let spider_neighbors: Vec<NodeId>
                = self.neighbors_of(spider).unwrap()
                .filter_map(|(id, _)| (id != state).then_some(id))
                .collect();
            self.remove_node(state);
            self.remove_node(spider);
            for sp_neighbor in spider_neighbors.into_iter() {
                let new = self.add_node(new_def);
                self.add_wire(new, sp_neighbor).ok();
            }
            true
        } else {
            false
        }
    }

    fn find_hfuse(&self) -> Option<HFuse> {
        let n0 = |_: &Diagram, _: &NodePath, _: &NodeId, n: &Node| {
            n.is_h()
        };
        let n1 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_h() && n.has_defarg() && dg.arity(id).unwrap() == 2
        };
        let n2 = |_: &Diagram, _: &NodePath, _: &NodeId, n: &Node| {
            n.is_h() && n.has_defarg()
        };
        self.find_path(self.nodes_inner(), &[n0, n1, n2], Vec::new())
            .map(|p| HFuse(p[0].0, p[1].0, p[2].0))
    }

    /// Simplify by applying the fuse rule for H-boxes.
    pub fn simplify_hfuse(&mut self) -> bool {
        if let Some(HFuse(h_arg, h_mid, h_end)) = self.find_hfuse() {
            let end_neighbors: Vec<NodeId>
                = self.neighbors_of(h_end).unwrap()
                .filter_map(|(id, _)| (id != h_mid).then_some(id))
                .collect();
            self.remove_io_gen(h_arg);
            self.remove_node(h_mid);
            self.remove_node(h_end);
            end_neighbors.into_iter()
                .for_each(|id| { self.add_wire(h_arg, id).ok(); });
            true
        } else {
            false
        }
    }

    fn find_habsorb(&self) -> Option<HAbsorb> {
        for (id, node) in self.nodes() {
            if node.is_h() && self.arity(id).unwrap() > 2 {
                let pi_xstates: Vec<NodeId>
                    = self.neighbors_of(id).unwrap()
                    .filter_map(|(id2, node2)| {
                        (
                            node2.is_x_and(|ph| phase_eq(ph, PI))
                            && self.arity(id2).unwrap() == 1
                        ).then_some(id2)
                    })
                    .collect();
                if !pi_xstates.is_empty() {
                    return Some(HAbsorb(pi_xstates));
                }
            }
        }
        None
    }

    /// Simplify by absorbing all X-states/effects of phase π into a neighboring
    /// H-box with arbitrary arity and argument.
    pub fn simplify_habsorb(&mut self) -> bool {
        if let Some(HAbsorb(pi_xstates)) = self.find_habsorb() {
            pi_xstates.into_iter()
                .for_each(|id| { self.remove_node(id); });
            true
        } else {
            false
        }
    }

    fn find_hexplode(&self) -> Option<HExplode> {
        let n0 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_x_and(|ph| phase_eq(ph, 0.0))
                && dg.arity(id).unwrap() == 1
        };
        let n1 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_h() && dg.arity(id).unwrap() > 2
        };
        self.find_path(self.nodes_inner(), &[n0, n1], Vec::new())
            .map(|path| HExplode(path[0].0, path[1].0))
    }

    /// Simplify by exploding a single phaseless X-state/effect through a
    /// neighboring H-box with arbitrary arity and argument.
    pub fn simplify_hexplode(&mut self) -> bool {
        if let Some(HExplode(state, h)) = self.find_hexplode() {
            let neighbors: Vec<NodeId>
                = self.neighbors_of(h).unwrap()
                .filter_map(|(id, _)| (id != state).then_some(id))
                .collect();
            self.remove_node(state);
            self.remove_node(h);
            for id in neighbors.into_iter() {
                let z = self.add_node(NodeDef::Z(0.0));
                self.add_wire(z, id).ok();
            }
            true
        } else {
            false
        }
    }

    fn find_hstate(&self) -> Option<HState> {
        self.nodes_inner()
            .find(|(id, node)| {
                node.is_h_and(|a| a == 1.0.into() || a == (-1.0).into())
                    && self.arity(id).unwrap() == 1
            })
            .map(|(id, _)| HState(id))
    }

    /// Simplify by converting an H-box with argument ±1 and one wire to a
    /// state represented by a Z-spider.
    pub fn simplify_hstate(&mut self) -> bool {
        if let Some(HState(h)) = self.find_hstate() {
            let (neighbor, nbr_node)
                = self.neighbors_of(h).unwrap()
                .next().unwrap();
            if nbr_node.is_h_and(|a| a == (-1.0).into()) {
                self.do_hstate_copy(HStateCopy(h, neighbor));
                return true;
            }
            let Some(Node::H(a)) = self.remove_node(h) else { unreachable!() };
            if a == 1.0.into() {
                let z = self.add_node(NodeDef::Z(0.0));
                self.add_wire(z, neighbor).ok();
            } else if a == (-1.0).into() {
                let z = self.add_node(NodeDef::Z(PI));
                self.add_wire(z, neighbor).ok();
            } else {
                unreachable!();
            }
            true
        } else {
            false
        }
    }

    fn find_hstate_copy(&self) -> Option<HStateCopy> {
        let n0 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            (
                n.is_h_and(|a| a == (-1.0).into())
                || n.is_z_and(|ph| phase_eq(ph, PI))
            ) && dg.arity(id).unwrap() == 1
        };
        let n1 = |_: &Diagram, _: &NodePath, _: &NodeId, n: &Node| {
            n.is_h_and(|a| a == (-1.0).into())
        };
        self.find_path(self.nodes_inner(), &[n0, n1], Vec::new())
            .map(|path| HStateCopy(path[0].0, path[1].0))
    }

    fn do_hstate_copy(&mut self, hstate_copy: HStateCopy) {
        let HStateCopy(state, hcopier) = hstate_copy;
        let neighbors: Vec<NodeId>
            = self.neighbors_of(hcopier).unwrap()
            .filter_map(|(id, _)| (id != state).then_some(id))
            .collect();
        self.remove_node(state);
        self.remove_node(hcopier);
        neighbors.into_iter()
            .for_each(|id| {
                let new = self.add_node(NodeDef::X(PI));
                self.add_wire(new, id).ok();
            });
    }

    /// Simplify by applying the all-H-box version of the state copy rule.
    pub fn simplify_hstate_copy(&mut self) -> bool {
        if let Some(hstate_copy) = self.find_hstate_copy() {
            self.do_hstate_copy(hstate_copy);
            true
        } else {
            false
        }
    }

    fn find_hmultistate(&self) -> Option<HMultiState> {
        self.nodes_inner()
            .find(|(_, node)| node.is_h_and(|a| a == 1.0.into()))
            .map(|(id, _)| HMultiState(id))
    }

    /// Simplify by expanding an H-box of arbitrary arity and argument 1 into a
    /// number of phaseless Z-states/effects equal to the original arity of the
    /// H-box.
    pub fn simplify_hmultistate(&mut self) -> bool {
        if let Some(HMultiState(h)) = self.find_hmultistate() {
            let neighbors: Vec<NodeId>
                = self.neighbors_of(h).unwrap()
                .map(|(id, _)| id)
                .collect();
            self.remove_node(h);
            neighbors.into_iter()
                .for_each(|id| {
                    let new = self.add_node(NodeDef::Z(0.0));
                    self.add_wire(new, id).ok();
                });
            true
        } else {
            false
        }
    }

    fn find_hhopf(&self) -> Option<HHopf> {
        let n0 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_z_and(|ph| phase_eq(ph, 0.0)) && dg.arity(id).unwrap() == 3
        };
        let n1 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            n.is_h_and(|a| a == (-1.0).into())
                && dg.arity(id).unwrap() == 3
                && dg.mutual_arity(id, p[0].0).unwrap() == 2
        };
        self.find_path(self.nodes_inner(), &[n0, n1], Vec::new())
            .map(|path| HHopf(path[0].0, path[1].0))
    }

    /// Simplify by applying the H-box version of the Hopf rule.
    pub fn simplify_hhopf(&mut self) -> bool {
        if let Some(HHopf(z, h)) = self.find_hhopf() {
            let z_neighbor
                = self.neighbors_of(z).unwrap()
                .filter_map(|(id, _)| (id != h).then_some(id))
                .next().unwrap();
            self.remove_node(z);
            self.add_wire(h, z_neighbor).ok();
            true
        } else {
            false
        }
    }

    fn find_hbialgebra(&self) -> Option<HBiAlgebra> {
        todo!()
    }

    /// Simplify by applying the H-box version of the bialgebra rule.
    #[allow(unused_variables)]
    pub fn simplify_hbialgebra(&mut self) -> bool {
        if let Some(hbialgebra) = self.find_hbialgebra() {
            todo!()
        } else {
            false
        }
    }

    fn find_havg(&self) -> Option<HAvg> {
        let n0 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_x_and(|ph| phase_eq(ph, PI))
                && dg.arity(id).unwrap() == 2
        };
        let n1 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_h() && dg.arity(id).unwrap() == 2
        };
        let n2 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_z_and(|ph| phase_eq(ph, 0.0))
                && dg.arity(id).unwrap() == 3
        };
        let n3 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            n.is_h()
                && dg.arity(id).unwrap() == 2
                && dg.mutual_arity(id, p[0].0).unwrap() == 1
        };
        self.find_path(self.nodes_inner(), &[n0, n1, n2, n3], Vec::new())
            .map(|p| HAvg(p[0].0, p[1].0, p[3].0, p[2].0))
    }

    /// Simplify by applying the averaging rule.
    pub fn simplify_havg(&mut self) -> bool {
        if let Some(HAvg(x, h1, h2, z)) = self.find_havg() {
            let neighbor: NodeId
                = self.neighbors_of(z).unwrap()
                .filter_map(|(id, _)| (id != h1 && id != h2).then_some(id))
                .next().unwrap();
            self.remove_node(x);
            self.remove_node(z);
            let Some(Node::H(a1))
                = self.remove_node(h1) else { unreachable!() };
            let Some(Node::H(a2))
                = self.remove_node(h2) else { unreachable!() };
            let new = self.add_node(NodeDef::H((a1 + a2) / 2.0));
            self.add_wire(new, neighbor).ok();
            true
        } else {
            false
        }
    }

    fn find_hmul(&self) -> Option<HMul> {
        todo!()
    }

    /// Simplify by applying the generalized multiplying rule.
    #[allow(unused_variables)]
    pub fn simplify_hmul(&mut self) -> bool {
        if let Some(hmul) = self.find_hmul() {
            todo!()
        } else {
            false
        }
    }

    fn find_hstate_mul(&self) -> Option<HStateMul> {
        for (id, node) in self.nodes_inner() {
            if node.is_z_and(|ph| phase_eq(ph, 0.0)) {
                let hstates: Vec<NodeId>
                    = self.neighbors_of(id).unwrap()
                    .filter_map(|(id2, node2)| {
                        (node2.is_h() && self.arity(id2).unwrap() == 1)
                            .then_some(id2)
                    })
                    .collect();
                if self.arity(id).unwrap() - hstates.len() == 1 {
                    return Some(HStateMul(hstates, id));
                }
            }
        }
        None
    }

    /// Simplify by applying the multiplying rule for H-box states.
    pub fn simplify_hstate_mul(&mut self) -> bool {
        if let Some(HStateMul(hstates, spider)) = self.find_hstate_mul() {
            let neighbor: NodeId
                = self.neighbors_of(spider).unwrap()
                .filter_map(|(id, _)| (!hstates.contains(&id)).then_some(id))
                .next().unwrap();
            let a: C64
                = hstates.into_iter()
                .map(|id| {
                    let Some(Node::H(a))
                        = self.remove_node(id) else { unreachable!() };
                    a
                })
                .product();
            self.remove_node(spider);
            let new = self.add_node(NodeDef::H(a));
            self.add_wire(new, neighbor).ok();
            true
        } else {
            false
        }
    }

    fn find_hintro(&self) -> Option<HIntro> {
        let n0 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_z_and(|ph| phase_eq(ph, 0.0))
                && dg.arity(id).unwrap() == 3
        };
        let n1 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_x_and(|ph| phase_eq(ph, PI))
                && dg.arity(id).unwrap() == 2
        };
        let n2 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_h() && dg.arity(id).unwrap() == 2
        };
        let n3 = |dg: &Diagram, _: &NodePath, id: &NodeId, n: &Node| {
            n.is_z_and(|ph| phase_eq(ph, 0.0))
                && dg.arity(id).unwrap() == 3
        };
        let n4 = |dg: &Diagram, p: &NodePath, id: &NodeId, n: &Node| {
            n == &p[2].1
                && dg.arity(id).unwrap() == 2
                && dg.mutual_arity(id, p[0].0).unwrap() == 1
        };
        self.find_path(self.nodes_inner(), &[n0, n1, n2, n3, n4], Vec::new())
            .map(|p| HIntro(p[0].0, p[1].0, p[2].0, p[4].0, p[3].0))
    }

    /// Simplify by applying the unary (reversed) version of the introduction
    /// rule.
    pub fn simplify_hintro(&mut self) -> bool {
        if let Some(HIntro(_z1, x, h1, h2, z2)) = self.find_hintro() {
            let z2_neighbor: NodeId
                = self.neighbors_of(z2).unwrap()
                .filter_map(|(id, _)| (id != h1 && id != h2).then_some(id))
                .next().unwrap();
            let Some(Node::H(a)) = self.remove_node(h1) else { unreachable!() };
            self.remove_node(h2);
            self.remove_node(x);
            self.remove_node(z2);
            let new = self.add_node(NodeDef::H(a));
            self.add_wire(new, z2_neighbor).ok();
            true
        } else {
            false
        }
    }

    /// Simplify a diagram by applying all simplification rules until no further
    /// simplifications can be made.
    pub fn simplify(&mut self) {
        use std::convert::identity;
        loop {
            if [
                self.simplify_istate(),
                self.simplify_hstate(),
                self.simplify_hmultistate(),
                self.simplify_hstate_copy(),
                self.simplify_state_copy(),
                self.simplify_color_flip(),
                self.simplify_fuse(),
                self.simplify_hfuse(),
                self.simplify_identity(),
                self.simplify_hloop(),
                self.simplify_hstate_mul(),
                self.simplify_habsorb(),
                self.simplify_hexplode(),
            ].into_iter().any(identity)
            { continue; }

            if [
                self.simplify_hopf(),
                self.simplify_h2hopf(),
                self.simplify_hhopf(),
            ].into_iter().any(identity)
            { continue; }

            if [
                self.simplify_h_euler(),
                self.simplify_phase_neg(),
            ].into_iter().any(identity)
            { continue; }

            if [
                // self.simplify_bialgebra(),
                // self.simplify_hbialgebra(),
                self.simplify_bit_bialgebra(),
            ].into_iter().any(identity)
            { continue; }

            if [
                // self.simplify_hmul(),
                self.simplify_havg(),
                self.simplify_hintro(),
            ].into_iter().any(identity)
            { continue; }

            if [
                self.reduce_h(),
                self.reduce_pi(),
            ].into_iter().any(identity)
            { continue; }

            if [
                self.simplify_hmove(),
                self.simplify_pi_commute(),
            ].into_iter().any(identity)
            { continue; }

            break;
        }
    }

    /// Find all nodes that are part of a scalar subgraph.
    ///
    /// A node is part of a scalar subgraph if there is no path from it to any
    /// `Input` or `Output` node. Note that the returned nodes may comprise more
    /// than one scalar.
    pub fn find_scalar_nodes(&self) -> HashMap<NodeId, Node> {
        // a node is part of a scalar subgraph if there is no path from it to
        // any Input or Output node
        //
        // all scalar subgraphs are found by starting with the set of all nodes
        // and removing those seen by BFS explorations starting at each of the
        // Input and Output nodes; everything that's left is part of a scalar
        // subgraph
        //
        // the returned value may contain multiple disconnected scalar subgraphs

        let mut nodes: HashMap<NodeId, Node> = self.nodes.clone();
        if !nodes.iter().any(|(_, node)| !node.is_generator()) {
            return nodes;
        }
        let mut to_visit: VecDeque<NodeId>;
        for (io, _) in self.inputs().chain(self.outputs()) {
            to_visit = vec![io].into();
            while let Some(id) = to_visit.pop_back() {
                for (id2, _) in self.neighbors_of(id).unwrap() {
                    if nodes.contains_key(&id2) {
                        to_visit.push_front(id2);
                    }
                }
                nodes.remove(&id);
            }
        }
        nodes
    }

    fn compute_scalar(&self, nodes: &HashMap<NodeId, Node>) -> C64 {
        // compute the value of the scalar by converting the subgraph found by
        // `find_scalar_nodes` to a ketbra diagram and contracting
        //
        // if `find_scalar_nodes` did its job right, the result of the
        // contraction is guaranteed to be an Element of a single term with no
        // inputs or outputs, the amplitude of which is the scalar
        //
        // conversion to a ketbra diagram is done by iterating over nodes and
        // analyzing input/output wires relative to what's already been seen
        //
        // this subgraph can be deformed arbitrarily, so no need to care about
        // the order of iteration
        //
        // self-wires are dealt with by adding two extra outgoing wires
        // immediately coupled to a spiderless Bell effect

        let mut elements: Vec<ketbra::Element>
            = Vec::with_capacity(nodes.len());
        let mut visited: HashSet<NodeId> = HashSet::default();
            // = HashSet::with_capacity(nodes.len());
        let qcount = self.input_counter.max(self.output_counter);
        let mut bell_wire = self.edge_id + qcount + 1..;
        let mut ins: Vec<usize>;
        let mut outs: Vec<usize>;
        let mut bell_wires: Vec<(usize, usize)>;
        for (&id, node) in nodes.iter() {
            visited.insert(id);

            ins = Vec::new();
            outs = Vec::new();
            bell_wires = Vec::new();
            for (id2, _) in self.neighbors_of(id).unwrap() {
                if id2 == id {
                    let b1 = bell_wire.next().unwrap();
                    let b2 = bell_wire.next().unwrap();
                    outs.push(b1);
                    outs.push(b2);
                    bell_wires.push((b1, b2));
                    continue;
                }
                if visited.contains(&id2) {
                    for wid in self.wires_between(id, id2).unwrap() {
                        ins.push(wid.0);
                    }
                } else {
                    for wid in self.wires_between(id, id2).unwrap() {
                        outs.push(wid.0);
                    }
                }
            }

            elements.push(node.as_element(ins, outs));
            if !bell_wires.is_empty() {
                for (b1, b2) in bell_wires.into_iter() {
                    elements.push(ketbra::Element::cap([b1, b2], None));
                }
            }
        }

        if elements.is_empty() {
            C64::new(1.0, 0.0)
        } else {
            ketbra::Diagram::new(elements)
                .contract().unwrap_or_else(|_| unreachable!())
                .as_scalar().unwrap_or_else(|| unreachable!())
        }
    }

    /// Find all nodes that are part of a scalar subgraph and compute their
    /// total product.
    ///
    /// See also [`find_scalar_nodes`][Self::find_scalar_nodes].
    pub fn get_scalar(&self) -> C64 {
        let nodes = self.find_scalar_nodes();
        self.compute_scalar(&nodes)
    }

    /// Find, compute, and remove all scalars from `self`, returning their total
    /// product.
    ///
    /// See also [`find_scalar_nodes`][Self::find_scalar_nodes] and
    /// [`remove_scalar_nodes`][Self::remove_scalar_nodes].
    pub fn remove_scalars(&mut self) -> C64 {
        let nodes = self.find_scalar_nodes();
        let scalar = self.compute_scalar(&nodes);
        nodes.into_iter().for_each(|(id, _)| { self.remove_node(id); });
        scalar
    }

    /// Like [`remove_scalars`][Self::remove_scalars], but do not actually
    /// compute any scalar values, only remove their nodes.
    ///
    /// See also [`find_scalar_nodes`][Self::find_scalar_nodes].
    pub fn remove_scalar_nodes(&mut self) {
        self.find_scalar_nodes()
            .into_iter()
            .for_each(|(id, _)| { self.remove_node(id); });
    }

    /// Convert `self` to a [`ketbra::Diagram`] representation.
    pub fn as_ketbra(&self) -> ketbra::Diagram {
        // assemble the ketbra diagram by iterating over nodes and analyzing
        // input/output wires relative to what's already been seen
        //
        // have to BFS explore starting from the input nodes in order to ensure
        // that input-adjacent nodes are placed first in the ketbra diagram
        //
        // we also want to have wire numbers in the ketbra diagram line up with
        // qubit indices for convenience, so if a given wire id (normally used
        // as-is for a ketbra wire index) coincides with a possible qubit index
        // (bounded from above by max{input_counter, output_counter}), shift it
        // by the maximum wire id in the (graph) diagram
        //
        // do this in two steps because nodes with paths to inputs/outputs need
        // to be visited in a special order, but everything else (i.e. part of a
        // scalar) doesn't
        //
        // self-wires are dealt with by adding two extra outgoing wires
        // immediately coupled to a spiderless Bell effect

        fn as_ketbra_inner(
            dg: &Diagram,
            visited: &mut HashSet<NodeId>,
            to_visit: &mut VecDeque<(NodeId, Node)>,
            elements: &mut Vec<ketbra::Element>,
        ) {
            let qcount = dg.input_counter.max(dg.output_counter);
            let nonq_wire = |id: usize| -> usize {
                if id < qcount { dg.edge_id + id } else { id }
            };
            let mut bell_wire = dg.edge_id + qcount + 1..;
            let mut ins: Vec<usize>;
            let mut outs: Vec<usize>;
            let mut bell_wires: Vec<(usize, usize)>;
            let mut empty_wires: Vec<(QubitId, QubitId)> = Vec::new();
            while let Some((id, node)) = to_visit.pop_back() {
                visited.insert(id);

                ins = Vec::new();
                outs = Vec::new();
                bell_wires = Vec::new();
                for (id2, node2) in dg.neighbors_of(id).unwrap() {
                    if id2 == id {
                        let b1 = bell_wire.next().unwrap();
                        let b2 = bell_wire.next().unwrap();
                        outs.push(b1);
                        outs.push(b2);
                        bell_wires.push((b1, b2));
                        continue;
                    }
                    if let (Node::Input(qin), Node::Output(qout))
                        = (node, node2)
                    {
                        visited.insert(id2);
                        empty_wires.push((qin, qout));
                        continue;
                    }
                    if
                        !visited.contains(&id2)
                            && !to_visit.contains(&(id2, node2))
                    {
                        to_visit.push_front((id2, node2));
                    }
                    if !node.is_generator() { break; }
                    if
                        node2.is_input()
                            || (!node2.is_output() && visited.contains(&id2))
                    {
                        if let Node::Input(qid) = node2 {
                            ins.push(qid.0);
                        } else {
                            for wid in dg.wires_between(id, id2).unwrap() {
                                ins.push(nonq_wire(wid.0));
                            }
                        }
                    } else if
                        node2.is_output()
                            || (!node2.is_input() && !visited.contains(&id2))
                    {
                        if let Node::Output(qid) = node2 {
                            outs.push(qid.0);
                        } else {
                            for wid in dg.wires_between(id, id2).unwrap() {
                                outs.push(nonq_wire(wid.0));
                            }
                        }
                    }
                }

                if node.is_generator() {
                    elements.push(node.as_element(ins, outs));
                }
                if !bell_wires.is_empty() {
                    for (b1, b2) in bell_wires.into_iter() {
                        elements.push(ketbra::Element::cap([b1, b2], None));
                    }
                }
            }

            if !empty_wires.is_empty() {
                empty_wires.sort_by(|(qin_l, _), (qin_r, _)| qin_l.cmp(qin_r));
                let (idents_in, mut idents_out): (Vec<QubitId>, Vec<QubitId>)
                    = empty_wires.into_iter().unzip();
                let mut swaps = Vec::<ketbra::Element>::new();
                let mut maybe_mismatch: Option<(usize, (&QubitId, &QubitId))>;
                loop {
                    maybe_mismatch
                        = idents_in.iter().zip(idents_out.iter()).enumerate()
                        .find(|(_, (qin, qout))| qin != qout);
                    if let Some((k, (qin, qswap))) = maybe_mismatch {
                        let Some((kswap, _))
                            = idents_out.iter().enumerate()
                            .find(|(_, qout)| qin == *qout)
                            else { unreachable!() };
                        swaps.push(ketbra::Element::swap([qin.0, qswap.0]));
                        idents_out.swap(k, kswap);
                    } else {
                        break;
                    }
                }
                swaps.reverse();
                elements.append(&mut swaps);
            }
        }

        // first step: all non-scalar nodes
        let mut elements: Vec<ketbra::Element>
            = Vec::with_capacity(self.nodes_inner().count());
        let mut visited: HashSet<NodeId> = HashSet::default();
            // = HashSet::with_capacity(self.nodes_inner().count());
        // init with input nodes to make sure they're seen first
        let mut to_visit: VecDeque<(NodeId, Node)> = self.inputs().collect();
        as_ketbra_inner(self, &mut visited, &mut to_visit, &mut elements);

        // second step: all nodes that aren't part of a scalar
        // reset `to_visit` with everything not already visited
        to_visit
            = self.nodes_inner()
            .filter(|(id, _)| !visited.contains(id))
            .collect();
        as_ketbra_inner(self, &mut visited, &mut to_visit, &mut elements);

        ketbra::Diagram::new(elements)
    }

    /// Return an object containing an encoding of `self` in the [dot
    /// language][dot-lang].
    ///
    /// Rendering this object using the default formatter will result in a full
    /// dot string representation of the diagram.
    ///
    /// [dot-lang]: https://en.wikipedia.org/wiki/DOT_(graph_description_language)
    pub fn to_graphviz(&self, name: &str) -> GraphResult<tabbycat::Graph> {
        use tabbycat::*;
        use tabbycat::attributes::*;
        use crate::vizdefs::*;
        // initial declarations
        let mut statements
            = StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rankdir(RankDir::LR)),
            )
            .add_attr(
                AttrType::Node,
                AttrList::new()
                    .add_pair(fontname(FONT))
                    .add_pair(fontsize(FONTSIZE))
                    .add_pair(margin(NODE_MARGIN))
                    ,
            );
        // ensure all inputs are in a subgraph at the same rank
        let mut inputs_subgraph_stmt
            = StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rank(RankType::Source)),
            );
        let mut inputs: Vec<(NodeId, Node)>
            = self.inputs()
            .chain(
                self.nodes_inner()
                .filter(|(id, _)| self.input_gens.contains_key(id))
            )
            .collect();
        inputs.sort_by(|(l, _), (r, _)| l.cmp(r));
        let mut prev: Option<usize> = None;
        for (NodeId(id), node) in inputs.into_iter() {
            let index
                = match node {
                    Node::Input(QubitId(index)) => index,
                    _ => self.input_gens[&NodeId(id)].0,
                };
            let attrs
                = if node.is_generator() {
                    node.graph_attrs()
                        .add_pair(xlabel(format!("In {}", index)))
                } else {
                    node.graph_attrs()
                };
            inputs_subgraph_stmt
                = inputs_subgraph_stmt.add_node(id.into(), None, Some(attrs));
            if let Some(prev_id) = prev {
                inputs_subgraph_stmt
                    = inputs_subgraph_stmt.add_edge(
                        Edge::head_node(
                            prev_id.into(), Some(Port::compass(Compass::South)),
                        )
                        .line_to_node(
                            id.into(), Some(Port::compass(Compass::North)),
                        )
                        .add_attrpair(style(Style::Invisible))
                    );
            }
            prev = Some(id);
        }
        statements
            = statements.add_subgraph(SubGraph::cluster(inputs_subgraph_stmt));
        // ensure all outputs are in a subgraph at the same rank
        let mut outputs_subgraph_stmt
            = StmtList::new()
            .add_attr(
                AttrType::Graph,
                AttrList::new().add_pair(rank(RankType::Sink)),
            );
        let mut outputs: Vec<(NodeId, Node)>
            = self.outputs()
            .chain(
                self.nodes_inner()
                .filter(|(id, _)| self.output_gens.contains_key(id))
            )
            .collect();
        outputs.sort_by(|(l, _), (r, _)| l.cmp(r));
        let mut prev: Option<usize> = None;
        for (NodeId(id), node) in outputs.into_iter() {
            let index
                = match node {
                    Node::Output(QubitId(index)) => index,
                    _ => self.output_gens[&NodeId(id)].0,
                };
            let attrs
                = if node.is_generator() {
                    node.graph_attrs()
                        .add_pair(xlabel(format!("Out {}", index)))
                } else {
                    node.graph_attrs()
                };
            outputs_subgraph_stmt
                = outputs_subgraph_stmt.add_node(id.into(), None, Some(attrs));
            if let Some(prev_id) = prev {
                outputs_subgraph_stmt
                    = outputs_subgraph_stmt.add_edge(
                        Edge::head_node(
                            prev_id.into(), Some(Port::compass(Compass::South)),
                        )
                        .line_to_node(
                            id.into(), Some(Port::compass(Compass::North)),
                        )
                        .add_attrpair(style(Style::Invisible))
                    );
            }
            prev = Some(id);
        }
        statements
            = statements.add_subgraph(SubGraph::cluster(outputs_subgraph_stmt));
        // add interior nodes
        for (NodeId(id), node) in self.nodes_inner() {
            if self.input_gens.contains_key(&NodeId(id))
                || self.output_gens.contains_key(&NodeId(id))
            {
                continue;
            }
            let attrs = node.graph_attrs();
            statements = statements.add_node(id.into(), None, Some(attrs));
        }
        // add wires
        for (_, Wire(NodeId(a), NodeId(b))) in self.wires() {
            statements
                = statements.add_edge(
                    Edge::head_node(a.into(), None)
                        .line_to_node(b.into(), None)
                );
        }
        GraphBuilder::default()
            .graph_type(GraphType::Graph)
            .strict(false)
            .id(Identity::quoted(name))
            .stmts(statements)
            .build()
            .map_err(GraphError::GraphVizError)
    }

    /// Like [`to_graphviz`][Self::to_graphviz], but render directly to a string
    /// and write it to `path`.
    pub fn save_graphviz<P>(&self, name: &str, path: P) -> GraphResult<()>
    where P: AsRef<Path>
    {
        let graphviz = self.to_graphviz(name)?;
        fs::OpenOptions::new()
            .write(true)
            .append(false)
            .create(true)
            .truncate(true)
            .open(path)?
            .write_all(format!("{}", graphviz).as_bytes())?;
        Ok(())
    }
}

/// Iterator type over the nodes in a diagram.
///
/// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
pub struct Nodes<'a> {
    iter: std::collections::hash_map::Iter<'a, NodeId, Node>
}

impl<'a> Iterator for Nodes<'a> {
    type Item = (NodeId, Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|(id, node)| (*id, *node))
    }
}

/// Iterator type over the interior (non-input/output) nodes in a diagram.
///
/// The iterator item type is `(`[`Node`]`, `[`Node`]`)`.
pub struct NodesInner<'a> {
    iter: Nodes<'a>
}

impl<'a> Iterator for NodesInner<'a> {
    type Item = (NodeId, Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, node)| node.is_generator())
    }
}

/// Iterator type over the inputs to a diagram.
///
/// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
pub struct Inputs<'a> {
    iter: std::collections::hash_map::Iter<'a, NodeId, Node>
}

impl<'a> Iterator for Inputs<'a> {
    type Item = (NodeId, Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, node)| node.is_input())
            .map(|(id, node)| (*id, *node))
    }
}

/// Iterator type over the outputs of a diagram.
///
/// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
pub struct Outputs<'a> {
    iter: std::collections::hash_map::Iter<'a, NodeId, Node>
}

impl<'a> Iterator for Outputs<'a> {
    type Item = (NodeId, Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, node)| node.is_output())
            .map(|(id, node)| (*id, *node))
    }
}

/// Iterator type over the wires in a diagram.
///
/// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
pub struct Wires<'a> {
    iter: std::collections::hash_map::Iter<'a, WireId, Wire>
}

impl<'a> Iterator for Wires<'a> {
    type Item = (WireId, Wire);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|(id, wire)| (*id, *wire))
    }
}

/// Iterator type over the interior (not connected to an input/output node)
/// wires in a diagram.
///
/// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
pub struct WiresInner<'a> {
    diagram: &'a Diagram,
    iter: Wires<'a>,
}

impl<'a> Iterator for WiresInner<'a> {
    type Item = (WireId, Wire);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, Wire(a, b))| {
            self.diagram.nodes[a].is_generator()
                && self.diagram.nodes[b].is_generator()
        })
    }
}


/// Iterator type over the endpoint data associated with the wires in a diagram.
///
/// The iterator item type is `(`[`WireId`]`, `[`WireEnds`]`)`.
pub struct WiresEnds<'a> {
    diagram: &'a Diagram,
    iter: std::collections::hash_map::Iter<'a, WireId, Wire>,
}

impl<'a> Iterator for WiresEnds<'a> {
    type Item = (WireId, WireEnds);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
            .map(|(id, Wire(a, b))| {
                (*id, WireEnds(self.diagram.nodes[a], self.diagram.nodes[b]))
            })
    }
}

/// Iterator type over the endpoint data associated with the interior (not
/// connected to an input/output node) wires in a diagram.
///
/// The iterator item type is `(`[`WireId`]`, `[`WireEnds`]`)`.
pub struct WiresEndsInner<'a> {
    iter: WiresEnds<'a>
}

impl<'a> Iterator for WiresEndsInner<'a> {
    type Item = (WireId, WireEnds);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, WireEnds(a, b))| {
            a.is_generator() && b.is_generator()
        })
    }
}

/// Iterator type over the wires attached to a particular node in a diagram.
///
/// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
pub struct WiresFrom<'a> {
    diagram: &'a Diagram,
    iter: std::collections::hash_set::Iter<'a, WireId>,
}

impl<'a> Iterator for WiresFrom<'a> {
    type Item = (WireId, Wire);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
            .map(|wire_id| (*wire_id, self.diagram.wires[wire_id]))
    }
}

/// Iterator type over the interior (not connected to an input/output node)
/// wires attached to a particular node in a diagram.
///
/// The iterator item type is `(`[`WireId`]`, `[`Wire`]`)`.
pub struct WiresFromInner<'a> {
    iter: WiresFrom<'a>
}

impl<'a> Iterator for WiresFromInner<'a> {
    type Item = (WireId, Wire);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, Wire(a, b))| {
            self.iter.diagram.nodes[a].is_generator()
                && self.iter.diagram.nodes[b].is_generator()
        })
    }
}

/// Iterator type over the wires attached to two particular nodes in a diagram.
///
/// The iterator item type is [`WireId`].
pub struct WiresBetween<'a> {
    iter: std::collections::hash_set::Iter<'a, WireId>,
}

impl<'a> Iterator for WiresBetween<'a> {
    type Item = WireId;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().copied()
    }
}

/// Iterator type over the interior (not connected to an input/output node)
/// wires attached to two particular nodes in a diagram.
///
/// The iterator item type is [`WireId`].
pub struct WiresBetweenInner<'a> {
    diagram: &'a Diagram,
    iter: WiresBetween<'a>,
}

impl<'a> Iterator for WiresBetweenInner<'a> {
    type Item = WireId;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|id| {
            let Wire(a, b) = &self.diagram.wires[id];
            self.diagram.nodes[a].is_generator()
                && self.diagram.nodes[b].is_generator()
        })
    }
}

/// Iterator type over the nodes neighboring a particular node in a diagram.
///
/// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
pub struct Neighbors<'a> {
    diagram: &'a Diagram,
    source: NodeId,
    iter: std::collections::hash_set::Iter<'a, WireId>,
}

impl<'a> Iterator for Neighbors<'a> {
    type Item = (NodeId, Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
            .map(|wire_id| {
                let other
                    = self.diagram.wires[wire_id].other_end_of(self.source)
                    .unwrap();
                (other, self.diagram.nodes[&other])
            })
    }
}

/// Iterator type over the interior (non-input/output) nodes neighboring a
/// particular node in a diagram.
///
/// The iterator item type is `(`[`NodeId`]`, `[`Node`]`)`.
pub struct NeighborsInner<'a> {
    iter: Neighbors<'a>
}

impl<'a> Iterator for NeighborsInner<'a> {
    type Item = (NodeId, Node);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find(|(_, node)| node.is_generator())
    }
}

/// Create a diagram using abbreviated syntax.
///
/// The first block defines the nodes with the syntax
/// ```text
/// <node_label> = <node_type> [ [arg] ]
/// ```
/// where `<node_type>` is a variant or method of [`NodeDef`]. Labels given to
/// nodes in this block are then used to define wire connections in the second
/// block with the syntax
/// ```text
/// <node1> -- <node2>
/// ```
/// All nodes' labels and corresponding IDs are returned in a
/// [`HashMap`]`<&'static `[`str`]`, `[`NodeId`]`>`.
///
/// The total return type is [`Result`]`<(`[`Diagram`]`, `[`HashMap`]`<&'static
/// `[`str`]`, `[`NodeId`]`>), `[`GraphError`]`>`.
///
/// The normal usage
/// ```
/// # use zx_calc::graph::*;
/// use std::f64::consts::PI;
///
/// # fn main() -> Result<(), GraphError> {
/// let mut diagram = Diagram::new();
///
/// let i = diagram.add_node(NodeDef::Input);
/// let o = diagram.add_node(NodeDef::Output);
/// let n1 = diagram.add_node(NodeDef::Z(PI));
/// let n2 = diagram.add_node(NodeDef::x());
/// let n3 = diagram.add_node(NodeDef::H((-1.0).into()));
/// diagram.add_wire(i, n1)?;
/// diagram.add_wire(o, n3)?;
/// diagram.add_wire(n1, n2)?;
/// diagram.add_wire(n1, n2)?;
/// diagram.add_wire(n2, n3)?;
/// # Ok(())
/// # }
/// ```
/// is equivalent to
/// ```
/// # use zx_calc::graph::*;
/// use zx_calc::diagram;
/// use std::f64::consts::PI;
///
/// diagram!(
///     nodes: {
///         i = input [],
///         o = output [],
///         n1 = Z [PI],
///         n2 = x [],
///         n3 = H [-1.0],
///     },
///     wires: {
///         i -- n1,
///         o -- n3,
///         n1 -- n2,
///         n1 -- n2,
///         n2 -- n3,
///     },
/// );
/// ```
/// where in the macro usage, the diagram is returned alongside the hashmap
/// ```text
/// {
///     "i" => NodeId(0),
///     "o" => NodeId(1),
///     "n1" => NodeId(2),
///     "n2" => NodeId(3),
///     "n3" => NodeId(4),
/// }
/// ```
#[macro_export]
macro_rules! diagram {
    (
        nodes: {
            $( $node_name:ident = $node_def:ident [$( $arg:expr ),* $(,)?] ),*
            $(,)?
        },
        wires: {
            $( $node1_name:ident -- $node2_name:ident ),* $(,)?
        } $(,)?
    ) => {
        {
            let mut diagram = $crate::graph::Diagram::new();
            $(
                let $node_name
                    = diagram.add_node(
                        NodeDef::$node_def($( ($arg).into() ),*)
                    );
            )*
            [$(
                diagram.add_wire($node1_name, $node2_name).map(|_| ()),
            )*]
            .into_iter()
            .collect::<Result<Vec<()>, $crate::graph::GraphError>>()
            .map(|_| {
                let nodes:
                    std::collections::HashMap<
                        &'static str,
                        $crate::graph::NodeId
                    >
                    = [$( (stringify!($node_name), $node_name) ),*]
                    .into_iter()
                    .collect();
                (diagram, nodes)
            })
        }
    }
}

