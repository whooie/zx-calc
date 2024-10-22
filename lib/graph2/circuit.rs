//! Alternative interface to [`Diagram`] from a circuit-based context.
//!
//! The default gate set is defined via [`Gate`], but can be extended using the
//! [`GateDiagram`] trait. See [`circuit_diagram!`][crate::circuit_diagram] for
//! example usage and abbreviated syntax.

use thiserror::Error;
use crate::{
    graph2::{ Diagram, Node, NodeId, QubitId },
    phase::Phase,
};

#[derive(Debug, Error)]
pub enum CircuitError {
    #[error("error in gate diagram: targeted {0} qubit(s), but the circuit only has {1}")]
    GateDiagramTooManyTargets(usize, usize),

    #[error("error in gate diagram: unequal or mismatching inputs and outputs: targeted {0} qubit(s), but diagram has {1} input(s) and {2} output(s)")]
    GateDiagramIO(usize, usize, usize),

    #[error("error in composition: cannot match {0} free output(s) with {1} free input(s)")]
    NonMatchingIO(usize, usize),

    #[error("diagram error: {0}")]
    GraphError(#[from] crate::graph2::GraphError),
}
pub type CircuitResult<T> = Result<T, CircuitError>;
#[allow(unused_imports)]
use CircuitError::*;

/// A unitary gate to apply in a quantum circuit.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Gate {
    /// Identity.
    I,
    /// Hadamard.
    H(QubitId),
    /// π-rotation about *x*.
    X(QubitId),
    /// Rotation about *x*.
    XRot(QubitId, Phase),
    /// π-rotation about *z*.
    Z(QubitId),
    /// Rotation about *z*.
    ZRot(QubitId, Phase),
    /// π-rotation about *x* on the second qubit, controlled by the first.
    CX(QubitId, QubitId),
    /// Rotation about *x* on the second qubit, controlled by the first.
    CXRot(QubitId, QubitId, Phase),
    /// π-rotation about *z* on the second qubit, controlled by the first.
    CZ(QubitId, QubitId),
    /// Rotation about *z* on the second qubit, controlled by the first.
    CZRot(QubitId, QubitId, Phase),
    /// Mølmer-Sørensen gate. This gate uses the *xx* definition.
    MS(QubitId, QubitId),
    /// Swap gate.
    Swap(QubitId, QubitId),
    /// Root-swap gate.
    SqrtSwap(QubitId, QubitId),
    /// Toffoli gate: π-rotation about *x* on the third qubit, controlled by the
    /// first and second.
    Toff(QubitId, QubitId, QubitId),
}

impl Gate {
    /// Return `true` if `self` is `I`.
    pub fn is_i(&self) -> bool { matches!(self, Self::I) }

    /// Return `true` if `self` is `H`.
    pub fn is_h(&self) -> bool { matches!(self, Self::H(..)) }

    /// Return `true` if `self` is `X`.
    pub fn is_x(&self) -> bool { matches!(self, Self::X(..)) }

    /// Return `true` if `self` is `XRot`.
    pub fn is_xrot(&self) -> bool { matches!(self, Self::XRot(..)) }

    /// Return `true` if `self` is `Z`.
    pub fn is_z(&self) -> bool { matches!(self, Self::Z(..)) }

    /// Return `true` if `self` is `ZRot`.
    pub fn is_zrot(&self) -> bool { matches!(self, Self::ZRot(..)) }

    /// Return `true` if `self` is `CX`.
    pub fn is_cx(&self) -> bool { matches!(self, Self::CX(..)) }

    /// Return `true` if `self` is `CXRot`.
    pub fn is_cxrot(&self) -> bool { matches!(self, Self::CXRot(..)) }

    /// Return `true` if `self` is `CZ`.
    pub fn is_cz(&self) -> bool { matches!(self, Self::CZ(..)) }

    /// Return `true` if `self` is `CZRot`.
    pub fn is_czrot(&self) -> bool { matches!(self, Self::CZRot(..)) }

    /// Return `true` if `self` is `MS`.
    pub fn is_molsor(&self) -> bool { matches!(self, Self::MS(..)) }

    /// Return `true` if `self` is `Swap`.
    pub fn is_swap(&self) -> bool { matches!(self, Self::Swap(..)) }

    /// Return `true` if `self` is `SqrtSwap`.
    pub fn is_sqrtswap(&self) -> bool { matches!(self, Self::SqrtSwap(..)) }

    /// Return `true` if `self` is `Toff`.
    pub fn is_toff(&self) -> bool { matches!(self, Self::Toff(..)) }

    /// Create a new X-rotation gate on qubit `k` with phase `(a / b) × 2π`.
    pub fn new_xrot(k: usize, a: i64, b: i64) -> Self {
        Self::XRot(k, Phase::new(a, b))
    }

    /// Create a new X-rotation gate on qubit `k` with π phase.
    pub fn xrot_pi(k: usize) -> Self { Self::XRot(k, Phase::pi()) }

    /// Create a new X-rotation gate on qubit `k` with π/2 phase.
    pub fn xrot_pi2(k: usize) -> Self { Self::XRot(k, Phase::pi2()) }

    /// Create a new X-rotation gate on qubit `k` with π/4 phase.
    pub fn xrot_pi4(k: usize) -> Self { Self::XRot(k, Phase::pi4()) }

    /// Create a new X-rotation gate on qubit `k` with π/8 phase.
    pub fn xrot_pi8(k: usize) -> Self { Self::XRot(k, Phase::pi8()) }

    /// Create a new X-rotation gate on qubit `k` with phase `2π / n`.
    pub fn xrot_frac(k: usize, n: i64) -> Self { Self::XRot(k, Phase::frac(n)) }

    /// Create a new Z-rotation gate on qubit `k` with phase `(a / b) × 2π`.
    pub fn new_zrot(k: usize, a: i64, b: i64) -> Self {
        Self::ZRot(k, Phase::new(a, b))
    }

    /// Create a new Z-rotation gate on qubit `k` with π phase.
    pub fn zrot_pi(k: usize) -> Self { Self::ZRot(k, Phase::pi()) }

    /// Create a new Z-rotation gate on qubit `k` with π/2 phase.
    pub fn zrot_pi2(k: usize) -> Self { Self::ZRot(k, Phase::pi2()) }

    /// Create a new Z-rotation gate on qubit `k` with π/4 phase.
    pub fn zrot_pi4(k: usize) -> Self { Self::ZRot(k, Phase::pi4()) }

    /// Create a new Z-rotation gate on qubit `k` with π/8 phase.
    pub fn zrot_pi8(k: usize) -> Self { Self::ZRot(k, Phase::pi8()) }

    /// Create a new Z-rotation gate on qubit `k` with phase `2π / n`.
    pub fn zrot_frac(k: usize, n: i64) -> Self { Self::ZRot(k, Phase::frac(n)) }

    /// Create a new CX-rotation gate on qubit `t`, controlled by qubit `c`,
    /// with phase `(a / b) × 2π`.
    pub fn new_cxrot(c: usize, t: usize, a: i64, b: i64) -> Self {
        Self::CXRot(c, t, Phase::new(a, b))
    }

    /// Create a new CX-rotation gate on qubit `t`, controlled by qubit `c`,
    /// with π phase.
    pub fn cxrot_pi(c: usize, t: usize) -> Self {
        Self::CXRot(c, t, Phase::pi()) }

    /// Create a new CX-rotation gate on qubit `t`, controlled by qubit `c`,
    /// with π/2 phase.
    pub fn cxrot_pi2(c: usize, t: usize) -> Self {
        Self::CXRot(c, t, Phase::pi2())
    }

    /// Create a new CX-rotation gate on qubit `t`, controlled by qubit `c`,
    /// with π/4 phase.
    pub fn cxrot_pi4(c: usize, t: usize) -> Self {
        Self::CXRot(c, t, Phase::pi4())
    }

    /// Create a new CX-rotation gate on qubit `t`, controlled by qubit `c`,
    /// with π/8 phase.
    pub fn cxrot_pi8(c: usize, t: usize) -> Self {
        Self::CXRot(c, t, Phase::pi8())
    }

    /// Create a new CX-rotation gate on qubit `t`, controlled by qubit `c`,
    /// with phase `2π / n`.
    pub fn cxrot_frac(c: usize, t: usize, n: i64) -> Self {
        Self::CXRot(c, t, Phase::frac(n))
    }

    /// Create a new CZ-rotation gate on qubit `t`, controlled by qubit `c`,
    /// with phase `(a / b) × 2π`.
    pub fn new_czrot(c: usize, t: usize, a: i64, b: i64) -> Self {
        Self::CZRot(c, t, Phase::new(a, b))
    }

    /// Create a new CZ-rotation gate on qubit `t`, controlled by qubit `c`,
    /// with π phase.
    pub fn czrot_pi(c: usize, t: usize) -> Self {
        Self::CZRot(c, t, Phase::pi()) }

    /// Create a new CZ-rotation gate on qubit `t`, controlled by qubit `c`,
    /// with π/2 phase.
    pub fn czrot_pi2(c: usize, t: usize) -> Self {
        Self::CZRot(c, t, Phase::pi2())
    }

    /// Create a new CZ-rotation gate on qubit `t`, controlled by qubit `c`,
    /// with π/4 phase.
    pub fn czrot_pi4(c: usize, t: usize) -> Self {
        Self::CZRot(c, t, Phase::pi4())
    }

    /// Create a new CZ-rotation gate on qubit `t`, controlled by qubit `c`,
    /// with π/8 phase.
    pub fn czrot_pi8(c: usize, t: usize) -> Self {
        Self::CZRot(c, t, Phase::pi8())
    }

    /// Create a new CZ-rotation gate on qubit `t`, controlled by qubit `c`,
    /// with phase `2π / n`.
    pub fn czrot_frac(c: usize, t: usize, n: i64) -> Self {
        Self::CZRot(c, t, Phase::frac(n))
    }

    /// Return `true` if `other` is the inverse of `self`.
    pub fn is_inv(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::I, Self::I) => true,
            (Self::H(a), Self::H(b)) => a == b,
            (Self::X(a), Self::X(b)) => a == b,
            (Self::XRot(a, ang_a), Self::XRot(b, ang_b))
                => a == b && *ang_a + *ang_b == Phase::zero(),
            (Self::Z(a), Self::Z(b)) => a == b,
            (Self::ZRot(a, ang_a), Self::ZRot(b, ang_b))
                => a == b && *ang_a + *ang_b == Phase::zero(),
            (Self::CX(c_a, t_a), Self::CX(c_b, t_b))
                => c_a == c_b && t_a == t_b,
            (Self::CXRot(c_a, t_a, ang_a), Self::CXRot(c_b, t_b, ang_b))
                => c_a == c_b && t_a == t_b && *ang_a + *ang_b == Phase::zero(),
            (Self::CZ(a_a, b_a), Self::CZ(a_b, b_b))
                => (a_a == a_b && b_a == b_b) || (a_a == b_b && a_b == b_a),
            (Self::CZRot(a_a, b_a, ang_a), Self::CZRot(a_b, b_b, ang_b))
                => (a_a == a_b && b_a == b_b && *ang_a + *ang_b == Phase::zero())
                || (a_a == b_b && a_b == b_a && *ang_a + *ang_b == Phase::zero()),
            (Self::MS(a_a, b_a), Self::MS(a_b, b_b))
                => (a_a == a_b && b_a == b_b) || (a_a == b_b && a_b == b_a),
            (Self::Swap(a_a, b_a), Self::Swap(a_b, b_b))
                => (a_a == a_b && b_a == b_b) || (a_a == b_b && a_b == b_a),
            (Self::SqrtSwap(a_a, b_a), Self::SqrtSwap(a_b, b_b))
                => (a_a == a_b && b_a == b_b) || (a_a == b_b && a_b == b_a),
            (Self::Toff(c0_a, c1_a, t_a), Self::Toff(c0_b, c1_b, t_b))
                => (c0_a == c0_b && c1_a == c1_b && t_a == t_b)
                || (c0_a == c1_b && c1_a == c0_b && t_a == t_b),
            _ => false,
        }
    }

    /// Return the inverse of `self`.
    pub fn inv(&self) -> Self {
        match *self {
            Self::I => Self::I,
            Self::H(k) => Self::H(k),
            Self::X(k) => Self::X(k),
            Self::XRot(k, ang) => Self::XRot(k, -ang),
            Self::Z(k) => Self::Z(k),
            Self::ZRot(k, ang) => Self::ZRot(k, -ang),
            Self::CX(c, t) => Self::CX(c, t),
            Self::CXRot(c, t, ang) => Self::CXRot(c, t, -ang),
            Self::CZ(a, b) => Self::CZ(a, b),
            Self::CZRot(a, b, ang) => Self::CZRot(a, b, -ang),
            Self::MS(a, b) => Self::MS(a, b),
            Self::Swap(a, b) => Self::Swap(a, b),
            Self::SqrtSwap(a, b) => Self::SqrtSwap(a, b),
            Self::Toff(c0, c1, t) => Self::Toff(c0, c1, t),
        }
    }
}

/// Interface to describe an abitrary operation on a subset of qubits using a
/// sub-diagram.
///
/// # Example
pub trait GateDiagram {
    /// Generate a sub-diagram describing the operator to apply to some number
    /// of qubits.
    fn generate(&self) -> Diagram;
}

/// The state of a single qubit.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum State {
    /// ∣0⟩ = ∣+⟩ + ∣–⟩
    Zero,
    /// ∣1⟩ = ∣+⟩ – ∣–⟩
    One,
    /// ∣+⟩ = ∣0⟩ + ∣1⟩
    Plus,
    /// ∣–⟩ = ∣0⟩ - ∣1⟩
    Minus,
    /// Empty, undefined state.
    Undef,
}

/// A post-selected measurement performed on a given single qubit.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Meas(pub QubitId, pub State);

/// A generic circuit operation.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum CircuitOp {
    /// A unitary gate.
    Gate(Gate),
    /// A post-selected, single-qubit measurement.
    Meas(Meas),
}

impl CircuitOp {
    /// Return `true` if `self` is `Gate`.
    pub fn is_gate(&self) -> bool { matches!(self, Self::Gate(_)) }

    /// Return `true` if `self` is `Meas`.
    pub fn is_meas(&self) -> bool { matches!(self, Self::Meas(_)) }
}

impl From<Gate> for CircuitOp {
    fn from(gate: Gate) -> Self { Self::Gate(gate) }
}

impl From<Meas> for CircuitOp {
    fn from(meas: Meas) -> Self { Self::Meas(meas) }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) struct OutWire {
    last: NodeId,
    out: NodeId,
}

/// Circuit-builder interface for a [`Diagram`].
///
/// The definitions used for supported [`Gate`]s are all chosen to use only
/// ZX-calculus primitives, except for [`Toff`][Gate::Toff], which uses a single
/// trinary H-box.
#[derive(Clone, Debug)]
pub struct CircuitDiagram {
    pub(crate) n: usize,
    pub(crate) diagram: Diagram,
    pub(crate) ins: Vec<NodeId>,
    pub(crate) outs: Vec<OutWire>,
}

macro_rules! insert_node {
    ($self:expr, $idx:expr, $nodedef:expr) => {
        if let Some(OutWire { last, out }) = $self.outs.get_mut($idx) {
            $self.diagram.remove_wires(*last, *out, None).unwrap();
            let new = $self.diagram.add_node($nodedef);
            $self.diagram.add_wire(*last, new).unwrap();
            $self.diagram.add_wire(new, *out).unwrap();
            *last = new;
            new
        } else {
            unreachable!()
        }
    }
}

impl AsRef<Diagram> for CircuitDiagram {
    fn as_ref(&self) -> &Diagram { &self.diagram }
}

impl CircuitDiagram {
    /// Create a new `CircuitDiagram` corresponding to the identity on `n`
    /// qubits.
    pub fn new(n: usize) -> Self {
        let mut diagram = Diagram::new();
        let (ins, outs): (Vec<NodeId>, Vec<OutWire>) =
            (0..n).map(|_| {
                let last = diagram.add_node(Node::Input);
                let out = diagram.add_node(Node::Output);
                (last, OutWire { last, out })
            })
            .unzip();
        Self { n, diagram, ins, outs }
    }

    /// Return the number of qubits.
    pub fn n(&self) -> usize { self.n }

    /// Return a reference to the underlying [`Diagram`].
    pub fn as_diagram(&self) -> &Diagram { &self.diagram }

    /// Unpack `self` into a bare [`Diagram`].
    pub fn into_diagram(self) -> Diagram { self.diagram }

    /// Set the input state for a given qubit.
    pub fn set_input(&mut self, k: usize, state: State) -> &mut Self {
        if let Some(id) = self.ins.get(k) {
            if let Some(node) = self.diagram.get_node_mut(*id) {
                match state {
                    State::Zero  => { *node = Node::x();    },
                    State::One   => { *node = Node::x_pi(); },
                    State::Plus  => { *node = Node::z();    },
                    State::Minus => { *node = Node::z_pi(); },
                    State::Undef => { *node = Node::Input;  },
                }
            }
        }
        self
    }

    /// Apply a mid-circuit measurement post-selected to a given outcome.
    ///
    /// The measured qubit is projected into the outcome state and remains alive
    /// for later operations.
    ///
    /// Passing [`State::Undef`] does nothing.
    pub fn measure_postsel(&mut self, k: usize, outcome: State) -> &mut Self {
        if let Some(OutWire { last, out }) = self.outs.get_mut(k) {
            self.diagram.remove_wires(*last, *out, None).unwrap();
            let node =
                match outcome {
                    State::Zero  => { Node::x()    },
                    State::One   => { Node::x_pi() },
                    State::Plus  => { Node::z()    },
                    State::Minus => { Node::z_pi() },
                    State::Undef => { return self; },
                };
            let effect = self.diagram.add_node(node);
            let state = self.diagram.add_node(node);
            self.diagram.add_wire(*last, effect).unwrap();
            self.diagram.add_wire(state, *out).unwrap();
            *last = state;
        }
        self
    }

    /// Apply a Hadamard gate to the `k`-th qubit.
    pub fn apply_h(&mut self, k: usize) -> &mut Self {
        if k >= self.n { return self; }
        insert_node!(self, k, Node::h());
        self
    }

    /// Apply an *X* gate to the `k`-th qubit.
    pub fn apply_x(&mut self, k: usize) -> &mut Self {
        if k >= self.n { return self; }
        insert_node!(self, k, Node::x_pi());
        self
    }

    /// Apply an *X* rotation gate to the `k`-th qubit.
    pub fn apply_xrot(&mut self, k: usize, ang: Phase) -> &mut Self {
        if k >= self.n { return self; }
        insert_node!(self, k, Node::X(ang));
        self
    }

    /// Apply a *Z* gate to the `k`-th qubit.
    pub fn apply_z(&mut self, k: usize) -> &mut Self {
        if k >= self.n { return self; }
        insert_node!(self, k, Node::z_pi());
        self
    }

    /// Apply a *Z* rotation gate to the `k`-th qubit.
    pub fn apply_zrot(&mut self, k: usize, ang: Phase) -> &mut Self {
        if k >= self.n { return self; }
        insert_node!(self, k, Node::Z(ang));
        self
    }

    /// Apply an *X* gate to the `t`-th qubit, controlled in the *z* basis by
    /// the `c`-th qubit.
    pub fn apply_cx(&mut self, c: usize, t: usize) -> &mut Self {
        if c >= self.n || t >= self.n || c == t { return self; }
        let z = insert_node!(self, c, Node::z());
        let x = insert_node!(self, t, Node::x());
        self.diagram.add_wire(z, x).unwrap();
        self
    }

    /// Apply an *X* rotation gate to the `t`-th qubit, controlled in the *z*
    /// basis by the `c`-th qubit.
    pub fn apply_cxrot(&mut self, c: usize, t: usize, ang: Phase) -> &mut Self {
        if c >= self.n || t >= self.n || c == t { return self; }
        self.apply_h(t) // TODO: see if we can do better than this
            .apply_czrot(c, t, ang)
            .apply_h(t)
    }

    /// Apply a *Z* gate controlled in the *z* basis to qubits `a` and `b`.
    ///
    /// This gate is symmetric with respect to its inputs.
    pub fn apply_cz(&mut self, a: usize, b: usize) -> &mut Self {
        if a >= self.n || b >= self.n || a == b { return self; }
        let za = insert_node!(self, a, Node::z());
        let zb = insert_node!(self, b, Node::z());
        let h = self.diagram.add_node(Node::h());
        self.diagram.add_wire(za, h).unwrap();
        self.diagram.add_wire(h, zb).unwrap();
        self
    }

    /// Apply a *Z* rotation gate controlled in the *z* basis to qubits `a` and
    /// `b`.
    ///
    /// This gate is symmetric with respect to its inputs.
    pub fn apply_czrot(&mut self, a: usize, b: usize, ang: Phase) -> &mut Self {
        if a >= self.n || b >= self.n || a == b { return self; }
        let za = insert_node!(self, a, Node::Z(ang / 2));
        let zb = insert_node!(self, b, Node::Z(ang / 2));
        let x = self.diagram.add_node(Node::x());
        let z = self.diagram.add_node(Node::Z(-(ang / 2)));
        self.diagram.add_wire(za, x).unwrap();
        self.diagram.add_wire(x, zb).unwrap();
        self.diagram.add_wire(x, z).unwrap();
        self
    }

    /// Apply a Mølmer-Sørensen gate to the `a`-th and `b`-th qubits.
    ///
    /// This gate is defined as a π/2 rotation about the *xx* axis in the
    /// appropriate two-qubit space and is symmetric with respect to its inputs.
    pub fn apply_ms(&mut self, a: usize, b: usize) -> &mut Self {
        if a >= self.n || b >= self.n || a == b { return self; }
        let xa = insert_node!(self, a, Node::x_pi2());
        let xb = insert_node!(self, b, Node::x_pi2());
        let h = self.diagram.add_node(Node::h());
        self.diagram.add_wire(xa, h).unwrap();
        self.diagram.add_wire(h, xb).unwrap();
        self
    }

    /// Apply a swap gate to the `a`-th and `b`-th qubits.
    ///
    /// This gate is symmetric with respect to its inputs.
    pub fn apply_swap(&mut self, a: usize, b: usize) -> &mut Self {
        if a >= self.n || b >= self.n || a == b { return self; }
        self.diagram.remove_wires(self.outs[a].last, self.outs[a].out, None)
            .unwrap();
        self.diagram.remove_wires(self.outs[b].last, self.outs[b].out, None)
            .unwrap();
        let mut tmp: NodeId = 0;
        std::mem::swap(&mut tmp, &mut self.outs[a].last);
        std::mem::swap(&mut tmp, &mut self.outs[b].last);
        std::mem::swap(&mut tmp, &mut self.outs[a].last);
        self.diagram.add_wire(self.outs[a].last, self.outs[a].out).unwrap();
        self.diagram.add_wire(self.outs[b].last, self.outs[b].out).unwrap();
        self
    }

    /// Apply a √swap gate to the `a`-th and `b`-th qubits.
    ///
    /// This gate is symmetric with respect to its inputs.
    pub fn apply_sqrt_swap(&mut self, a: usize, b: usize) -> &mut Self {
        self.apply_cx(a, b)
            .apply_h(a)
            .apply_czrot(b, a, Phase::pi2())
            .apply_h(a)
            .apply_cx(a, b)
    }

    /// Apply an *X* gate to the `t`-th qubit, controlled by the `c0`-th and
    /// `c1`-th qubits in the *z* basis.
    ///
    /// **Note:** Unlike all other supported gates, this gate uses the
    /// ZH-calculus-based definition, which requires use of a single trinary
    /// H-box.
    pub fn apply_toff(&mut self, c0: usize, c1: usize, t: usize) -> &mut Self {
        let z0 = insert_node!(self, c0, Node::z());
        let z1 = insert_node!(self, c1, Node::z());
        let xt = insert_node!(self, t,  Node::x());
        let h0 = self.diagram.add_node(Node::h());
        let h1 = self.diagram.add_node(Node::h());
        self.diagram.add_wire(z0, h0).unwrap();
        self.diagram.add_wire(z1, h0).unwrap();
        self.diagram.add_wire(h0, h1).unwrap();
        self.diagram.add_wire(h1, xt).unwrap();
        self
    }

    /// Apply a gate.
    pub fn apply_gate(&mut self, gate: Gate) -> &mut Self {
        match gate {
            Gate::I => self,
            Gate::H(k) => self.apply_h(k),
            Gate::X(k) => self.apply_x(k),
            Gate::XRot(k, ang) => self.apply_xrot(k, ang),
            Gate::Z(k) => self.apply_z(k),
            Gate::ZRot(k, ang) => self.apply_zrot(k, ang),
            Gate::CX(c, t) => self.apply_cx(c, t),
            Gate::CXRot(c, t, ang) => self.apply_cxrot(c, t, ang),
            Gate::CZ(a, b) => self.apply_cz(a, b),
            Gate::CZRot(a, b, ang) => self.apply_czrot(a, b, ang),
            Gate::MS(a, b) => self.apply_ms(a, b),
            Gate::Swap(a, b) => self.apply_swap(a, b),
            Gate::SqrtSwap(a, b) => self.apply_sqrt_swap(a, b),
            Gate::Toff(c0, c1, t) => self.apply_toff(c0, c1, t),
        }
    }

    /// Apply a sequence of gates.
    pub fn apply_circuit<'a, I>(&mut self, gates: I) -> &mut Self
    where I: IntoIterator<Item = &'a Gate>
    {
        gates.into_iter().copied().for_each(|g| { self.apply_gate(g); });
        self
    }

    /// Apply an arbitrary operation to a number of qubits, where the operation
    /// is described by some [diagram][GateDiagram].
    ///
    /// The generated diagram must have equal numbers of inputs and outputs
    /// matching the number of unique qubit indices specified in `qubits`, and
    /// must not target more than the number of qubits in the circuit. Does
    /// nothing if any target qubit indices are out of bounds.
    #[allow(unused_variables, unused_mut)]
    pub fn apply_gate_diagram<I, G>(&mut self, qubits: I, gate: &G)
        -> CircuitResult<&mut Self>
    where
        I: IntoIterator<Item = usize>,
        G: GateDiagram,
    {
        let mut targets: Vec<QubitId> = Vec::new();
        qubits.into_iter()
            .for_each(|q| { if !targets.contains(&q) { targets.push(q); } });
        let ntarg = targets.len();
        if ntarg > self.n {
            return Err(GateDiagramTooManyTargets(ntarg, self.n));
        }
        targets.sort_unstable();
        let subdg = gate.generate();
        let ninput = subdg.count_inputs();
        let noutput = subdg.count_outputs();
        if ntarg != ninput || ntarg != noutput {
            return Err(GateDiagramIO(ntarg, ninput, noutput));
        }
        todo!()
    }

    /// Create a copy of `self` with swapped inputs and outputs, the signs of
    /// all spiders' phases flipped, and all H-boxes' arguments conjugated.
    pub fn adjoint(&self) -> Self { self.clone().into_adjoint() }

    /// Swap inputs and outputs, flip the signs of all spiders' phases, and
    /// conjugate all H-boxes' arguments, consuming `self`.
    pub fn into_adjoint(mut self) -> Self {
        self.adjoint_mut();
        self
    }

    /// Swap inputs and outputs, flip the signs of all spiders' phases, and
    /// conjugate all H-boxes' arguments, modifying `self` in place.
    pub fn adjoint_mut(&mut self) -> &mut Self {
        let outs_new: Vec<OutWire> =
            self.ins.iter().copied()
            .map(|in_id| {
                let inner_neighbor =
                    self.diagram.neighbors_of(in_id).unwrap()
                    .next().unwrap();
                OutWire { last: inner_neighbor.0, out: in_id }
            })
            .collect();
        let ins_new: Vec<NodeId> =
            self.outs.iter()
            .map(|OutWire { last: _, out }| *out)
            .collect();
        self.diagram.adjoint_mut();
        self.ins = ins_new;
        self.outs = outs_new;
        self
    }

    /// Return the tensor product of `self` and `other`, copying both.
    pub fn tensor(&self, other: &Self) -> Self {
        self.clone().into_tensor(other.clone())
    }

    /// Return the tensor product of `self` and `other`, consuming both.
    pub fn into_tensor(mut self, other: Self) -> Self {
        self.tensor_with(other);
        self
    }

    /// Compute the tensor product of `self` and `other`, consuming `other` and
    /// modifying `self` in place.
    pub fn tensor_with(&mut self, other: Self) -> &mut Self {
        let node_shift = self.diagram.nodes.len();
        let Self { n, diagram, mut ins, mut outs } = other;
        self.n += n;
        self.diagram.tensor_with(diagram);
        ins.iter_mut().for_each(|n| { *n += node_shift; });
        outs.iter_mut()
            .for_each(|o| {
                o.last += node_shift;
                o.out += node_shift;
            });
        self.ins.append(&mut ins);
        self.outs.append(&mut outs);
        self
    }

    /// Return the composition `self ∘ other`, attempting to match the outputs
    /// of `other` to the inputs of `self`, copying both.
    ///
    /// This operation will fail if `self` and `other` contain different numbers
    /// of qubits.
    ///
    /// See also [`compose_rev`][Self::compose_rev].
    pub fn compose(&self, other: &Self) -> CircuitResult<Self> {
        if self.n != other.n { return Err(NonMatchingIO(other.n, self.n)); }
        self.clone().into_compose(other.clone())
    }

    /// Return the composition `self ∘ other`, attempting to match the outputs
    /// of `other` to the inputs of `self`, consuming both.
    ///
    /// This operation will fail if `self` and `other` contain different numbers
    /// of qubits.
    ///
    /// See also [`into_compose_rev`][Self::into_compose_rev].
    pub fn into_compose(mut self, other: Self) -> CircuitResult<Self> {
        self.compose_with(other)?;
        Ok(self)
    }

    /// Compute the composition `self ∘ other`, attempting to match the outputs
    /// of `other` to the inputs of `self`, consuming `other` and modifying
    /// `self` in place.
    ///
    /// This operation will fail if `self` and `other` contain different numbers
    /// of qubits.
    ///
    /// See also [`compose_with_rev`][Self::compose_with_rev].
    pub fn compose_with(&mut self, other: Self) -> CircuitResult<&mut Self> {
        if self.n != other.n { return Err(NonMatchingIO(other.n, self.n)); }
        let Self { n: _, diagram, ins, outs: _ } = other;
        self.diagram.compose_with(diagram)
            .expect("unexpected composition error");
        self.ins = ins;
        Ok(self)
    }

    /// Return the composition `other ∘ self`, attempting to match the outputs
    /// of `self` to the inputs of `other`, copying both.
    ///
    /// This operation will fail if `self` and `other` contain different numbers
    /// of qubits.
    ///
    /// See also [`compose`][Self::compose].
    pub fn compose_rev(&self, other: &Self) -> CircuitResult<Self> {
        if self.n != other.n { return Err(NonMatchingIO(self.n, other.n)); }
        self.clone().into_compose_rev(other.clone())
    }

    /// Return the composition `other ∘ self`, attempting to match the outputs
    /// of `self` to the inputs of `other`, consuming both.
    ///
    /// This operation will fail if `self` and `other` contain sdifferent
    /// numbers of qubits.
    ///
    /// Seee also [`into_compose`][Self::into_compose].
    pub fn into_compose_rev(mut self, other: Self) -> CircuitResult<Self> {
        self.compose_with_rev(other)?;
        Ok(self)
    }

    /// Compute the composition `other ∘ self`, attempting to match the outputs
    /// of `self` to the inputs of `other`, consuming `other` and modifying
    /// `self` in place.
    ///
    /// This operation will fail if `self` and `other` contain different numbers
    /// of qubits.
    ///
    /// See also [`compose_with`][Self::compose_with].
    pub fn compose_with_rev(&mut self, other: Self) -> CircuitResult<&mut Self> {
        if self.n != other.n { return Err(NonMatchingIO(self.n, other.n)); }
        let Self { n: _, diagram, ins: _, outs } = other;
        self.diagram.compose_with_rev(diagram)
            .expect("unexpected composition error");
        self.outs = outs;
        Ok(self)
    }

    /// Apply a general circuit operation.
    pub fn apply_op<O>(&mut self, op: O) -> &mut Self
    where O: Into<CircuitOp>
    {
        match op.into() {
            CircuitOp::Gate(gate) =>
                self.apply_gate(gate),
            CircuitOp::Meas(Meas(k, outcome)) =>
                self.measure_postsel(k, outcome),
        }
    }

}

/// Create a [`CircuitDiagram`] using abbreviated syntax.
///
/// The first argument defines the number of qubits in the circuit, and the
/// input states for any or all of them are defined in the following block.
/// Circuit operations are subsequently applied using the
/// [`CircuitDiagram::apply_op`] interface. This macro returns a single value of
/// type `CircuitDiagram`.
///
/// # Arguments
/// - `$n`: `usize`,
/// - `$qubit`: [`QubitId`]
/// - `$state`: `impl Into<`[`State`]`>`
/// - `$op`: `impl Into<`[`CircuitOp`]`>`
///
/// # Example
/// The usage
/// ```
/// # use zx_calc::graph2::circuit::*;
/// use zx_calc::circuit_diagram;
///
/// let outcome_a = State::Zero;
/// let outcome_b = State::One;
///
/// circuit_diagram!(
///     qubits: 3,
///     inputs: { 1 => State::Zero, 2 => State::Zero },
///     ops: {
///         // prepare a Bell state on qubits 1 and 2
///         Gate::H(1);
///         Gate::CX(1, 2);
///
///         // teleportation circuit
///         Gate::CX(0, 1);
///         Gate::H(0);
///         Meas(0, outcome_a);
///         Meas(1, outcome_b);
///         if outcome_b == State::One { Gate::X(2) } else { Gate::I };
///         if outcome_a == State::One { Gate::Z(2) } else { Gate::I };
///     }
/// );
/// ```
/// is equivalent to
/// ```
/// # use zx_calc::graph2::circuit::*;
/// let outcome_a = State::Zero;
/// let outcome_b = State::One;
///
/// let mut diagram = CircuitDiagram::new(3);
/// diagram.set_input(1, State::Zero);
/// diagram.set_input(2, State::One);
///
/// // prepare a Bell state on qubits 1 and 2
/// diagram.apply_op(Gate::H(1));
/// diagram.apply_op(Gate::CX(1, 2));
///
/// // teleportation circuit
/// diagram.apply_op(Gate::CX(0, 1));
/// diagram.apply_op(Gate::H(0));
/// diagram.apply_op(Meas(0, outcome_a));
/// diagram.apply_op(Meas(1, outcome_b));
/// diagram.apply_op(if outcome_b == State::One { Gate::X(2) } else { Gate::I });
/// diagram.apply_op(if outcome_a == State::One { Gate::Z(2) } else { Gate::I });
/// ```
#[macro_export]
macro_rules! circuit_diagram {
    (
        qubits: $n:expr,
        inputs: { $( $qubit:expr => $state:expr ),* $(,)? },
        ops: { $( $op:expr );* $(;)? } $(,)?
    ) => {
        {
            let mut diagram = $crate::graph2::circuit::CircuitDiagram::new($n);
            $( diagram.set_input($qubit, $state); )*
            $( diagram.apply_op($op); )*
            diagram
        }
    }
}

pub use crate::circuit_diagram;

