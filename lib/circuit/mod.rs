use thiserror::Error;
use crate::phase::Phase;

#[derive(Debug, Error)]
pub enum CircuitError {
    #[error("error in gate diagram: targeted {0} qubit(s), but the circuit only has {1}")]
    GateDiagramTooManyTargets(usize, usize),

    #[error("error in gate diagram: unequal or mismatching inputs and outputs: targeted {0} qubit(s), but diagram has {1} input(s) and {2} output(s)")]
    GateDiagramIO(usize, usize, usize),

    #[error("error in composition: cannot match {0} free output(s) with {1} free input(s)")]
    NonMatchingIO(usize, usize),

    #[error("diagram error: {0}")]
    GraphError(#[from] crate::graph::GraphError),
}
pub type CircuitResult<T> = Result<T, CircuitError>;
#[allow(unused_imports)]
use CircuitError::*;

pub(crate) mod graph;
pub(crate) mod tensor;

/// A unitary gate to apply in a quantum circuit.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Gate {
    /// Identity.
    I,
    /// Hadamard.
    H(usize),
    /// π-rotation about *x*.
    X(usize),
    /// Rotation about *x*.
    XRot(usize, Phase),
    /// π-rotation about *z*.
    Z(usize),
    /// Rotation about *z*.
    ZRot(usize, Phase),
    /// π-rotation about *x* on the second qubit, controlled by the first.
    CX(usize, usize),
    /// Rotation about *x* on the second qubit, controlled by the first.
    CXRot(usize, usize, Phase),
    /// π-rotation about *z* on the second qubit, controlled by the first.
    CZ(usize, usize),
    /// Rotation about *z* on the second qubit, controlled by the first.
    CZRot(usize, usize, Phase),
    /// Mølmer-Sørensen gate. This gate uses the *xx* definition.
    MS(usize, usize),
    /// Swap gate.
    Swap(usize, usize),
    /// Root-swap gate.
    SqrtSwap(usize, usize),
    /// Toffoli gate: π-rotation about *x* on the third qubit, controlled by the
    /// first and second.
    Toff(usize, usize, usize),
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

    /// Create a new X-rotation gate on qubit `k` with π phase.
    pub fn xrot_pi(k: usize) -> Self { Self::XRot(k, Phase::pi()) }

    /// Create a new X-rotation gate on qubit `k` with π/2 phase.
    pub fn xrot_pi2(k: usize) -> Self { Self::XRot(k, Phase::pi2()) }

    /// Create a new X-rotation gate on qubit `k` with π/4 phase.
    pub fn xrot_pi4(k: usize) -> Self { Self::XRot(k, Phase::pi4()) }

    /// Create a new X-rotation gate on qubit `k` with π/8 phase.
    pub fn xrot_pi8(k: usize) -> Self { Self::XRot(k, Phase::pi8()) }

    /// Create a new X-rotation gate on qubit `k` with phase `(a / b) × 2π`.
    pub fn xrot_frac(k: usize, a: i64, b: i64) -> Self {
        Self::XRot(k, Phase::new(a, b))
    }

    /// Create a new Z-rotation gate on qubit `k` with π phase.
    pub fn zrot_pi(k: usize) -> Self { Self::ZRot(k, Phase::pi()) }

    /// Create a new Z-rotation gate on qubit `k` with π/2 phase.
    pub fn zrot_pi2(k: usize) -> Self { Self::ZRot(k, Phase::pi2()) }

    /// Create a new Z-rotation gate on qubit `k` with π/4 phase.
    pub fn zrot_pi4(k: usize) -> Self { Self::ZRot(k, Phase::pi4()) }

    /// Create a new Z-rotation gate on qubit `k` with π/8 phase.
    pub fn zrot_pi8(k: usize) -> Self { Self::ZRot(k, Phase::pi8()) }

    /// Create a new Z-rotation gate on qubit `k` with phase `(a / b) × 2π`.
    pub fn zrot_frac(k: usize, a: i64, b: i64) -> Self {
        Self::ZRot(k, Phase::new(a, b))
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
    /// with phase `(a / b) × 2π`.
    pub fn cxrot_frac(c: usize, t: usize, a: i64, b: i64) -> Self {
        Self::CXRot(c, t, Phase::new(a, b))
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
    /// with phase `(a / b) × 2π`.
    pub fn czrot_frac(c: usize, t: usize, a: i64, b: i64) -> Self {
        Self::CZRot(c, t, Phase::new(a, b))
    }

    /// Return `true` if `other` is the inverse of `self`.
    pub fn is_inv(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::I, Self::I) => true,
            (Self::H(a), Self::H(b)) => a == b,
            (Self::X(a), Self::X(b)) => a == b,
            (Self::XRot(a, ang_a), Self::XRot(b, ang_b)) =>
                a == b && *ang_a + *ang_b == Phase::zero(),
            (Self::Z(a), Self::Z(b)) => a == b,
            (Self::ZRot(a, ang_a), Self::ZRot(b, ang_b)) =>
                a == b && *ang_a + *ang_b == Phase::zero(),
            (Self::CX(c_a, t_a), Self::CX(c_b, t_b)) =>
                c_a == c_b && t_a == t_b,
            (Self::CXRot(c_a, t_a, ang_a), Self::CXRot(c_b, t_b, ang_b)) =>
                c_a == c_b && t_a == t_b && *ang_a + *ang_b == Phase::zero(),
            (Self::CZ(a_a, b_a), Self::CZ(a_b, b_b)) =>
                (a_a == a_b && b_a == b_b) || (a_a == b_b && a_b == b_a),
            (Self::CZRot(a_a, b_a, ang_a), Self::CZRot(a_b, b_b, ang_b)) =>
                (a_a == a_b && b_a == b_b && *ang_a + *ang_b == Phase::zero())
                || (a_a == b_b && a_b == b_a && *ang_a + *ang_b == Phase::zero()),
            (Self::MS(a_a, b_a), Self::MS(a_b, b_b)) =>
                (a_a == a_b && b_a == b_b) || (a_a == b_b && a_b == b_a),
            (Self::Swap(a_a, b_a), Self::Swap(a_b, b_b)) =>
                (a_a == a_b && b_a == b_b) || (a_a == b_b && a_b == b_a),
            (Self::SqrtSwap(a_a, b_a), Self::SqrtSwap(a_b, b_b)) =>
                (a_a == a_b && b_a == b_b) || (a_a == b_b && a_b == b_a),
            (Self::Toff(c0_a, c1_a, t_a), Self::Toff(c0_b, c1_b, t_b)) =>
                (c0_a == c0_b && c1_a == c1_b && t_a == t_b)
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
pub trait GateDiagram<D> {
    /// Generate a sub-diagram describing the operator to apply to some number
    /// of qubits.
    fn generate(&self) -> D;
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
pub struct Meas(pub usize, pub State);

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

/// `CircuitData` is a sealed trait to govern diagram composition -- we want
/// `Circuit` to be generic over diagram implementations (tensor- and
/// graph-based), but ensure proper behavior under gate application and diagram
/// composition.
pub(crate) mod private { pub trait CircuitDataSeal { } }
pub(crate) use private::CircuitDataSeal;

/// Trait for diagram types to back [`CircuitDiagram`].
///
/// This trait is sealed to ensure proper behavior under gate application and
/// diagram composition.
pub trait CircuitData: Sized + CircuitDataSeal {
    /// Associated error type for general operations.
    type Error;

    /// Return the number of free input qubit wires.
    fn num_inputs(&self) -> usize;

    /// Return the number of free output qubit wires.
    fn num_outputs(&self) -> usize;

    /// Apply a [`Gate`].
    ///
    /// Applying the gate to a non-existant qubit does nothing.
    fn apply_gate(&mut self, gate: Gate);

    /// Apply an arbitrary operation described by a (generally complex)
    /// [diagram][GateDiagram].
    ///
    /// The generated diagram must have equal numbers of inputs and outputs
    /// matching the number of unique qubit indices specified in `qubits`, and
    /// must not target more than the number of qubits in the circuit. Does
    /// nothing if any target qubit indices are out of bounds.
    fn apply_gate_diagram<G>(&mut self, gate: &G) -> Result<(), Self::Error>
    where G: GateDiagram<Self>;

    /// Apply a mid-circuit measurement post-selected to a specific outcome.
    ///
    /// The measured qubit is projected into the outcome state and remains alive
    /// for later operations.
    ///
    /// Passing [`State::Undef`] or applying to a non-existant qubit does
    /// nothing.
    fn measure_postsel(&mut self, qubit: usize, outcome: State);

    /// Conjugate `self` in place, swapping inputs and outputs and conjugating
    /// all complex numbers.
    fn conjugate(&mut self);

    /// Return the tensor product of `self ⊗ other`, consuming both.
    fn tensor(self, other: Self) -> Self;

    /// Return the composition `self ∘ rhs`, attempting to match the outputs of
    /// `other` to the inputs of `self`, consuming both.
    ///
    /// This operation will fail if `self` and `other` contain different numbers
    /// of qubits.
    fn compose(self, rhs: Self) -> Result<Self, Self::Error>;
}

/// A quantum circuit backed by an implementation `A`.
#[derive(Copy, Clone, Debug)]
pub struct CircuitDiagram<A>(pub(crate) A);

