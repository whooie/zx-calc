use thiserror::Error;
#[allow(unused_imports)]
use crate::{
    graph::{
        DiagramData,
        Diagram,
        NodeId,
        TPhase,
        Spider,
        ZX, ZXNode, ZXWire,
        ZH, ZHNode,
        CT, CTNode, CTWire,
    },
    phase::Phase,
};
use super::{ CircuitData, Gate, State };

/// Errors for fallible operations on graph-backed circuit diagrams.
#[derive(Debug, Error)]
pub enum GrError {
    /// Returned when composition is attempted between two diagrams that have
    /// non-matching inputs and outputs.
    #[error("cannot match {0} free output(s) with {1} free input(s)")]
    ComposeIO(usize, usize),

    /// Returned by operations on the underlying [`Diagram`].
    #[error("{0}")]
    GraphError(#[from] crate::graph::GraphError),
}
/// Results with `GrError` errors.
pub type GrResult<T> = Result<T, GrError>;
use GrError::*;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) struct IOWire<A>
where A: DiagramData
{
    io: NodeId,
    inner: A::Wire,
}

/// Graph-based circuit diagram type for use with
/// [`CircuitDiagram`][super::CircuitDiagram].
///
/// The definitions used for supported gates are chosen to use only ZX-calculus
/// primitives except for the Toffoli gate, which uses a single trinary H-box in
/// the ZH implementation.
#[derive(Clone, Debug)]
pub struct Gr<A>
where A: DiagramData
{
    pub(crate) diagram: Diagram<A>,
    pub(crate) ins: Vec<IOWire<A>>,
    pub(crate) outs: Vec<IOWire<A>>,
}

impl Gr<ZX> {
}

#[allow(unused_variables, unused_mut)]
impl CircuitData for Gr<ZX> {
    type Gate = Gate<Phase>;
    type Diagram = Diagram<ZX>;
    type Error = GrError;

    fn num_inputs(&self) -> usize { self.ins.len() }

    fn num_outputs(&self) -> usize { self.outs.len() }

    fn apply_gate(&mut self, gate: &Self::Gate) {
        todo!()
    }

    fn measure_postsel(&mut self, qubit: usize, outcome: State) {
        todo!()
    }

    fn conjugate(&mut self) {
        todo!()
    }

    fn tensor(self, other: Self) -> Self {
        todo!()
    }

    fn compose(self, rhs: Self) -> Result<Self, Self::Error> {
        todo!()
    }

    fn unwrap(self) -> Self::Diagram { self.diagram }
}

impl Gr<ZH> {
    pub(crate) fn insert_node(&mut self, idx: usize, node: ZHNode) -> NodeId {
        let io_id: usize;
        let new_inner: usize;
        if let Some(IOWire { io, inner }) = self.outs.get_mut(idx) {
            io_id = *io;
            self.diagram.remove_wires(*io, *inner, None).unwrap();
            let new = self.diagram.add_node(node);
            self.diagram.add_wire(*inner, new).unwrap();
            self.diagram.add_wire(new, *io).unwrap();
            *inner = new;
            new_inner = new;
        } else {
            panic!()
        }
        let mb_in = self.ins.iter_mut().find(|iow| iow.inner == io_id);
        if let Some(IOWire { io: _, inner }) = mb_in { *inner = new_inner; }
        new_inner
    }

    /// Create a new `Gr<ZH>` corresponding to the identity on `n` qubits.
    pub fn new(n: usize) -> Self {
        let mut diagram = Diagram::<ZH>::new();
        let (ins, outs): (Vec<IOWire<ZH>>, Vec<IOWire<ZH>>) =
            (0..n).map(|_| {
                let input = diagram.add_node(ZHNode::Input);
                let output = diagram.add_node(ZHNode::Output);
                (
                    IOWire { io: input, inner: output },
                    IOWire { io: output, inner: input },
                )
            })
            .unzip();
        Self { diagram, ins, outs }
    }

    /// Set the input state for a given qubit.
    ///
    /// No-op if the qubit doesn't exist.
    pub fn set_input(&mut self, k: usize, state: State) {
        match state {
            State::Zero => {
                self.diagram.apply_state_input(k, Spider::x()).ok();
            },
            State::One => {
                self.diagram.apply_state_input(k, Spider::x_pi()).ok();
            },
            State::Plus => {
                self.diagram.apply_state_input(k, Spider::z()).ok();
            },
            State::Minus => {
                self.diagram.apply_state_input(k, Spider::z_pi()).ok();
            },
            State::Undef => {
                self.diagram.remove_state_input(k).ok();
            },
        }
    }

    /// Apply a Hadamard gate to the `k`-th qubit.
    pub fn apply_h(&mut self, k: usize) {
        if k >= self.outs.len() { return; }
        self.insert_node(k, ZHNode::h());
    }

    /// Apply an *X* gate to the `k`-th qubit.
    pub fn apply_x(&mut self, k: usize) {
        if k >= self.outs.len() { return; }
        self.insert_node(k, ZHNode::x_pi());
    }

    /// Apply an *X* rotation gate to the `k`-th qubit.
    pub fn apply_xrot(&mut self, k: usize, ang: Phase) {
        if k >= self.outs.len() { return; }
        self.insert_node(k, ZHNode::X(ang));
    }

    /// Apply a *Z* gate to the `k`-th qubit.
    pub fn apply_z(&mut self, k: usize) {
        if k >= self.outs.len() { return; }
        self.insert_node(k, ZHNode::z_pi());
    }

    /// Apply a *Z* rotation gate to the `k`-th qubit.
    pub fn apply_zrot(&mut self, k: usize, ang: Phase) {
        if k >= self.outs.len() { return; }
        self.insert_node(k, ZHNode::Z(ang));
    }

    /// Apply an *X* gate to the `t`-th qubit, controlled in the *z* basis by
    /// the `c`-th qubit.
    pub fn apply_cx(&mut self, c: usize, t: usize) {
        if c >= self.outs.len() || t >= self.outs.len() || c == t { return; }
        let z = self.insert_node(c, ZHNode::z());
        let x = self.insert_node(t, ZHNode::x());
        self.diagram.add_wire(z, x).unwrap();
    }

    /// Apply an *X* rotation gate to the `t`-th qubit, controlled in the *z*
    /// basis by the `c`-th qubit.
    pub fn apply_cxrot(&mut self, c: usize, t: usize, ang: Phase) {
        if c >= self.outs.len() || t >= self.outs.len() || c == t { return; }
        self.apply_h(t); // TODO: see if we can do better than this
        self.apply_czrot(c, t, ang);
        self.apply_h(t);
    }

    /// Apply a *Z* gate controlled in the *z* basis to qubits `a` and `b`.
    ///
    /// This gate is symmetric with respect to its inputs.
    pub fn apply_cz(&mut self, a: usize, b: usize) {
        if a >= self.outs.len() || b >= self.outs.len() || a == b { return; }
        let za = self.insert_node(a, ZHNode::z());
        let zb = self.insert_node(b, ZHNode::z());
        let h = self.diagram.add_node(ZHNode::h());
        self.diagram.add_wire(za, h).unwrap();
        self.diagram.add_wire(h, zb).unwrap();
    }

    /// Apply a *Z* rotation gate controlled in the *z* basis to qubits `a` and
    /// `b`.
    ///
    /// This gate is symmetric with respect to its inputs.
    pub fn apply_czrot(&mut self, a: usize, b: usize, ang: Phase) {
        if a >= self.outs.len() || b >= self.outs.len() || a == b { return; }
        let za = self.insert_node(a, ZHNode::Z(ang / 2));
        let zb = self.insert_node(b, ZHNode::Z(ang / 2));
        let x = self.diagram.add_node(ZHNode::x());
        let z = self.diagram.add_node(ZHNode::Z(-(ang / 2)));
        self.diagram.add_wire(za, x).unwrap();
        self.diagram.add_wire(x, zb).unwrap();
        self.diagram.add_wire(x, z).unwrap();
    }

    /// Apply a Mølmer-Sørensen gate to the `a`-th and `b`-th qubits.
    ///
    /// This gate is defined as a π/2 rotation about the *xx* axis in the
    /// appropriate two-qubit space and is symmetric with respect to its inputs.
    pub fn apply_ms(&mut self, a: usize, b: usize) {
        if a >= self.outs.len() || b >= self.outs.len() || a == b { return; }
        let xa = self.insert_node(a, ZHNode::x_pi2());
        let xb = self.insert_node(b, ZHNode::x_pi2());
        let h = self.diagram.add_node(ZHNode::h());
        self.diagram.add_wire(xa, h).unwrap();
        self.diagram.add_wire(h, xb).unwrap();
    }

    /// Apply a swap gate to the `a`-th and `b`-th qubits.
    ///
    /// This gate is symmetric with respect to its inputs.
    pub fn apply_swap(&mut self, a: usize, b: usize) {
        if a >= self.outs.len() || b >= self.outs.len() || a == b { return; }
        self.diagram.remove_wires(self.outs[a].inner, self.outs[a].io, None)
            .unwrap();
        self.diagram.remove_wires(self.outs[b].inner, self.outs[b].io, None)
            .unwrap();
        let mut tmp: NodeId = 0;
        std::mem::swap(&mut tmp, &mut self.outs[a].inner);
        std::mem::swap(&mut tmp, &mut self.outs[b].inner);
        std::mem::swap(&mut tmp, &mut self.outs[a].inner);
        self.diagram.add_wire(self.outs[a].inner, self.outs[a].io).unwrap();
        self.diagram.add_wire(self.outs[b].inner, self.outs[b].io).unwrap();
    }

    /// Apply a √swap gate to the `a`-th and `b`-th qubits.
    ///
    /// This gate is symmetric with respect to its inputs.
    pub fn apply_sqrt_swap(&mut self, a: usize, b: usize) {
        self.apply_cx(a, b);
        self.apply_h(a);
        self.apply_czrot(b, a, Phase::pi2());
        self.apply_h(a);
        self.apply_cx(a, b);
    }

    /// Apply an *X* gate to the `t`-th qubit, controlled by the `c0`-th and
    /// `c1`-th qubits in the *z* basis.
    ///
    /// **Note:** Unlike all other supported gates, this gate uses the
    /// ZH-calculus-based definition, which requires use of a single trinary
    /// H-box.
    pub fn apply_toff(&mut self, c0: usize, c1: usize, t: usize) {
        if c0 >= self.outs.len()
            || c1 >= self.outs.len()
            || t >= self.outs.len()
        {
            return;
        }
        let z0 = self.insert_node(c0, ZHNode::z());
        let z1 = self.insert_node(c1, ZHNode::z());
        let xt = self.insert_node(t,  ZHNode::x());
        let h0 = self.diagram.add_node(ZHNode::h());
        let h1 = self.diagram.add_node(ZHNode::h());
        self.diagram.add_wire(z0, h0).unwrap();
        self.diagram.add_wire(z1, h0).unwrap();
        self.diagram.add_wire(h0, h1).unwrap();
        self.diagram.add_wire(h1, xt).unwrap();
    }

}

#[allow(unused_variables, unused_mut)]
impl CircuitData for Gr<ZH> {
    type Gate = Gate<Phase>;
    type Diagram = Diagram<ZH>;
    type Error = GrError;

    fn num_inputs(&self) -> usize { self.ins.len() }

    fn num_outputs(&self) -> usize { self.outs.len() }

    fn apply_gate(&mut self, gate: &Self::Gate) {
        match *gate {
            Gate::Id => { },
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

    fn measure_postsel(&mut self, qubit: usize, outcome: State) {
        let node =
            match outcome {
                State::Zero  => { ZHNode::x()    },
                State::One   => { ZHNode::x_pi() },
                State::Plus  => { ZHNode::z()    },
                State::Minus => { ZHNode::z_pi() },
                State::Undef => { return;        },
            };
        let io_id: usize;
        let new_inner: usize;
        if let Some(IOWire { io, inner }) = self.outs.get_mut(qubit) {
            io_id = *io;
            self.diagram.remove_wires(*io, *inner, None).unwrap();
            let effect = self.diagram.add_node(node);
            let state = self.diagram.add_node(node);
            self.diagram.add_wire(*inner, effect).unwrap();
            self.diagram.add_wire(state, *io).unwrap();
            *inner = state;
            new_inner = effect;
        } else {
            return;
        }
        let mb_in = self.ins.iter_mut().find(|iow| iow.inner == io_id);
        if let Some(IOWire { io: _, inner }) = mb_in { *inner = new_inner; }
    }

    fn conjugate(&mut self) {
        self.diagram.adjoint_mut();
        std::mem::swap(&mut self.ins, &mut self.outs);
    }

    fn tensor(mut self, other: Self) -> Self {
        let Self { diagram, mut ins, mut outs } = other;
        let shift = self.diagram.tensor_with(diagram);
        ins.iter_mut()
            .for_each(|iow| { iow.io += shift; iow.inner += shift; });
        outs.iter_mut()
            .for_each(|iow| { iow.io += shift; iow.inner += shift; });
        self.ins.append(&mut ins);
        self.outs.append(&mut outs);
        self
    }

    fn compose(mut self, rhs: Self) -> Result<Self, Self::Error> {
        if self.ins.len() != rhs.outs.len() {
            return Err(ComposeIO(rhs.outs.len(), self.ins.len()));
        }
        let Self { diagram, mut ins, outs: _ } = rhs;
        let shift = self.diagram.compose_with(diagram)?;
        ins.iter_mut()
            .for_each(|iow| { iow.io += shift; iow.inner += shift; });
        self.ins = ins;
        Ok(self)
    }

    fn unwrap(self) -> Self::Diagram { self.diagram }
}

impl Gr<CT> {
}

#[allow(unused_variables, unused_mut)]
impl CircuitData for Gr<CT> {
    type Gate = Gate<TPhase>;
    type Diagram = Diagram<CT>;
    type Error = GrError;

    fn num_inputs(&self) -> usize { self.ins.len() }

    fn num_outputs(&self) -> usize { self.outs.len() }

    fn apply_gate(&mut self, gate: &Self::Gate) {
        todo!()
    }

    fn measure_postsel(&mut self, qubit: usize, outcome: State) {
        todo!()
    }

    fn conjugate(&mut self) {
        todo!()
    }

    fn tensor(self, other: Self) -> Self {
        todo!()
    }

    fn compose(self, rhs: Self) -> Result<Self, Self::Error> {
        todo!()
    }

    fn unwrap(self) -> Self::Diagram { self.diagram }
}

