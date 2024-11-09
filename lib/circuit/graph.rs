// //! Alternative interface to [`Diagram`] from a circuit-based context.
// //!
// //! The default gate set is defined via [`Gate`], but can be extended using the
// //! [`GateDiagram`] trait. See [`graph_circuit!`][crate::graph_circuit] for
// //! example usage and abbreviated syntax.
//
// use thiserror::Error;
// use crate::{
//     circuit::{ CircuitOp, Gate, GateDiagram, Meas, State },
//     graph::{ Diagram, Node, NodeId, QubitId },
//     phase::Phase,
// };
//
// #[derive(Debug, Error)]
// pub enum GrError {
//     #[error("targeted {0} qubit(s), but the circuit only has {1}")]
//     GateDiagramTooManyTargets(usize, usize),
//
//     #[error("unequal or mismatching inputs and outputs: targeted {0} qubit(s), but diagram has {1} input(s) and {2} output(s)")]
//     GateDiagramIO(usize, usize, usize),
//
//     #[error("cannot match {0} free output(s) with {1} free input(s)")]
//     NonMatchingIO(usize, usize),
//
//     #[error("{0}")]
//     GraphError(#[from] crate::graph::GraphError),
// }
// pub type GrResult<T> = Result<T, GrError>;
// use GrError::*;
//
// #[derive(Copy, Clone, Debug, PartialEq, Eq)]
// pub(crate) struct OutWire {
//     last: NodeId,
//     out: NodeId,
// }
//
// /// Graph-based circuit diagram type for use with [`Circuit`][super::Circuit].
// ///
// /// The definitions used for supported [`Gate`]s are all chosen to use only
// /// ZX-calculus primitives, except for [`Toff`][Gate::Toff], which uses a single
// /// trinary H-box.
// #[derive(Clone, Debug)]
// pub struct Gr {
//     pub(crate) n: usize,
//     pub(crate) diagram: Diagram,
//     pub(crate) ins: Vec<NodeId>,
//     pub(crate) outs: Vec<OutWire>,
// }
//
// impl AsRef<Diagram> for Gr {
//     fn as_ref(&self) -> &Diagram { &self.diagram }
// }
//
// impl Gr {
//     pub(crate) fn insert_node(&mut self, idx: QubitId, node: Node) -> NodeId {
//         if let Some(OutWire { last, out }) = self.outs.get_mut(idx) {
//             self.diagram.remove_wires(*last, *out, None).unwrap();
//             let new = self.diagram.add_node(node);
//             self.diagram.add_wire(*last, new).unwrap();
//             self.diagram.add_wire(new, *out).unwrap();
//             *last = new;
//             new
//         } else {
//             unreachable!()
//         }
//     }
//
//     /// Create a new `Gr` corresponding to the identity on `n`
//     /// qubits.
//     pub fn new(n: usize) -> Self {
//         let mut diagram = Diagram::new();
//         let (ins, outs): (Vec<NodeId>, Vec<OutWire>) =
//             (0..n).map(|_| {
//                 let last = diagram.add_node(Node::Input);
//                 let out = diagram.add_node(Node::Output);
//                 (last, OutWire { last, out })
//             })
//             .unzip();
//         Self { n, diagram, ins, outs }
//     }
//
//     /// Return the number of qubits.
//     pub fn n(&self) -> usize { self.n }
//
//     /// Return a reference to the underlying [`Diagram`].
//     pub fn as_diagram(&self) -> &Diagram { &self.diagram }
//
//     /// Unpack `self` into a bare [`Diagram`].
//     pub fn into_diagram(self) -> Diagram { self.diagram }
//
//     /// Set the input state for a given qubit.
//     pub fn set_input(&mut self, k: usize, state: State) -> &mut Self {
//         if let Some(id) = self.ins.get(k) {
//             if let Some(node) = self.diagram.get_node_mut(*id) {
//                 match state {
//                     State::Zero  => { *node = Node::x();    },
//                     State::One   => { *node = Node::x_pi(); },
//                     State::Plus  => { *node = Node::z();    },
//                     State::Minus => { *node = Node::z_pi(); },
//                     State::Undef => { *node = Node::Input;  },
//                 }
//             }
//         }
//         self
//     }
//
//     /// Apply a mid-circuit measurement post-selected to a given outcome.
//     ///
//     /// The measured qubit is projected into the outcome state and remains alive
//     /// for later operations.
//     ///
//     /// Passing [`State::Undef`] does nothing.
//     pub fn measure_postsel(&mut self, k: usize, outcome: State) -> &mut Self {
//         if let Some(OutWire { last, out }) = self.outs.get_mut(k) {
//             self.diagram.remove_wires(*last, *out, None).unwrap();
//             let node =
//                 match outcome {
//                     State::Zero  => { Node::x()    },
//                     State::One   => { Node::x_pi() },
//                     State::Plus  => { Node::z()    },
//                     State::Minus => { Node::z_pi() },
//                     State::Undef => { return self; },
//                 };
//             let effect = self.diagram.add_node(node);
//             let state = self.diagram.add_node(node);
//             self.diagram.add_wire(*last, effect).unwrap();
//             self.diagram.add_wire(state, *out).unwrap();
//             *last = state;
//         }
//         self
//     }
//
//     /// Apply a Hadamard gate to the `k`-th qubit.
//     pub fn apply_h(&mut self, k: usize) -> &mut Self {
//         if k >= self.n { return self; }
//         self.insert_node(k, Node::h());
//         self
//     }
//
//     /// Apply an *X* gate to the `k`-th qubit.
//     pub fn apply_x(&mut self, k: usize) -> &mut Self {
//         if k >= self.n { return self; }
//         self.insert_node(k, Node::x_pi());
//         self
//     }
//
//     /// Apply an *X* rotation gate to the `k`-th qubit.
//     pub fn apply_xrot(&mut self, k: usize, ang: Phase) -> &mut Self {
//         if k >= self.n { return self; }
//         self.insert_node(k, Node::X(ang));
//         self
//     }
//
//     /// Apply a *Z* gate to the `k`-th qubit.
//     pub fn apply_z(&mut self, k: usize) -> &mut Self {
//         if k >= self.n { return self; }
//         self.insert_node(k, Node::z_pi());
//         self
//     }
//
//     /// Apply a *Z* rotation gate to the `k`-th qubit.
//     pub fn apply_zrot(&mut self, k: usize, ang: Phase) -> &mut Self {
//         if k >= self.n { return self; }
//         self.insert_node(k, Node::Z(ang));
//         self
//     }
//
//     /// Apply an *X* gate to the `t`-th qubit, controlled in the *z* basis by
//     /// the `c`-th qubit.
//     pub fn apply_cx(&mut self, c: usize, t: usize) -> &mut Self {
//         if c >= self.n || t >= self.n || c == t { return self; }
//         let z = self.insert_node(c, Node::z());
//         let x = self.insert_node(t, Node::x());
//         self.diagram.add_wire(z, x).unwrap();
//         self
//     }
//
//     /// Apply an *X* rotation gate to the `t`-th qubit, controlled in the *z*
//     /// basis by the `c`-th qubit.
//     pub fn apply_cxrot(&mut self, c: usize, t: usize, ang: Phase) -> &mut Self {
//         if c >= self.n || t >= self.n || c == t { return self; }
//         self.apply_h(t) // TODO: see if we can do better than this
//             .apply_czrot(c, t, ang)
//             .apply_h(t)
//     }
//
//     /// Apply a *Z* gate controlled in the *z* basis to qubits `a` and `b`.
//     ///
//     /// This gate is symmetric with respect to its inputs.
//     pub fn apply_cz(&mut self, a: usize, b: usize) -> &mut Self {
//         if a >= self.n || b >= self.n || a == b { return self; }
//         let za = self.insert_node(a, Node::z());
//         let zb = self.insert_node(b, Node::z());
//         let h = self.diagram.add_node(Node::h());
//         self.diagram.add_wire(za, h).unwrap();
//         self.diagram.add_wire(h, zb).unwrap();
//         self
//     }
//
//     /// Apply a *Z* rotation gate controlled in the *z* basis to qubits `a` and
//     /// `b`.
//     ///
//     /// This gate is symmetric with respect to its inputs.
//     pub fn apply_czrot(&mut self, a: usize, b: usize, ang: Phase) -> &mut Self {
//         if a >= self.n || b >= self.n || a == b { return self; }
//         let za = self.insert_node(a, Node::Z(ang / 2));
//         let zb = self.insert_node(b, Node::Z(ang / 2));
//         let x = self.diagram.add_node(Node::x());
//         let z = self.diagram.add_node(Node::Z(-(ang / 2)));
//         self.diagram.add_wire(za, x).unwrap();
//         self.diagram.add_wire(x, zb).unwrap();
//         self.diagram.add_wire(x, z).unwrap();
//         self
//     }
//
//     /// Apply a Mølmer-Sørensen gate to the `a`-th and `b`-th qubits.
//     ///
//     /// This gate is defined as a π/2 rotation about the *xx* axis in the
//     /// appropriate two-qubit space and is symmetric with respect to its inputs.
//     pub fn apply_ms(&mut self, a: usize, b: usize) -> &mut Self {
//         if a >= self.n || b >= self.n || a == b { return self; }
//         let xa = self.insert_node(a, Node::x_pi2());
//         let xb = self.insert_node(b, Node::x_pi2());
//         let h = self.diagram.add_node(Node::h());
//         self.diagram.add_wire(xa, h).unwrap();
//         self.diagram.add_wire(h, xb).unwrap();
//         self
//     }
//
//     /// Apply a swap gate to the `a`-th and `b`-th qubits.
//     ///
//     /// This gate is symmetric with respect to its inputs.
//     pub fn apply_swap(&mut self, a: usize, b: usize) -> &mut Self {
//         if a >= self.n || b >= self.n || a == b { return self; }
//         self.diagram.remove_wires(self.outs[a].last, self.outs[a].out, None)
//             .unwrap();
//         self.diagram.remove_wires(self.outs[b].last, self.outs[b].out, None)
//             .unwrap();
//         let mut tmp: NodeId = 0;
//         std::mem::swap(&mut tmp, &mut self.outs[a].last);
//         std::mem::swap(&mut tmp, &mut self.outs[b].last);
//         std::mem::swap(&mut tmp, &mut self.outs[a].last);
//         self.diagram.add_wire(self.outs[a].last, self.outs[a].out).unwrap();
//         self.diagram.add_wire(self.outs[b].last, self.outs[b].out).unwrap();
//         self
//     }
//
//     /// Apply a √swap gate to the `a`-th and `b`-th qubits.
//     ///
//     /// This gate is symmetric with respect to its inputs.
//     pub fn apply_sqrt_swap(&mut self, a: usize, b: usize) -> &mut Self {
//         self.apply_cx(a, b)
//             .apply_h(a)
//             .apply_czrot(b, a, Phase::pi2())
//             .apply_h(a)
//             .apply_cx(a, b)
//     }
//
//     /// Apply an *X* gate to the `t`-th qubit, controlled by the `c0`-th and
//     /// `c1`-th qubits in the *z* basis.
//     ///
//     /// **Note:** Unlike all other supported gates, this gate uses the
//     /// ZH-calculus-based definition, which requires use of a single trinary
//     /// H-box.
//     pub fn apply_toff(&mut self, c0: usize, c1: usize, t: usize) -> &mut Self {
//         let z0 = self.insert_node(c0, Node::z());
//         let z1 = self.insert_node(c1, Node::z());
//         let xt = self.insert_node(t,  Node::x());
//         let h0 = self.diagram.add_node(Node::h());
//         let h1 = self.diagram.add_node(Node::h());
//         self.diagram.add_wire(z0, h0).unwrap();
//         self.diagram.add_wire(z1, h0).unwrap();
//         self.diagram.add_wire(h0, h1).unwrap();
//         self.diagram.add_wire(h1, xt).unwrap();
//         self
//     }
//
//     /// Apply a gate.
//     pub fn apply_gate(&mut self, gate: Gate) -> &mut Self {
//         match gate {
//             Gate::I => self,
//             Gate::H(k) => self.apply_h(k),
//             Gate::X(k) => self.apply_x(k),
//             Gate::XRot(k, ang) => self.apply_xrot(k, ang),
//             Gate::Z(k) => self.apply_z(k),
//             Gate::ZRot(k, ang) => self.apply_zrot(k, ang),
//             Gate::CX(c, t) => self.apply_cx(c, t),
//             Gate::CXRot(c, t, ang) => self.apply_cxrot(c, t, ang),
//             Gate::CZ(a, b) => self.apply_cz(a, b),
//             Gate::CZRot(a, b, ang) => self.apply_czrot(a, b, ang),
//             Gate::MS(a, b) => self.apply_ms(a, b),
//             Gate::Swap(a, b) => self.apply_swap(a, b),
//             Gate::SqrtSwap(a, b) => self.apply_sqrt_swap(a, b),
//             Gate::Toff(c0, c1, t) => self.apply_toff(c0, c1, t),
//         }
//     }
//
//     /// Apply a sequence of gates.
//     pub fn apply_circuit<'a, I>(&mut self, gates: I) -> &mut Self
//     where I: IntoIterator<Item = &'a Gate>
//     {
//         gates.into_iter().copied().for_each(|g| { self.apply_gate(g); });
//         self
//     }
//
//     /// Apply an arbitrary operation to a number of qubits, where the operation
//     /// is described by some [diagram][GateDiagram].
//     ///
//     /// The generated diagram must have equal numbers of inputs and outputs
//     /// matching the number of unique qubit indices specified in `qubits`, and
//     /// must not target more than the number of qubits in the circuit. Does
//     /// nothing if any target qubit indices are out of bounds.
//     #[allow(unused_variables, unused_mut)]
//     pub fn apply_gate_diagram<I, G>(&mut self, qubits: I, gate: &G)
//         -> GrResult<&mut Self>
//     where
//         I: IntoIterator<Item = usize>,
//         G: GateDiagram<Self>,
//     {
//         let mut targets: Vec<QubitId> = Vec::new();
//         qubits.into_iter()
//             .for_each(|q| { if !targets.contains(&q) { targets.push(q); } });
//         let ntarg = targets.len();
//         if ntarg > self.n {
//             return Err(GateDiagramTooManyTargets(ntarg, self.n));
//         }
//         targets.sort_unstable();
//         let subdg = gate.generate();
//         let ninput = subdg.diagram.count_inputs();
//         let noutput = subdg.diagram.count_outputs();
//         if ntarg != ninput || ntarg != noutput {
//             return Err(GateDiagramIO(ntarg, ninput, noutput));
//         }
//         todo!()
//     }
//
//     /// Swap inputs and outputs, flip the signs of all spiders' phases, and
//     /// conjugate all H-boxes' arguments as well as the global scalar, consuming
//     /// `self`.
//     pub fn adjoint(mut self) -> Self {
//         self.adjoint_mut();
//         self
//     }
//
//     /// Swap inputs and outputs, flip the signs of all spiders' phases, and
//     /// conjugate all H-boxes' arguments as well as the global scalar, modifying
//     /// `self` in place.
//     pub fn adjoint_mut(&mut self) -> &mut Self {
//         let outs_new: Vec<OutWire> =
//             self.ins.iter().copied()
//             .map(|in_id| {
//                 let inner_neighbor =
//                     self.diagram.neighbors_of(in_id).unwrap()
//                     .next().unwrap();
//                 OutWire { last: inner_neighbor.0, out: in_id }
//             })
//             .collect();
//         let ins_new: Vec<NodeId> =
//             self.outs.iter()
//             .map(|OutWire { last: _, out }| *out)
//             .collect();
//         self.diagram.adjoint_mut();
//         self.ins = ins_new;
//         self.outs = outs_new;
//         self
//     }
//
//     /// Return the tensor product of `self` and `other`, consuming both.
//     pub fn tensor(mut self, other: Self) -> Self {
//         self.tensor_with(other);
//         self
//     }
//
//     /// Compute the tensor product of `self` and `other`, consuming `other` and
//     /// modifying `self` in place.
//     pub fn tensor_with(&mut self, other: Self) -> &mut Self {
//         let node_shift = self.diagram.nodes.len();
//         let Self { n, diagram, mut ins, mut outs } = other;
//         self.n += n;
//         self.diagram.tensor_with(diagram);
//         ins.iter_mut().for_each(|n| { *n += node_shift; });
//         outs.iter_mut()
//             .for_each(|o| {
//                 o.last += node_shift;
//                 o.out += node_shift;
//             });
//         self.ins.append(&mut ins);
//         self.outs.append(&mut outs);
//         self
//     }
//
//     /// Return the composition `self ∘ other`, attempting to match the outputs
//     /// of `other` to the inputs of `self`, consuming both.
//     ///
//     /// This operation will fail if `self` and `other` contain different numbers
//     /// of qubits.
//     ///
//     /// See also [`compose_rev`][Self::compose_rev].
//     pub fn compose(mut self, other: Self) -> GrResult<Self> {
//         self.compose_with(other)?;
//         Ok(self)
//     }
//
//     /// Compute the composition `self ∘ other`, attempting to match the outputs
//     /// of `other` to the inputs of `self`, consuming `other` and modifying
//     /// `self` in place.
//     ///
//     /// This operation will fail if `self` and `other` contain different numbers
//     /// of qubits.
//     ///
//     /// See also [`compose_with_rev`][Self::compose_with_rev].
//     pub fn compose_with(&mut self, other: Self) -> GrResult<&mut Self> {
//         if self.n != other.n { return Err(NonMatchingIO(other.n, self.n)); }
//         let Self { n: _, diagram, ins, outs: _ } = other;
//         self.diagram.compose_with(diagram)
//             .expect("unexpected composition error");
//         self.ins = ins;
//         Ok(self)
//     }
//
//     /// Return the composition `other ∘ self`, attempting to match the outputs
//     /// of `self` to the inputs of `other`, consuming both.
//     ///
//     /// This operation will fail if `self` and `other` contain sdifferent
//     /// numbers of qubits.
//     ///
//     /// Seee also [`compose`][Self::compose].
//     pub fn compose_rev(mut self, other: Self) -> GrResult<Self> {
//         self.compose_with_rev(other)?;
//         Ok(self)
//     }
//
//     /// Compute the composition `other ∘ self`, attempting to match the outputs
//     /// of `self` to the inputs of `other`, consuming `other` and modifying
//     /// `self` in place.
//     ///
//     /// This operation will fail if `self` and `other` contain different numbers
//     /// of qubits.
//     ///
//     /// See also [`compose_with`][Self::compose_with].
//     pub fn compose_with_rev(&mut self, other: Self) -> GrResult<&mut Self> {
//         if self.n != other.n { return Err(NonMatchingIO(self.n, other.n)); }
//         let Self { n: _, diagram, ins: _, outs } = other;
//         self.diagram.compose_with_rev(diagram)
//             .expect("unexpected composition error");
//         self.outs = outs;
//         Ok(self)
//     }
//
//     /// Apply a general circuit operation.
//     pub fn apply_op<O>(&mut self, op: O) -> &mut Self
//     where O: Into<CircuitOp>
//     {
//         match op.into() {
//             CircuitOp::Gate(gate) =>
//                 self.apply_gate(gate),
//             CircuitOp::Meas(Meas(k, outcome)) =>
//                 self.measure_postsel(k, outcome),
//         }
//     }
//
// }

