//! Provides tools for conventional circuit notation.

pub enum Gate {
    Z(f64),
    X(f64),
    Y(f64),
    H(usize),
    CNOT(usize, usize),
    Toff(usize, usize, usize),
}

pub struct Circuit {
    qubits: usize,
    gates: Vec<Gate>,
}

