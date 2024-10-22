//! Provides tools for conventional circuit notation.

pub enum Gate {
    X(usize, f64),
    Y(usize, f64),
    Z(usize, f64),
    H(usize),
    CX(usize, usize, f64),
    CY(usize, usize, f64),
    CZ(usize, usize, f64),
    Toff(usize, usize, usize),
    Swap(usize, usize),
}

pub struct Circuit {
    qubits: usize,
    gates: Vec<Gate>,
}

