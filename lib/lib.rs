#![allow(dead_code, non_snake_case, non_upper_case_globals)]

//! This package contains various tools for working with quantum circuits of
//! qubits and interfacing with the [ZX- and ZH-calculi][wiki].
//!
//! - [`circuit`] provides constructs for ordinary descriptions of quantum
//!   circuits using the conventional [circuit notation][qcircuits].
//! - [`sliced`] computes diagram end products using naive tensor
//!   multiplication, with dense and sparse backing.
//! - [`graph`] implements the full set of rewrite rules of the ZX- and
//!   ZH-calculus using proper graph-based implementations.
//!
//! [wiki]: https://en.wikipedia.org/wiki/ZX-calculus
//! [qcircuits]: https://en.wikipedia.org/wiki/Quantum_circuit
//! [dirac-bra-ket]: https://en.wikipedia.org/wiki/Bra%E2%80%93ket_notation
//!
//! # See also
//! - [PyZX](https://github.com/Quantomatic/pyzx): a Python implementation of
//!   the ZX-calculus and its rewrite rules.
//! - [QuiZX](https://github.com/Quantomatic/quizx/tree/master): a Rust
//!   implementation of the above.
//!
//! # Further reading
//! - J. van de Wetering, "ZX-calculus for the working quantum computer
//!   scientist." [arXiv:2012.13966](https://arxiv.org/abs/2012.13966)
//! - B. Coecke, "Basic ZX-calculus for students and professionals."
//!   [arXiv:2303.03163](https://arxiv.org/abs/2303.03163)
//! - R. Moyard, "Introduction to the ZX-calculus."
//!   [Pennylane](https://pennylane.ai/qml/demos/tutorial_zx_calculus/)
//! - H. Bombin *et al.*, "Unifying flavors of fault tolerance with the ZX
//!   calculus." [arXiv:2303.08829](https://arxiv.org/abs/2303.08829)
//!

pub mod circuit;
pub mod phase;
pub mod complex;
pub mod graph;
pub mod sliced;
pub mod indexmap;
pub(crate) mod vizdefs;

/// Re-exported for use with the [`c!`] macro.
pub extern crate num_complex;

/// Handy macro to create `num_complex::Complex64`s from more natural and
/// succinct syntax.
///
/// ```
/// use std::f64::consts::PI;
/// use num_complex::Complex64;
/// use zx_calc::c;
///
/// assert_eq!( c!(i (-1.0)),    Complex64::new(0.0, -1.0)      );
/// assert_eq!( c!(e PI),        Complex64::cis(PI)             );
/// assert_eq!( c!(1.0),         Complex64::new(1.0, 0.0)       );
/// assert_eq!( c!(1.0 + i 1.0), Complex64::new(1.0, 1.0)       );
/// assert_eq!( c!(1.0 - i 1.0), Complex64::new(1.0, -1.0)      );
/// assert_eq!( c!(1.0 + 1.0 i), Complex64::new(1.0, 1.0)       );
/// assert_eq!( c!(1.0 - 1.0 i), Complex64::new(1.0, -1.0)      );
/// assert_eq!( c!(1.0, 1.0),    Complex64::new(1.0, 1.0)       );
/// assert_eq!( c!(1.0, e PI),   Complex64::from_polar(1.0, PI) );
/// ```
#[macro_export]
macro_rules! c {
    ( i $im:expr )
        => { $crate::num_complex::Complex64::new(0.0, $im) };
    ( e $ph:expr )
        => { $crate::num_complex::Complex64::cis($ph) };
    ( $re:expr )
        => { $crate::num_complex::Complex64::new($re, 0.0) };
    ( $re:literal + i $im:literal )
        => { $crate::num_complex::Complex64::new($re, $im) };
    ( $re:literal - i $im:literal )
        => { $crate::num_complex::Complex64::new($re, -$im) };
    ( $re:literal + $im:literal i )
        => { $crate::num_complex::Complex64::new($re, $im) };
    ( $re:literal - $im:literal i )
        => { $crate::num_complex::Complex64::new($re, -$im) };
    ( $r:expr, e $ph:expr )
        => { $crate::num_complex::Complex64::from_polar($r, $ph) };
    ( $re:expr, $im:expr )
        => { $crate::num_complex::Complex64::new($re, $im) };
}

