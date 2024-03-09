#![allow(dead_code, non_snake_case, non_upper_case_globals)]

//! This package contains various tools for working with quantum circuits of
//! qubits and interfacing with the [ZX- and ZH-calculi][wiki].
//!
//! - [`circuit`] provides constructs for ordinary descriptions of quantum
//! circuits using the conventional [circuit notation][qcircuits].
//! - [`ketbra`] computes the complete linear maps or end products represented
//! by diagrams in the ZX(H)-calculus based on naive manipulation of [kets and
//! bras][dirac-bra-ket].
//! - [`graph`] implements many rewrite rules in the calculus using proper
//! graph-based representations, enabling automated diagram simplification.
//!
//! [wiki]: https://en.wikipedia.org/wiki/ZX-calculus
//! [qcircuits]: https://en.wikipedia.org/wiki/Quantum_circuit
//! [dirac-bra-ket]: https://en.wikipedia.org/wiki/Bra%E2%80%93ket_notation
//!
//! # See also
//! - [PyZX](https://github.com/Quantomatic/pyzx): a Python implementation of
//! the ZX-calculus and its rewrite rules.
//! - [QuiZX](https://github.com/Quantomatic/quizx/tree/master): a Rust
//! implementation of the above.
//!
//! # Further reading
//! - B. Coecke, "Basic ZX-calculus for students and professionals."
//! [arXiv:2303.03163](https://arxiv.org/abs/2303.03163)
//! - J. van de Wetering, "ZX-calculus for the working quantum computer
//! scientist." [arXiv:2012.13966](https://arxiv.org/abs/2012.13966)
//! - R. Moyard, "Introduction to the ZX-calculus."
//! [Pennylane](https://pennylane.ai/qml/demos/tutorial_zx_calculus/)
//! - H. Bombin *et al.*, "Unifying flavors of fault tolerance with the ZX
//! calculus." [arXiv:2303.08829](https://arxiv.org/abs/2303.08829)
//!

pub mod circuit;
pub mod ketbra;
pub mod graph;
pub(crate) mod vizdefs;

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

// pub extern crate num_rational;
// /// Create [`num_rational::Rational64`]s with succinct syntax.
// ///
// /// *Panics if rational approximation is impossible, e.g. for NaN or ±inf.*
// ///
// /// ```
// /// use std::f64::consts::PI;
// /// use num_rational::Rational64;
// /// use zx_calc::r;
// ///
// /// assert_eq!(
// ///     r!(PI),
// ///     Rational64::approximate_float(PI)
// ///         .expect("r!: can't approximate floating-point value")
// /// );
// /// assert_eq!( r!(42, 11), Rational64::new(42, 11) );
// /// ```
// #[macro_export]
// macro_rules! r {
//     ( $x:expr ) => {
//         $crate::num_rational::Rational64::approximate_float($x)
//             .expect("r!: can't approximate floating-point value")
//     };
//     ( $numer:expr, $denom:expr ) => {
//         $crate::num_rational::Rational64::new($numer, $denom)
//     };
// }

// /// Create [`num_complex::Complex`]`<`[`num_rational::Rational64`]`>`s with
// /// succinct syntax.
// ///
// /// *Panics if rational approximation is impossible, e.g. for NaN for ±inf.*
// ///
// /// ```
// /// use std::f64::consts::PI;
// /// use num_complex::Complex;
// /// use zx_calc::{ cc, r };
// ///
// /// assert_eq!( cc!(i (-1.0)),    Complex::new(r!(0.0), r!(-1.0))      );
// /// assert_eq!(
// ///     cc!(e PI),
// ///     {
// ///         let z = num_complex::Complex64::cis(PI);
// ///         Complex::new(r!(z.re), r!(z.im))
// ///     }
// /// );
// /// assert_eq!( cc!(1.0),         Complex::new(r!(1.0), r!(0.0))       );
// /// assert_eq!( cc!(1.0 + i 1.0), Complex::new(r!(1.0), r!(1.0))       );
// /// assert_eq!( cc!(1.0 - i 1.0), Complex::new(r!(1.0), r!(-1.0))      );
// /// assert_eq!( cc!(1.0 + 1.0 i), Complex::new(r!(1.0), r!(1.0))       );
// /// assert_eq!( cc!(1.0 - 1.0 i), Complex::new(r!(1.0), r!(-1.0))      );
// /// assert_eq!(
// ///     cc!(2.0, e PI),
// ///     {
// ///         let z = num_complex::Complex64::from_polar(2.0, PI);
// ///         Complex::new(r!(z.re), r!(z.im))
// ///     }
// /// );
// /// assert_eq!( cc!(1.0, 1.0),    Complex::new(r!(1.0), r!(1.0))       );
// /// ```
// ///
// /// See also [`r!`].
// #[macro_export]
// macro_rules! cc {
//     ( i $im:expr )
//         => { $crate::num_complex::Complex::new(r!(0.0), r!($im)) };
//     ( e $ph:expr ) => {
//         {
//             let z = num_complex::Complex64::cis($ph);
//             Complex::new(r!(z.re), r!(z.im))
//         }
//     };
//     ( $re:expr )
//         => { $crate::num_complex::Complex::new(r!($re), r!(0.0)) };
//     ( $re:literal + i $im:literal )
//         => { $crate::num_complex::Complex::new(r!($re), r!($im)) };
//     ( $re:literal - i $im:literal )
//         => { $crate::num_complex::Complex::new(r!($re), r!(-$im)) };
//     ( $re:literal + $im:literal i )
//         => { $crate::num_complex::Complex::new(r!($re), r!($im)) };
//     ( $re:literal - $im:literal i )
//         => { $crate::num_complex::Complex::new(r!($re), r!(-$im)) };
//     ( $r:expr, e $ph:expr ) => {
//         {
//             let z = num_complex::Complex64::from_polar($r, $ph);
//             Complex::new(r!(z.re), r!(z.im))
//         }
//     };
//     ( $re:expr, $im:expr )
//         => { $crate::num_complex::Complex::new(r!($re), r!($im)) };
// }

