#![allow(dead_code, non_snake_case, non_upper_case_globals)]

pub mod graph;
pub mod ketbra;

pub extern crate num_complex;
/// Handy macro to create `num_complex::Complex64`s from more natural and
/// succinct syntax.
///
/// ```
/// # use zx_calc::c;
/// use num_complex::Complex64;
/// use std::f64::consts::PI;
///
/// assert_eq!( c!(1.0),         Complex64::new(1.0, 0.0)       );
/// assert_eq!( c!(i 1.0),       Complex64::new(0.0, 1.0)       );
/// assert_eq!( c!(e PI),        Complex64::cis(PI)             );
/// assert_eq!( c!(1.0 + i 1.0), Complex64::new(1.0, 1.0)       );
/// assert_eq!( c!(1.0 - i 1.0), Complex64::new(1.0, -1.0)      );
/// assert_eq!( c!(1.0 + 1.0 i), Complex64::new(1.0, 1.0)       );
/// assert_eq!( c!(1.0 - 1.0 i), Complex64::new(1.0, -1.0)      );
/// assert_eq!( c!(1.0, 1.0),    Complex64::new(1.0, 1.0)       );
/// assert_eq!( c!(1.0, e PI),   Complex64::from_polar(1.0, PI) );
/// ```
#[macro_export]
macro_rules! c {
    ( $re:expr )
        => { $crate::num_complex::Complex64::new($re, 0.0) };
    ( i $im:expr )
        => { $crate::num_complex::Complex64::new(0.0, $im) };
    ( e $ph:expr )
        => { $crate::num_complex::Complex64::cis($ph) };
    ( $re:literal + i $im:literal )
        => { $crate::num_complex::Complex64::new($re, $im) };
    ( $re:literal - i $im:literal )
        => { $crate::num_complex::Complex64::new($re, -$im) };
    ( $re:literal + $im:literal i )
        => { $crate::num_complex::Complex64::new($re, $im) };
    ( $re:literal - $im:literal i )
        => { $crate::num_complex::Complex64::new($re, -$im) };
    ( $re:expr, $im:expr )
        => { $crate::num_complex::Complex64::new($re, $im) };
    ( $r:expr, e $ph:expr )
        => { $crate::num_complex::Complex64::from_polar($r, $ph) };
}

