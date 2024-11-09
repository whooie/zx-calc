// //! Alternative interface to [`Diagram`] from a circuit-based context.
// //!
// //! The default gate set is defined via [`Gate`], but can be extended using the
// //! [`GateDiagram`] trait. See [`tensor_circuit!`][crate::tensor_circuit] for
// //! example usage and abbreviated syntax.
//
// use thiserror::Error;
// use crate::{ tensor::Diagram, phase::Phase };
//
// #[derive(Debug, Error)]
// pub enum CircuitError {
// }
