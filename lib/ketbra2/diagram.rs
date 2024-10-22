use crate::ketbra2::Element;

/// Represents a series of [`Element`]s as "slices" of a ZX(H)-diagram, each of
/// whose inputs and outputs can, in principle, locally span the entire lateral
/// space of wires.
#[derive(Clone, Debug)]
pub struct Diagram {
    slices: Vec<Element>
}


