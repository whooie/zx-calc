use num_complex::Complex64 as C64;
use crate::ketbra2::{ Element, KBResult, KetBra };

/// Represents a series of [`Element`]s as "slices" of a ZX(H)-diagram, each of
/// whose inputs and outputs can, in principle, locally span the entire lateral
/// space of wires.
#[derive(Clone, Debug)]
pub struct Diagram {
    elems: Vec<Element>,
    // scalar: C64,
}

impl FromIterator<Element> for Diagram {
    fn from_iter<I>(iter: I) -> Self
    where I: IntoIterator<Item = Element>
    {
        Self { elems: iter.into_iter().collect() /*, scalar: 1.0_f64.into()*/ }
    }
}

impl From<Element> for Diagram {
    fn from(elem: Element) -> Self {
        Self { elems: vec![elem] /*, scalar: 1.0_f64.into()*/ }
    }
}

impl Diagram {
    /// Create a new diagram from a list of [`Element`]s.
    ///
    /// `Element`s should follow the order in which they would be drawn in a
    /// real ZX-diagram; i.e. the referse of the order in which dot-products are
    /// taken.
    pub fn new<I>(elems: I) -> Self
    where I: IntoIterator<Item = Element>
    {
        elems.into_iter().collect()
    }

    // /// Set an overall scalar factor.
    // pub fn with_scalar(mut self, z: C64) -> Self {
    //     self.scalar = z;
    //     self
    // }
    //
    // /// Set an overall scalar factor.
    // pub fn set_scalar(&mut self, z: C64) -> &mut Self {
    //     self.scalar = z;
    //     self
    // }

    /// Fold over all slices of the diagram with the [dot
    /// product][Element::dot] and return the result as a single [`Element`].
    pub fn contract(self) -> KBResult<Element> {
        self.elems.into_iter()
            .try_fold(
                None,
                |mb, elem| {
                    mb.map(|acc: Element| acc.into_then(elem)).transpose()
                },
            )
            .map(|mb_res| {
                mb_res.unwrap_or_else(|| {
                    KetBra::new(C64::from(1.0), [], []).into()
                })
            })
    }

    /// Compose `self` with `rhs` by attaching the outputs of `rhs` to the
    /// inputs of `self`.
    ///
    /// See also [`compose_rev`][Self::compose_rev] for composition
    /// following a left-to-right diagram placement.
    pub fn compose(self, rhs: Self) -> Self {
        rhs.elems.into_iter()
            .chain(self.elems)
            .collect()
    }

    /// Compose `self` with `rhs` by attaching the outputs of `self` to the
    /// inputs of `rhs`.
    ///
    /// See also [`compose`][Self::compose`] for composition following the usual
    /// `self ∘ rhs` operation.
    pub fn compose_rev(self, rhs: Self) -> Self {
        self.elems.into_iter()
            .chain(rhs.elems)
            .collect()
    }

    /// Return all input and output wire indices.
    pub fn ins_outs(&self) -> (Vec<usize>, Vec<usize>) {
        let mut ins: Vec<usize> = Vec::new();
        let mut outs: Vec<usize> = Vec::new();
        for element in self.elems.iter() {
            element.ins().into_iter()
                .for_each(|k| {
                    if outs.contains(&k) {
                        vec_remove_elem(&mut outs, &k);
                    } else {
                        ins.push(k);
                    }
                });
            element.outs().into_iter()
                .for_each(|k| { outs.push(k); });
        }
        (ins, outs)
    }

}

fn vec_remove_elem<T>(v: &mut Vec<T>, target: &T) -> Option<T>
where T: PartialEq
{
    v.iter().enumerate()
        .find_map(|(k, t)| (t == target).then_some(k))
        .map(|k0| v.swap_remove(k0))
}

/// Iterator type for [`Diagram`]s over references to each [`Element`] slice.
///
/// The iterator item type is `&`[`Element`].
#[derive(Clone, Debug)]
pub struct DiagramIter<'a> {
    iter: std::slice::Iter<'a, Element>
}

impl<'a> Iterator for DiagramIter<'a> {
    type Item = &'a Element;

    fn next(&mut self) -> Option<Self::Item> { self.iter.next() }
}

impl<'a> DoubleEndedIterator for DiagramIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> { self.iter.next_back() }
}

impl<'a> ExactSizeIterator for DiagramIter<'a> {
    fn len(&self) -> usize { self.iter.len() }
}

impl<'a> std::iter::FusedIterator for DiagramIter<'a> { }

/// Iterator type for [`Diagram`]s over mutable references to each [`Element`]
/// slice.
///
/// The iterator item type is `&mut `[`Element`].
#[derive(Debug)]
pub struct DiagramIterMut<'a> {
    iter: std::slice::IterMut<'a, Element>
}

impl<'a> Iterator for DiagramIterMut<'a> {
    type Item = &'a mut Element;

    fn next(&mut self) -> Option<Self::Item> { self.iter.next() }
}

impl<'a> DoubleEndedIterator for DiagramIterMut<'a> {
    fn next_back(&mut self) -> Option<Self::Item> { self.iter.next_back() }
}

impl<'a> ExactSizeIterator for DiagramIterMut<'a> {
    fn len(&self) -> usize { self.iter.len() }
}

impl<'a> std::iter::FusedIterator for DiagramIterMut<'a> { }

/// Iterator type for [`Diagram`]s over each [`Element`] slice.
///
/// The iterator item type is [`Element`].
#[derive(Clone, Debug)]
pub struct DiagramIntoIter {
    iter: std::vec::IntoIter<Element>
}

impl Iterator for DiagramIntoIter {
    type Item = Element;

    fn next(&mut self) -> Option<Self::Item> { self.iter.next() }
}

impl DoubleEndedIterator for DiagramIntoIter {
    fn next_back(&mut self) -> Option<Self::Item> { self.iter.next_back() }
}

impl ExactSizeIterator for DiagramIntoIter {
    fn len(&self) -> usize { self.iter.len() }
}

impl std::iter::FusedIterator for DiagramIntoIter { }

