use num_complex::Complex64 as C64;
use crate::phase::Phase;
use super::*;

/// Convert all unary H-boxes into unary Z-spiders.
///
/// This is the comprehensive version of [`HState`], searching more efficiently
/// for all unary H-boxes.
///
/// ![h_state][h_state]
#[embed_doc_image::embed_doc_image("h_state", "assets/rules/HState.svg")]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct HStateAll;

/// Output of [`HStateAll::find`].
#[derive(Debug)]
pub struct HStateAllData<'a> {
    pub(crate) dg: &'a mut Diagram,
    pub(crate) hh: Vec<(NodeId, C64)>,
}

impl<'a> HStateAllData<'a> {
    /// Return the number of unary H-boxes found with argument ±1.
    pub fn len(&self) -> usize { self.hh.len() }

    /// Return `true` if the number of H-boxes found is zero.
    pub fn is_empty(&self) -> bool { self.hh.is_empty() }

    /// Return a reference to all found H-boxes with their arguments.
    pub fn groups(&self) -> &Vec<(NodeId, C64)> { &self.hh }
}

impl RuleSeal for HStateAll { }
impl RuleFinder for HStateAll {
    type Output<'a> = HStateAllData<'a>;

    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        const EPSILON: f64 = 1e-12;
        const I: C64 = C64 { re: 0.0, im: 1.0 };
        let mut hh: Vec<(NodeId, C64)> = Vec::new();
        for (id, node) in dg.nodes_inner() {
            if node.is_h() && dg.arity(id).unwrap() == 1 {
                let arg = node.arg().unwrap();
                if (arg + 1.0).norm() < EPSILON
                    || (arg - 1.0).norm() < EPSILON
                    || (arg + I).norm() < EPSILON
                    || (arg - I).norm() < EPSILON
                {
                    hh.push((id, arg));
                }
            }
        }
        if hh.is_empty() {
            None
        } else {
            Some(HStateAllData { dg, hh })
        }
    }
}

impl<'a> RuleSeal for HStateAllData<'a> { }
impl<'a> Rule for HStateAllData<'a> {
    fn simplify(self) {
        const EPSILON: f64 = 1e-12;
        const I: C64 = C64 { re: 0.0, im: 1.0 };
        let Self { dg, hh } = self;
        for (h, arg) in hh.into_iter() {
            let ph =
                if (arg + 1.0).norm() < EPSILON {
                    Phase::pi()
                } else if (arg - 1.0).norm() < EPSILON {
                    Phase::zero()
                } else if (arg + I).norm() < EPSILON {
                    -Phase::pi2()
                } else {
                    Phase::pi()
                };
            let n = dg.get_node_mut(h).unwrap();
            *n = Node::Z(ph);
        }
    }
}

