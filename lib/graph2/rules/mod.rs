#![allow(unused_imports)]

use crate::graph2::{ NodeId, Diagram, node::Node };

pub(crate) fn fst<T, U>(pair: (T, U)) -> T { pair.0 }

pub(crate) fn snd<T, U>(pair: (T, U)) -> U { pair.1 }

pub(crate) fn rev<T, U>(pair: (T, U)) -> (U, T) { (pair.1, pair.0) }

/// `Rule` and `RuleFinder` are sealed traits to govern "rule types" -- we want
/// to expose a modular, genericizable per-rule API, but restrict creation and
/// use of the rules so that no invalid operations are performed on diagrams.
pub(crate) mod private { pub trait RuleSeal { } }
pub(crate) use private::RuleSeal;

/// A trait for types that can explore a [`Diagram`] and find a particular
/// structure to simplify.
///
/// Types implementing this trait (usually unit structs) have the only purpose
/// of denoting the existence of a particular rewrite rule. This is to allow for
/// programmatic manipulation of rewrite rules.
///
/// This trait is sealed to ensure that no invalid operations are performed on
/// diagrams; toward this purpose, all types created by this trait via
/// [`find`][RuleFinder::find] (i.e. [`Self::Output`]) store a mutable reference
/// to the diagram they will simplify.
pub trait RuleFinder: private::RuleSeal {
    /// The type representing the instantiated (but not executed) rewrite rule.
    type Output<'a>: Rule;

    /// Explore a [`Diagram`] to find a particular structure. Returns `None` if
    /// the structure does not exist.
    fn find(self, diagram: &mut Diagram) -> Option<Self::Output<'_>>;
}

/// A trait representing an unexecuted rewrite rule on a [`Diagram`].
///
/// This trait is sealed to ensure that no invalid operations are performed on
/// diagrams; toward this purpose, all types implementing this trait are not
/// constructable as such (see [`RuleFinder`]) and store a mutable reference to
/// the diagram they will simplify. (Further, no type implementing this trait
/// implements [`Clone`] or [`Copy`].)
pub trait Rule: private::RuleSeal {
    /// Execute the rewrite rule, consuming self and releasing the inner hold on
    /// the diagram.
    fn simplify(self);
}

/// Master list of available rewrite rules.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Simplify {
    // ZX
    Bialgebra,
    BitBialgebra,
    ColorFlip,
    ColorFlipAll,
    Fuse,
    FuseAll,
    H2Hopf,
    HEuler,
    HEulerColorFlip,
    HLoop,
    HMove,
    HSelfLoop,
    HSelfLoopAll,
    Hopf,
    Identity,
    IdentityAll,
    Istate,
    PhaseNeg,
    PhaseNegAll,
    PiBookend,
    PiCommute,
    SpiderSelfLoop,
    SpiderSelfLoopAll,
    StateCopy,

    // ZH
    HAbsorb,
    HAbsorbAll,
    HAvg,
    HBialgebra,
    HExplode,
    HExplodeAll,
    HFuse,
    HHopf,
    HIntro,
    HMultiState,
    HMultiStateAll,
    HState,
    HStateAll,
    HStateCopy,
    HStateCopyAll,
    HStateMul,
    HStateMulAll,
}

impl private::RuleSeal for Simplify { }
impl RuleFinder for Simplify {
    type Output<'a> = SimplifyData<'a>;

    #[allow(unused_variables, unused_mut)]
    fn find(self, diagram: &mut Diagram) -> Option<Self::Output<'_>> {
        todo!()
    }
}

/// Output [`Simplify::find`].
#[derive(Debug)]
pub enum SimplifyData<'a> {
    Identity(IdentityData<'a>),
}

impl<'a> private::RuleSeal for SimplifyData<'a> { }
impl<'a> Rule for SimplifyData<'a> {
    fn simplify(self) {
        todo!()
    }
}

// ZX
mod bialgebra;
pub use bialgebra::*;
mod bit_bialgebra;
pub use bit_bialgebra::*;
mod color_flip;
pub use color_flip::*;
mod color_flip_all;
pub use color_flip_all::*;
mod fuse;
pub use fuse::*;
mod fuse_all;
pub use fuse_all::*;
mod h2_hopf;
pub use h2_hopf::*;
mod h_euler;
pub use h_euler::*;
mod h_euler_color_flip;
pub use h_euler_color_flip::*;
mod h_loop;
pub use h_loop::*;
mod h_move;
pub use h_move::*;
mod h_self_loop;
pub use h_self_loop::*;
mod h_self_loop_all;
pub use h_self_loop_all::*;
mod hopf;
pub use hopf::*;
mod identity;
pub use identity::*;
mod identity_all;
pub use identity_all::*;
mod istate;
pub use istate::*;
mod phase_neg;
pub use phase_neg::*;
mod phase_neg_all;
pub use phase_neg_all::*;
mod pi_bookend;
pub use pi_bookend::*;
mod pi_commute;
pub use pi_commute::*;
mod spider_self_loop;
pub use spider_self_loop::*;
mod spider_self_loop_all;
pub use spider_self_loop_all::*;
mod state_copy;
pub use state_copy::*;

// ZH
mod h_absorb;
pub use h_absorb::*;
mod h_absorb_all;
pub use h_absorb_all::*;
mod h_avg;
pub use h_avg::*;
mod h_bialgebra;
pub use h_bialgebra::*;
mod h_explode;
pub use h_explode::*;
mod h_explode_all;
pub use h_explode_all::*;
mod h_fuse;
pub use h_fuse::*;
mod h_hopf;
pub use h_hopf::*;
mod h_intro;
pub use h_intro::*;
mod h_multi_state;
pub use h_multi_state::*;
mod h_multi_state_all;
pub use h_multi_state_all::*;
mod h_state;
pub use h_state::*;
mod h_state_all;
pub use h_state_all::*;
mod h_state_copy;
pub use h_state_copy::*;
mod h_state_copy_all;
pub use h_state_copy_all::*;
mod h_state_mul;
pub use h_state_mul::*;
mod h_state_mul_all;
pub use h_state_mul_all::*;

type NodePath<'a> = [(NodeId, &'a Node)];
type PathNodeSpec<'a> = fn(&Diagram, &NodePath<'a>, NodeId, &Node) -> bool;

impl Diagram {
    // Recursively find a finite sequence of nodes in the diagram that satisfy a
    // sequence of predicates, given a locally defined set of nodes to search
    // over at each point in the path. The initial call to this function accepts
    // any iterator as the set of nodes to search, but recursive calls will only
    // use the set of neighbors of the last matching path node. Nodes with IDs
    // matching any of this already in the path are automatically skipped. Any
    // matching path can only have as many nodes as the number of predicates
    // provided. Nodes matching a predicate are pushed onto the accumulator,
    // which is returned as `Some` if all predicates are matched, otherwise
    // `None` is returned instead.
    pub(crate) fn find_path<'a, I>(
        &'a self,
        nodes: I,
        specs: &[PathNodeSpec<'a>],
    ) -> Option<Vec<(NodeId, &'a Node)>>
    where I: Iterator<Item = (NodeId, &'a Node)>
    {
        pub(crate) fn find_path_inner<'b, I>(
            dg: &'b Diagram,
            nodes: I,
            specs: &[PathNodeSpec<'b>],
            acc: &mut Vec<(NodeId, &'b Node)>,
        ) -> bool
        where I: Iterator<Item = (NodeId, &'b Node)>
        {
            if let Some(spec) = specs.first() {
                for (id, node) in nodes {
                    if acc.iter().any(|(path_id, _)| *path_id == id) {
                        continue;
                    }
                    if spec(dg, acc.as_ref(), id, node) {
                        acc.push((id, node));
                        let save_n = acc.len();
                        let rec =
                            find_path_inner(
                                dg,
                                dg.neighbors_of(id).unwrap(),
                                &specs[1..],
                                acc,
                            );
                        if rec {
                            return true;
                        } else {
                            acc.truncate(save_n);
                            continue;
                        }
                    }
                }
                false
            } else {
                true
            }
        }

        let mut acc: Vec<(NodeId, &Node)> = Vec::with_capacity(specs.len());
        find_path_inner(self, nodes, specs, &mut acc)
            .then_some(acc)
    }

    /// Shortcut to calling [`RuleFinder::find`] on `self`, returning any output.
    pub fn find_rule<R>(&mut self, rule: R) -> Option<R::Output<'_>>
    where R: RuleFinder
    {
        rule.find(self)
    }

    /// Find and immediately apply a rewrite rule, returning `true` if the rule
    /// was successfully applied.
    pub fn simplify_rule<R>(&mut self, rule: R) -> bool
    where R: RuleFinder
    {
        rule.find(self)
            .map(|rule_data| { rule_data.simplify(); })
            .is_some()
    }

    /// Repeatedly apply a rewrite rule a maximum of `max` times, returning the
    /// number of times the rule was actually applied.
    pub fn simplify_rule_n<R>(
        &mut self,
        rule: R,
        max: Option<usize>,
    ) -> usize
    where R: RuleFinder + Clone
    {
        if let Some(n) = max {
            (0..n)
                .take_while(|_| self.simplify_rule(rule.clone()))
                .count()
        } else {
            (0..)
                .take_while(|_| self.simplify_rule(rule.clone()))
                .count()
        }
    }
}

