//! Rewrite rules for diagram simplification.
//!
//! Application of particular rules is facilitated through the [`RuleFinder`]
//! and [`Rule`] traits. Usually their methods are not called as such; a
//! particular rule is instead represented by a type implementing `RuleFinder`,
//! which is passed to [`Diagram::find_rule`], [`Diagram::simplify_rule`], or
//! [`Diagram::simplify_rule_n`].
//!
//! Rules are chosen to prefer the elimination of wires and H-boxes, and convert
//! X-spiders to Z-spiders where possible.
//!
//! See [`Simplify`] for a master list of all simplify rules combined into a
//! single type.
//!
//! # Example
//!
//! TODO

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
    /// [`Bialgebra`]
    Bialgebra,
    /// [`ColorFlip`]
    ColorFlip,
    /// [`ColorFlipAll`]
    ColorFlipAll,
    /// [`Fuse`]
    Fuse,
    /// [`FuseAll`]
    FuseAll,
    /// [`FuseMulti`]
    FuseMulti,
    /// [`H2Hopf`]
    H2Hopf,
    /// [`HEuler`]
    HEuler,
    /// [`HEulerColorFlip`]
    HEulerColorFlip,
    /// [`HLoop`]
    HLoop,
    /// [`HMove`]
    HMove,
    /// [`HSelfLoop`]
    HSelfLoop,
    /// [`HSelfLoopAll`]
    HSelfLoopAll,
    /// [`Hopf`]
    Hopf,
    /// [`IState`]
    IState,
    /// [`IStateAll`]
    IStateAll,
    /// [`Identity`]
    Identity,
    /// [`IdentityAll`]
    IdentityAll,
    /// [`PhaseNeg`]
    PhaseNeg,
    /// [`PhaseNegAll`]
    PhaseNegAll,
    /// [`PiCommute`]
    PiCommute,
    /// [`ScalarPair`]
    ScalarPair,
    /// [`ScalarPairAll`]
    ScalarPairAll,
    /// [`SpiderSelfLoop`]
    SpiderSelfLoop,
    /// [`SpiderSelfLoopAll`]
    SpiderSelfLoopAll,
    /// [`StateCopy`]
    StateCopy,
    /// [`StateCopyAll`]
    StateCopyAll,

    // ZH
    /// [`HAbsorb`]
    HAbsorb,
    /// [`HAbsorbAll`]
    HAbsorbAll,
    /// [`HAvg`]
    HAvg,
    /// [`HBialgebra`]
    HBialgebra,
    /// [`HExplode`]
    HExplode,
    /// [`HExplodeAll`]
    HExplodeAll,
    /// [`HFuse`]
    HFuse,
    /// [`HHopf`]
    HHopf,
    /// [`HIntro`]
    HIntro,
    /// [`HMultiState`]
    HMultiState,
    /// [`HMultiStateAll`]
    HMultiStateAll,
    /// [`HState`]
    HState,
    /// [`HStateAll`]
    HStateAll,
    /// [`HStateCopy`]
    HStateCopy,
    /// [`HStateCopyAll`]
    HStateCopyAll,
    /// [`HStateMul`]
    HStateMul,
    /// [`HStateMulAll`]
    HStateMulAll,
}

impl private::RuleSeal for Simplify { }
impl RuleFinder for Simplify {
    type Output<'a> = SimplifyData<'a>;

    #[allow(unused_variables, unused_mut)]
    fn find(self, dg: &mut Diagram) -> Option<Self::Output<'_>> {
        match self {
            Self::Bialgebra         => Some(SimplifyData::Bialgebra(Bialgebra.find(dg)?)),
            Self::ColorFlip         => Some(SimplifyData::ColorFlip(ColorFlip.find(dg)?)),
            Self::ColorFlipAll      => Some(SimplifyData::ColorFlipAll(ColorFlipAll.find(dg)?)),
            Self::Fuse              => Some(SimplifyData::Fuse(Fuse.find(dg)?)),
            Self::FuseAll           => Some(SimplifyData::FuseAll(FuseAll.find(dg)?)),
            Self::FuseMulti         => Some(SimplifyData::FuseMulti(FuseMulti.find(dg)?)),
            Self::H2Hopf            => Some(SimplifyData::H2Hopf(H2Hopf.find(dg)?)),
            Self::HEuler            => Some(SimplifyData::HEuler(HEuler.find(dg)?)),
            Self::HEulerColorFlip   => Some(SimplifyData::HEulerColorFlip(HEulerColorFlip.find(dg)?)),
            Self::HLoop             => Some(SimplifyData::HLoop(HLoop.find(dg)?)),
            Self::HMove             => Some(SimplifyData::HMove(HMove.find(dg)?)),
            Self::HSelfLoop         => Some(SimplifyData::HSelfLoop(HSelfLoop.find(dg)?)),
            Self::HSelfLoopAll      => Some(SimplifyData::HSelfLoopAll(HSelfLoopAll.find(dg)?)),
            Self::Hopf              => Some(SimplifyData::Hopf(Hopf.find(dg)?)),
            Self::IState            => Some(SimplifyData::IState(IState.find(dg)?)),
            Self::IStateAll         => Some(SimplifyData::IStateAll(IStateAll.find(dg)?)),
            Self::Identity          => Some(SimplifyData::Identity(Identity.find(dg)?)),
            Self::IdentityAll       => Some(SimplifyData::IdentityAll(IdentityAll.find(dg)?)),
            Self::PhaseNeg          => Some(SimplifyData::PhaseNeg(PhaseNeg.find(dg)?)),
            Self::PhaseNegAll       => Some(SimplifyData::PhaseNegAll(PhaseNegAll.find(dg)?)),
            Self::PiCommute         => Some(SimplifyData::PiCommute(PiCommute.find(dg)?)),
            Self::ScalarPair        => Some(SimplifyData::ScalarPair(ScalarPair.find(dg)?)),
            Self::ScalarPairAll     => Some(SimplifyData::ScalarPairAll(ScalarPairAll.find(dg)?)),
            Self::SpiderSelfLoop    => Some(SimplifyData::SpiderSelfLoop(SpiderSelfLoop.find(dg)?)),
            Self::SpiderSelfLoopAll => Some(SimplifyData::SpiderSelfLoopAll(SpiderSelfLoopAll.find(dg)?)),
            Self::StateCopy         => Some(SimplifyData::StateCopy(StateCopy.find(dg)?)),
            Self::StateCopyAll      => Some(SimplifyData::StateCopyAll(StateCopyAll.find(dg)?)),
            Self::HAbsorb           => Some(SimplifyData::HAbsorb(HAbsorb.find(dg)?)),
            Self::HAbsorbAll        => Some(SimplifyData::HAbsorbAll(HAbsorbAll.find(dg)?)),
            Self::HAvg              => Some(SimplifyData::HAvg(HAvg.find(dg)?)),
            Self::HBialgebra        => Some(SimplifyData::HBialgebra(HBialgebra.find(dg)?)),
            Self::HExplode          => Some(SimplifyData::HExplode(HExplode.find(dg)?)),
            Self::HExplodeAll       => Some(SimplifyData::HExplodeAll(HExplodeAll.find(dg)?)),
            Self::HFuse             => Some(SimplifyData::HFuse(HFuse.find(dg)?)),
            Self::HHopf             => Some(SimplifyData::HHopf(HHopf.find(dg)?)),
            Self::HIntro            => Some(SimplifyData::HIntro(HIntro.find(dg)?)),
            Self::HMultiState       => Some(SimplifyData::HMultiState(HMultiState.find(dg)?)),
            Self::HMultiStateAll    => Some(SimplifyData::HMultiStateAll(HMultiStateAll.find(dg)?)),
            Self::HState            => Some(SimplifyData::HState(HState.find(dg)?)),
            Self::HStateAll         => Some(SimplifyData::HStateAll(HStateAll.find(dg)?)),
            Self::HStateCopy        => Some(SimplifyData::HStateCopy(HStateCopy.find(dg)?)),
            Self::HStateCopyAll     => Some(SimplifyData::HStateCopyAll(HStateCopyAll.find(dg)?)),
            Self::HStateMul         => Some(SimplifyData::HStateMul(HStateMul.find(dg)?)),
            Self::HStateMulAll      => Some(SimplifyData::HStateMulAll(HStateMulAll.find(dg)?)),
        }
    }
}

/// Output [`Simplify::find`].
#[derive(Debug)]
pub enum SimplifyData<'a> {
    // ZX
    /// [`BialgebraData`]
    Bialgebra(BialgebraData<'a>),
    /// [`ColorFlipData`]
    ColorFlip(ColorFlipData<'a>),
    /// [`ColorFlipAllData`]
    ColorFlipAll(ColorFlipAllData<'a>),
    /// [`FuseData`]
    Fuse(FuseData<'a>),
    /// [`FuseAllData`]
    FuseAll(FuseAllData<'a>),
    /// [`FuseMultiData`]
    FuseMulti(FuseMultiData<'a>),
    /// [`H2HopfData`]
    H2Hopf(H2HopfData<'a>),
    /// [`HEulerData`]
    HEuler(HEulerData<'a>),
    /// [`HEulerColorFlipData`]
    HEulerColorFlip(HEulerColorFlipData<'a>),
    /// [`HLoopData`]
    HLoop(HLoopData<'a>),
    /// [`HMoveData`]
    HMove(HMoveData<'a>),
    /// [`HSelfLoopData`]
    HSelfLoop(HSelfLoopData<'a>),
    /// [`HSelfLoopAllData`]
    HSelfLoopAll(HSelfLoopAllData<'a>),
    /// [`HopfData`]
    Hopf(HopfData<'a>),
    /// [`IStateData`]
    IState(IStateData<'a>),
    /// [`IStateAllData`]
    IStateAll(IStateAllData<'a>),
    /// [`IdentityData`]
    Identity(IdentityData<'a>),
    /// [`IdentityAllData`]
    IdentityAll(IdentityAllData<'a>),
    /// [`PhaseNegData`]
    PhaseNeg(PhaseNegData<'a>),
    /// [`PhaseNegAllData`]
    PhaseNegAll(PhaseNegAllData<'a>),
    /// [`PiCommuteData`]
    PiCommute(PiCommuteData<'a>),
    /// [`ScalarPairData`]
    ScalarPair(ScalarPairData<'a>),
    /// [`ScalarPairAllData`]
    ScalarPairAll(ScalarPairAllData<'a>),
    /// [`SpiderSelfLoopData`]
    SpiderSelfLoop(SpiderSelfLoopData<'a>),
    /// [`SpiderSelfLoopAllData`]
    SpiderSelfLoopAll(SpiderSelfLoopAllData<'a>),
    /// [`StateCopyData`]
    StateCopy(StateCopyData<'a>),
    /// [`StateCopyAllData`]
    StateCopyAll(StateCopyAllData<'a>),

    // ZH
    /// [`HAbsorbData`]
    HAbsorb(HAbsorbData<'a>),
    /// [`HAbsorbAllData`]
    HAbsorbAll(HAbsorbAllData<'a>),
    /// [`HAvgData`]
    HAvg(HAvgData<'a>),
    /// [`HBialgebraData`]
    HBialgebra(HBialgebraData<'a>),
    /// [`HExplodeData`]
    HExplode(HExplodeData<'a>),
    /// [`HExplodeAllData`]
    HExplodeAll(HExplodeAllData<'a>),
    /// [`HFuseData`]
    HFuse(HFuseData<'a>),
    /// [`HHopfData`]
    HHopf(HHopfData<'a>),
    /// [`HIntroData`]
    HIntro(HIntroData<'a>),
    /// [`HMultiStateData`]
    HMultiState(HMultiStateData<'a>),
    /// [`HMultiStateAllData`]
    HMultiStateAll(HMultiStateAllData<'a>),
    /// [`HStateData`]
    HState(HStateData<'a>),
    /// [`HStateAllData`]
    HStateAll(HStateAllData<'a>),
    /// [`HStateCopyData`]
    HStateCopy(HStateCopyData<'a>),
    /// [`HStateCopyAllData`]
    HStateCopyAll(HStateCopyAllData<'a>),
    /// [`HStateMulData`]
    HStateMul(HStateMulData<'a>),
    /// [`HStateMulAllData`]
    HStateMulAll(HStateMulAllData<'a>),
}

impl<'a> private::RuleSeal for SimplifyData<'a> { }
impl<'a> Rule for SimplifyData<'a> {
    fn simplify(self) {
        match self {
            Self::Bialgebra(data)         => data.simplify(),
            Self::ColorFlip(data)         => data.simplify(),
            Self::ColorFlipAll(data)      => data.simplify(),
            Self::Fuse(data)              => data.simplify(),
            Self::FuseAll(data)           => data.simplify(),
            Self::FuseMulti(data)         => data.simplify(),
            Self::H2Hopf(data)            => data.simplify(),
            Self::HEuler(data)            => data.simplify(),
            Self::HEulerColorFlip(data)   => data.simplify(),
            Self::HLoop(data)             => data.simplify(),
            Self::HMove(data)             => data.simplify(),
            Self::HSelfLoop(data)         => data.simplify(),
            Self::HSelfLoopAll(data)      => data.simplify(),
            Self::Hopf(data)              => data.simplify(),
            Self::IState(data)            => data.simplify(),
            Self::IStateAll(data)         => data.simplify(),
            Self::Identity(data)          => data.simplify(),
            Self::IdentityAll(data)       => data.simplify(),
            Self::PhaseNeg(data)          => data.simplify(),
            Self::PhaseNegAll(data)       => data.simplify(),
            Self::PiCommute(data)         => data.simplify(),
            Self::ScalarPair(data)        => data.simplify(),
            Self::ScalarPairAll(data)     => data.simplify(),
            Self::SpiderSelfLoop(data)    => data.simplify(),
            Self::SpiderSelfLoopAll(data) => data.simplify(),
            Self::StateCopy(data)         => data.simplify(),
            Self::StateCopyAll(data)      => data.simplify(),
            Self::HAbsorb(data)           => data.simplify(),
            Self::HAbsorbAll(data)        => data.simplify(),
            Self::HAvg(data)              => data.simplify(),
            Self::HBialgebra(data)        => data.simplify(),
            Self::HExplode(data)          => data.simplify(),
            Self::HExplodeAll(data)       => data.simplify(),
            Self::HFuse(data)             => data.simplify(),
            Self::HHopf(data)             => data.simplify(),
            Self::HIntro(data)            => data.simplify(),
            Self::HMultiState(data)       => data.simplify(),
            Self::HMultiStateAll(data)    => data.simplify(),
            Self::HState(data)            => data.simplify(),
            Self::HStateAll(data)         => data.simplify(),
            Self::HStateCopy(data)        => data.simplify(),
            Self::HStateCopyAll(data)     => data.simplify(),
            Self::HStateMul(data)         => data.simplify(),
            Self::HStateMulAll(data)      => data.simplify(),
        }
    }
}

// ZX
mod bialgebra;
pub use bialgebra::*;
mod color_flip;
pub use color_flip::*;
mod color_flip_all;
pub use color_flip_all::*;
mod fuse;
pub use fuse::*;
mod fuse_all;
pub use fuse_all::*;
mod fuse_multi;
pub use fuse_multi::*;
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
mod i_state;
pub use i_state::*;
mod i_state_all;
pub use i_state_all::*;
mod identity;
pub use identity::*;
mod identity_all;
pub use identity_all::*;
mod phase_neg;
pub use phase_neg::*;
mod phase_neg_all;
pub use phase_neg_all::*;
mod pi_commute;
pub use pi_commute::*;
mod scalar_pair;
pub use scalar_pair::*;
mod scalar_pair_all;
pub use scalar_pair_all::*;
mod spider_self_loop;
pub use spider_self_loop::*;
mod spider_self_loop_all;
pub use spider_self_loop_all::*;
mod state_copy;
pub use state_copy::*;
mod state_copy_all;
pub use state_copy_all::*;

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

