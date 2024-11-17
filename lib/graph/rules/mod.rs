//! Rewrite rules for diagram simplification.
//!
//! Application of particular rules is facilitated through the [`RuleFinder`]
//! and [`Rule`] traits. Usually their methods are not called as such; a
//! particular rule is instead represented by a type implementing `RuleFinder`,
//! which is passed to [`Diagram::find_rule`], [`Diagram::simplify_rule`], or
//! [`Diagram::simplify_rule_n`].
//!
//! Rules are chosen to prefer the elimination of wires and H-boxes, and convert
//! X-spiders to Z-spiders where possible. H-boxes carry an extra scalar factor
//! of 1/√2 when they are binary with argument –1 so that they coincide with the
//! usual definition of the Hadamard gate.
//!
//! See [`SimplifyZX`], [`SimplifyZH`], and [`SimplifyCT`] for master lists of
//! applicable rewrite rules.
//!
//! # Example
//!
//! TODO

use crate::graph::{ NodeId, Diagram, DiagramData, WireData };
pub(crate) use crate::graph::{
    ZX, ZXNode, ZXWire,
    ZH, ZHNode,
    CT, CTNode, CTWire,
};

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
pub trait RuleFinder<A>
where A: DiagramData
{
    /// The type representing the instantiated (but not executed) rewrite rule.
    type Output<'a>: Rule<A> where A: 'a;

    /// Explore a [`Diagram`] to find a particular structure. Returns `None` if
    /// the structure does not exist.
    fn find(self, diagram: &mut Diagram<A>) -> Option<Self::Output<'_>>;
}

/// A trait representing an unexecuted rewrite rule on a [`Diagram`].
///
/// This trait is sealed to ensure that no invalid operations are performed on
/// diagrams; toward this purpose, all types implementing this trait are not
/// constructable as such (see [`RuleFinder`]) and store a mutable reference to
/// the diagram they will simplify. (Further, no type implementing this trait
/// implements [`Clone`] or [`Copy`].)
pub trait Rule<A> {
    /// Execute the rewrite rule, consuming self and releasing the inner hold on
    /// the diagram.
    fn simplify(self);
}

pub(crate) type NodePath<'a, N> = [(NodeId, &'a N)];
pub(crate) type NodeSpec<'a, A, W, N> =
    fn(&Diagram<A>, &NodePath<'a, N>, &'a W, &'a N) -> bool;
pub(crate) type NodeSpec0<'a, A, N> =
    fn(&Diagram<A>, NodeId, &'a N) -> bool;

pub(crate) type ZXNodeSpec0<'a> = NodeSpec0<'a, ZX, ZXNode>;
pub(crate) type ZXNodeSpec<'a> = NodeSpec<'a, ZX, ZXWire, ZXNode>;
pub(crate) type ZHNodeSpec0<'a> = NodeSpec0<'a, ZH, ZHNode>;
pub(crate) type ZHNodeSpec<'a> = NodeSpec<'a, ZH, NodeId, ZHNode>;
pub(crate) type CTNodeSpec0<'a> = NodeSpec0<'a, CT, CTNode>;
pub(crate) type CTNodeSpec<'a> = NodeSpec<'a, CT, CTWire, CTNode>;

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
mod h_absorb;
pub use h_absorb::*;
mod h_absorb_all;
pub use h_absorb_all::*;
mod h_avg;
pub use h_avg::*;
mod h_bialgebra;
pub use h_bialgebra::*;
mod h_euler;
pub use h_euler::*;
mod h_euler_color_flip;
pub use h_euler_color_flip::*;
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
mod h_loop;
pub use h_loop::*;
mod h_loop_all;
pub use h_loop_all::*;
mod h_move;
pub use h_move::*;
mod h_multi_state;
pub use h_multi_state::*;
mod h_multi_state_all;
pub use h_multi_state_all::*;
mod h_self_loop;
pub use h_self_loop::*;
mod h_self_loop_all;
pub use h_self_loop_all::*;
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

/// Master list of rules applicable to [`Diagram`]`<`[`ZX`]`>`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum SimplifyZX {
    Bialgebra,
    ColorFlip,
    ColorFlipAll,
    Fuse,
    FuseAll,
    FuseMulti,
    H2Hopf,
    HEuler,
    HEulerColorFlip,
    HLoop,
    HLoopAll,
    HMove,
    Hopf,
    IState,
    IStateAll,
    Identity,
    IdentityAll,
    PhaseNeg,
    PhaseNegAll,
    PiCommute,
    ScalarPair,
    ScalarPairAll,
    SpiderSelfLoop,
    SpiderSelfLoopAll,
    StateCopy,
    StateCopyAll,
}

/// Output of [`SimplifyZX::find`].
#[derive(Debug)]
pub enum SimplifyZXData<'a> {
    BialgebraData(BialgebraData<'a, ZX>),
    ColorFlipData(ColorFlipData<'a, ZX>),
    ColorFlipAllData(ColorFlipAllData<'a, ZX>),
    FuseData(FuseData<'a, ZX>),
    FuseAllData(FuseAllData<'a, ZX>),
    FuseMultiData(FuseMultiData<'a, ZX>),
    H2HopfData(H2HopfData<'a, ZX>),
    HEulerData(HEulerData<'a, ZX>),
    HEulerColorFlipData(HEulerColorFlipData<'a, ZX>),
    HLoopData(HLoopData<'a, ZX>),
    HLoopAllData(HLoopAllData<'a, ZX>),
    HMoveData(HMoveData<'a, ZX>),
    HopfData(HopfData<'a, ZX>),
    IStateData(IStateData<'a, ZX>),
    IStateAllData(IStateAllData<'a, ZX>),
    IdentityData(IdentityData<'a, ZX>),
    IdentityAllData(IdentityAllData<'a, ZX>),
    PhaseNegData(PhaseNegData<'a, ZX>),
    PhaseNegAllData(PhaseNegAllData<'a, ZX>),
    PiCommuteData(PiCommuteData<'a, ZX>),
    ScalarPairData(ScalarPairData<'a, ZX>),
    ScalarPairAllData(ScalarPairAllData<'a, ZX>),
    SpiderSelfLoopData(SpiderSelfLoopData<'a, ZX>),
    SpiderSelfLoopAllData(SpiderSelfLoopAllData<'a, ZX>),
    StateCopyData(StateCopyData<'a, ZX>),
    StateCopyAllData(StateCopyAllData<'a, ZX>),
}

impl RuleFinder<ZX> for SimplifyZX {
    type Output<'a> = SimplifyZXData<'a>;

    fn find(self, dg: &mut Diagram<ZX>) -> Option<Self::Output<'_>> {
        use SimplifyZXData::*;
        match self {
            Self::Bialgebra => Bialgebra.find(dg).map(BialgebraData),
            Self::ColorFlip => ColorFlip.find(dg).map(ColorFlipData),
            Self::ColorFlipAll => ColorFlipAll.find(dg).map(ColorFlipAllData),
            Self::Fuse => Fuse.find(dg).map(FuseData),
            Self::FuseAll => FuseAll.find(dg).map(FuseAllData),
            Self::FuseMulti => FuseMulti.find(dg).map(FuseMultiData),
            Self::H2Hopf => H2Hopf.find(dg).map(H2HopfData),
            Self::HEuler => HEuler.find(dg).map(HEulerData),
            Self::HEulerColorFlip => HEulerColorFlip.find(dg).map(HEulerColorFlipData),
            Self::HLoop => HLoop.find(dg).map(HLoopData),
            Self::HLoopAll => HLoopAll.find(dg).map(HLoopAllData),
            Self::HMove => HMove.find(dg).map(HMoveData),
            Self::Hopf => Hopf.find(dg).map(HopfData),
            Self::IState => IState.find(dg).map(IStateData),
            Self::IStateAll => IStateAll.find(dg).map(IStateAllData),
            Self::Identity => Identity.find(dg).map(IdentityData),
            Self::IdentityAll => IdentityAll.find(dg).map(IdentityAllData),
            Self::PhaseNeg => PhaseNeg.find(dg).map(PhaseNegData),
            Self::PhaseNegAll => PhaseNegAll.find(dg).map(PhaseNegAllData),
            Self::PiCommute => PiCommute.find(dg).map(PiCommuteData),
            Self::ScalarPair => ScalarPair.find(dg).map(ScalarPairData),
            Self::ScalarPairAll => ScalarPairAll.find(dg).map(ScalarPairAllData),
            Self::SpiderSelfLoop => SpiderSelfLoop.find(dg).map(SpiderSelfLoopData),
            Self::SpiderSelfLoopAll => SpiderSelfLoopAll.find(dg).map(SpiderSelfLoopAllData),
            Self::StateCopy => StateCopy.find(dg).map(StateCopyData),
            Self::StateCopyAll => StateCopyAll.find(dg).map(StateCopyAllData),
        }
    }
}

impl<'a> Rule<ZX> for SimplifyZXData<'a> {
    fn simplify(self) {
        match self {
            Self::BialgebraData(data) => data.simplify(),
            Self::ColorFlipData(data) => data.simplify(),
            Self::ColorFlipAllData(data) => data.simplify(),
            Self::FuseData(data) => data.simplify(),
            Self::FuseAllData(data) => data.simplify(),
            Self::FuseMultiData(data) => data.simplify(),
            Self::H2HopfData(data) => data.simplify(),
            Self::HEulerData(data) => data.simplify(),
            Self::HEulerColorFlipData(data) => data.simplify(),
            Self::HLoopData(data) => data.simplify(),
            Self::HLoopAllData(data) => data.simplify(),
            Self::HMoveData(data) => data.simplify(),
            Self::HopfData(data) => data.simplify(),
            Self::IStateData(data) => data.simplify(),
            Self::IStateAllData(data) => data.simplify(),
            Self::IdentityData(data) => data.simplify(),
            Self::IdentityAllData(data) => data.simplify(),
            Self::PhaseNegData(data) => data.simplify(),
            Self::PhaseNegAllData(data) => data.simplify(),
            Self::PiCommuteData(data) => data.simplify(),
            Self::ScalarPairData(data) => data.simplify(),
            Self::ScalarPairAllData(data) => data.simplify(),
            Self::SpiderSelfLoopData(data) => data.simplify(),
            Self::SpiderSelfLoopAllData(data) => data.simplify(),
            Self::StateCopyData(data) => data.simplify(),
            Self::StateCopyAllData(data) => data.simplify(),
        }
    }
}

/// Master list of rules applicable to [`Diagram`]`<`[`ZH`]`>`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum SimplifyZH {
    Bialgebra,
    ColorFlip,
    ColorFlipAll,
    Fuse,
    FuseAll,
    FuseMulti,
    H2Hopf,
    HEuler,
    HEulerColorFlip,
    HLoop,
    HLoopAll,
    HMove,
    Hopf,
    IState,
    IStateAll,
    Identity,
    IdentityAll,
    PhaseNeg,
    PhaseNegAll,
    PiCommute,
    ScalarPair,
    ScalarPairAll,
    SpiderSelfLoop,
    SpiderSelfLoopAll,
    StateCopy,
    StateCopyAll,
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
    HSelfLoop,
    HSelfLoopAll,
    HState,
    HStateAll,
    HStateCopy,
    HStateCopyAll,
    HStateMul,
    HStateMulAll,
}

/// Output of [`SimplifyZH::find`].
#[derive(Debug)]
pub enum SimplifyZHData<'a> {
    BialgebraData(BialgebraData<'a, ZH>),
    ColorFlipData(ColorFlipData<'a, ZH>),
    ColorFlipAllData(ColorFlipAllData<'a, ZH>),
    FuseData(FuseData<'a, ZH>),
    FuseAllData(FuseAllData<'a, ZH>),
    FuseMultiData(FuseMultiData<'a, ZH>),
    H2HopfData(H2HopfData<'a, ZH>),
    HEulerData(HEulerData<'a, ZH>),
    HEulerColorFlipData(HEulerColorFlipData<'a, ZH>),
    HLoopData(HLoopData<'a, ZH>),
    HLoopAllData(HLoopAllData<'a, ZH>),
    HMoveData(HMoveData<'a, ZH>),
    HopfData(HopfData<'a, ZH>),
    IStateData(IStateData<'a, ZH>),
    IStateAllData(IStateAllData<'a, ZH>),
    IdentityData(IdentityData<'a, ZH>),
    IdentityAllData(IdentityAllData<'a, ZH>),
    PhaseNegData(PhaseNegData<'a, ZH>),
    PhaseNegAllData(PhaseNegAllData<'a, ZH>),
    PiCommuteData(PiCommuteData<'a, ZH>),
    ScalarPairData(ScalarPairData<'a, ZH>),
    ScalarPairAllData(ScalarPairAllData<'a, ZH>),
    SpiderSelfLoopData(SpiderSelfLoopData<'a, ZH>),
    SpiderSelfLoopAllData(SpiderSelfLoopAllData<'a, ZH>),
    StateCopyData(StateCopyData<'a, ZH>),
    StateCopyAllData(StateCopyAllData<'a, ZH>),
    HAbsorbData(HAbsorbData<'a, ZH>),
    HAbsorbAllData(HAbsorbAllData<'a, ZH>),
    HAvgData(HAvgData<'a, ZH>),
    HBialgebraData(HBialgebraData<'a, ZH>),
    HExplodeData(HExplodeData<'a, ZH>),
    HExplodeAllData(HExplodeAllData<'a, ZH>),
    HFuseData(HFuseData<'a, ZH>),
    HHopfData(HHopfData<'a, ZH>),
    HIntroData(HIntroData<'a, ZH>),
    HMultiStateData(HMultiStateData<'a, ZH>),
    HMultiStateAllData(HMultiStateAllData<'a, ZH>),
    HSelfLoopData(HSelfLoopData<'a, ZH>),
    HSelfLoopAllData(HSelfLoopAllData<'a, ZH>),
    HStateData(HStateData<'a, ZH>),
    HStateAllData(HStateAllData<'a, ZH>),
    HStateCopyData(HStateCopyData<'a, ZH>),
    HStateCopyAllData(HStateCopyAllData<'a, ZH>),
    HStateMulData(HStateMulData<'a, ZH>),
    HStateMulAllData(HStateMulAllData<'a, ZH>),
}

impl RuleFinder<ZH> for SimplifyZH {
    type Output<'a> = SimplifyZHData<'a>;

    fn find(self, dg: &mut Diagram<ZH>) -> Option<Self::Output<'_>> {
        use SimplifyZHData::*;
        match self {
            Self::Bialgebra => Bialgebra.find(dg).map(BialgebraData),
            Self::ColorFlip => ColorFlip.find(dg).map(ColorFlipData),
            Self::ColorFlipAll => ColorFlipAll.find(dg).map(ColorFlipAllData),
            Self::Fuse => Fuse.find(dg).map(FuseData),
            Self::FuseAll => FuseAll.find(dg).map(FuseAllData),
            Self::FuseMulti => FuseMulti.find(dg).map(FuseMultiData),
            Self::H2Hopf => H2Hopf.find(dg).map(H2HopfData),
            Self::HEuler => HEuler.find(dg).map(HEulerData),
            Self::HEulerColorFlip => HEulerColorFlip.find(dg).map(HEulerColorFlipData),
            Self::HLoop => HLoop.find(dg).map(HLoopData),
            Self::HLoopAll => HLoopAll.find(dg).map(HLoopAllData),
            Self::HMove => HMove.find(dg).map(HMoveData),
            Self::Hopf => Hopf.find(dg).map(HopfData),
            Self::IState => IState.find(dg).map(IStateData),
            Self::IStateAll => IStateAll.find(dg).map(IStateAllData),
            Self::Identity => Identity.find(dg).map(IdentityData),
            Self::IdentityAll => IdentityAll.find(dg).map(IdentityAllData),
            Self::PhaseNeg => PhaseNeg.find(dg).map(PhaseNegData),
            Self::PhaseNegAll => PhaseNegAll.find(dg).map(PhaseNegAllData),
            Self::PiCommute => PiCommute.find(dg).map(PiCommuteData),
            Self::ScalarPair => ScalarPair.find(dg).map(ScalarPairData),
            Self::ScalarPairAll => ScalarPairAll.find(dg).map(ScalarPairAllData),
            Self::SpiderSelfLoop => SpiderSelfLoop.find(dg).map(SpiderSelfLoopData),
            Self::SpiderSelfLoopAll => SpiderSelfLoopAll.find(dg).map(SpiderSelfLoopAllData),
            Self::StateCopy => StateCopy.find(dg).map(StateCopyData),
            Self::StateCopyAll => StateCopyAll.find(dg).map(StateCopyAllData),
            Self::HAbsorb => HAbsorb.find(dg).map(HAbsorbData),
            Self::HAbsorbAll => HAbsorbAll.find(dg).map(HAbsorbAllData),
            Self::HAvg => HAvg.find(dg).map(HAvgData),
            Self::HBialgebra => HBialgebra.find(dg).map(HBialgebraData),
            Self::HExplode => HExplode.find(dg).map(HExplodeData),
            Self::HExplodeAll => HExplodeAll.find(dg).map(HExplodeAllData),
            Self::HFuse => HFuse.find(dg).map(HFuseData),
            Self::HHopf => HHopf.find(dg).map(HHopfData),
            Self::HIntro => HIntro.find(dg).map(HIntroData),
            Self::HMultiState => HMultiState.find(dg).map(HMultiStateData),
            Self::HMultiStateAll => HMultiStateAll.find(dg).map(HMultiStateAllData),
            Self::HSelfLoop => HSelfLoop.find(dg).map(HSelfLoopData),
            Self::HSelfLoopAll => HSelfLoopAll.find(dg).map(HSelfLoopAllData),
            Self::HState => HState.find(dg).map(HStateData),
            Self::HStateAll => HStateAll.find(dg).map(HStateAllData),
            Self::HStateCopy => HStateCopy.find(dg).map(HStateCopyData),
            Self::HStateCopyAll => HStateCopyAll.find(dg).map(HStateCopyAllData),
            Self::HStateMul => HStateMul.find(dg).map(HStateMulData),
            Self::HStateMulAll => HStateMulAll.find(dg).map(HStateMulAllData),
        }
    }
}

impl<'a> Rule<ZH> for SimplifyZHData<'a> {
    fn simplify(self) {
        match self {
            Self::BialgebraData(data) => data.simplify(),
            Self::ColorFlipData(data) => data.simplify(),
            Self::ColorFlipAllData(data) => data.simplify(),
            Self::FuseData(data) => data.simplify(),
            Self::FuseAllData(data) => data.simplify(),
            Self::FuseMultiData(data) => data.simplify(),
            Self::H2HopfData(data) => data.simplify(),
            Self::HEulerData(data) => data.simplify(),
            Self::HEulerColorFlipData(data) => data.simplify(),
            Self::HLoopData(data) => data.simplify(),
            Self::HLoopAllData(data) => data.simplify(),
            Self::HMoveData(data) => data.simplify(),
            Self::HopfData(data) => data.simplify(),
            Self::IStateData(data) => data.simplify(),
            Self::IStateAllData(data) => data.simplify(),
            Self::IdentityData(data) => data.simplify(),
            Self::IdentityAllData(data) => data.simplify(),
            Self::PhaseNegData(data) => data.simplify(),
            Self::PhaseNegAllData(data) => data.simplify(),
            Self::PiCommuteData(data) => data.simplify(),
            Self::ScalarPairData(data) => data.simplify(),
            Self::ScalarPairAllData(data) => data.simplify(),
            Self::SpiderSelfLoopData(data) => data.simplify(),
            Self::SpiderSelfLoopAllData(data) => data.simplify(),
            Self::StateCopyData(data) => data.simplify(),
            Self::StateCopyAllData(data) => data.simplify(),
            Self::HAbsorbData(data) => data.simplify(),
            Self::HAbsorbAllData(data) => data.simplify(),
            Self::HAvgData(data) => data.simplify(),
            Self::HBialgebraData(data) => data.simplify(),
            Self::HExplodeData(data) => data.simplify(),
            Self::HExplodeAllData(data) => data.simplify(),
            Self::HFuseData(data) => data.simplify(),
            Self::HHopfData(data) => data.simplify(),
            Self::HIntroData(data) => data.simplify(),
            Self::HMultiStateData(data) => data.simplify(),
            Self::HMultiStateAllData(data) => data.simplify(),
            Self::HSelfLoopData(data) => data.simplify(),
            Self::HSelfLoopAllData(data) => data.simplify(),
            Self::HStateData(data) => data.simplify(),
            Self::HStateAllData(data) => data.simplify(),
            Self::HStateCopyData(data) => data.simplify(),
            Self::HStateCopyAllData(data) => data.simplify(),
            Self::HStateMulData(data) => data.simplify(),
            Self::HStateMulAllData(data) => data.simplify(),
        }
    }
}

/// Master list of rules applicable to [`Diagram`]`<`[`CT`]`>`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum SimplifyCT {
}

/// Output of [`SimplifyCT::find`].
#[derive(Debug)]
pub enum SimplifyCTData<'a> {
    PlaceHolder(std::marker::PhantomData<&'a ()>)
}

#[allow(unused_variables, unused_mut)]
impl RuleFinder<CT> for SimplifyCT {
    type Output<'a> = SimplifyCTData<'a>;

    fn find(self, dg: &mut Diagram<CT>) -> Option<Self::Output<'_>> {
        todo!()
        // use SimplifyZHData::*;
        // match self {
        // }
    }
}

impl<'a> Rule<CT> for SimplifyCTData<'a> {
    fn simplify(self) {
        todo!()
        // match self {
        // }
    }
}

impl<A> Diagram<A>
where A: DiagramData
{
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
        init: NodeSpec0<'a, A, A::Node>,
        rest: &[NodeSpec<'a, A, A::Wire, A::Node>],
    ) -> Option<Vec<(NodeId, &'a A::Node)>>
    where I: Iterator<Item = (NodeId, &'a A::Node)>
    {
        fn find_path_inner<'b, B, I>(
            dg: &'b Diagram<B>,
            nodes: I,
            specs: &[NodeSpec<'b, B, B::Wire, B::Node>],
            acc: &mut Vec<(NodeId, &'b B::Node)>,
        ) -> bool
        where
            B: DiagramData,
            I: Iterator<Item = (&'b B::Wire, &'b B::Node)>
        {
            if let Some(spec) = specs.first() {
                for (wire, node) in nodes.into_iter() {
                    if acc.iter().any(|(path_id, _)| *path_id == wire.id()) {
                        continue;
                    }
                    if spec(dg, acc.as_ref(), wire, node) {
                        acc.push((wire.id(), node));
                        let save_n = acc.len();
                        let nnb = dg.neighbors(wire.id()).unwrap();
                        let rec = find_path_inner(dg, nnb, &specs[1..], acc);
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

        let mut acc: Vec<(NodeId, &A::Node)> =
            Vec::with_capacity(rest.len() + 1);
        for (id, node) in nodes.into_iter() {
            if init(self, id, node) {
                acc.push((id, node));
                let nnb = self.neighbors(id).unwrap();
                if find_path_inner(self, nnb, rest, &mut acc) {
                    return Some(acc);
                }
                acc.clear();
            }
        }
        None
    }

    /// Shortcut to calling [`RuleFinder::find`] on `self`, returning any 
    /// output.
    pub fn find_rule<R>(&mut self, rule: R) -> Option<R::Output<'_>>
    where R: RuleFinder<A>
    {
        rule.find(self)
    }

    /// Find and immediately apply a rewrite rule, returning `true` if the rule
    /// was successfully applied.
    pub fn simplify_rule<R>(&mut self, rule: R) -> bool
    where R: RuleFinder<A>
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
    where R: RuleFinder<A> + Clone
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

