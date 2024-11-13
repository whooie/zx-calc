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
//! # Example
//!
//! TODO

use crate::graph2::{ NodeId, Diagram, DiagramData, WireData };
pub(crate) use crate::graph2::{
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
pub(crate) type NodeSpec<'a, A, W, N> = fn(&Diagram<A>, &NodePath<'a, N>, &'a W, &'a N) -> bool;
pub(crate) type NodeSpec0<'a, A, N> = fn(&Diagram<A>, NodeId, &'a N) -> bool;
pub(crate) type ZXNodeSpec0<'a> = NodeSpec0<'a, ZX, ZXNode>;
pub(crate) type ZXNodeSpec<'a> = NodeSpec<'a, ZX, ZXWire, ZXNode>;
pub(crate) type ZHNodeSpec0<'a> = NodeSpec0<'a, ZH, ZHNode>;
pub(crate) type ZHNodeSpec<'a> = NodeSpec<'a, ZH, NodeId, ZHNode>;
pub(crate) type CTNodeSpec0<'a> = NodeSpec0<'a, CT, CTNode>;
pub(crate) type CTNodeSpec<'a> = NodeSpec<'a, CT, CTWire, CTNode>;

mod bialgebra;
pub use bialgebra::*;

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

