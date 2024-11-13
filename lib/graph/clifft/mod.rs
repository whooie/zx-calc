//! Specialization of items in [`graph`][crate::graph] to Clifford+*T* quantum
//! circuits.
//!
#![cfg_attr(
    feature = "doc-images",
    cfg_attr(
        all(),
        doc = ::embed_doc_image::embed_image!("hadamard_wire", "assets/hadamard_wire.svg"),
        doc = ::embed_doc_image::embed_image!("star_def", "assets/star_def.svg"),
        doc = ::embed_doc_image::embed_image!("star_wire", "assets/star_wire.svg"),
        doc = ::embed_doc_image::embed_image!("auto_wire_rules", "assets/auto_wire_rules.svg"),
    )
)]
#![cfg_attr(
    not(feature = "doc-images"),
    doc = "**Doc images are not enabled**. Compile with feature `doc-images` and Rust version >= 1.54 to enable."
)]
//!
//! Here, we restrict to only the "graph-like" ZX-calculus (i.e. only Z-spiders
//! and empty, binary H-boxes) with phases equal to an integer multiple of
//! *π*/4. Additionally, this module's tools are implemented with a greater
//! focus on performance than explanation, as with the rest of this crate.
//! Namely, the H-box is eliminated as an explicit generator, being replaced by
//! the "Hadamard wire", which is drawn as a dashed blue line.
//!
//! ![hadamard_wire][hadamard_wire]
//!
//! Further, the star generator is introduced to enable concise representations
//! of multi-control gates
//!
//! ![star_def][star_def]
//!
//! which is, like H-boxes, then replaced with a "star wire" as a dashed orange
//! line.
//!
//! ![star_wire][star_wire]
//!
//! Finally, this diagram automatically implements Hopf- and self loop-like
//! rules at graph construction time, rather than as explicit rewrite rules:
//!
//! ![auto_wire_rules][auto_wire_rules]
//!
//! One will likely want to deal with this diagram implementation mostly through
//! the high-level [`Circuit`][crate::circuit::CircuitDiagram] interface, since
//! many common circuit patterns become obfuscated when expressed under these
//! constraints.
//!
//! For more information on the above, see the paper whose methods are
//! implemented in this module, [arXiv:2307.01803][speedy-contraction].
//!
//! [speedy-contraction]: https://arxiv.org/abs/2307.01803

pub mod complex;
pub mod tphase;
pub mod diagram;

mod tnode;
pub(crate) use tnode::*;

