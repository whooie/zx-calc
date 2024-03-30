# zx-calc

*WARNING: This project is in a very early state! Main tools are usable (see
the list of features below), but may be buggy and/or error-prone.*

Contains tools for working in the ZX-calculus, a diagrammatic language for
reasoning about linear maps between two-level quantum states.

## Features
- Multiple diagram representations
    - [x] Ranked, tensorless ket-bra forms
    - [x] Graph-based forms
    - [x] Render to graphviz
    - [x] Conversion between representations
- Diagram tools
    - [x] Scalar computation
    - [x] Evaluation (based on conversion to ketbra form)
    - [ ] Diagram composition/tensoring
- Diagram simplification (graph-based only)
    - [ ] ZX rewrite rules
        - [x] Z/X/HH identity removal
        - [x] Spider fusion
        - [x] Hopf rule
        - [x] "H2-Hopf" rule (identical colors with intervening H-boxes)
        - [x] Spider color change
        - [x] Remove self-loops with single H-boxes
        - [x] Euler-angle expansion of an H-box between two spiders of the same
          color
        - [x] Color-change and fuse for Z-H-X sandwiches
        - [x] π-commute and fuse for Z-X-Z/X-Z-X sandwiches
        - [x] Heuristic-based general reduction of H-boxes via color change
        - [x] Heuristic-based general reduction of π-spiders via π-commute
        - [ ] Bialgebra rule (partial: 2x2 case with nπ phases)
        - [x] State/effect copying
        - [x] Equivalence of spider-states with ±π/2 phases
    - [ ] ZH rewrite rules
        - [x] H-box fusion
        - [x] π-state/effect absorption
        - [x] State/effect explosion through an H-box
        - [x] H-state/effect conversion to Z-state/effect
        - [x] π Z-state/effect copying through H-boxes
        - [x] State/effect expansion of H-boxes with label 1
        - [x] H-box version of the Hopf rule
        - [ ] H-box version of the bialgebra rule
        - [x] H-box averaging rule
        - [x] H-state multiplication rule
        - [ ] Generalized H-box multiplication rule
        - [x] H-box wire introduction rule (reversed)
    - [ ] Single, exhaustive "simplify" (partial: waiting on final rewrite rules
      and testing with different strategies)
    - [x] Scalar computation and removal
    - [ ] Equality testing via simplification
- Circuits
    - [ ] Basic ranked circuit representation
    - [ ] Read/write OpenQASM
    - [ ] Render to Graphviz
    - [ ] Conversion to diagram representations
    - [ ] Circuit extraction from diagrams
- [ ] Time/space optimization

## Helpful resources
* B. Coecke, "Basic ZX-calculus for students and professionals."
  [arXiv:2303.03163](https://arxiv.org/abs/2303.03163)
* J. van de Wetering, "ZX-calculus for the working quantum computer scientist."
  [arXiv:2012.13966](https://arxiv.org/abs/2012.13966)
* R. Moyard, "Introduction to the ZX-calculus."
  [Pennylane](https://pennylane.ai/qml/demos/tutorial_zx_calculus/)
* H. Bombin *et al.*, "Unifying flavors of fault tolerance with the ZX
  calculus." [arXiv:2303.08829](https://arxiv.org/abs/2303.08829)

## See also
* [PyZX](https://github.com/Quantomatic/pyzx): a Python implementation of the
  ZX-calculus and its rewrite rules
* [QuiZX](https://github.com/Quantomatic/quizx/tree/master): a Rust
  implementation of the above.

