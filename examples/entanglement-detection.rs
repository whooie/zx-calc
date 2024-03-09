use std::{
    f64::consts::PI,
    time::Instant,
};
use zx_calc::{ diagram, graph::* };

fn timeit<F, T>(mut f: F) -> (T, f64)
where F: FnMut() -> T
{
    let t0 = Instant::now();
    let out: T = f();
    (out, (Instant::now() - t0).as_secs_f64())
}

// entanglement detection circuit
//
// ∣0⟩ -----------*--------X--------------*-----*--
//                |        |              |     |
// ∣0⟩ -----------|--*--X--*--H--X-----X--|--S--X--
//                |  |  |        |     |  |  
// ∣0⟩ --H--S--*--X--|--|--------*--X--|--X--------
//             |     |  |           |  |   
// ∣0⟩ --------X-----X--*-----------*--*-----------
//
fn main() -> anyhow::Result<()> {
    // create the diagram
    // the dropped value is a hashmap of all node IDs
    let (mut diagram, nodes) = diagram!(
        nodes: {
            i0 = input [],
            i1 = input [],
            i2 = input [],
            i3 = input [],
            h0 = h [],
            h1 = h [],
            s0 = Z [PI / 2.0],
            s1 = Z [PI / 2.0],
            cnot0z = z [],
            cnot0x = x [],
            cnot1z = z [],
            cnot1x = x [],
            cnot2z = z [],
            cnot2x = x [],
            cnot3z = z [],
            cnot3x = x [],
            cnot4z = z [],
            cnot4x = x [],
            cnot5z = z [],
            cnot5x = x [],
            cnot6z = z [],
            cnot6x = x [],
            cnot7z = z [],
            cnot7x = x [],
            cnot8z = z [],
            cnot8x = x [],
            cnot9z = z [],
            cnot9x = x [],
            o0 = output [],
            o1 = output [],
            o2 = output [],
            o3 = output [],
        },
        wires: {
            i0 -- cnot1z,
            cnot1z -- cnot4x,
            cnot4x -- cnot8z,
            cnot8z -- cnot9z,
            cnot9z -- o0,

            i1 -- cnot2z,
            cnot2z -- cnot3x,
            cnot3x -- cnot4z,
            cnot4z -- h1,
            h1 -- cnot5x,
            cnot5x -- cnot7x,
            cnot7x -- s1,
            s1 -- cnot9x,
            cnot9x -- o1,

            i2 -- h0,
            h0 -- s0,
            s0 -- cnot0z,
            cnot0z -- cnot1x,
            cnot1x -- cnot5z,
            cnot5z -- cnot6x,
            cnot6x -- cnot8x,
            cnot8x -- o2,

            i3 -- cnot0x,
            cnot0x -- cnot2x,
            cnot2x -- cnot3z,
            cnot3z -- cnot6z,
            cnot6z -- cnot7z,
            cnot7z -- o3,

            cnot0z -- cnot0x,
            cnot1z -- cnot1x,
            cnot2z -- cnot2x,
            cnot3z -- cnot3x,
            cnot4z -- cnot4x,
            cnot5z -- cnot5x,
            cnot6z -- cnot6x,
            cnot7z -- cnot7x,
            cnot8z -- cnot8x,
            cnot9z -- cnot9x,
        },
    )?;
    diagram.apply_state(nodes["i0"], Spider::x())?;
    diagram.apply_state(nodes["i1"], Spider::x())?;
    diagram.apply_state(nodes["i2"], Spider::x())?;
    diagram.apply_state(nodes["i3"], Spider::x())?;

    // both representations can be rendered as DOT language
    print!("save init graph to graphviz ... ");
    let (res, t) = timeit(|| diagram.save_graphviz("graph1", "graph1.gv"));
    println!("{:.3e} secs", t);
    res?;

    // scalar subgraphs can be removed
    print!("remove graph scalars ... ");
    let (scalar, t) = timeit(|| diagram.remove_scalars());
    println!("{:.3e} secs", t);
    println!("scalar = {:.3e}", scalar);

    // simplify the diagram using a naive, exhaustive strategy to continually
    // apply all available rewrite rules
    print!("graph simplifications ... ");
    let (_, t) = timeit(|| { diagram.simplify(); });
    println!("{:.3e} secs", t);

    print!("save simplified graph to graphviz ... ");
    let (res, t) = timeit(|| diagram.save_graphviz("graph2", "graph2.gv"));
    println!("{:.3e} secs", t);
    res?;

    Ok(())
}

