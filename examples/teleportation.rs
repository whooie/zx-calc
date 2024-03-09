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

// demonstrate that the standard state teleportation circuit works
//
//         ∣ψ⟩ ------*--H-- measurement A
//                   |
//              / ---X----- measurement B
// ∣00⟩ + ∣11⟩ | 
//              \ ---------(X if B)---(Z if A)--- ∣ψ⟩
//
fn main() -> anyhow::Result<()> {
    // results from measurements on Alice's wires
    const A: bool = false;
    const B: bool = true;

    // create the diagram
    // the dropped value is a hashmap of all node IDs
    let (mut diagram, _) = diagram!(
        nodes: {
            i = input [],
            z = z [], // CNOT
            x = x [], // |
            h = h [],
            eff0 = X [f64::from(A) * PI],
            eff1 = X [f64::from(B) * PI],
            corrx = X [f64::from(B) * PI],
            corrz = Z [f64::from(A) * PI],
            o = output [],
        },
        wires: {
            i -- z,
            z -- h,
            h -- eff0,
            z -- x,
            x -- eff1,
            x -- corrx,
            corrx -- corrz,
            corrz -- o,
        },
    )?;

    // convert the default graph representation to a ket-bra representation
    print!("convert to ketbra ... ");
    let (ketbras, t) = timeit(|| diagram.as_ketbra());
    println!("{:.3e} secs", t);
    println!("ketbras =\n{:.3}", ketbras);

    // both representations can be rendered as DOT language
    print!("save init graph to graphviz ... ");
    let (res, t) = timeit(|| diagram.save_graphviz("graph1", "graph1.gv"));
    println!("{:.3e} secs", t);
    res?;

    print!("save ketbra to graphviz ... ");
    let (res, t)
        = timeit(|| ketbras.save_graphviz("ketbra_conv1", "ketbra_conv1.gv"));
    println!("{:.3e} secs", t);
    res?;

    // convert the ket-bra representation back to a graph
    print!("convert ketbra to graph ... ");
    let (res, t) = timeit(|| ketbras.as_graph());
    println!("{:.3e} secs", t);

    let graph = res?;
    print!("save ketbra graph to graphviz ... ");
    let (res, t)
        = timeit(|| graph.save_graphviz("ketbra_conv2", "ketbra_conv2.gv"));
    println!("{:.3e} secs", t);
    res?;

    // brute-force computation of the diagram with the ket-bras
    print!("contract ketbra ... ");
    let (contr, t) = timeit(|| ketbras.contract().unwrap());
    println!("{:.3e} secs", t);
    println!("output =\n{}", format!("{:.3}", contr).replace(" + ", "\n+ "));

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

    // // alternatively, choose which rules to apply and when; returned bools
    // // indicate whether the rule was successfully applied
    // let (simps, t)
    //     = timeit(|| {
    //         vec![
    //             diagram.simplify_bit_bialgebra(),
    //             diagram.simplify_fuse(),
    //             diagram.simplify_fuse(),
    //             diagram.simplify_hopf(),
    //             diagram.simplify_identity(),
    //             diagram.simplify_identity(),
    //         ]
    //     });
    // println!("{:.3e}", t);
    // println!("simplifications = {:?}", simps);

    print!("save simplified graph to graphviz ... ");
    let (res, t) = timeit(|| diagram.save_graphviz("graph2", "graph2.gv"));
    println!("{:.3e} secs", t);
    res?;

    Ok(())
}

