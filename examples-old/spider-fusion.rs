use std::{
    f64::consts::TAU,
    time::Instant,
};
use rand::{ thread_rng, Rng };
use zx_calc::graph::*;

const N: usize = 10_000;

fn timeit<F, T>(mut f: F) -> (T, f64)
where F: FnMut() -> T
{
    let t0 = Instant::now();
    let out: T = f();
    (out, (Instant::now() - t0).as_secs_f64())
}

// Generate ten thousand spiders with randomly chosen phases and ten thousand
// randomly chosen wires, and time how long it takes to fuse them all together.
fn main() -> anyhow::Result<()> {
    let mut rng = thread_rng();
    let mut diagram
        = Diagram::from_nodes(
            (0..N).map(|_| NodeDef::Z(TAU * rng.gen::<f64>())));
    for _ in 0..N {
        diagram.add_wire(
            rng.gen_range(0..N), rng.gen_range(0..N))?;
    }
    diagram.add_input_wire(rng.gen_range(0..N))?;
    diagram.add_output_wire(rng.gen_range(0..N))?;

    println!("# spiders = # wires = {}", N);
    print!("fuse all spiders ... ");
    let (_, t) = timeit(|| { while diagram.simplify_fuse() { } });
    println!("{:.3} secs", t);

    Ok(())
}

