# pdf-set

Provide PModel that is a set of PDF (probability density function) for range coder

## installation

`Cargo.toml`

```toml
[dependencies]
# this pmodel
pdf_set = {git="https://github.com/diegodox/pdf-set.git"}
# also we need rangecoder
range_coder = {git="https://github.com/diegodox/range_coder_rust.git", branch="carry_less_without_freq_table"}
```
## example

```rust
struct GaussianDist {
    h: f64,
    w: f64,
    m: u8,
}
impl PDF for GaussianDist {
    fn freq(&self, v: usize) -> f64 {
        self.h
            * self.w
            * (-1.0
                * self.w
                * self.w
                * (std::cmp::max(v, self.m as usize) - std::cmp::min(v, self.m as usize)) as f64
                * (std::cmp::max(v, self.m as usize) - std::cmp::min(v, self.m as usize)) as f64)
                .exp()
    }
}
fn test() {
    let ansewr = vec![34, 45, 128, 255, 0];
    let g1 = GaussianDist {
        h: std::f64::MAX,
        w: std::f64::MIN_POSITIVE,
        m: 128,
    };
    let g2 = GaussianDist {
        h: 10.0,
        w: 2.0,
        m: 30,
    };
    let g3 = GaussianDist {
        h: 2.0,
        w: 5.0,
        m: 70,
    };
    let set = PDFSet::new(vec![g1, g2, g3]);
    let pm = set.finalize();
    let mut encoder = Encoder::new();
    for i in &ansewr {
        encoder.encode(&pm, *i);
    }
    encoder.finish();
    let data = encoder.data().clone();
    let mut decoder = Decoder::new();
    decoder.set_data(data);
    decoder.decode_start();
    let mut decoded = Vec::new();
    for _ in 0..ansewr.len() {
        decoded.push(decoder.decode_one_alphabet(&pm));
    }
    assert_eq!(ansewr, decoded);
}
```
