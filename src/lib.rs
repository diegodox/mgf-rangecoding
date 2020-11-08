//! 任意の確率分布の和から，レンジコーダの確率モデルを生成する  
//! 確率密度関数を表すトレイト: PDF  
//! トレイトPDFの集合: PDFSet  
//! PDFSetを量子化した確率密度関数: QuantizedPDFSet  
//! QuantizedPDFSetはRangeCoderのPModelを実装  

pub use range_coder;
use range_coder::decoder::Decoder;
use range_coder::pmodel::PModel;
/// a set of probability density functions.
pub struct PDFSet<T: PDF> {
    pdf_list: Vec<T>,
}
impl<T: PDF> PDFSet<T> {
    pub fn new(vec: Vec<T>) -> Self {
        Self { pdf_list: vec }
    }
    pub fn add_pdf(&mut self, pdf: T) {
        self.pdf_list.push(pdf);
    }
    pub fn finalize(self) -> QuantizedPDFSet {
        const RANGE_MAX: usize = std::u8::MAX as usize;
        const RANGE_SIZE: usize = RANGE_MAX + 1;
        const RANGE: std::ops::RangeInclusive<usize> = 0..=RANGE_MAX;
        let (freq_src, tot_freq_src) = {
            let mut freq_src = Vec::with_capacity(RANGE_SIZE);
            let tot_freq = RANGE
                .into_iter()
                // 確率質量関数の確率の合計を計算する
                .map(|x| {
                    self.pdf_list
                        .iter()
                        .map(|p| p.freq(x as usize))
                        .sum::<f64>()
                })
                // 累積確率を計算する
                .fold(0f64, |cum, freq| {
                    // 頻度表に登録する
                    freq_src.push(freq);
                    cum + freq
                });
            (freq_src, tot_freq)
        };
        // 量子化
        let (freq, cum_freq) = {
            /// 各値に底上げとして1ずつ割り振るので，maxから引いておく
            const MAX_TOT_FREQ: u32 = std::u32::MAX - (std::u8::MAX as u32 + 1);
            let mut freq = Vec::with_capacity(RANGE_SIZE);
            let mut cum_freq = Vec::with_capacity(RANGE_SIZE);
            RANGE
                .into_iter()
                // 整数へ丸めた頻度を計算（1の底上げもする）
                .map(|x| (MAX_TOT_FREQ as f64 * (freq_src[x as usize] / tot_freq_src)) as u32 + 1)
                // 累積頻度の計算
                .scan(0, |cum, freq| {
                    let cum_clone = cum.clone();
                    *cum += freq;
                    Some((freq, cum_clone))
                })
                // 頻度表に登録
                .for_each(|(f, cum)| {
                    freq.push(f);
                    cum_freq.push(cum);
                });
            (freq, cum_freq)
        };
        QuantizedPDFSet { freq, cum_freq }
    }
}
/// probability density function
pub trait PDF {
    fn freq(&self, v: usize) -> f64;
}
pub struct QuantizedPDFSet {
    freq: Vec<u32>,
    cum_freq: Vec<u32>,
}
impl PModel for QuantizedPDFSet {
    fn c_freq(&self, index: usize) -> u32 {
        self.freq[index]
    }
    fn cum_freq(&self, index: usize) -> u32 {
        self.cum_freq[index]
    }
    fn total_freq(&self) -> u32 {
        *self.cum_freq.last().unwrap() + *self.freq.last().unwrap()
    }
    fn find_index(&self, decoder: &Decoder) -> usize {
        let mut left = 0;
        let mut right = std::u8::MAX as usize;
        let rfreq = (decoder.data() - decoder.range_coder().lower_bound())
            / decoder.range_coder().range_par_total(self.total_freq());
        while left < right {
            let mid = (left + right) / 2;
            let mid_cum = self.cum_freq(mid + 1);
            if mid_cum as u64 <= rfreq {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        left
    }
}
impl std::fmt::Debug for QuantizedPDFSet {
    fn fmt(&self, _f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..=std::u8::MAX {
            println!("{:03}: {}", i, self.c_freq(i as usize));
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use crate::PDFSet;
    use crate::QuantizedPDFSet;
    use crate::PDF;
    use range_coder::{decoder::Decoder, encoder::Encoder, pmodel::PModel};
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
                    * (std::cmp::max(v, self.m as usize) - std::cmp::min(v, self.m as usize))
                        as f64
                    * (std::cmp::max(v, self.m as usize) - std::cmp::min(v, self.m as usize))
                        as f64)
                    .exp()
        }
    }
    fn simple_pmodel() -> QuantizedPDFSet {
        let g1 = GaussianDist {
            h: 10.0,
            w: 5.0,
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
        let set = PDFSet {
            pdf_list: vec![g1, g2, g3],
        };
        set.finalize()
    }
    fn large_pmodel() -> QuantizedPDFSet {
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
        let set = PDFSet {
            pdf_list: vec![g1, g2, g3],
        };
        set.finalize()
    }
    #[test]
    fn it_works() {
        let ansewr = vec![34, 45, 128, 255, 0];
        let pm = simple_pmodel();
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
        println!("{:?}", pm);
    }
    #[test]
    fn large_pdf() {
        let ansewr = vec![34, 45, 128, 255, 0];
        let pm = large_pmodel();
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
    #[test]
    fn use_full_range() {
        let pm = large_pmodel();
        assert!((std::u32::MAX as f64)
            .mul_add(-0.9999999, pm.total_freq() as f64)
            .is_sign_positive());
    }
    #[test]
    fn use_full_range2() {
        let pm = simple_pmodel();
        assert!((std::u32::MAX as f64)
            .mul_add(-0.9999999, pm.total_freq() as f64)
            .is_sign_positive());
    }
}
