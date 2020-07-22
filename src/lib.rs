//! 任意の確率分布の和から，レンジコーダの確率モデルを生成する  
//! 確率密度関数を表すトレイト: PDF  
//! トレイトPDFの集合: PDFSet  
//! PDFSetを量子化した確率密度関数: QuantizedPDFSet  
//! QuantizedPDFSetはRangeCoderのPModelを実装  
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
        let mut freq_src = Vec::new();
        let mut cum_freq_src = Vec::new();
        {
            // freq,cum_freqをf64で計算
            let mut cum_freq_tmp = 0.0;
            for x in 0..=std::u8::MAX {
                let mut freq = 0.0;
                for p in &self.pdf_list {
                    freq += p.freq(x as usize);
                }
                freq_src.push(freq);
                cum_freq_src.push(cum_freq_tmp);
                cum_freq_tmp += freq;
            }
        }
        // 量子化
        let mut freq = Vec::new();
        let mut cum_freq = Vec::new();
        {
            // 各値に底上げとして1ずつ割り振るので，maxから引いておく
            let max = std::u32::MAX - (std::u8::MAX as u32 + 1);
            {
                let tot_freq_src = cum_freq_src.last().unwrap() + freq_src.last().unwrap();
                let mut cum_freq_tmp = 0;
                for x in 0..=std::u8::MAX {
                    let f = max * (freq_src[x as usize] / tot_freq_src) as u32 + 1;
                    freq.push(f);
                    cum_freq.push(cum_freq_tmp);
                    cum_freq_tmp += f;
                }
            }
        }
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

#[cfg(test)]
mod tests {
    use crate::PDFSet;
    use crate::PDF;
    use range_coder::decoder::Decoder;
    use range_coder::encoder::Encoder;
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
                    * (std::cmp::max(v, self.m as usize) - std::cmp::min(v, self.m as usize)) as f64
                   ).exp()
        }
    }
    #[test]
    fn it_works() {
        let ansewr = vec![34, 45, 128, 255, 0];
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
    #[test]
    fn large_pdf() {
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
        let set = PDFSet {
            pdf_list: vec![g1, g2, g3],
        };
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
}
