use crate::constants;
use ark_bn254::Fr;
use ark_ff::fields::Field;
use ark_ff::Zero;
use std::str::FromStr;

pub trait Sponge {
    fn absorb();
    fn squeeze();
    fn hash();
}

pub struct PoseidonParams {
    t: usize,
    alpha: usize,
    num_f: usize,
    num_p: usize,
    pub mds_matrix: Vec<Vec<Fr>>,
    pub round_constants: Vec<Fr>,
}

pub struct Poseidon {
    state: Vec<Fr>,
    params: PoseidonParams,
}

fn load_constants(t: usize, num_f: usize, num_p: usize) -> (Vec<Vec<Fr>>, Vec<Fr>) {
    let (c_strij, m_strij) = constants::constants();

    let c_str = &c_strij[t];
    let m_str = &m_strij[t];

    let mut round_constants: Vec<Fr> = Vec::with_capacity(t * (num_f + num_p));
    for i in 0..t {
        let a: Fr = Fr::from_str(c_str[i]).expect("a");
        round_constants.push(a);
    }

    let mut mds: Vec<Vec<Fr>> = Vec::new();
    for i in 0..t {
        let mut mds_i: Vec<Fr> = Vec::with_capacity(t);
        for j in 0..t {
            mds_i.push(Fr::from_str(m_str[i][j]).expect("s"));
        }
        mds.push(mds_i);
    }

    (mds, round_constants)
}

impl PoseidonParams {
    fn new(t: usize) -> Self {
        let mut params = PoseidonParams {
            t,
            alpha: 5,
            num_f: 57,
            num_p: 8,
            mds_matrix: Vec::new(),
            round_constants: Vec::new(),
        };
        (params.mds_matrix, params.round_constants) =
            load_constants(params.t, params.num_f, params.num_p);
        params
    }

    fn alpha(&self) -> usize {
        self.alpha
    }

    fn width(&self) -> usize {
        self.t
    }

    fn num_f(&self) -> usize {
        self.num_f
    }

    fn num_p(&self) -> usize {
        self.num_p
    }
}

impl Poseidon {
    fn new(t: usize, state: Vec<Fr>) -> Self {
        Poseidon {
            state,
            params: PoseidonParams::new(t),
        }
    }

    fn sbox_full(&mut self) {
        for i in 0..self.params.width() {
            let temp = self.state[0];
            self.state[i] = self.state[0].square();
            self.state[i] = self.state[0].square();
            self.state[i] *= temp;
        }
    }
    fn sbox_partial(&mut self) {
        let temp = self.state[0];
        self.state[0] = self.state[0].square();
        self.state[0] = self.state[0].square();
        self.state[0] *= temp;
    }

    fn sbox(&mut self, round_i: usize) {
        if round_i < self.params.num_p / 2 || round_i > self.params.num_p / 2 + self.params.num_f {
            self.sbox_full()
        } else {
            self.sbox_partial()
        }
    }

    fn mix(&mut self, round_i: usize) {
        let mut new_state: Vec<Fr> = Vec::with_capacity(self.params.width());
        for (i, new_state_i) in new_state.iter_mut().enumerate().take(self.params.width()) {
            *new_state_i = Fr::zero();
            for j in 0..self.params.width() {
                *new_state_i += self.state[j] * self.params.mds_matrix[i][j];
            }
        }
        self.state = new_state
    }

    fn ark(&mut self, ith: usize) {
        for i in 0..self.params.width() {
            self.state[i] += &self.params.round_constants[ith * self.params.width() + i];
        }
    }

    pub fn hash(&mut self) -> Result<Fr, String> {
        let num_rounds = self.params.num_f + self.params.num_p;
        for i in 0..num_rounds {
            self.ark(i);
            self.sbox(i);
            self.mix(i);
        }

        Ok(self.state[0])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_works() {
        let state: Vec<Fr> = vec![Fr::from(1), Fr::from(2)];
        let mut pos = Poseidon::new(2, state);

        println!("state: {:?}", pos.state);
        let out = pos.hash().unwrap_or(Fr::zero());
        println!("out: {}", out);
    }
}
