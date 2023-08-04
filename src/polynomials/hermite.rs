use polynomial::Polynomial;  
//use std::cmp::PartialEq;
//use std::ops::{Mul, Div};
use std::cmp::min;

// Following the numpy way of doing things, we implement in this file some
// conversions between Hermite polynomials and standard basis polynomials

pub fn hermemulx(p: &Polynomial<f32>) -> Polynomial<f32> {
    let data = p.data();
    //let mut new_data = Vec::<T>::new();
    if data.len() == 1 && data[0] == 0.0 {
        return Polynomial::new(vec![0.0]);
    }

    let mut new_data = vec![0.0; data.len()+1];
    new_data[0] = 0.0;
    new_data[1] = data[0];

    for i in 1..data.len() {
        new_data[i + 1] = data[i];
        new_data[i - 1] += data[i] * i as f32;
    }
    Polynomial::new(new_data)
}

pub fn hermadd(p: &Polynomial<f32>, q: &Polynomial<f32>) -> Polynomial<f32>{
    //let mut used_coeffics: Vec<f32>;
    let relevant_length = min(p.data().len(), q.data().len());
    let larger = if p.data().len() > q.data().len() {
        p
    } else {
        q
    };

    let p_coeffics = p
        .data()
        .iter()
        .take(relevant_length);

    let q_coeffics = q
        .data()
        .iter()
        .take(relevant_length);

    let unused_coeffics = larger
        .data()
        .iter()
        .skip(relevant_length)
        .collect::<Vec<&f32>>();

    let mut used_coeffics = p_coeffics
        .zip(q_coeffics)
        .map(|(x, y)| x + y)
        .collect::<Vec<f32>>();

    used_coeffics
        .extend(unused_coeffics);

    Polynomial::new(used_coeffics)
}

pub fn poly2herme(p: &Polynomial<f32>) -> Polynomial<f32> {
    let data = p.data();
    let deg = data.len()-1;

    let mut res = Polynomial::new(vec![data[deg]]);
    let range = (0..deg).rev();
    for index in range {
        res = hermadd(&hermemulx(&res), &Polynomial::new(vec![data[index]]));
    }
    res
}

pub fn herme2poly(p: &Polynomial<f32>) -> Polynomial<f32> {
    let data = p.data();
    let n = data.len();

    // at n = 1 or n = 2, the polynomial is already in standard basis
    if n == 1 || n == 2 {
        return Polynomial::from(p.to_owned());
    } else {
        let c0: f32 = data[n - 2];
        let c1: f32 = data[n - 1];
        let range = (2..n).rev(); // it's inclusive though...
        let mut p0  = Polynomial::new(vec![c0]);
        let mut p1  = Polynomial::new(vec![c1]);
        for i in range {
            //println!("i: {:?} p0 at beginning: {:?}; p1 at beginning: {:?}", i, p0, p1);
            let tmp = p0; // tmp is a number
            p0 = Polynomial::new(vec![data[i-2]]) - &p1 * Polynomial::new(vec![(i-1) as f32]);
            p1 = tmp + polymulx(&p1);
        }
        //println!("p0: {:?}; p1:  {:?}", p0, p1);
        hermadd(&p0, &polymulx(&p1))
    }
}

/// Multiply a polynomial by x
fn polymulx(p: &Polynomial<f32>) -> Polynomial<f32> {
    p * Polynomial::new(vec![0.0, 1.0])
}
