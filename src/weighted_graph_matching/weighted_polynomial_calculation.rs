use std::mem::size_of;
use itertools::Itertools;
use std::collections::HashSet;
use std::iter;
use crate::traits::Graph;



const MAX_NODES: usize = size_of::<usize>()* 8;
const NUM_SIZE: usize = size_of::<usize>() * 8;
const POLY_SIZE: usize = size_of::<usize>() * 8;




/// THis calculates the weighted matching polynomial from a matrix of weights,
/// via the method of calculating the sum of products of the weights.
/// This is the function endpoint for the permutation-based method described below. Too slow!
pub fn weighted_polynomial_calculation(weights: &[f32; MAX_NODES * MAX_NODES], graph_size: usize) -> [f32; POLY_SIZE] {
    let mut poly = [0.0; POLY_SIZE];
    let range = 0..=graph_size;
    for (i, coeff_index) in range.enumerate() {
        // coefficient
        let coefficient = weighted_coefficient_calculation(weights, graph_size, coeff_index);

        // add to the polynomial
        poly[i] = coefficient;
    }
    poly 
}

/// This uses a permutation-based mechanism which is too slow!
pub fn weighted_coefficient_calculation(weights: &[f32; MAX_NODES * MAX_NODES], graph_size: usize, coeffic: usize) -> f32 {
    // giben the polynomial poly, and the weights, calculate the weighted 
    // coefficient as marked by coeffic
    if coeffic == 0 {
        return 1.0
    }
    // the following iterator is now an iterator over the indices of non-zero weights.
    let non_zero_weight_indices = weights
        .iter()
        .enumerate()
        .filter(|(_, x)| **x != 0.0) // filter to positive
        .filter(|(i, _)| {  // filter to upper diagonal
            let res = (i / graph_size, i % graph_size);
            res.0 <= res.1
         })
        .map(|(i, _)| i)
        .collect::<Vec<usize>>();

    // We want to produce an iterator whose element is a tuple/iterator of 
    // indices. The product of each of the indices is one of the weights;
    // the sum of these products is the coefficient
    let my_iterators = iter::repeat(non_zero_weight_indices.clone())
        .take(coeffic/2)
        .multi_cartesian_product()
        .map(|v| {
            let mut v = v;
            v.sort();
            v
        })
        .unique();


    // for each of the k index sets, keep the pairs such that the column and row of each of the elements is not shared by
    // any of the other indices.
    let full_hashsets = my_iterators
        .clone()
        .map(|v| // looks like 1, 3, 11
                {
                 let hashsets = v
                  .iter() // turn into an iterator
                  //.map(|x| {
                      //println!("x: {}", x);
                      //x
                  //})
                  .map(|i| HashSet::from([i / graph_size, i % graph_size])) // convert to hashsets of the column and row
                  .fold(HashSet::new(), |acc, hs| {
                      acc
                          .union(&hs)
                          .cloned()
                          .collect()
                  }); // fold into a single hashset
                     hashsets
                });
    let iterator = my_iterators.zip(full_hashsets)
                         .filter(|(_, hs)| hs.len() == coeffic) // filter out the hashsets that have a length of 0;
                         .map(|(set, _)| set); 

    // now get the coefficient using this big iterator
    let coeffic = iterator
        .map(|v| 
             {
                v.iter()
                .map(|i| weights[*i])
                .product::<f32>()
             })
        .sum::<f32>();
    coeffic
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_weighted_coefficient_calculation() {
        let true_weights: [f32; 16] = [0.0, 2.0, 0.0, 2.0, 
                                       2.0, 0.0, 2.0, 0.0,
                                       0.0, 2.0, 0.0, 2.0, 
                                       2.0, 0.0, 2.0, 0.0];
        let true_poly = [3, 0, 6, 0, 1];

        //let poly = Polynomial::new(vec![8.0, 0.0, 8.0, 0.0, 1.0]);
        let mut poly = [0 as u64; size_of::<u64>() * 8];
        let mut weights: [f32; 4096] = [0.0; 4096];
        weights[..16].copy_from_slice(&true_weights);
        poly[..5].copy_from_slice(&true_poly);
        let weighted_coefficient = weighted_coefficient_calculation(&weights, 4, 2);
        assert_eq!(weighted_coefficient, 8.0);
    }

}
