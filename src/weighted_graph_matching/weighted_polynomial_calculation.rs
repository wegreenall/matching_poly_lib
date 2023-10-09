use std::mem::size_of;
use itertools::{iproduct, izip};
use itertools::Itertools;
use permutator::{Permutation, KPermutationIterator, k_permutation};
use std::collections::HashSet;
use std::iter;
use crate::{BinaryGraph, calculate_matching_polynomial_pointer_addresses};
use crate::traits::Graph;



const MAX_NODES: usize = size_of::<usize>()* 8;
const NUM_SIZE: usize = size_of::<usize>() * 8;
const POLY_SIZE: usize = size_of::<usize>() * 8;

/// The following two functions together calculate the weighted matching polynomial by
/// building the address and then in situ calculating the weight.
pub fn weighted_matching_polynomial_addresses(graph: BinaryGraph, weights: &[f32; MAX_NODES * MAX_NODES]) -> [f32; size_of::<usize>()*8] {
        let poly: &mut [f32; POLY_SIZE] = &mut [0.0; POLY_SIZE];
        let address = 1;

        for i in 0..POLY_SIZE {
            poly[i] = 0.0;
        }
        _calculate_weighted_matching_polynomial_static_addresses(graph, weights, poly, address);
        *poly
}

fn _calculate_weighted_matching_polynomial_static_addresses<T: Graph>(graph: T, weights: &[f32; MAX_NODES * MAX_NODES], poly: &mut [f32; POLY_SIZE], current_address: usize) {
    if graph.edgeless() {
        let node_count = graph.edgeless_node_count();
        poly[node_count] += weight_from_address(current_address, weights, graph.initial_graph_size());
    } else {
        let (graph_prime, graph_prime_prime) = graph.get_graph_primes();
        _calculate_weighted_matching_polynomial_static_addresses(graph_prime_prime, weights, poly, current_address << 1);
        _calculate_weighted_matching_polynomial_static_addresses(graph_prime, weights, poly, (current_address << 1) + 1); // put this at the end as I think
    }
}


/// The following function will accept a sequence of addresses, and weights,
/// and then build the weighted matching polynomial from them.
pub fn weighted_matching_polynomial_from_addresses(addresses: Vec<usize>, weights: &[f32; MAX_NODES * MAX_NODES], graph_size: usize) -> [f32; size_of::<usize>()*8] {
    let poly: &mut [f32; POLY_SIZE] = &mut [0.0; POLY_SIZE];
    for i in 0..POLY_SIZE {
        poly[i] = 0.0;
    }
    addresses.iter().for_each(|address| {
        // get the relevant coefficient
        let coefficient = graph_size - 2 * (address.count_zeros() - address.leading_zeros()) as usize;
        poly[coefficient] += weight_from_address(*address, weights, graph_size);
    });
    *poly
}

/// Given an address, and a set of weights, 
/// calculates the relevant weight for the path that is implied by the address.
/// graph_size should be the current graph_size
pub fn weight_from_address(address: usize, weights: &[f32; MAX_NODES * MAX_NODES], graph_size: usize) -> f32 {
    // start with a 1-weight
    let mut weight = 1.0;
    // the length of the address
    let address_len = NUM_SIZE - address.leading_zeros() as usize;
    // subract 1 from the address length to avoid the "cap" bit
    let mut mask = 1 << (address_len-2); // this is the mask

    // an all-node mask. As we process the address, we will zero out removed
    // nodes as bits in this mask
    let mut valid_nodes = (1 << (graph_size + 1)) - 1;
    let mut weights_to_skip = 0;

    while mask > 0 {
        let bit = address & mask;
        if bit != mask { // i.e. it's a right-leg in the binary tree
            weight *= get_next_weight(weights,
                                      weights_to_skip,
                                      graph_size,
                                      &mut valid_nodes);
            weights_to_skip = 0;
        } else {
            weights_to_skip += 1;
        }
        mask >>= 1;
    }
    // when we have exhausted the mask, _weight_ is the product of the weights 
    // down the tree
    weight
}

/// Given the current address, graph size, and a marker containing the valid 
/// nodes. 
/// The algorithm works as follows. A given address signifies a path through 
/// the graph. That path has a weight attached, and we want to alter the 
/// appropriate coefficient in the polynomial. To do this we need to calculate
/// the weight of the path.
/// 
/// First, we start at node 0.
/// An address is a binary number like 
///            11100100
/// The first bit is the cap bit; thereon a 0 refers to a right-leg, and a 1
/// refers to a left-leg. If it's a right leg, it marks the deletion  of a pair
/// of nodes and postpending of a graph weight to the path weight. The aim of
/// the function below is to calculate the postpending graph weight.
///
/// To do this, we travel through the graph adjacency skipping over edges as
/// recorded as left-legs. When we get to the corresponding right-leg, we alter
/// the set of relevant nodes that we are aiming to work with and then we
/// return the calculated. 
///
/// 
fn get_next_weight(//address: usize,
                       weights: &[f32; MAX_NODES * MAX_NODES],
                       //address_len: usize,
                       //mask: usize,
                       weights_to_skip: usize,
                       graph_size: usize, 
                       valid_nodes: &mut usize) -> f32 {

    // capture the fact that the weights are to be skipped
    let mut skipped_weights = 0;

    let mut current_node =  valid_nodes.trailing_zeros() as usize;

    let mut index: usize =  0; // put index into the outer scope

    // build out the index
    'nodes: loop {
        let range = ((current_node * graph_size)..((current_node + 1) * graph_size)).rev();
        for i in range {
            if weights[i] != 0.0 {
                if skipped_weights == weights_to_skip {
                    index = i;
                    break 'nodes;
                }
                skipped_weights += 1;
            }
        }
        current_node += (*valid_nodes >> (current_node+1)).trailing_zeros() as usize + 1;
        if current_node > graph_size{
            break 'nodes
        };
    }

    // if the weight is found to connect between current_node and say node '3'  
    // then we need to convert: 11111111 -> 11110110 
    // This is because the 4rd bit from the right is the 4th node in the graph
    // current node is then always the number of trailing zeroes
    let valid_nodes_mask = !((1 << (index % graph_size))) & !(1 << current_node);
    *valid_nodes &= valid_nodes_mask;

    // now we have the index of the weight we want
    weights[index]
}

/// THis calculates the weighted matching polynomial from a matrix of weights,
/// via the method of calculating the sum of products of the weights.
/// This is the function endpoint for the permutation-based method described below. Too slow!
pub fn weighted_polynomial_calculation(weights: &[f32; MAX_NODES * MAX_NODES], graph_size: usize) -> [f32; POLY_SIZE] {
    let mut poly = [0.0; POLY_SIZE];
    let range = (0..=graph_size).rev();
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

    #[test]
    fn test_weighted_polynomial_static_addresses() {
        let data = [
            0b1101, 0b110, 0b11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let true_weights: [f32; 16] = [0.0, 2.0, 0.0, 2.0, 
                                       2.0, 0.0, 2.0, 0.0,
                                       0.0, 2.0, 0.0, 2.0, 
                                       2.0, 0.0, 2.0, 0.0];
        let graph = BinaryGraph::from(data); // the fully connected graph
 
        //let poly = Polynomial::new(vec![8.0, 0.0, 8.0, 0.0, 1.0]);
        //let mut poly = [0.0 as f32; size_of::<f32>() * 8];
        let mut weights: [f32; 4096] = [0.0; 4096];
        weights[..16].copy_from_slice(&true_weights);
        //poly[..5].copy_from_slice(&true_poly);
        let weighted_polynomial = weighted_matching_polynomial_addresses(graph, &weights);
        assert_eq!(weighted_polynomial[..5], [8.0, 0.0, 8.0, 0.0, 1.0]);
    }

    #[test]
    fn test_weighted_matching_polynomial_from_addresses() {
        let data = [
            0b1101, 0b110, 0b11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let true_weights: [f32; 16] = [0.0, 2.0, 0.0, 2.0, 
                                       2.0, 0.0, 2.0, 0.0,
                                       0.0, 2.0, 0.0, 2.0, 
                                       2.0, 0.0, 2.0, 0.0];
        let graph = BinaryGraph::from(data); // the fully connected graph
 
        // get the addresses from the graph:
        let (_, addresses) = calculate_matching_polynomial_pointer_addresses(graph);
        
        let mut weights: [f32; 4096] = [0.0; 4096];
        weights[..16].copy_from_slice(&true_weights);
        let weighted_polynomial = weighted_matching_polynomial_from_addresses(addresses, &weights, graph.initial_graph_size());
        assert_eq!(weighted_polynomial[..5], [8.0, 0.0, 8.0, 0.0, 1.0]);
    }

}
