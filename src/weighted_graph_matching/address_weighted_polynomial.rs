use std::mem::size_of;
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
/// and then  uild the weighted matching polynomial from them.
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

#[cfg(test)]
mod tests {
    use super::*;
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
