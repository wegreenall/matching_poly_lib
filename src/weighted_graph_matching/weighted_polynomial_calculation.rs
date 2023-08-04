use std::mem::size_of;
use itertools::iproduct;
use itertools::Itertools;
use permutator::{Permutation, KPermutationIterator, k_permutation};


const MAX_NODES: usize = size_of::<usize>()* 8;
const NUM_SIZE: usize = size_of::<usize>() * 8;
/// Given an address, and a set of weights, 
/// calculates the relevant weight for the path that is implied by the address.
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
pub fn get_next_weight(//address: usize,
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

pub fn weighted_coefficient_calculation(poly: [u64; size_of::<usize>()*8], weights: &[f32; MAX_NODES * MAX_NODES], graph_size: usize, coeffic: usize) -> f64 {
    // giben the polynomial poly, and the weights, calculate the weighted 
    // coefficient as marked by coeffic
    let mut weighted_coefficient = 1.0;    

    // can probably just pass the value
    let missing_nodes = graph_size - 2 * coeffic;

    // get all the combinations of these missing node numbers
    let mut non_zero_weights = weights
        .iter()
        .enumerate()
        .filter(|(_, x)| **x != 0.0) // filter to positive
        .filter(|(i, _)| {  // filter to upper diagonal
            let res = (i / graph_size, i % graph_size);
            res.0 <= res.1
        })
        .map(|(i, x)| (i, *x))
        .collect::<Vec<(usize, f32)>>();

    for stuff in non_zero_weights.clone() {
        println!("{:?}", stuff);
    }
    //let permutator = HeapPermutationIterator::new(&mut non_zero_weights);
    let mut counter = 0;
    let mut perms = Vec::<Vec<&(usize, f32)>>::new();

    let k_permutation_iterator = KPermutationIterator::new(&mut non_zero_weights, coeffic);
    for perm in k_permutation_iterator {
        println!("{:?}", perm);
    }

    // now take the product?
    //let mut multi_prod = missing_nodes;
    //let mut prod = iproduct!(non_zero_weights.clone(), non_zero_weights.clone())
        //.filter(|((i, _),(j,_))| {
            //let res = (i / graph_size, i % graph_size);
            //let res2 = (j / graph_size, j % graph_size);
            //// not same row
            //// not same column
            //// not column and row
            //// not row and column
            //res.0 != res2.0 && res.1 != res2.1 && res.1 != res2.0 && res.0 != res2.1
        //});
    //for weight in prod {
        //println!("{:?}", weight);
    //}        
    

    // we now have the indices of the non-zero weights
    weighted_coefficient
    
}
