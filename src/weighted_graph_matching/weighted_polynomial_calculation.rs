use std::mem::size_of;

const MAX_NODES: usize = size_of::<usize>()* 8;
const NUM_SIZE: usize = size_of::<usize>() * 8;
/// Given an address, and a set of weights, 
/// calculates the relevant weight for the path that is implied by the address.
///
/// Remember, we calculate a Left as a 1
pub fn weight_from_address(address: usize, weights: &[f32; MAX_NODES * MAX_NODES], graph_size: usize) -> f32 {
    // start with a 1-weight
    let mut weight = 1.0;
    // the length of the address
    let address_len = NUM_SIZE - address.leading_zeros() as usize;
    // subract 1 from the address length to avoid the "cap" bit
    let mut mask = 1 << (address_len-2); // this is the mask
    println!("address: {:b}", address);
    println!("mask: {:b}", mask);

    // an all-node mask. As we process the address, we will zero out removed
    // nodes as bits in this mask
    let mut valid_nodes = (1 << (graph_size + 1)) - 1;
    let mut weights_to_skip = 0;
    println!("initial valid_nodes: {:b}", valid_nodes);
    while mask > 0 {
        let bit = address & mask;
        println!("bit: {:b}", bit);
        if bit != mask { // i.e. it's a right-leg in the binary tree
            println!("Now changing the weight...");
            weight *= get_next_weight( address, weights,  address_len, mask, weights_to_skip, graph_size, &mut valid_nodes);
            println!("valid_nodes: {:b}", valid_nodes);
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

/// Given the current address, and the mask, can we get the "next weight"?
pub fn get_next_weight(address: usize,
                       weights: &[f32; MAX_NODES * MAX_NODES],
                       address_len: usize,
                       mask: usize,
                       weights_to_skip: usize,
                       graph_size: usize, 
                       valid_nodes: &mut usize) -> f32 {

    // capture the fact that the weights are to be skipped
    let mut skipped_weights = 0;
    let mut skipped_zeros = 0;

    let mut current_node =  valid_nodes.trailing_zeros() as usize;
    println!("current_node: {}", current_node);

    // build the range
    //println!("starting_index: {}", starting_index);
    //println!("range: {:?}", range);

    let mut index: usize =  0; // put index into the outer scope
    println!("\n");    
    println!("weights_to_skip: {}", weights_to_skip);

    // build out the index
    'nodes: loop {
        println!("current_node: {}", current_node);
        let range = ((current_node * graph_size)..((current_node + 1) * graph_size)).rev();
        println!("range: {:?}", range);
        println!("starting the range loop again");
        for i in range {
            println!("i: {}", i);
            println!("weights[i]: {}", weights[i]);
            if weights[i] != 0.0 {
                println!("skipped_weights: {}", skipped_weights);
                if skipped_weights == weights_to_skip {
                    index = i;
                    println!("index: {}", index);
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
    println!("index at end of process: {}", index);

    // if the weight is found to connect between current_node and say node '3'  
    // then we need to convert: 11111111 -> 11110110 
    // This is because the 4rd bit from the right is the 4th node in the graph
    // current node is then always the number of trailing zeroes
    println!("index {} graph_size {}", index, graph_size);
    println!("index % graph_size {}", index % graph_size);
    println!("current_node {}", current_node);
    let valid_nodes_mask = !((1 << (index % graph_size))) & !(1 << current_node);
    println!("valid_nodes_mask: {:b}", valid_nodes_mask);
    *valid_nodes &= valid_nodes_mask;
    println!("valid_nodes having just changed them: {:b}", valid_nodes);
    // now we have the index of the weight we want
    weights[index]
}
