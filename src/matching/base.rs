#[allow(unreachable_code)]
use crate::matching::Graph;
use std::mem::size_of;


/// This file contains functions that
/// are used to calculate the matching polynomials for graphs.
/// This method will use static memory for the polynomial.
/// Since the computation graph for the polynomial requires
/// no reads until returning to the caller, and the 
/// write operations are done in any order (since we 
/// calculate the signless matching polynomial), we can
/// avoid large numbers of allocations by using static memory.
///

const POLY_SIZE: usize = size_of::<usize>()*8;
//static mut POLY: [u64; size_of::<usize>()*8] = [0; size_of::<usize>()*8]; 
pub fn calculate_matching_polynomial_static(graph: Graph) -> [u64; size_of::<usize>()*8] {
   let poly: &mut [u64; POLY_SIZE] = &mut [0; POLY_SIZE];
   unsafe {
       // clear the memory
       for i in 0..size_of::<usize>()*8 {
           poly[i] = 0;
       }
       _calculate_matching_polynomial_static(graph, poly);
       return *poly
   }
}

/*
 * The caching algorithm will work as follows:
 *  - after the return from each branch, we will store the current POLY
 *  - from these, calculate the change to the polynomial that results
 *    from that branch.
 *  - Since each update to the polynomial is just a positive sum,
 *    we can cache these diffs to get the addition corresponding to that
 *    graph.
 *  - store this diff alongside a hashed version of the graph
 *  Ideally, we could calculate a canonical graph form that allow us to 
 *  compare graphs without having to do a full graph isomorphism check.
 *   So that we can always 
 *  get two graphs to be comparable
 *  */

unsafe fn _calculate_matching_polynomial_static(graph: Graph, poly: &mut [u64; POLY_SIZE]) {
    if graph.edgeless() {
        let node_count = graph.edgeless_node_count();
        poly[node_count] += 1;
    } else {
        let (graph_prime, graph_prime_prime) = graph.get_graph_primes();
        _calculate_matching_polynomial_static(graph_prime_prime, poly);
        _calculate_matching_polynomial_static(graph_prime, poly); // put this at the end as I think
    }
}

//fn _get_graph_primes(graph: Graph) -> (Graph, Graph) {
    
//}
