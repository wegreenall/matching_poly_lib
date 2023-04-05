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

static mut POLY: [u64; size_of::<usize>()*8] = [0; size_of::<usize>()*8]; 
pub fn calculate_matching_polynomial_static(graph: Graph) -> [u64; size_of::<usize>()*8] {
   unsafe {
       // clear the memory
       POLY = [0; size_of::<usize>()*8];
       _calculate_matching_polynomial_static(graph);
       return POLY;
   }
}

unsafe fn _calculate_matching_polynomial_static(graph: Graph) {
    if graph.edgeless() {
        let node_count = graph.edgeless_node_count();
        POLY[node_count] += 1;
    } else {
        let (graph_prime, graph_prime_prime) = graph.get_graph_primes();
        _calculate_matching_polynomial_static(graph_prime);
        _calculate_matching_polynomial_static(graph_prime_prime);
    }
}
