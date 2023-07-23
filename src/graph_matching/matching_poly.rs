#[allow(unreachable_code)]
use polynomial::Polynomial;
//use crate::graph_matching::BinaryGraph;
use crate::binary_graph_matching::BinaryGraph;
use std::mem::size_of;
use crate::polynomials::{herme2poly, poly2herme};
use crate::traits::Graph;

/// This file contains functions that are used to calculate the matching
/// polynomials for graphs. This method will use static memory for the
/// POLYNOMIAL as opposed to the GRAPH. Since the computation graph for the
/// polynomial requires no reads until returning to the caller, and the write
/// operations are done in any order (since we calculate the signless matching
/// polynomial), we can avoid large numbers of allocations by using static memory.
     
const POLY_SIZE: usize = size_of::<usize>()*8;

pub fn calculate_matching_polynomial_pointer(graph: BinaryGraph) -> [u64; size_of::<usize>()*8] {
   let poly: &mut [u64; POLY_SIZE] = &mut [0; POLY_SIZE];
   //unsafe {
       // clear the memory
       for i in 0..size_of::<usize>()*8 {
           poly[i] = 0;
       }
       _calculate_matching_polynomial_static(Box::new(graph), poly);
       return *poly
   //}
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

/// the following function assumes that the graph is to have its polynomial calculated
/// in the standard basis.
fn _calculate_matching_polynomial_static(graph: Box<BinaryGraph>, poly: &mut [u64; POLY_SIZE]) {
    if graph.edgeless() {
        let node_count = graph.edgeless_node_count();
        poly[node_count] += 1;
    } else {
        let (graph_prime, graph_prime_prime) = graph.get_graph_primes();
        _calculate_matching_polynomial_static(graph_prime_prime, poly);
        _calculate_matching_polynomial_static(graph_prime, poly); // put this at the end as I think
    }
}

/// the following function, aimed at being a drop-in replacement for the one above,
/// assumes that the graph is to have its polynomial calculated adaptively as the 
/// density of the relevant subraph changes over the course of the algorithm.
fn _calculate_matching_polynomial_static_adaptive(graph: Box<BinaryGraph>, poly: &mut [u64; POLY_SIZE], complement: bool) {
    if graph.edgeless() {
        let node_count = graph.edgeless_node_count();
        if complement {
            // run the update as if in the Hermite basis
            poly[node_count] += 1;
        } else {
            // run the update as if in the standard basis
            poly[node_count] += 1;
        }
        poly[node_count] += 1;
    } else {
        let (graph_prime, graph_prime_prime) = graph.get_graph_primes();

        if graph_prime_prime.density() >= 0.5 {
            let complement = true;
            _calculate_matching_polynomial_static_adaptive(graph_prime_prime, poly, complement);
        } else {
            _calculate_matching_polynomial_static_adaptive(graph_prime, poly, complement); // put this at the end as I think
        }
    }
}

//pub fn _calculate_matching_polynomial_binary(graph: BinaryGraph) -> Polynomial<u64> {
    //// the base case for the process is that the graph is edgeless.
    //// This means that, of the remaining nodes, each of their integer
    //// representations is a power of two.
    //if graph.edgeless() { // i.e. we're at the base case.
        //// produce a sequence of coefficients the same length as the number of vertices
        //let mut coeffics = vec![0; graph.edgeless_node_count()];
        //coeffics.push(1);
        //let poly = Polynomial::new(coeffics);
        //return poly
    //} else {
        //// G' = G - {an edge e}
        //// G'' = G - {the nodes connected to the edge e}
        //let (graph_prime, graph_prime_prime) = graph.get_graph_primes();

        //let poly_1 = _calculate_matching_polynomial_binary(graph_prime);
        //let poly_2 = _calculate_matching_polynomial_binary(graph_prime_prime);
        //let poly = poly_1 + poly_2;
        //return poly
    //}
//} 
