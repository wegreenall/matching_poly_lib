use crate::{binary_graph_matching::BinaryGraph, polynomials::poly2herme};
use std::mem::size_of;
use crate::traits::Graph;
use polynomial::Polynomial;

/// This file contains functions that are used to calculate the matching
/// polynomials for graphs. This method will use static memory for the
/// POLYNOMIAL as opposed to the GRAPH. Since the computation graph for the
/// polynomial requires no reads until returning to the caller, and the write
/// operations are done in any order (since we calculate the signless matching
/// polynomial), we can avoid large numbers of allocations by using static
/// memory.
const POLY_SIZE: usize = size_of::<usize>()*8;

pub fn calculate_matching_polynomial_pointer(graph: BinaryGraph) -> [u64; size_of::<usize>()*8] {
   let poly: &mut [u64; POLY_SIZE] = &mut [0; POLY_SIZE];
   for i in 0..size_of::<usize>()*8 {
       poly[i] = 0;
   }
   _calculate_matching_polynomial_static(graph, poly);
   return *poly
}

pub fn calculate_matching_polynomial_pointer_addresses(graph: BinaryGraph) -> ([u64; size_of::<usize>()*8], Vec<usize>) {
        let poly: &mut [u64; POLY_SIZE] = &mut [0; POLY_SIZE];
        let mut addresses: Vec<usize> = Vec::new(); 
        let address = 1;

        for i in 0..size_of::<usize>()*8 {
            poly[i] = 0;
        }
        _calculate_matching_polynomial_static_addresses(graph, poly, address, &mut addresses);
        return (*poly, addresses)
}

pub fn calculate_matching_polynomial_adaptive(graph: BinaryGraph) -> ([u64; size_of::<usize>()*8]) {
    let poly: &mut [u64; POLY_SIZE] = &mut [0; POLY_SIZE];
    let hermites: &mut [f32; POLY_SIZE * POLY_SIZE] = &mut [0.0; POLY_SIZE * POLY_SIZE];

    // cache the Hermite polynomials
    for i in 0..graph.initial_graph_size() + 1 {
        // make a monomial that is of order i
        let data = [0.0; POLY_SIZE];

        let mut coeffics = data[..i].to_vec();
        coeffics.push(1.0);
        let monomial = Polynomial::new(coeffics);
        let hermite = poly2herme(&monomial);
        hermites[i*POLY_SIZE..i*POLY_SIZE + i + 1].copy_from_slice(&hermite.data());
    }
    // set the blank polynomial
    for i in 0..size_of::<usize>()*8 {
        poly[i] = 0;
    }
    
    if graph.density() >= 0.5 {
        _calculate_matching_polynomial_static_adaptive(graph, poly, true, hermites);
    } else {
        _calculate_matching_polynomial_static_adaptive(graph, poly, false, hermites);

    }
    return *poly
}

// From here are the recursive functions called by the functions above. 

/// the following function assumes that the graph is to have its polynomial
/// calculated in the standard basis.
fn _calculate_matching_polynomial_static<T: Graph>(graph: T, poly: &mut [u64; POLY_SIZE]) {
    if graph.edgeless() {
        let node_count = graph.edgeless_node_count();
        poly[node_count] += 1;
        //println!("edge reached!");
    } else {
        let (graph_prime, graph_prime_prime) = graph.get_graph_primes();
        _calculate_matching_polynomial_static(graph_prime_prime, poly);
        println!("poly: {:?}", poly);
        _calculate_matching_polynomial_static(graph_prime, poly); // put this at the end as I think
    }
}

fn _calculate_matching_polynomial_static_addresses<T: Graph>(graph: T, poly: &mut [u64; POLY_SIZE], current_address: usize, addresses: &mut Vec<usize>) {
    if graph.edgeless() {
        let node_count = graph.edgeless_node_count();
        poly[node_count] += 1;
        addresses.push(current_address);
    } else {
        let (graph_prime, graph_prime_prime) = graph.get_graph_primes();
        _calculate_matching_polynomial_static_addresses(graph_prime_prime, poly, current_address << 1, addresses);
        _calculate_matching_polynomial_static_addresses(graph_prime, poly, (current_address << 1) + 1, addresses); // put this at the end as I think
    }
}

/// the following function, aimed at being a drop-in replacement for
/// calculate_matching_polynomial_static, assumes that the graph is to have its
/// polynomial calculated adaptively as the density of the relevant subraph
/// changes over the course of the algorithm.
fn _calculate_matching_polynomial_static_adaptive<T: Graph>(graph: T, poly: &mut [u64; POLY_SIZE], complement: bool, hermites: &[f32; POLY_SIZE * POLY_SIZE]) {
    if graph.edgeless() {
        let node_count = graph.edgeless_node_count();
        if complement {
            // run the update as if in the Hermite basis
            let hermite_coeffics = hermites[node_count*POLY_SIZE..node_count*POLY_SIZE + node_count + 1].to_vec();
            for i in 0..node_count + 1 {
                poly[i] += hermite_coeffics[i] as u64;
            }
        } else {
            // run the update as if in the standard basis
            poly[node_count] += 1;
        }
    } else {
        let (graph_prime, graph_prime_prime) = graph.get_graph_primes();

        // Handle graph_prime_prime
        if graph_prime_prime.density() >= 0.5 {
            let graph_prime_prime = graph_prime_prime.complement();
            _calculate_matching_polynomial_static_adaptive(graph_prime_prime, poly, false, hermites);
        } else {
            _calculate_matching_polynomial_static_adaptive(graph_prime_prime, poly, false, hermites);
        }
        // Handle graph_prime
        if graph_prime.density() >= 0.5 {
            let graph_prime = graph_prime.complement();
            _calculate_matching_polynomial_static_adaptive(graph_prime, poly, false, hermites);
        } else {
            _calculate_matching_polynomial_static_adaptive(graph_prime, poly, false, hermites);
        }
    }
}

pub fn _calculate_matching_polynomial_binary<T: Graph>(graph: T) -> Polynomial<u64> {
    // the base case for the process is that the graph is edgeless.
    // This means that, of the remaining nodes, each of their integer
    // representations is a power of two.
    if graph.edgeless() { // i.e. we're at the base case.
        // produce a sequence of coefficients the same length as the number of vertices
        //println!("Hit edgeless graph! with {} nodes", graph.edgeless_node_count());
        let mut coeffics = vec![0; graph.edgeless_node_count()];
        coeffics.push(1);
        let poly = Polynomial::new(coeffics);
        //println!("Polynomial: {:?}", poly);
        //println!("graph {:?}", graph.data);
        return poly
    } else {
        // get G' and G''
        // G' = G - an edge
        // G'' = G - the nodes connected to the edge removed to get G'
        let (graph_prime, graph_prime_prime) = graph.get_graph_primes();

        let poly_1 = _calculate_matching_polynomial_binary(graph_prime);
        let poly_2 = _calculate_matching_polynomial_binary(graph_prime_prime);
        let poly = poly_1 + poly_2;
        return poly
    }
} 
