use crate::{binary_graph_matching::BinaryGraph, polynomials::poly2herme, polynomials::herme2poly};
use std::mem::size_of;
use crate::traits::Graph;
use polynomial::Polynomial;

fn modulo(x: i64, y: i64) -> i64 {
    ((x % y) + y) % y
}
fn ceil(x: i64, y: i64) -> i64 {
    (x + y - 1) / y
}

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
        for i in 0..POLY_SIZE {
            poly[i] = 0;
        }
        _calculate_matching_polynomial_static(graph, poly);
        return *poly
}

pub fn calculate_matching_polynomial_pointer_addresses(graph: BinaryGraph) -> ([u64; size_of::<usize>()*8], Vec<usize>) {
        let poly: &mut [u64; POLY_SIZE] = &mut [0; POLY_SIZE];
        let mut addresses: Vec<usize> = Vec::new(); 
        let address = 1;

        for i in 0..POLY_SIZE {
            poly[i] = 0;
        }
        _calculate_matching_polynomial_static_addresses(graph, poly, address, &mut addresses);
        return (*poly, addresses)
}

pub fn calculate_matching_polynomial_adaptive(graph: BinaryGraph) -> [i64; POLY_SIZE] {
        let poly: &mut [i64; POLY_SIZE] = &mut [0; POLY_SIZE];
        let hermites: &mut [f32; POLY_SIZE * POLY_SIZE] = &mut [0.0; POLY_SIZE * POLY_SIZE];

        // cache the Hermite polynomials
        for i in 0..graph.initial_graph_size() + 1 {
            // make a monomial that is of order i
            let data = [0.0; POLY_SIZE];
            let mut coeffics = data[..i].to_vec();
            coeffics.push(1.0);
            let monomial = Polynomial::new(coeffics);

            // and treat it like it is a Hermite polynomial
            let hermite = herme2poly(&monomial);
            hermites[i*POLY_SIZE..i*POLY_SIZE + i + 1].copy_from_slice(&hermite.data());
        }

        // set the blank polynomial
        for i in 0..POLY_SIZE {
            poly[i] = 0;
        }
        
        // now run the recursive function that does it
        _calculate_matching_polynomial_static_adaptive(graph, poly, false, hermites, 1);
        return *poly
}

// From here are the recursive functions called by the functions above. 

/// the following function assumes that the graph is to have its polynomial
/// calculated in the standard basis.
fn _calculate_matching_polynomial_static<T: Graph>(graph: T, poly: &mut [u64; POLY_SIZE]) {
    if graph.edgeless() {
        let node_count = graph.edgeless_node_count();
        poly[node_count] += 1;
    } else {
        let (graph_prime, graph_prime_prime) = graph.get_graph_primes();
        _calculate_matching_polynomial_static(graph_prime_prime, poly);
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
fn _calculate_matching_polynomial_static_adaptive<T: Graph>(mut graph: T, poly: &mut [i64; POLY_SIZE], mut complement: bool, hermites: &[f32; POLY_SIZE * POLY_SIZE], mut sign_coeffic: i64) {
    let graph_density = graph.density();
    if graph_density >= 0.5 && !complement {
        complement = true;
        sign_coeffic = (-1 as i64).pow((graph.initial_graph_size() as u32 - graph.graph_size() as u32)/2) as i64;
        graph = graph.complement();
    }
    if graph.edgeless() {
        let graph_size = graph.graph_size();
        if complement {
            // run the update as if in the Hermite basis
            let hermite_coeffics = hermites[(graph_size) * POLY_SIZE..(graph_size) * POLY_SIZE + graph_size + 1].to_vec();
            for i in 0..graph_size + 1 {
                poly[i] += sign_coeffic * hermite_coeffics[i] as i64;
            }
        } else {
                // run the update as if in the standard basis
                poly[graph_size] += (-1 as i64).pow((graph.initial_graph_size() as u32 - graph_size as u32)/2);
        }
    } else {
        let (graph_prime, graph_prime_prime) = graph.get_graph_primes();
        _calculate_matching_polynomial_static_adaptive(graph_prime_prime,
                                                           poly,
                                                           complement,
                                                           hermites, 
                                                           sign_coeffic);
        // Handle graph_prime
        _calculate_matching_polynomial_static_adaptive(graph_prime,
                                                       poly,
                                                       complement,
                                                       hermites, 
                                                       sign_coeffic);
    }
}

pub fn _calculate_matching_polynomial_binary<T: Graph>(graph: T) -> Polynomial<u64> {
    // the base case for the process is that the graph is edgeless.
    // This means that, of the remaining nodes, each of their integer
    // representations is a power of two.
    if graph.edgeless() { // i.e. we're at the base case.
        // produce a sequence of coefficients the same length as the number of vertices
        let mut coeffics = vec![0; graph.edgeless_node_count()];
        coeffics.push(1);
        let poly = Polynomial::new(coeffics);
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
