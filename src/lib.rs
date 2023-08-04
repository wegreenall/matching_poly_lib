//pub mod matching;
pub mod graph_matching;
pub mod matching_raw_memory;
pub mod weighted_graph_matching;
pub mod petgraph;
pub mod binary_graph_matching;
pub mod polynomials;

use traits::Graph;
use polynomial::Polynomial;
use polynomials::{poly2herme,  hermadd, hermemulx , herme2poly};
use crate::weighted_graph_matching::{WeightedGraph, _calculate_weighted_matching_polynomial_binary, get_weighted_deck, weight_from_address, weighted_coefficient_calculation};

pub use binary_graph_matching::{calculate_matching_polynomial_pointer, calculate_matching_polynomial_pointer_addresses, calculate_matching_polynomial_adaptive};
pub use binary_graph_matching::{BinaryGraph, _calculate_matching_polynomial_binary};
use matching_raw_memory::{calculate_matching_polynomial_raw, GraphProperties};
use std::mem::size_of;

pub mod traits;
#[cfg(test)]
mod tests {
    use super::*;

    //#[test]
    fn static_polynomial_calculation_hard() {
        let data = [
            83397, 39080, 19209, 8859, 4503, 3121, 1130, 734, 364, 147, 103, 41, 28, 8, 6, 3, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let graph = BinaryGraph::from(data); // the fully connected graph
        let graph_size = graph.graph_size();
        let matching_poly = calculate_matching_polynomial_pointer(graph);
        let graph2 = BinaryGraph::from(data); // the fully connected graph
        let matching_poly_2 = _calculate_matching_polynomial_binary(graph2);

        assert_eq!(&matching_poly[..=graph_size], matching_poly_2.data());
    }

    #[test]
    fn static_polynomial_calculation_fully_connected() {
        let data = [
            512, 256, 128, 64, 32, 16, 8, 4, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let graph = BinaryGraph::from(data); // the fully connected graph
        let graph_size = graph.graph_size();
        let matching_poly = calculate_matching_polynomial_pointer(graph);

        let graph2 = BinaryGraph::from(data); // the fully connected graph
        let matching_poly_2 = _calculate_matching_polynomial_binary(graph2);

        assert_eq!(&matching_poly[..=graph_size], matching_poly_2.data());
    }

    //#[test]
    fn raw_polynomial_calculation_hard() {
        let data = [
            83397, 39080, 19209, 8859, 4503, 3121, 1130, 734, 364, 147, 103, 41, 28, 8, 6, 3, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let graph = data; // the fully connected graph
        let graph_size = graph.graph_size();
        let matching_poly = calculate_matching_polynomial_raw(graph);

        let graph2 = BinaryGraph::from(data); // the fully connected graph
        let matching_poly_2 = _calculate_matching_polynomial_binary(graph2);

        assert_eq!(&matching_poly[..=graph_size], matching_poly_2.data());
    }

    #[test]
    fn raw_polynomial_calculation_fully_connected() {
        let data = [
            512, 256, 128, 64, 32, 16, 8, 4, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];

        let graph = data; // the fully connected graph
        let graph_size = graph.graph_size();
        let matching_poly = calculate_matching_polynomial_raw(graph);

        let graph2 = BinaryGraph::from(data); // the fully connected graph
        let matching_poly_2 = _calculate_matching_polynomial_binary(graph2);

        assert_eq!(&matching_poly[..=graph_size], matching_poly_2.data());
    }

    #[test]
    /// Tests the weighted matching polynomial setup
    /// Remember, the weights need to be passed as the whole matrix fit in to
    /// the full 4096 set of weights.
    fn weighted_polynomial_test() {
        let data = [
            10, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let true_weights_1: [f32; 16] = [1.0, 1.0, 1.0, 1.0, 
                              1.0, 1.0, 1.0, 1.0,
                              1.0, 1.0, 1.0, 1.0, 
                              1.0, 1.0, 1.0, 1.0];
        let true_weights_2: [f32; 16] = [2.0, 2.0, 2.0, 2.0, 
                              2.0, 2.0, 2.0, 2.0,
                              2.0, 2.0, 2.0, 2.0, 
                              2.0, 2.0, 2.0, 2.0];
        let mut weights_1: [f32; 4096] = [0.0; 4096];
        let mut weights_2: [f32; 4096] = [0.0; 4096];
        weights_1[..16].copy_from_slice(&true_weights_1);
        weights_2[..16].copy_from_slice(&true_weights_2);
        let weighted_graph_1 = WeightedGraph::from(data, weights_1); // the fully connected graph
        let weighted_graph_2 = WeightedGraph::from(data, weights_2); // the fully connected graph
        let graph_1 = BinaryGraph::from(data); // the fully connected graph
        let graph_2 = BinaryGraph::from(data); // the fully connected graph
                                                        //
        let weighted_matching_polynomial_1 = _calculate_weighted_matching_polynomial_binary(weighted_graph_1);
        let weighted_matching_polynomial_2 = _calculate_weighted_matching_polynomial_binary(weighted_graph_2);
        let matching_polynomial_1 = _calculate_matching_polynomial_binary(graph_1);
        let matching_polynomial_2 = _calculate_matching_polynomial_binary(graph_2);
        assert_eq!(matching_polynomial_1.data(), &[1, 0, 3, 0, 1]);
        assert_eq!(matching_polynomial_2.data(), &[1, 0, 3, 0, 1]);
        assert_eq!(weighted_matching_polynomial_1.data(), &[1.0, 0.0, 3.0, 0.0, 1.0]);
        assert_eq!(weighted_matching_polynomial_2.data(), &[16.0, 0.0, 12.0, 0.0, 1.0]);
    }

    #[test]
    fn weighted_polynomial_test_2() {
        let data = [
            9, 6, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let true_weights_1: [f32; 16] = [1.0, 1.0, 1.0, 1.0, 
                              1.0, 1.0, 1.0, 1.0,
                              1.0, 1.0, 1.0, 1.0, 
                              1.0, 1.0, 1.0, 1.0];
        let true_weights_2: [f32; 16] = [2.0, 2.0, 2.0, 2.0, 
                              2.0, 2.0, 2.0, 2.0,
                              2.0, 2.0, 2.0, 2.0, 
                              2.0, 2.0, 2.0, 2.0];
        let mut weights_1: [f32; 4096] = [0.0; 4096];
        let mut weights_2: [f32; 4096] = [0.0; 4096];
        weights_1[..16].copy_from_slice(&true_weights_1);
        weights_2[..16].copy_from_slice(&true_weights_2);
        let weighted_graph_1 = WeightedGraph::from(data, weights_1); // the fully connected graph
        let weighted_graph_2 = WeightedGraph::from(data, weights_2); // the fully connected graph
        let graph_1 = BinaryGraph::from(data); // the fully connected graph
        let graph_2 = BinaryGraph::from(data); // the fully connected graph
                                                        //
        let weighted_matching_polynomial_1 = _calculate_weighted_matching_polynomial_binary(weighted_graph_1);
        let weighted_matching_polynomial_2 = _calculate_weighted_matching_polynomial_binary(weighted_graph_2);
        let matching_polynomial_1 = _calculate_matching_polynomial_binary(graph_1);
        let matching_polynomial_2 = _calculate_matching_polynomial_binary(graph_2);
        println!("matching polynomial 2: {:?}", weighted_matching_polynomial_2.data());          
        assert_eq!(matching_polynomial_1.data(), &[1, 0, 3, 0, 1]);
        assert_eq!(matching_polynomial_2.data(), &[1, 0, 3, 0, 1]);
        assert_eq!(weighted_matching_polynomial_1.data(), &[1.0, 0.0, 3.0, 0.0, 1.0]);
        assert_eq!(weighted_matching_polynomial_2.data(), &[16.0, 0.0, 12.0, 0.0, 1.0]);
    }

    #[test]
    fn weighted_polynomial_test_3() {
        let data = [
            10, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let true_weights_1: [f32; 16] = [1.0, 1.0, 1.0, 1.0, 
                                         1.0, 1.0, 1.0, 1.0,
                                         1.0, 1.0, 1.0, 1.0, 
                                         1.0, 1.0, 1.0, 1.0];

        let true_weights_2: [f32; 16] = [1.0, 1.0, -2.0, 1.0, 
                                         1.0, 1.0, 1.0, -2.0,
                                         -2.0, 1.0, 1.0, -2.0, 
                                         1.0, -2.0, -2.0, 1.0];

        let mut weights_1: [f32; 4096] = [0.0; 4096];
        let mut weights_2: [f32; 4096] = [0.0; 4096];
        weights_1[..16].copy_from_slice(&true_weights_1);
        weights_2[..16].copy_from_slice(&true_weights_2);
        let weighted_graph_1 = WeightedGraph::from(data, weights_1); // the fully connected graph
        let weighted_graph_2 = WeightedGraph::from(data, weights_2); // the fully connected graph
        let graph_1 = BinaryGraph::from(data); // the fully connected graph
        let graph_2 = BinaryGraph::from(data); // the fully connected graph
                                                        //
        let weighted_matching_polynomial_1 = _calculate_weighted_matching_polynomial_binary(weighted_graph_1);
        let weighted_matching_polynomial_2 = _calculate_weighted_matching_polynomial_binary(weighted_graph_2);
        let matching_polynomial_1 = _calculate_matching_polynomial_binary(graph_1);
        let matching_polynomial_2 = _calculate_matching_polynomial_binary(graph_2);

        println!("weightd matching polynomial 1: {:?}", weighted_matching_polynomial_1.data());          
        println!("weighted matching polynomial 2: {:?}", weighted_matching_polynomial_2.data());          
          
        assert_eq!(matching_polynomial_1.data(), &[1, 0, 3, 0, 1]);
        assert_eq!(matching_polynomial_2.data(), &[1, 0, 3, 0, 1]);
        assert_eq!(weighted_matching_polynomial_1.data(), &[1.0, 0.0, 3.0, 0.0, 1.0]);
        assert_eq!(weighted_matching_polynomial_2.data(), &[16.0, 0.0, 12.0, 0.0, 1.0]);
    }

    #[test]
    fn poly2herme_test_1() {
        let poly = Polynomial::new(vec![0.0, 1.0, 2.0, 3.0]);
        let herm = poly2herme(&poly);
        assert_eq!(herm.data(), &[2.0, 10.0, 2.0, 3.0]);
    }

    #[test]
    fn poly2herme_test_2() {
        let poly = Polynomial::new(vec![2.0, 0.0, -4.0, 0.0, 1.0]);
        let herm = poly2herme(&poly);
        assert_eq!(herm.data(), &[1.0, 0.0, 2.0, 0.0, 1.0]);
    }

    #[test]
    fn herme2poly_test() {
        let herm = Polynomial::new(vec![2.0, 10.0, 2.0, 3.0]);
        let poly = herme2poly(&herm);
        assert_eq!(poly.data(), &[0.0, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn herme2poly_test_2() {
        let herm = Polynomial::new(vec![1.0, 0.0, 2.0, 0.0, 1.0]);
        let poly = herme2poly(&herm);
        assert_eq!(poly.data(), &[2.0, 0.0, -4.0, 0.0, 1.0]);
    }

    #[test]
    fn herme2poly_test_3() {
        // n = 0
        let herm = Polynomial::new(vec![1.0]);
        let poly = herme2poly(&herm);
        assert_eq!(poly.data(), &[1.0]);

        // n = 1
        let herm = Polynomial::new(vec![0.0, 1.0]);
        let poly = herme2poly(&herm);
        assert_eq!(poly.data(), &[0.0, 1.0]);

        // n = 2
        let herm = Polynomial::new(vec![0.0, 0.0, 1.0]);
        let poly = herme2poly(&herm);
        assert_eq!(poly.data(), &[-1.0, 0.0, 1.0]);

        // n = 3
        let herm = Polynomial::new(vec![0.0, 0.0, 0.0, 1.0]);
        let poly = herme2poly(&herm);
        assert_eq!(poly.data(), &[0.0, -3.0, 0.0, 1.0]);

        // n = 4
        let herm = Polynomial::new(vec![0.0, 0.0, 0.0, 0.0, 1.0]);
        let poly = herme2poly(&herm);
        assert_eq!(poly.data(), &[3.0, 0.0, -6.0, 0.0, 1.0]);
    }

    #[test]
    fn hermadd_test() {
        let poly_1 = Polynomial::new(vec![1.0, 2.0, 3.0]);
        let poly_2 = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
        let result = hermadd(&poly_1, &poly_2);
        assert_eq!(result.data(), &[2.0, 4.0, 6.0, 4.0]);
    }

    #[test]
    fn hermemul_test() {
        let poly_1 = Polynomial::new(vec![1.0, 2.0, 3.0]);
        let mul = hermemulx(&poly_1);
        let test_data = vec![2.0, 7.0, 2.0, 3.0];
        assert_eq!(mul.data(), &test_data);
    }

    #[test]
    fn weight_from_address_test() {
        // set up the basic, initial data
        let true_weights_1: [f32; 16] = [0.0, 0.0, 0.0, 3.0, 
                                         0.0, 0.0, 2.0, 0.0,
                                         0.0, 2.0, 0.0, 4.0, 
                                         3.0, 0.0, 4.0, 0.0];
        let mut weights_1: [f32; 4096] = [0.0; 4096];
        //let data = [
            //10, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            //0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            //0, 0, 0, 0, 0, 0, 0, 0, 0,
        //];
        // set up the weights for the graph
        weights_1[..16].copy_from_slice(&true_weights_1);
        //let weighted_graph_1 = WeightedGraph::from(data, weights_1); // the fully connected graph
        //let weighted_matching_polynomial = _calculate_weighted_matching_polynomial_binary(weighted_graph_1);
        //assert_eq!(weighted_matching_polynomial.data(), &[12.0, 0.0, 5.0, 0.0, 1.0]);

        // one weight
        let address = 0b1110; // lllr
        let weight = weight_from_address(address, &weights_1, 4);
        assert_eq!(weight, 4.0);

        // two weights
        let address = 0b100; // lllr
        let weight = weight_from_address(address, &weights_1, 4);
        assert_eq!(weight, 6.0);

        let address = 0b110;
        let weight = weight_from_address(address, &weights_1, 4);
        assert_eq!(weight, 2.0);
    }

    #[test]
    fn weight_from_address_test_2() {
        // set up the basic, initial data
        let true_weights_1: [f32; 16] = [0.0, 0.0, 3.0, 0.0, 
                                         0.0, 0.0, 0.0, 2.0,
                                         3.0, 0.0, 0.0, 4.0, 
                                         0.0, 2.0, 4.0, 0.0];
        let mut weights_1: [f32; 4096] = [0.0; 4096];

        // set up the weights for the graph
        weights_1[..16].copy_from_slice(&true_weights_1);
        //let weighted_graph_1 = WeightedGraph::from(data, weights_1); // the fully connected graph
        //let weighted_matching_polynomial = _calculate_weighted_matching_polynomial_binary(weighted_graph_1);
        //assert_eq!(weighted_matching_polynomial.data(), &[12.0, 0.0, 5.0, 0.0, 1.0]);

        // one weight
        let address = 0b1110; // lllr
        let weight = weight_from_address(address, &weights_1, 4);
        assert_eq!(weight, 4.0);

        // two weights
        let address = 0b100; // lrr
        let weight = weight_from_address(address, &weights_1, 4);
        assert_eq!(weight, 6.0);

        let address = 0b110; // llr
        let weight = weight_from_address(address, &weights_1, 4);
        assert_eq!(weight, 2.0);
    }

    #[test]
    fn get_addresses() {
        let data = [
            10, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let true_weights_1: [f32; 16] = [0.0, 0.0, 3.0, 0.0, 
                                         0.0, 0.0, 0.0, 2.0,
                                         3.0, 0.0, 0.0, 4.0, 
                                         0.0, 2.0, 4.0, 0.0];
        let mut weights_1: [f32; 4096] = [0.0; 4096];
        weights_1[..16].copy_from_slice(&true_weights_1);

        let graph = BinaryGraph::from(data);
        let mut addresses: Vec<usize>; 
        let (_, address_vec) = calculate_matching_polynomial_pointer_addresses(graph);
        addresses = address_vec;
        addresses.sort();
        let mut true_addresses = vec![0b1111, 0b1110, 0b110, 0b101, 0b0100];
        true_addresses.sort();
        assert_eq!(addresses, true_addresses);
    }

    #[test]
    fn get_polynomial_with_addresses() {
        let true_weights_1: [f32; 16] = [0.0, 0.0, 3.0, 0.0, 
                                         0.0, 0.0, 0.0, 2.0,
                                         3.0, 0.0, 0.0, 4.0, 
                                         0.0, 2.0, 4.0, 0.0];
        let mut weights_1: [f32; 4096] = [0.0; 4096];
        weights_1[..16].copy_from_slice(&true_weights_1);

        let true_addresses = vec![0b1111, 0b1110, 0b110, 0b101, 0b0100];
        let mut polynomial_coefficients = vec![1.0, 0.0, 3.0, 0.0, 1.0];
        true_addresses.iter().for_each(|x| {
            let address = *x;
            let weight = weight_from_address(address, &weights_1, 4);

            // get the number of zeroes
            polynomial_coefficients[4 - 2 * (address.count_zeros() - address.leading_zeros()) as usize] += weight - 1.0;
        });
        assert_eq!(vec![6.0, 0.0, 9.0, 0.0, 1.0], polynomial_coefficients);
    }

    #[test]
    fn get_polynomial_with_addresses_2() {
        let data = [
            0b1101, 0b110, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let true_weights: [f32; 16] = [0.0, 2.0, 0.0, 2.0, 
                                       2.0, 0.0, 2.0, 0.0,
                                       0.0, 2.0, 0.0, 2.0, 
                                       2.0, 0.0, 2.0, 0.0];
        let mut weights_1: [f32; 4096] = [0.0; 4096];
        weights_1[..16].copy_from_slice(&true_weights);

        let graph = BinaryGraph::from(data);
        let (_, addresses) = calculate_matching_polynomial_pointer_addresses(graph);
        let mut polynomial_coefficients: Vec<f32> = vec![0.0, 0.0, 0.0, 0.0, 0.0];

        addresses.iter().for_each(|x| {
            let address = *x;
            let weight = weight_from_address(address, &weights_1, 4);

            // get the number of zeroes
            let index = 4 - 2 * (address.count_zeros() - address.leading_zeros()) as usize;
            polynomial_coefficients[index] = polynomial_coefficients[index] as f32 + weight;
        });
        assert_eq!(vec![8.0, 0.0, 8.0, 0.0, 1.0], polynomial_coefficients);
    }

    #[test]
    fn test_complement_fc() {
        // fully_connected data
        let fc_data = [
            0b11111, 0b1111, 0b111, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ];
        let graph = BinaryGraph::from(fc_data);

        let complement = graph.complement();
        assert_eq!([0b10000, 0b1000, 0b100, 0b10, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0], complement.data());
    }

    #[test]
    fn test_complement_standard() {
        // now, standard data
        let standard_data = [
            0b11101, 0b1001, 0b110, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let graph_2 = BinaryGraph::from(standard_data);
        let complement_2 = graph_2.complement();
        assert_eq!([0b10010, 0b1110, 0b101, 0b10, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0], complement_2.data());
    }

    #[test]
    fn test_complement_standard_missing() {
        // now, standard data, but with a missing node or 2
        let standard_data_missing = [
            0b11001, 0b1001, 0, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let graph_3 = BinaryGraph::from_graph_subset(standard_data_missing, 5);
        // print the graph
        for element in graph_3.data().iter().take(5) {
            println!("element {:b}", element);
        }

        // print the complement
        let complement_3 = graph_3.complement();
        for element in complement_3.data().iter().take(5) {
            println!("element {:b}", element);
        }
        assert_eq!([0b10010, 0b1010, 0, 0b10, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0], complement_3.data());
    }

    #[test]
    fn matching_polynomial_static() {
        let fc_data = [
            0b1111, 0b111, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ];
        let graph = BinaryGraph::from(fc_data);

        let matching_polynomial = calculate_matching_polynomial_pointer(graph);
        assert_eq!(matching_polynomial[..graph.graph_size()+1], [3, 0, 6, 0, 1]);

    }

    #[test]
    fn matching_polynomial_static_2() {
        let standard_data = [
            0b11001, 0b1001, 0b110, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let graph = BinaryGraph::from(standard_data);

        let matching_polynomial = calculate_matching_polynomial_pointer(graph);
        println!("standard data matching polynomial {:?}", matching_polynomial);
        assert_eq!(matching_polynomial[..graph.graph_size()+1], [0, 4, 0, 5, 0, 1]);
    }

    #[test]
    fn matching_polynomial_static_3() {
        let standard_data_missing = [
            0b11001, 0b1001, 0, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let graph = BinaryGraph::from(standard_data_missing);
        let matching_polynomial = calculate_matching_polynomial_pointer(graph);
        assert_eq!(matching_polynomial[..graph.graph_size()+1], [1, 0, 4, 0, 1]);
    }

    #[test]
    fn matching_polynomial_adaptive() {
        // grpah size is 4
        let fc_data = [
            0b1111, 0b111, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0
        ];
        let graph = BinaryGraph::from(fc_data);

        let matching_polynomial = calculate_matching_polynomial_adaptive(graph);
        println!("We got this far...");
        assert_eq!(matching_polynomial[..graph.graph_size()+1], [3, 0, -6, 0, 1]);
        // grpah size is 5
        let fc_data = [
            0b11111, 0b1111, 0b111, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0
        ];
        let graph = BinaryGraph::from(fc_data);

        let matching_polynomial = calculate_matching_polynomial_adaptive(graph);
        println!("We got this far...");
        assert_eq!(matching_polynomial[..graph.graph_size()+1], [0, 15, 0, -10, 0, 1]);
    }

    #[test]
    fn matching_polynomial_adaptive_2() {
        // grpah size is 4
        // now, standard data
        let standard_data = [
            0b11001, 0b1001, 0b110, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
        ];
        let graph = BinaryGraph::from(standard_data);

        let matching_polynomial = calculate_matching_polynomial_adaptive(graph);
        assert_eq!(matching_polynomial[..graph.graph_size()+1], [0, 4, 0, -5, 0, 1]);
    }

    #[test]
    fn matching_polynomial_adaptive_3() {
        let standard_data_without_missing = [
            0b1101, 0b101, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0,
        ];
        let graph = BinaryGraph::from(standard_data_without_missing);

        let matching_polynomial = calculate_matching_polynomial_adaptive(graph);
        println!("About to test the standard data without missing");
        assert_eq!(matching_polynomial[..graph.graph_size()+1], [1, 0, -4, 0, 1]);
    }


    #[test]
    fn matching_polynomial_adaptive_4() {
        let standard_data_missing = [
            0b11001, 0b1001, 0, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0,0
        ];
        //let graph = BinaryGraph::from(standard_data_missing);
        let graph = BinaryGraph::from_graph_subset(standard_data_missing, 5);
        println!("graph complement:{:?}", graph.complement());
        //let graph = BinaryGraph::from_graph_subset(standard_data_missing, 5);

        let matching_polynomial = calculate_matching_polynomial_adaptive(graph);
        println!("About to test the standard data with missing");
        assert_eq!(matching_polynomial[..graph.graph_size()+1], [1, 0, -4, 0, 1]);
    }
    #[test]
    fn matching_polynomial_adaptive_5() {
        let chain_data = [
            0b11000, 0b1100, 0b110, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ];
        let graph = BinaryGraph::from(chain_data);

        let matching_polynomial = calculate_matching_polynomial_adaptive(graph);
        println!("Number of nodes: {}", graph.graph_size());
        assert_eq!(matching_polynomial[..graph.graph_size()+1], [0, 3, 0, -4, 0, 1]);
    }

    #[test]
    fn matching_polynomial_adaptive_6() {
        let chain_data = [
            0b110000, 0b11000, 0b1100, 0b110, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ];
        let graph = BinaryGraph::from(chain_data);

        let matching_polynomial = calculate_matching_polynomial_adaptive(graph);
        println!("Number of nodes: {}", graph.graph_size());
        assert_eq!(matching_polynomial[..graph.graph_size()+1], [-1, 0, 6, 0, -5, 0, 1]);
    }
    #[test]
    fn matching_polynomial_adaptive_7() {
        let chain_data = [
            0b1100000, 0b110000, 0b11000, 0b1100, 0b110, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ];
        let graph = BinaryGraph::from(chain_data);

        let matching_polynomial = calculate_matching_polynomial_adaptive(graph);
        //println!("");
        assert_eq!(matching_polynomial[..graph.graph_size()+1], [0, -4, 0, 10, 0, -6, 0, 1]);
    }
    #[test]
    fn matching_polynomial_adaptive_8() {
        let chain_data = [
            0b11000000, 0b1100000, 0b110000, 0b11000, 0b1100, 0b110, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ];
        let graph = BinaryGraph::from(chain_data);

        let matching_polynomial = calculate_matching_polynomial_adaptive(graph);
        //println!("About to test the standard data without missing");
        assert_eq!(matching_polynomial[..graph.graph_size()+1], [1, 0, -10, 0, 15, 0, -7, 0, 1]);
    }


    #[test]
    fn test_density() {
        let fc_data = [
            0b11111, 0b1111, 0b111, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0
        ];
        let graph = BinaryGraph::from(fc_data);
        assert_eq!(graph.density(), 1.0);
    }

    #[test]
    fn test_edge_count() {
        let fc_data = [
            0b11111, 0b1111, 0b111, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0
        ];
        let graph = BinaryGraph::from(fc_data);
        assert_eq!(graph.edge_count(), 10);
    }

    #[test]
    fn test_edge_count_2() {
        let standard_data_missing = [
            0b11001, 0b1001, 0, 0b11, 0b1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0,
        ];
        let graph = BinaryGraph::from(standard_data_missing);
        assert_eq!(graph.edge_count(), 4);
    }

    #[test]
    fn test_coefficient_calculation() {
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
        let weighted_coefficients = weighted_coefficient_calculation(poly, &weights, 4, 2);
        //assert_eq!(weighted_coefficients[..5], [8.0, 0.0, 8.0, 0.0, 1.0]);
    }


}
