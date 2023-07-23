//pub mod matching;
pub mod graph_matching;
pub mod matching_raw_memory;
pub mod weighted_graph_matching;
pub mod petgraph;
pub mod binary_graph_matching;
pub mod polynomials;


pub use graph_matching::calculate_matching_polynomial_pointer;
use traits::Graph;
use polynomial::Polynomial;
use polynomials::{poly2herme,  hermadd, hermemulx , herme2poly};
pub use binary_graph_matching::{BinaryGraph, _calculate_matching_polynomial_binary};
use crate::weighted_graph_matching::{WeightedGraph, _calculate_weighted_matching_polynomial_binary, get_weighted_deck, weight_from_address};

use matching_raw_memory::{calculate_matching_polynomial_raw, GraphData, GraphProperties};

mod traits;
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
        let matching_poly_2 = _calculate_matching_polynomial_binary(Box::new(graph2));

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
        let matching_poly_2 = _calculate_matching_polynomial_binary(Box::new(graph2));

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
        let matching_poly_2 = _calculate_matching_polynomial_binary(Box::new(graph2));

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
        let matching_poly_2 = _calculate_matching_polynomial_binary(Box::new(graph2));

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
        let matching_polynomial_1 = _calculate_matching_polynomial_binary(Box::new(graph_1));
        let matching_polynomial_2 = _calculate_matching_polynomial_binary(Box::new(graph_2));
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
        let matching_polynomial_1 = _calculate_matching_polynomial_binary(Box::new(graph_1));
        let matching_polynomial_2 = _calculate_matching_polynomial_binary(Box::new(graph_2));
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
        let matching_polynomial_1 = _calculate_matching_polynomial_binary(Box::new(graph_1));
        let matching_polynomial_2 = _calculate_matching_polynomial_binary(Box::new(graph_2));

        println!("weightd matching polynomial 1: {:?}", weighted_matching_polynomial_1.data());          
        println!("weighted matching polynomial 2: {:?}", weighted_matching_polynomial_2.data());          
          
        assert_eq!(matching_polynomial_1.data(), &[1, 0, 3, 0, 1]);
        assert_eq!(matching_polynomial_2.data(), &[1, 0, 3, 0, 1]);
        assert_eq!(weighted_matching_polynomial_1.data(), &[1.0, 0.0, 3.0, 0.0, 1.0]);
        assert_eq!(weighted_matching_polynomial_2.data(), &[16.0, 0.0, 12.0, 0.0, 1.0]);
    }

    #[test]
    fn poly2herme_test() {
        let poly = Polynomial::new(vec![0.0, 1.0, 2.0, 3.0]);
        let herm = poly2herme(&poly);
        assert_eq!(herm.data(), &[2.0, 10.0, 2.0, 3.0]);
    }

    #[test]
    fn herme2poly_test() {
        let herm = Polynomial::new(vec![2.0, 10.0, 2.0, 3.0]);
        let poly = herme2poly(&herm);
        assert_eq!(poly.data(), &[0.0, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn hermadd_test() {
        let poly_1 = Polynomial::new(vec![1.0, 2.0, 3.0]);
        let poly_2 = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
        let result = hermadd(&poly_1, &poly_2);
        //println!("result: {:?}", result.data());
        assert_eq!(result.data(), &[2.0, 4.0, 6.0, 4.0]);
    }

    #[test]
    fn hermemul_test() {
        let poly_1 = Polynomial::new(vec![1.0, 2.0, 3.0]);
        let mul = hermemulx(&poly_1);
        let test_data = vec![2.0, 7.0, 2.0, 3.0];
        println!("mul: {:?}", mul.data());
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
        println!("true weight for address {:b}: {}", address, weight);
        assert_eq!(weight, 4.0);
        println!("\n");

        // two weights
        let address = 0b100; // lllr
        let weight = weight_from_address(address, &weights_1, 4);
        println!("true weight for address {:b}: {}", address, weight);
        assert_eq!(weight, 6.0);
    }

}
