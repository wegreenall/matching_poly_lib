//pub mod matching;
pub mod graph_matching;
pub mod matching_raw_memory;
pub mod weighted_graph_matching;
pub mod petgraph;
pub mod binary_graph_matching;

pub use graph_matching::{Graph, calculate_matching_polynomial_pointer, _calculate_matching_polynomial_binary};
//use matching::calculate_matching_polynomial_pointer;
//use matching::{
    //get_deck, get_matching_polies_stable_graph,
//};
use crate::weighted_graph_matching::{WeightedGraph, _calculate_weighted_matching_polynomial_binary, get_weighted_deck};

use matching_raw_memory::{calculate_matching_polynomial_raw, GraphData, GraphProperties};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn static_polynomial_calculation_hard() {
        let data = [
            83397, 39080, 19209, 8859, 4503, 3121, 1130, 734, 364, 147, 103, 41, 28, 8, 6, 3, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let graph = Graph::from(data); // the fully connected graph
        let graph_size = graph.graph_size();
        let matching_poly = calculate_matching_polynomial_pointer(graph);
        //println!("matching poly: {:?}", matching_poly);
        let graph2 = Graph::from(data); // the fully connected graph
        let matching_poly_2 = _calculate_matching_polynomial_binary(graph2);

        //println!("matching poly 2: {:?}", matching_poly_2);
        assert_eq!(&matching_poly[..=graph_size], matching_poly_2.data());
    }

    #[test]
    fn static_polynomial_calculation_fully_connected() {
        let data = [
            512, 256, 128, 64, 32, 16, 8, 4, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let graph = Graph::from(data); // the fully connected graph
        let graph_size = graph.graph_size();
        let matching_poly = calculate_matching_polynomial_pointer(graph);

        let graph2 = Graph::from(data); // the fully connected graph
        let matching_poly_2 = _calculate_matching_polynomial_binary(graph2);

        assert_eq!(&matching_poly[..=graph_size], matching_poly_2.data());
    }

    #[test]
    fn raw_polynomial_calculation_hard() {
        let data = [
            83397, 39080, 19209, 8859, 4503, 3121, 1130, 734, 364, 147, 103, 41, 28, 8, 6, 3, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let graph = data; // the fully connected graph
        let graph_size = graph.graph_size();
        let matching_poly = calculate_matching_polynomial_raw(graph);

        let graph2 = Graph::from(data); // the fully connected graph
        let matching_poly_2 = _calculate_matching_polynomial_binary(graph2);

        //println!("matching poly 2: {:?}", matching_poly_2);
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

        let graph2 = Graph::from(data); // the fully connected graph
        let matching_poly_2 = _calculate_matching_polynomial_binary(graph2);

        //println!("matching poly 2: {:?}", matching_poly_2);
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
        let graph_1 = Graph::from(data); // the fully connected graph
        let graph_2 = Graph::from(data); // the fully connected graph
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
        let graph_1 = Graph::from(data); // the fully connected graph
        let graph_2 = Graph::from(data); // the fully connected graph
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
}
