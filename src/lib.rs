pub mod matching;
pub mod matching_raw_memory;

use matching::calculate_matching_polynomial_static;
use matching::{
    get_deck, get_matching_polies_stable_graph, Graph, _calculate_matching_polynomial_binary,
};

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
        let matching_poly = calculate_matching_polynomial_static(graph);
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
        let matching_poly = calculate_matching_polynomial_static(graph);

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
}
