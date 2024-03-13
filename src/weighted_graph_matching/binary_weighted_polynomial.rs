use polynomial::Polynomial;
use crate::WeightedGraph;
use crate::traits::Graph;

/// The weighted matching polynomial is defined via the edge-deletion recurrence:
/// Q(G, x) = Q(G - e, x) + w(e)^2 * Q(G - N(e), x) 
/// where e is an edge in G, N(e) is the pair of nodes connected to e,
/// and w(e) is the weight associated with e. 
pub fn _calculate_weighted_matching_polynomial_binary(weighted_graph: WeightedGraph) -> Polynomial<f32> {
    // the base case for the process is that the graph is edgeless.
    // This means that, of the remaining nodes, each of their integer
    // representations is a power of two.
    if weighted_graph.graph.edgeless() { // i.e. we're at the base case.
        // produce a sequence of coefficients the same length as the number of vertices
        let mut coeffics = vec![0.0; weighted_graph.graph.edgeless_node_count()];
        coeffics.push(1.0);
        let poly = Polynomial::new(coeffics);
        return poly

    } else {
        // get G' and G''
        // G' = G - an edge
        // G'' = G - the nodes connected to the edge removed to get G'
        let (graph_prime, graph_prime_prime, weight) = weighted_graph.get_graph_primes();
        //println!("{}", weight);

        let poly_1 = _calculate_weighted_matching_polynomial_binary(graph_prime);
        let poly_2 = _calculate_weighted_matching_polynomial_binary(graph_prime_prime);
         
        // convert the weight to a 1d polynomial to make it multiplicable
        let new_poly = Polynomial::new(vec![weight]);
        let poly = poly_1 + new_poly * poly_2;
        return poly
    }
} 
