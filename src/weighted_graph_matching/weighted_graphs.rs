use crate::binary_graph_matching::BinaryGraph;
use crate::traits::Graph;
use polynomial::Polynomial;
use std::mem;

const MAX_NODES: usize = mem::size_of::<usize>()*8;

#[derive(Debug, Clone, Copy)]
pub struct WeightedGraph{
    pub graph: BinaryGraph,
    pub weights: [f32; MAX_NODES*MAX_NODES],
}

impl WeightedGraph {
    pub const fn new() -> WeightedGraph {
         let blank_weights = [0.0 as f32; MAX_NODES*MAX_NODES];
         WeightedGraph {
             graph: BinaryGraph::new(),
             weights: blank_weights,
         }
   }

   pub fn from(data: [usize; mem::size_of::<usize>()*8], weights: [f32; MAX_NODES*MAX_NODES]) -> WeightedGraph {
        WeightedGraph {
            graph: BinaryGraph::from(data),
            weights,
        }
    }

    pub fn graph_size(&self) -> usize {
        self.graph.graph_size()
    }

   pub fn get_graph_primes(self) -> (WeightedGraph, WeightedGraph, f32) {
       let mut new_graph = self.graph; // Should copy
       let mut new_graph2 = self.graph; // Should copy

       // get the relevant edge
       let (start_node, end_node, graph_size) = self.graph.get_relevant_edge();
       
       // G' = G - e
       new_graph.remove_edge(start_node, end_node, graph_size);
       
       // G'' = G - {v, w} where {w, v} are the nodes connected to e
       new_graph2.remove_node(start_node, graph_size);
       new_graph2.remove_node(end_node, graph_size);
       let weight = self.weights[start_node * self.graph.initial_graph_size() + end_node];

       // clean the weights for each of the new ones
       (WeightedGraph{graph: new_graph, weights: self.weights},
        WeightedGraph{graph: new_graph2, weights: self.weights},
        weight) // pull the right weight from the weights
   }
}

pub fn get_weighted_deck(weighted_graph: &WeightedGraph) -> Vec<WeightedGraph> {
    let mut deck = Vec::<WeightedGraph>::new();
    let graph_size = weighted_graph.graph.graph_size();
    for i in 0..graph_size {
        //println!("current graph: {}", current_graph);
        let mut current_graph = weighted_graph.graph.clone();
        current_graph.remove_node(i, graph_size); 
        let current_weighted_graph = WeightedGraph{graph: current_graph, weights: weighted_graph.weights};
        deck.push(current_weighted_graph.clone());
    }
    deck
}
/// The weighted matching polynomial is defined via the edge-deletion recurrence:
/// Q(G, x) = Q(G - e, x) + w(e)^2 * Q(G - N(e), x) 
/// where e is an edge in G, and N(e) is the pair of nodes connected to e,
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

        let poly_1 = _calculate_weighted_matching_polynomial_binary(graph_prime);
        let poly_2 = _calculate_weighted_matching_polynomial_binary(graph_prime_prime);
         
        // convert the weight to a 1d polynomial to make it multiplicable
        //let squared_weight = weight.checked_pow(2);
        let squared_weight = weight.powi(2);
        //let squared_weight = if let x = squared_weight { x } else { return Polynomial::new(vec![0.0]) };
        let new_poly = Polynomial::new(vec![squared_weight]);
        let poly = poly_1 + new_poly * poly_2;
        return poly
    }
} 
