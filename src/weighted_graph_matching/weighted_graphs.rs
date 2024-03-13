use crate::binary_graph_matching::BinaryGraph;
use crate::traits::Graph;
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
