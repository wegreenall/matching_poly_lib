use std::fmt::Debug; 

pub trait Graph: Clone + Debug {
    fn remove_node(&mut self, node: usize, graph_size : usize);
    fn remove_edge(&mut self, node1: usize, node2: usize, graph_size: usize);
    fn edgeless_node_count(&self) -> usize;
    fn graph_size(&self) -> usize;
    fn edge_count(&self) -> usize;
    fn edgeless(&self) -> bool;
    fn initial_graph_size(&self) -> usize;
    fn get_graph_primes(self) -> (Self, Self);
    fn get_relevant_edge(&self) -> (usize, usize, usize);
    fn density(&self) -> f32;
    fn complement(&self) -> Self;
}

pub fn get_deck<T: Graph>(graph: T) -> Vec<T>{
    let mut deck = Vec::<T>::new();
    let graph_size = graph.graph_size();
    for i in 0..graph_size {
        //println!("current graph: {}", current_graph);
        let mut current_graph = graph.clone();
        current_graph.remove_node(i, graph_size); 
        deck.push(current_graph.clone());
    }
    deck
}
