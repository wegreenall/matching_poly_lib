pub trait Graph {
    fn remove_node(&mut self, node: usize, graph_size : usize);
    fn remove_edge(&mut self, node1: usize, node2: usize, graph_size: usize);
    fn edgeless_node_count(&self) -> usize;
    fn graph_size(&self) -> usize;
    fn edge_count(&self) -> usize;
    fn edgeless(&self) -> bool;
    fn initial_graph_size(&self) -> usize;
    fn get_graph_primes(self) -> (Box<Self>, Box<Self>);
    fn get_relevant_edge(&self) -> (usize, usize, usize);
    fn density(&self) -> f32;
}

