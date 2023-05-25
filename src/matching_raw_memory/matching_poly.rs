//#[allow(unreachable_code)]
use std::mem::size_of;

/// This file contains functions that
/// are used to calculate the matching polynomials for graphs.
/// This method will use static memory for the polynomial.
/// Since the computation graph for the polynomial requires
/// no reads until returning to the caller, and the 
/// write operations are done in any order (since we 
/// calculate the signless matching polynomial), we can
/// avoid large numbers of allocations by using static memory.
///

// Set up the size of a polynomial
const POLY_SIZE: usize = size_of::<usize>()*8;
const MAX_NODES: usize = size_of::<usize>()*8;

pub type GraphData = [usize; POLY_SIZE];

// set up the polynomial memory
static mut POLY: [u64; POLY_SIZE] = [0; POLY_SIZE];

// set up the graph memory
static mut GRAPH_MEMORY: [GraphData; POLY_SIZE * POLY_SIZE] = [[0; POLY_SIZE]; POLY_SIZE * POLY_SIZE];

pub fn calculate_matching_polynomial_raw(graph: GraphData) -> [u64; size_of::<usize>()*8] {
   //let graph_memory: &mut [GraphData; POLY_SIZE * POLY_SIZE] = &mut [[0; POLY_SIZE]; POLY_SIZE * POLY_SIZE];
   //let poly: &mut [u64; POLY_SIZE] = &mut [0; POLY_SIZE];
   unsafe {
       // clear the memory
       for i in 0..size_of::<usize>()*8 {
           POLY[i] = 0;
       }
       for i in 0..POLY_SIZE * POLY_SIZE {
           GRAPH_MEMORY[i] = graph;
       }
       // copy the graph into the first slot
       GRAPH_MEMORY[0] = graph;
       _calculate_matching_polynomial_raw(0);
       return POLY;
   }
}

unsafe fn _calculate_matching_polynomial_raw(depth: usize) {
    let graph = GRAPH_MEMORY[depth];

    if graph.edgeless() {
        let node_count = graph.edgeless_node_count();
        POLY[node_count] += 1;
    } else { // not in the base case
        //GRAPH_MEMORY[depth] = graph;
        GRAPH_MEMORY[depth+1] = graph;

        let (start_node, end_node, graph_size) = graph.get_relevant_edge();

        // G'
        GRAPH_MEMORY[depth].remove_edge(start_node, end_node, graph_size);

        // G''
        GRAPH_MEMORY[depth+1].remove_nodes(start_node, end_node, graph_size);

        //GRAPH_MEMORY[depth+1].remove_node(start_node, graph_size);
        //GRAPH_MEMORY[depth+1].remove_node(end_node, graph_size);

        _calculate_matching_polynomial_raw(depth+1); // put this at the end as I think
        _calculate_matching_polynomial_raw(depth);
    }
}

pub trait GraphProperties {
    fn get_relevant_edge(&self) -> (usize, usize, usize);
    fn edgeless(&self) -> bool;
    fn graph_size(&self) -> usize;
    fn remove_nodes(&mut self, start_node: usize, end_node:usize, graph_size : usize);
    fn remove_node(&mut self, node: usize, graph_size : usize);
    fn remove_edge(&mut self, node1: usize, node2: usize, graph_size: usize);
    fn edgeless_node_count(&self) -> usize;
}
impl GraphProperties for GraphData {
    fn edgeless(&self) -> bool {
        self.iter()
            .all(|x| x == &(0 as usize) || x.is_power_of_two())
    }

    fn graph_size(&self) -> usize{
        self.iter()
            .filter(|x| x> &&(0 as usize)) // i.e. get the ones that are valid
            .count()
    }

    /// Removes two nodes in s single pass, as we have to loop through the whole of
    /// the graph anyway
    fn remove_nodes(&mut self, start_node: usize, end_node: usize, graph_size : usize) {
        // remove node from adjacency list
        self[start_node] = 0;
        self[end_node] = 0;

        // Now remove its connected edges
        self
            .iter_mut()
            .for_each(|x| *x &= !((1<<graph_size.saturating_sub(start_node+1)) | (1<<graph_size.saturating_sub(end_node+1))));
    }

    fn remove_node(&mut self, node: usize, graph_size : usize) {
        // remove node from adjacency list
        self[node] = 0;

        // Now remove its connected edges
        self
            .iter_mut()
            .for_each(|x| *x &= !(1<<graph_size.saturating_sub(node+1)));
    }

    fn remove_edge(&mut self, start_node: usize, end_node: usize, graph_size: usize) {
        // remove edge from adjacency list
        let shift = graph_size.saturating_sub(end_node);
         
        // zero the ith from the right
        self[start_node] &= (!(1 << (shift-1))) as usize; 
    }

    fn edgeless_node_count(&self) -> usize {
        // count the number of nodes in the graph
        // To count the number of nodes in th graph,
        // count the number of the bits in the graph
        // Since, by assumption, the graph is edgeless,
        // we can just count the number of 1s in the graph
        self.iter()
            .sum::<usize>()
            .count_ones() as usize
    }

    fn get_relevant_edge(&self) -> (usize, usize, usize) {
        // since the nodes are ordered in INCREASING order of degree, we can
        // just drop the last(right-most in the binary representation)
        // edge we find, on the first still-relevant node.
        // starting_node: the index of the first node that still has edges from it

        // if drop_most_connected_edge is true, we delete the first connected edge
        // we find on the first relevant node. Otherwise, we delete the "last"
        // edge for the first relevant node.
        let drop_first_connected_edge: bool = false;

        let starting_node = self
            .iter()
            .enumerate()
            .filter(|(_, x)| (x > &&(0 as usize)))
            .filter(|(_, x)| !(x.is_power_of_two()))
            .next()
            .unwrap()
            .0;

        //let starting_node = match starting_node {
            //Some((i, _)) => i,
            //None => panic!("Graph Debug: {:?}", self),
        //};

        // the first relevant node is a number like: (1  0  0  1  1  0  1)
        // the next one would be e.g.:               (0  1  1  0  1  1  0) 
        // i.e. 1 on the diagonal.
        // The RELEVANT data, not including the diagonal one to represent its inclusion in the
        // graph, is then the node minus the 1 on the diagonal.
        // Later, we will calculate:
        //              node_data - (1<<node_index)
        let starting_node_data = self[starting_node];
        
        // now we have the relevant starting node, we can calculate the edge to drop
        // by finding the LAST bit that is set to 1. First, however, we need to
        // calculate the offset due to non-usize sized grpahs.
                                                                                
        //The first relevant node has some leading zeros up to its relevant diagonal,
        //a 1, and then a set of leading zeros up to the first edge.
        // Comparison point: number of zeros from start of adjacency until the graph information
        // starts.
        let comparison_point = starting_node_data.leading_zeros() as usize - starting_node;
        let graph_size = MAX_NODES.saturating_sub(comparison_point);
        // clean starting node data: the integer corresponding to the node with its first power of two
        // removed
        //  the edge to drop goes between the starting node and the end of the first edge
        let clean_starting_node_data = starting_node_data &!(1<<(graph_size - starting_node - 1));
        let end_node: usize;
        //if drop_first_connected_edge {
            //end_node = clean_starting_node_data.leading_zeros() as usize - comparison_point;
        //} else {
            // the edge to drop goes between the starting node and the end of its last edge
        let trailing_zeros = clean_starting_node_data.trailing_zeros() as usize + 1;
        end_node = graph_size.saturating_sub(trailing_zeros);
        //}
        (starting_node, end_node, graph_size)
    }
}

pub fn get_deck(graph: &GraphData) -> Vec<GraphData> {
    let mut deck: Vec<GraphData> = Vec::new();
    let graph_size = graph.graph_size();
    for i in 0..graph_size {
        let mut current_graph = graph.clone();
        current_graph.remove_node(i, graph_size);
        deck.push(current_graph.clone());
    }
    deck
}
