use core::fmt;
use std::mem::size_of;
use polynomial::Polynomial;
use std::mem;
use std::cmp::PartialEq;
//use crate::traits::Graph;

const MAX_NODES: usize = mem::size_of::<usize>()*8;

/// We represent graphs as a seequence of integers in which each bit represents
/// an edge. The first bit however represents whether that node is contained
/// in the graph; removing a node implies zero-ing this bit.
#[derive(Debug, Clone, Copy)]
pub struct Graph {
    data: [usize; size_of::<usize>()*8],
    initial_graph_size: usize,
}


impl Graph {
   pub const fn new() -> Graph {
        let blank_data = [0; size_of::<usize>()*8];
         Graph {
               data: blank_data,
               initial_graph_size: 0,
         }
   }

   pub fn from(data: [usize; size_of::<usize>()*8]) -> Graph {
        Graph {
            data,
            initial_graph_size: data
                .iter()
                .filter(|x| x> &&(0 as usize))
                .count(),
        }
    }

    pub fn data(self) -> [usize; size_of::<usize>()*8] {
        self.data
    }

    pub fn remove_node(&mut self, node: usize, graph_size : usize) {
        // remove node from adjacency list
        self.data[node] = 0;

        // Now remove its connected edges
        self.data
            .iter_mut()
            .for_each(|x| *x &= !(1<<graph_size.saturating_sub(node+1)));
    }

    pub fn remove_edge(&mut self, node1: usize, node2: usize, graph_size: usize) {
        // remove edge from adjacency list
        let shift = graph_size.saturating_sub(node2);
         
        // zero the ith from the right
        self.data[node1] &= (!(1 << (shift-1))) as usize; 
    }

    pub fn edgeless_node_count(&self) -> usize {
        // count the number of nodes in the graph
        // To count the number of nodes in th graph,
        // count the number of the bits in the graph
        // Since, by assumption, the graph is edgeless,
        // we can just count the number of 1s in the graph
        self.data
            .iter()
            .sum::<usize>()
            .count_ones() as usize
    }

    pub fn graph_size(&self) -> usize{
        self.data
            .iter()
            .filter(|x| x> &&(0 as usize)) // i.e. get the ones that are valid
            .count()
    }
    pub fn initial_graph_size(&self) -> usize{
        self.initial_graph_size   }

    /// checks whether the graph is edgeless, i.e. if each of the elements
    /// is a power of two or 0
    pub fn edgeless(&self) -> bool {
        self.data
            .iter()
            .all(|x| x == &(0 as usize) || x.is_power_of_two())
    } 

    pub fn get_graph_primes(self) -> (Graph, Graph) {
        let mut new_graph = self; // Should copy
        let mut new_graph2 = self; // Should copy

        // get the relevant edge
        let (start_node, end_node, graph_size) = self.get_relevant_edge();

        // G' = G - e
        new_graph.remove_edge(start_node, end_node, graph_size);
        
        // G'' = G - {v, w} where {w, v} are the nodes connected to e
        new_graph2.remove_node(start_node, graph_size);
        new_graph2.remove_node(end_node, graph_size);
        (new_graph, new_graph2)
    }

    /// To step through and calculate the matching polynomial, we use the edge
    /// remove recurrence:
    /// m(g, x) = m(G', x) - m(G'', x)
    ///
    /// G' = G - e
    /// G'' = G - {v, w} where {w, v} are the nodes
    /// at the ends of e.
    ///
    /// Thus, we get get the "relevant edge e" which is the first edge in the
    /// first remaining node. Since the nodes are ordered in decreasing order 
    /// of degree, dropping the first edge we find drops the most edges from 
    /// the graph, since the nodes at its ends will be the "most connected" 
    /// nodes.
    pub fn get_relevant_edge(&self) -> (usize, usize, usize) {
        // since the nodes are ordered in INCREASING order of degree, we can
        // just drop the last(right-most in the binary representation)
        // edge we find, on the first still-relevant node.
        // starting_node: the index of the first node that still has edges from it


        let starting_node = self.data
            .iter()
            .enumerate()
            .filter(|(_, x)| (x > &&(0 as usize)))
            .filter(|(_, x)| !(x.is_power_of_two()))
            .next()
            .unwrap()
            .0;

        // the first relevant node is a number like: (1  0  0  1  1  0  1)
        // the next one would be e.g.:               (0  1  1  0  1  1  0) 
        // i.e. 1 on the diagonal.
        // The RELEVANT data, not including the diagonal one to represent its inclusion in the
        // graph, is then the node minus the 1 on the diagonal.
        // Later, we will calculate:
        //              node_data - (1<<node_index)
        let starting_node_data = self.data[starting_node];
        
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

        // if drop_most_connected_edge is true, we delete the first connected edge
        // we find on the first relevant node. Otherwise, we delete the "last"
        // edge for the first relevant node.
        let drop_first_connected_edge: bool = false;
        if drop_first_connected_edge {
            end_node = clean_starting_node_data.leading_zeros() as usize - comparison_point;
        } else {
            // the edge to drop goes between the starting node and the end of its last edge
            let trailing_zeros = clean_starting_node_data.trailing_zeros() as usize + 1;
            end_node = graph_size.saturating_sub(trailing_zeros);
        }
        let print_stuff: bool = false;
        if print_stuff {
            println!("\n");
            println!("starting_node {}", starting_node);
            println!("end_node {}", end_node);
            println!("starting_node_data {:b}", starting_node_data);
            println!("MAX NODES: {}", MAX_NODES);
            println!("graph_size: {} ", graph_size);
            println!("comparison_point {}", comparison_point);
            println!("clean_starting_node_data {:b}", clean_starting_node_data);
            println!("clean_starting_node_data leading zeros: {}", clean_starting_node_data.leading_zeros());
        }
        (starting_node, end_node, graph_size)
    }
}


/// Equality of the graphs will have be done based on if they are isomorphic to 
/// each other. This problem is, in general, NP-complete. However, we can
/// use the fact that the graphs are undirected and unweighted to simplify the
/// problem. We can then use the fact that the graphs are ordered in decreasing
/// order of degree to simplify the problem further.
/// 
impl PartialEq for Graph {
    fn eq(&self, other: &Graph) -> bool {
        self.data == other.data
    }
}

impl std::fmt::Display for Graph {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut result = write!(f, "\n");
        for x in self.data.iter() {
            result = write!(f, "{}\n", x);
        }
        result
    }
}

pub fn get_deck(graph: &Graph) -> Vec<Graph>{
    let mut deck = Vec::<Graph>::new();
    let graph_size = graph.graph_size();
    for i in 0..graph_size {
        //println!("current graph: {}", current_graph);
        let mut current_graph = graph.clone();
        current_graph.remove_node(i, graph_size); 
        deck.push(current_graph.clone());
    }
    deck
}


