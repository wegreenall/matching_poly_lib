mod example_petgraphs;
mod graphs;
mod petgraph;
mod matching_poly;

pub use self::graphs::{_calculate_matching_polynomial_binary,_calculate_weighted_matching_polynomial_binary, get_deck, Graph};
pub use self::petgraph::{_calculate_matching_polynomial, get_matching_polies_stable_graph};
pub use self::matching_poly::calculate_matching_polynomial_pointer;
