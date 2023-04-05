mod matchings;
mod binary_representation;
mod petgraph;
mod base;

pub use self::binary_representation::{_calculate_matching_polynomial_binary, get_deck, Graph};
pub use self::petgraph::{_calculate_matching_polynomial, get_matching_polies_stable_graph};
pub use self::base::calculate_matching_polynomial_static;
