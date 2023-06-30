mod graphs;
mod matching_poly;

pub use self::graphs::{_calculate_matching_polynomial_binary, get_deck, Graph};
pub use self::matching_poly::calculate_matching_polynomial_pointer;
