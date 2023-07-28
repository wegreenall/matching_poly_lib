mod binary_representation;
mod matching_poly;

pub use binary_representation::BinaryGraph;
pub use self::matching_poly::{ calculate_matching_polynomial_pointer,
    calculate_matching_polynomial_pointer_addresses,
    calculate_matching_polynomial_adaptive,
    _calculate_matching_polynomial_binary};
