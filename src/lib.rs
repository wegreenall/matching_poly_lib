pub mod matching;

use matching::{get_matching_polies_stable_graph, get_deck, Graph, _calculate_matching_polynomial_binary};
use matching::{calculate_matching_polynomial_static};

pub fn rs_test_function() {
    println!("Printing from the test function in the library");
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
