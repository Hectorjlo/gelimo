//! Helper macros and floating point comparison utilities for integration tests.

#[allow(dead_code)]
pub fn assert_close(a: f64, b: f64, tol: f64) {
    let diff = (a - b).abs();
    assert!(
        diff <= tol,
        "Values differ by {}, expected tolerance {}",
        diff,
        tol
    );
}
