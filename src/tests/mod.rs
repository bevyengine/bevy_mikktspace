//! Tests that should be completed in-crate.
//! This allows using `#[cfg(test)]` within the crate to add additional checks.
//! Tests in a `tests/` folder at the repository root are treated as external
//! to the crate, and thus the crate is _not_ compiled with `test`.

mod obj;
pub(crate) mod shadow;
