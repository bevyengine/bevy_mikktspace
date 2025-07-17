# Fuzzing

## Setup

See the [Rust Fuzz Book](https://rust-fuzz.github.io/book/introduction.html) for details on getting started.
In short:

- Install nightly with `rustup install nightly`
- Install `cargo-fuzz` with `cargo install cargo-fuzz`
- Run one of the fuzzers with `cargo +nightly fuzz run <fuzz target>` where `<fuzz target>` is one of the binaries as specified in this [Cargo.toml](./Cargo.toml).
