# opera-ms: one binary, the Rust pipeline with its C++ compiled in. The
# external tools come from PATH (conda environment: environment.yml).

install: opera-ms

opera-ms:
	cargo build --release --locked
	cp target/release/opera-ms ./opera-ms

test:
	cargo fmt --check && cargo clippy --all-targets --locked -- -D warnings && cargo test --locked

clean:
	cargo clean
	rm -f ./opera-ms

.PHONY: install opera-ms test clean
