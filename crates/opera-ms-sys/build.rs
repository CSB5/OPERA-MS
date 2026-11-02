//! Compiles the C++ of `cpp/` (SIGMA, OPERA-LG, short_read_analysis) and the
//! cxx bridge into this crate.
//!
//! The C++ is always built with `-std=c++11 -O3 -ffp-contract=off`, whatever
//! the Cargo profile, so that debug and release builds run the same code.
//! `-ffp-contract=off` disables fused multiply-add: results then don't depend
//! on the compiler's target or on inlining. Each program's `main()` is
//! reused: `*_main.cpp` is compiled with `-Dmain=<name>_main` and called by
//! `cpp/bridge.cc`.

use std::path::{Path, PathBuf};

fn sources(dir: &Path, skip: &str) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = std::fs::read_dir(dir)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", dir.display()))
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|e| e == "cpp"))
        .filter(|p| p.file_name().is_some_and(|n| n != skip))
        .collect();
    files.sort();
    files
}

fn cpp_build() -> cc::Build {
    let mut build = cc::Build::new();
    build
        .cpp(true)
        .std("c++11")
        .opt_level(3)
        .debug(false)
        .warnings(false)
        .flag_if_supported("-w")
        .flag("-ffp-contract=off");
    build
}

fn main() {
    let cpp = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("cpp");
    let sigma = cpp.join("sigma");
    let opera = cpp.join("opera-lg");
    let sra = cpp.join("short_read_analysis");

    // The bridge first: it uses the libraries below, which must come after
    // it on the link line.
    cxx_build::bridge("src/cpp.rs")
        .file("cpp/bridge.cc")
        .include(&opera)
        .std("c++17")
        .opt_level(3)
        .flag("-ffp-contract=off")
        .compile("opera_ms_bridge");

    for (name, file, includes) in [
        (
            "sigma_main",
            sigma.join("sigma_main.cpp"),
            vec![sigma.clone()],
        ),
        (
            "opera_lg_main",
            opera.join("opera_main.cpp"),
            vec![opera.clone()],
        ),
        (
            "short_read_analysis_main",
            sra.join("short_read_analysis_main.cpp"),
            vec![sigma.clone(), opera.clone()],
        ),
    ] {
        let mut build = cpp_build();
        build.file(&file).define("main", name);
        for include in includes {
            build.include(include);
        }
        build.compile(&format!("opera_ms_{name}"));
        println!("cargo:rerun-if-changed={}", file.display());
    }

    let mut lib = cpp_build();
    lib.files(sources(&sigma, "sigma_main.cpp"))
        .files(sources(&opera, "opera_main.cpp"));
    lib.include(&sigma).include(&opera);
    lib.compile("opera_ms_cpp");

    for dir in [&sigma, &opera] {
        println!("cargo:rerun-if-changed={}", dir.display());
    }
    println!("cargo:rerun-if-changed=cpp/bridge.cc");
    println!("cargo:rerun-if-changed=cpp/bridge.h");
}
