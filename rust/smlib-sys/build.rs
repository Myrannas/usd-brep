//! Builds NVIDIA's SMLib kernel (`source/SMLib`, `source/SM_API`) and our C shim
//! with the `cc` crate, caching both static libraries outside `target/`.
//!
//! The cache lives in `$SMLIB_CACHE_DIR`, by default `~/Library/Caches/smlib` on
//! macOS and `$XDG_CACHE_HOME/smlib` (or `~/.cache/smlib`) elsewhere. It holds two
//! kinds of entry, each a directory named by a content hash:
//!
//! - `kernel-<hash>/` with `libsmlib_kernel.a` (SMLib) and `libsmlib_api.a`
//!   (SM_API): keyed by every kernel source and header, the compiler and its
//!   version, the full compiler command line, target and opt level.
//! - `shim-<hash>/libsmlib_shim.a`: keyed by the kernel key plus the shim sources,
//!   so editing the shim recompiles one file rather than the whole kernel.
//!
//! On a hit the script only emits link flags. On a miss it compiles into a temp
//! directory beside the entry and renames it into place while holding
//! `<entry>.lock`, so parallel builds from several worktrees share one result.

use std::env;
use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use sha2::{Digest, Sha256};

/// Bump to invalidate every cache entry when the build procedure changes in a
/// way the hashed inputs don't capture.
const CACHE_SCHEMA: &str = "smlib-sys cache v1";

const KERNEL_LIB: &str = "smlib_kernel";
const API_LIB: &str = "smlib_api";
const SHIM_LIB: &str = "smlib_shim";

/// Upstream's macOS flags (`macos_common_buildoptions` in premake5-utils.lua),
/// minus `-arch`, which cc derives from the Cargo target. `-std=c++17` follows
/// `-ansi`, as upstream, so it wins.
const CXX_FLAGS: &[&str] = &[
    "-ansi",
    "-Wall",
    "-Wno-missing-braces",
    "-std=c++17",
    "-Wno-inconsistent-missing-override",
    "-Wno-overloaded-virtual",
    "-Wno-unused-but-set-variable",
    "-Wno-unused-variable",
    "-Wno-parentheses-equality",
];

/// Upstream's macOS release defines (premake5.lua), plus `SM_API_EXPORTS` from
/// SM_API's project.
const DEFINES: &[(&str, Option<&str>)] = &[
    ("NDEBUG", None),
    ("NOMINMAX", None),
    ("SM_GRAPHICS_CALLBACKS", Some("1")),
    ("SM_API_EXPORTS", None),
];

/// Environment variables that change the compiler command or the cache location.
const TRACKED_ENV: &[&str] = &[
    "SMLIB_CACHE_DIR",
    "CXX",
    "CXXFLAGS",
    "CRATE_CC_NO_DEFAULTS",
    "MACOSX_DEPLOYMENT_TARGET",
];

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let repo_root = manifest_dir.join("../..").canonicalize().unwrap();
    let source = repo_root.join("source");
    let shim_dir = manifest_dir.join("shim");
    let bindings = manifest_dir.join("src/bindings.rs");
    let target = env::var("TARGET").unwrap();

    for dir in ["SMLib/src", "SMLib/inc", "SM_API/src", "SM_API/inc"] {
        println!("cargo:rerun-if-changed={}", source.join(dir).display());
    }
    println!("cargo:rerun-if-changed={}", shim_dir.display());
    println!("cargo:rerun-if-changed={}", bindings.display());
    for var in TRACKED_ENV {
        println!("cargo:rerun-if-env-changed={var}");
    }
    for var in ["CXX", "CXXFLAGS"] {
        println!("cargo:rerun-if-env-changed={var}_{target}");
        println!(
            "cargo:rerun-if-env-changed={var}_{}",
            target.replace('-', "_")
        );
    }

    let shim_header = shim_dir.join("shim.h");
    sync_bindings(&shim_header, &bindings);

    let kernel_inputs = {
        let mut files = Vec::new();
        for dir in ["SMLib/inc", "SM_API/inc"] {
            files.extend(list_files(&source.join(dir), None));
        }
        files.extend(list_files(&source.join("SMLib/src"), Some("cpp")));
        files.extend(list_files(&source.join("SM_API/src"), Some("cpp")));
        files
    };
    let smlib_sources = list_files(&source.join("SMLib/src"), Some("cpp"));
    let api_sources = list_files(&source.join("SM_API/src"), Some("cpp"));
    let shim_inputs = list_files(&shim_dir, None);
    let shim_sources: Vec<PathBuf> = shim_inputs
        .iter()
        .filter(|p| p.extension().is_some_and(|e| e == "cpp"))
        .cloned()
        .collect();

    let base = base_build();
    let compiler = base.try_get_compiler().expect("no C++ compiler found");

    let hash_start = Instant::now();
    let mut kernel_key = Sha256::new();
    kernel_key.update(CACHE_SCHEMA);
    hash_toolchain(&mut kernel_key, &compiler, &target, &repo_root);
    hash_files(&mut kernel_key, &repo_root, &kernel_inputs);
    let kernel_hash = hex(&kernel_key.finalize());

    let mut shim_key = Sha256::new();
    shim_key.update(&kernel_hash);
    hash_files(&mut shim_key, &repo_root, &shim_inputs);
    let shim_hash = hex(&shim_key.finalize());
    let hash_ms = hash_start.elapsed().as_millis();

    // Include order matters: SMLib/inc and SM_API/inc each have a StdAfx.h, and
    // upstream builds each project with its own directory first.
    let smlib_inc = source.join("SMLib/inc");
    let api_inc = source.join("SM_API/inc");
    let root = cache_root();
    let kernel_entry = ensure_entry(&root, &format!("kernel-{}", &kernel_hash[..32]), |out| {
        let mut smlib = base.clone();
        smlib.include(&smlib_inc);
        let mut api = base.clone();
        api.include(&api_inc).include(&smlib_inc);
        [
            compile(&smlib, out, &smlib_sources, KERNEL_LIB),
            compile(&api, out, &api_sources, API_LIB),
        ]
        .join("; ")
    });
    let shim_entry = ensure_entry(&root, &format!("shim-{}", &shim_hash[..32]), |out| {
        let mut shim = base.clone();
        shim.include(&shim_dir)
            .include(&api_inc)
            .include(&smlib_inc);
        compile(&shim, out, &shim_sources, SHIM_LIB)
    });

    println!("smlib-sys: hashed kernel and shim sources in {hash_ms} ms");
    for entry in [&shim_entry, &kernel_entry] {
        match &entry.built {
            Some(report) => println!(
                "cargo:warning=smlib-sys: built {}: {report}",
                entry.dir.display()
            ),
            None => println!("smlib-sys: cache hit at {}", entry.dir.display()),
        }
        println!("cargo:rustc-link-search=native={}", entry.dir.display());
    }
    // Static archives are searched in order: the shim calls SM_API, which
    // calls SMLib.
    for lib in [SHIM_LIB, API_LIB, KERNEL_LIB] {
        println!("cargo:rustc-link-lib=static={lib}");
    }
    if target.contains("apple") {
        println!("cargo:rustc-link-lib=c++");
    } else if !target.contains("msvc") {
        println!("cargo:rustc-link-lib=stdc++");
    }
}

/// The compiler configuration shared by SMLib, SM_API and the shim, minus
/// include directories.
fn base_build() -> cc::Build {
    let mut build = cc::Build::new();
    build
        .cpp(true)
        .cargo_metadata(false)
        .cargo_warnings(false)
        .emit_rerun_if_env_changed(false)
        // Debug info would split the cache between dev and release profiles
        // and triple the archive size; the kernel is always built optimised.
        .debug(false)
        // cc adds -Wall -Wextra; keep -Wall and upstream's suppressions, which
        // come after it. (`warnings(false)` would pass -w and hide everything.)
        .extra_warnings(false);
    for (name, value) in DEFINES {
        build.define(name, *value);
    }
    for flag in CXX_FLAGS {
        build.flag(flag);
    }
    build
}

/// Hashes everything about the toolchain that changes the objects: the
/// compiler binary and version, its full argument list (with this checkout's
/// path replaced, so worktrees share entries), the target and the opt level.
fn hash_toolchain(key: &mut Sha256, compiler: &cc::Tool, target: &str, repo_root: &Path) {
    let root = repo_root.to_string_lossy();
    let mut text = String::new();
    writeln!(text, "compiler {}", compiler.path().display()).unwrap();
    writeln!(
        text,
        "version {}",
        command_output(compiler.path(), &["--version"])
    )
    .unwrap();
    for arg in compiler.args() {
        writeln!(
            text,
            "arg {}",
            arg.to_string_lossy().replace(&*root, "<root>")
        )
        .unwrap();
    }
    for (name, value) in compiler.get_envs() {
        writeln!(
            text,
            "env {}={}",
            name.to_string_lossy(),
            value.to_string_lossy()
        )
        .unwrap();
    }
    writeln!(text, "target {target}").unwrap();
    writeln!(
        text,
        "opt-level {}",
        env::var("OPT_LEVEL").unwrap_or_default()
    )
    .unwrap();
    if target.contains("apple") {
        let sdk = command_output(
            Path::new("xcrun"),
            &["--sdk", "macosx", "--show-sdk-version"],
        );
        writeln!(text, "sdk {sdk}").unwrap();
    }
    key.update(text.as_bytes());
}

fn hash_files(key: &mut Sha256, repo_root: &Path, files: &[PathBuf]) {
    for file in files {
        let rel = file.strip_prefix(repo_root).unwrap_or(file);
        let bytes = fs::read(file).unwrap_or_else(|e| panic!("reading {}: {e}", file.display()));
        key.update(rel.to_string_lossy().as_bytes());
        key.update((bytes.len() as u64).to_le_bytes());
        key.update(&bytes);
    }
}

/// Sorted regular files in `dir`, optionally only those with `extension`.
fn list_files(dir: &Path, extension: Option<&str>) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = fs::read_dir(dir)
        .unwrap_or_else(|e| panic!("reading {}: {e}", dir.display()))
        .map(|entry| entry.unwrap().path())
        .filter(|p| p.is_file())
        .filter(|p| extension.is_none_or(|ext| p.extension().is_some_and(|e| e == ext)))
        .collect();
    files.sort();
    files
}

fn command_output(program: &Path, args: &[&str]) -> String {
    match Command::new(program).args(args).output() {
        Ok(out) => String::from_utf8_lossy(&out.stdout).trim().to_owned(),
        Err(_) => String::new(),
    }
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().fold(String::new(), |mut s, b| {
        write!(s, "{b:02x}").unwrap();
        s
    })
}

fn cache_root() -> PathBuf {
    if let Some(dir) = env::var_os("SMLIB_CACHE_DIR").filter(|d| !d.is_empty()) {
        return PathBuf::from(dir);
    }
    let home = env::var_os("HOME").map(PathBuf::from);
    if cfg!(target_os = "macos") {
        if let Some(home) = home {
            return home.join("Library/Caches/smlib");
        }
    } else if let Some(xdg) = env::var_os("XDG_CACHE_HOME").filter(|d| !d.is_empty()) {
        return PathBuf::from(xdg).join("smlib");
    } else if let Some(home) = home {
        return home.join(".cache/smlib");
    }
    PathBuf::from(env::var("OUT_DIR").unwrap()).join("smlib-cache")
}

struct Entry {
    dir: PathBuf,
    /// What the build did, when this run built the entry rather than hitting.
    built: Option<String>,
}

/// Returns the cache entry `root/name`, building it with `build` first if it
/// doesn't exist. Builds happen in a temp directory under an exclusive lock and
/// are renamed into place, so a present entry is always complete.
fn ensure_entry(root: &Path, name: &str, build: impl FnOnce(&Path) -> String) -> Entry {
    let dir = root.join(name);
    if dir.is_dir() {
        return Entry { dir, built: None };
    }
    fs::create_dir_all(root).unwrap_or_else(|e| panic!("creating {}: {e}", root.display()));
    let lock_path = root.join(format!("{name}.lock"));
    let lock = File::create(&lock_path)
        .unwrap_or_else(|e| panic!("creating {}: {e}", lock_path.display()));
    lock.lock()
        .unwrap_or_else(|e| panic!("locking {}: {e}", lock_path.display()));
    if dir.is_dir() {
        // Another build finished it while we waited for the lock.
        return Entry { dir, built: None };
    }
    // Holding the lock, any temp dir for this entry is a failed earlier build.
    let tmp_prefix = format!("{name}.tmp-");
    for stale in fs::read_dir(root).into_iter().flatten().flatten() {
        if stale.file_name().to_string_lossy().starts_with(&tmp_prefix) {
            let _ = fs::remove_dir_all(stale.path());
        }
    }
    let tmp = root.join(format!("{tmp_prefix}{}", std::process::id()));
    fs::create_dir_all(&tmp).unwrap_or_else(|e| panic!("creating {}: {e}", tmp.display()));
    let report = build(&tmp);
    fs::rename(&tmp, &dir).unwrap_or_else(|e| panic!("renaming {} into place: {e}", tmp.display()));
    drop(lock);
    let _ = fs::remove_file(&lock_path);
    Entry {
        dir,
        built: Some(report),
    }
}

/// Collects every line of compiler stderr, to count warnings and keep a log.
#[derive(Default)]
struct Log(Mutex<Vec<String>>);

impl cc::BuildMessageLogger for Log {
    fn log(
        &self,
        kind: cc::BuildMessageKind,
        msg: cc::BuildMessage<'_>,
        _extra: &dyn std::any::Any,
    ) {
        let line = match kind {
            cc::BuildMessageKind::StderrForwarding => msg.to_string(),
            _ => format!("[cc] {msg}"),
        };
        self.0.lock().unwrap().push(line);
    }
}

/// Compiles `sources` into `out/lib<lib>.a`, writes `out/build.log` with the
/// compiler command, timing and every diagnostic, and returns a summary.
fn compile(base: &cc::Build, out: &Path, sources: &[PathBuf], lib: &str) -> String {
    let log = Arc::new(Log::default());
    let compiler = base.try_get_compiler().expect("no C++ compiler found");
    let mut build = base.clone();
    let obj = out.join(format!("obj-{lib}"));
    build
        .out_dir(&obj)
        .files(sources)
        .message_logger(Some(log.clone()));

    let start = Instant::now();
    let result = build.try_compile(lib);
    let secs = start.elapsed().as_secs_f64();

    let lines = std::mem::take(&mut *log.0.lock().unwrap());
    let warnings = lines.iter().filter(|l| l.contains(": warning:")).count();
    let summary = format!(
        "lib{lib}.a from {} files in {secs:.1} s, {warnings} compiler warnings",
        sources.len()
    );

    let mut file = File::options()
        .create(true)
        .append(true)
        .open(out.join("build.log"))
        .unwrap();
    writeln!(file, "{summary}").unwrap();
    // The compiler and its arguments only: the full Command would also record
    // this process's environment, which can hold credentials.
    write!(file, "compiler: {}", compiler.path().display()).unwrap();
    for arg in compiler.args() {
        write!(file, " {}", arg.to_string_lossy()).unwrap();
    }
    writeln!(file).unwrap();
    for line in &lines {
        writeln!(file, "{line}").unwrap();
    }
    drop(file);

    if let Err(e) = result {
        for line in lines.iter().filter(|l| l.contains("error")).take(40) {
            println!("cargo:warning={line}");
        }
        panic!(
            "compiling {lib} failed after {secs:.1} s: {e}\nfull log: {}",
            out.join("build.log").display()
        );
    }
    // The archive is all we need; objects would double the entry's size.
    let archive = format!("lib{lib}.a");
    fs::rename(obj.join(&archive), out.join(&archive))
        .expect("moving the archive out of its obj dir");
    let _ = fs::remove_dir_all(&obj);
    summary
}

/// Keeps `src/bindings.rs` in step with `shim.h`. The bindings carry the hash
/// of the header they were generated from; normal builds only check it, so
/// they don't need libclang. Building with `--features regenerate-bindings`
/// reruns bindgen.
fn sync_bindings(header: &Path, bindings: &Path) {
    let header_hash = hex(&Sha256::digest(fs::read(header).expect("reading shim.h")));
    let stamp = format!("// shim.h sha256: {header_hash}");

    #[cfg(feature = "regenerate-bindings")]
    {
        let generated = bindgen::Builder::default()
            .header(header.to_string_lossy())
            .allowlist_function("sm_.*")
            .allowlist_type("sm_.*")
            .allowlist_var("SM_.*")
            .default_enum_style(bindgen::EnumVariation::Consts)
            .prepend_enum_name(false)
            .layout_tests(true)
            .use_core()
            .generate()
            .expect("bindgen failed on shim.h");
        let text = format!(
            "// @generated by bindgen from shim/shim.h. Do not edit; rebuild with\n\
             // `--features regenerate-bindings` after changing shim.h.\n{stamp}\n\n{generated}"
        );
        fs::write(bindings, text).expect("writing src/bindings.rs");
    }

    let current = fs::read_to_string(bindings).unwrap_or_default();
    if !current.lines().take(4).any(|l| l == stamp) {
        panic!(
            "src/bindings.rs is out of date with shim/shim.h. Regenerate it with \
             `cargo build -p smlib-sys --features regenerate-bindings` (needs libclang)."
        );
    }
}
