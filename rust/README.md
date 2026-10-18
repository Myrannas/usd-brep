# Rust crates for USD BRep's kernel

- `smlib-sys`: builds SMLib (`source/SMLib`) and SM_API (`source/SM_API`) with
  the `cc` crate, plus a C shim (`smlib-sys/shim`), and exposes the shim through
  committed bindgen bindings.
- `smlib`: the safe API. `Kernel`, `Brep` (primitives, booleans, fillets,
  chamfers, volume, manifold check) and `Mesh`, a tessellation that keeps its
  links to the B-rep: per-triangle face ids, per-vertex surface (u, v),
  per-face surface parameters and per-edge polylines.

## Fillets and chamfers

Edges are picked by geometric query, never by index, since indices change
with every operation. An `EdgeQuery` tests each edge's `EdgeGeometry`
(sample points, length, and the adjacent faces' outward normals at its
midpoint): `MidpointOnPlane`, `InPlane`, `OnCircle`, `BetweenFaces` (by
normals), `Sharp` (by the angle between normals), and `And` / `Or` / `Not`.

- `Brep::fillet` wraps `SmApiFilletEdgesPerEdge`: each `EdgeBlend` gives its
  edges a radius and a circular, linear or blend cross-section;
  `FilletOptions` sets G1/G2/G3 continuity for blends. A linear section is a
  chamfer set back by the radius on each face of a 90° edge.
- `Brep::chamfer` is `SmApiChamferFillet` restricted to the picked edges: a
  constant chamfer width with beveled corners. It leaves the kernel's global
  merge off, which `SmApiChamferFillet` turns on and which breaks chamfers on
  circles (see `PATCHES.md`).
- Both consume the brep and fail with `SHIM_ERR_BROKE_SOLID` when a manifold
  solid comes back non-manifold, which the kernel can do while reporting
  success.
- There is no undo: `SmApiRemoveFillet` is a no-op (`KNOWN_ISSUES.md`), so
  rebuild from the unfilleted solid. The default 15° surface angle avoids
  the inside-out triangles next to fillets listed there; the tests check
  every triangle's winding.

```sh
cargo test -p smlib
```

## Build cache

The kernel takes about 30 s to compile on a 10-core M-series Mac. `build.rs`
compiles it once per distinct input and caches the static libraries outside
`target/`, in `$SMLIB_CACHE_DIR` (default `~/Library/Caches/smlib` on macOS,
`$XDG_CACHE_HOME/smlib` or `~/.cache/smlib` elsewhere):

- `kernel-<hash>/`: `libsmlib_kernel.a` and `libsmlib_api.a`, keyed by every
  kernel source and header, the compiler and its version, the full compiler
  command line, target and opt level.
- `shim-<hash>/`: `libsmlib_shim.a`, keyed by the kernel key plus the shim, so a
  shim change recompiles one file.

A `cargo clean`, a new worktree or another branch with the same kernel and shim
hits the cache and only links. Entries are built in a temp directory and renamed
into place under a lock file, so parallel builds share one compile. Each entry's
`build.log` has the compile time, the compiler command and every diagnostic.
Entries are never evicted; delete old ones by hand.

Consumers should build `smlib-sys` optimised in every profile, since opt level
is part of the key and the kernel is very slow at `-O0`:

```toml
[profile.dev.package.smlib-sys]
opt-level = 3
```

## Regenerating the bindings

`smlib-sys/src/bindings.rs` is generated from `shim/shim.h` and committed, so
builds don't need libclang. Its header records the hash of the `shim.h` it came
from, and the build fails if they disagree. After changing `shim.h`:

```sh
cargo build -p smlib-sys --features regenerate-bindings
```

## Threading

SMLib keeps process-global state and is not thread-safe. `Kernel::init` binds
the kernel to the first thread that calls it; `Kernel` and `Brep` are neither
`Send` nor `Sync`, so every kernel call happens on that thread. `Mesh` is plain
data and can be sent anywhere.
