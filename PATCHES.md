# Patches to NVIDIA's kernel source

This fork adds Rust crates in `rust/` and leaves NVIDIA's kernel alone wherever
it can. Any change to files under `source/` is one commit per fix, listed here
against the issue it addresses (a `KNOWN_ISSUES.md` entry or a bug we found), so
upstream merges stay easy to review.

Base: upstream `v1.0.0` (`c9c5979`).

| Commit | Files | Issue | Change |
|---|---|---|---|
| | | | |

None so far. `source/SMLib` and `source/SM_API` build unmodified with Apple
clang 21 on Apple Silicon through `rust/smlib-sys/build.rs`, and the `smlib`
tests pass against them.

## Notes from the bring-up (no patch needed)

- `SmTopologySweep.cpp:1773` triggers `-Wmisleading-indentation`: an `if`
  around `SetSame(pE, pSweptEdge)` was commented out and the statement kept its
  indentation. Harmless, and the only compiler warning in the kernel build.
- `SmApiCreateCylinder` builds an `SmCone` with a zero semi-angle, not an
  `SmCylinder`; the shim reports it as a cylinder via `SmCone::IsCylinder()`.
- `SmApiCreateContext` is the only setup call needed; neither TBB nor USD is
  linked on macOS.

## Notes from the fillet work (no patch needed)

Fillets and chamfers needed no kernel changes either. What the shim works
around:

- `SmApiChamferFillet` breaks chamfers on circular edges. On the top rim of a
  10 mm hole through a 100 × 60 × 40 box, widths of about 3.5 to 5.7 (legs 2.5
  to 4) return `SM_SUCCESS` but leave an open shell: the chamfer face is
  added, the neighbouring faces aren't trimmed to it, and one rail comes back
  split into dozens of unshared edges. A 20 mm hole fails at leg 3 too. The
  cause is the executive set-up, not the constant-distance solver: it calls
  `SetDoGlobalMerge(TRUE)` while topology insertion stays on by default, and
  `SmFilletExecutive.h` says only one merge mode should be on. A constant-
  radius solver set up the same way fails the same way; either solver with
  global merge off chamfers every size correctly, and box chamfers with
  beveled corners come out identical. So the shim's `sm_brep_chamfer_edges`
  copies `SmApiChamferFillet` except that call, and is 22 to 56% faster for it
  in the three cases timed. Upstream's `SmApiChamferFillet` itself is left
  unpatched, since nothing here calls it.

  Global merge is the executive's fallback for large cuts, where a fillet
  removes whole pieces of a face. Two such chamfers tried (cutting past a 2 mm
  boss, overlapping chamfers on a 4 mm plate) fail cleanly with the brep
  untouched whether it's on or off. A corner edge that misses a side face
  errors out with it off rather than being left to the merge, and
  `SmFilletNx1Corner` turns it on itself for corners that need it.
- Fillets can also report success on a broken result: radius 25 on every
  edge of a 40 mm-high box returns `SM_SUCCESS` and a non-manifold brep. The
  shim checks manifoldness before and after each fillet or chamfer and
  returns `SM_SHIM_ERR_BROKE_SOLID` instead.
- `SmApiFilletEdgesPerEdge` calls `SetDoGlobalMerge(TRUE)` before
  `SetFilletParameters`, which resets it to `FALSE`, so global merge is off.
  Circular sections take the per-edge branch, which uses a 0.1 arc
  tolerance instead of the single-type branch's `dXSectTol`. Neither showed
  up in the results: box fillets come out as exact cylinders and sphere
  corners, with volumes within 1e-7 of analytic.

## Notes from the labels work (no patch needed)

Labels on faces and edges ride on SMLib's attribute system, and it delivered
every notification the tests look for: attributes survive booleans (including
faces split into pieces and faces brought over from the second operand),
fillets, chamfers, brep copies, and the face merges after a union. No kernel
change was needed.

- SMLib's default merge rule keeps the surviving face's attribute, and which
  face survives depends on operand order: `a ∪ b` and `b ∪ a` of two flush
  boxes keep different labels on the merged side faces. The shim's
  `LabelAttribute` overrides `SmAttribute::Merge` so the older label (the
  lower feature id) wins either way.
- Faces a fillet or chamfer creates, and edges a boolean or fillet creates,
  arrive with no attribute. That's by design, not a missing notification:
  they come from nothing. Callers label them after the operation.
- `SmPrimitiveCreation::CreateCurveSweepFromFaces` sweeps copies of the
  face's edge curves, so its side faces have no link back to the edges they
  came from. The shim's `sm_brep_sweep_face` finds each side's source edge by
  where its start edge lies.
- `SmApiFaceComputeCentroid` clamps its relative accuracy to 1e-4, which puts
  a hole placed at a face's centroid up to a few hundredths of a millimetre
  off on a 500 mm face. The shim calls `SmFace::ComputePreciseProperties`
  directly at 1e-7.
- `SmSurface::IsPlanar`'s normal comes from a bounding pseudo-box and can face
  either way, so the shim evaluates a planar face's normal from its surface
  and orients it out of the solid by region, as `outward_normal` does.
