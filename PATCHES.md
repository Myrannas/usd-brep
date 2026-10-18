# Patches to NVIDIA's kernel source

This fork adds Rust crates in `rust/` and leaves NVIDIA's kernel alone wherever
it can. Any change to files under `source/` is one commit per fix, listed here
against the issue it addresses (a `KNOWN_ISSUES.md` entry or a bug we found), so
upstream merges stay easy to review.

Base: upstream `v1.0.0` (`c9c5979`).

| Commit | Files | Issue | Change |
|---|---|---|---|
| `a81ed1b` | `SmOffsetGeometryCreation.cpp` | Found offsetting faces: in a local operation, `ExtendTrimEdges` lifts a tangent edge from whichever side pointer order picks, which for a face that kept its surface is the old edge; the solid comes out open. | For `SM_OO_LOCAL_OPERATION`, lift from the side whose lifted midpoint lies on the other side's new surface, measured to the whole surface (exact for analytic surfaces and planar B-splines), not the G1-extended one. |
| `fdc9c9b` | `SmOffsetGeometryCreation.cpp` | Found offsetting faces: when the curve replacing an open edge is closed (a plane meeting a sphere), `MateEdges` keeps the trim interval inside the curve's parameter range, which can be the far arc: an inside-out solid with no gaps. | Keep the arc nearer the original edge, re-anchoring the 3D and UV curves when it runs past the curve's seam, and only if they still agree. |
| `c282328` | `SmOffsetGeometryCreation.cpp` | Found offsetting faces: `TrimExtendedEdge` trims against a not-yet-mated tangent neighbour's old UV curve, which in a local operation can describe the old corner. | For `SM_OO_LOCAL_OPERATION`, also offer the base vertex as a trim point when it lies on the new curve; the shim moves vertices it can place there for the operation. |
| `4d13c5e` | `SmTess.h`, `SmTess.cpp` | Found keeping tessellated faces across rebuilds: `Phase1SetupBrep` builds every face's surface cache and splits every edge, so a small edit still pays for the whole solid's edge preparation. | Add `SmTess::Phase1SetupFaces`, which sets up only the faces given (with caches for the faces sharing their edges, which split only those edges) exactly as `Phase1SetupBrep` would, and reports refinement or anything `Phase1SetupBrep` would recover from; and a `SplitFaceEdgesToUVSize` overload that splits only edges shared with given faces. Existing calls are unchanged. |

`source/SM_API` is unmodified. Each offset patch changes only
`SM_OO_LOCAL_OPERATION` (`SmBrep::LocalOperation` and the shim's face
offset); solid offsets and shells take the code paths they always did. The
tessellation patch only adds methods: `SmApiTessellate` and
`Phase1SetupBrep` take the code paths they always did.

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

## Notes from the face-offset work

`Brep::offset_faces` moves any set of faces along their outward normals, the
way Parasolid's face offset does. It is `SmBrep::LocalOperation`: every face
gets a surface and `SmOffsetExecutive::DoExtendedOffset` with
`SM_OO_LOCAL_OPERATION` and extend-convex-edges re-mates every edge and merges
the result. Used as is, it moves a plain or chamfered box's top face exactly
but not much else. What goes wrong, and what the shim
(`rust/smlib-sys/shim/shim.cpp`) does about it:

- A fillet on a moved face must move too: left in place it can't meet the
  face (outward) or is cut through (inward). The shim gives each face a
  motion. A picked plane translates by d along its normal, and a picked cone
  along its axis (a cone's offset is that translation); other curved faces
  take `SmSurface::CreateOffsetSurface`. Blends joined to a moving face through
  smooth edges translate by the shortest move that keeps them tangent to every
  neighbour, a small least-squares solve: a fillet under a raised top rises,
  one between a raised top and a drafted wall slides along the wall, one
  between two pushed faces follows both. A blend the move leaves in place
  (a vertical fillet under a raised top) stretches instead.
- The executive's tangent path works by lifting a face's old UV trim curve
  onto its new surface. That suits a true offset, where every face moves with
  its parameterisation, not a wall that keeps its surface while its top edge
  moves; which side it lifts from, and whether a neighbouring tangent edge has
  been mated yet, follow pointer-keyed map order, so results varied between
  runs (ASLR) as well as distances. Patch `a81ed1b` fixes the side.
- Edges re-intersected when nothing about them changed shape lose precision,
  and fail where two equal fillets meet at a corner: their mitre ends at the
  corner of the face they round, where the cylinders touch tangentially, a
  singular point of their intersection. So an edge moves rigidly by v when
  every face at both of its ends lands where v would put it (it moves by v, or
  v moves it within its own surface). The shim moves every such edge itself,
  before the executive mates any, putting the moved edge on each face's copy
  and marking it mated, and the executive leaves it alone. A straight edge
  that stays on its own line while its ends move along it (a wall's edge with
  a vertical fillet under a raised top) is put in the same way, as the line
  between its ends' predicted places: mated, it came from one face's lifted
  trim curve, and which face followed pointer order, so the same part rebuilt
  to different last digits, and different triangles, from run to run.
- The executive picks each extended edge's end among candidate points by
  nearness to the original vertex, which after a long move can be the wrong
  candidate, and through an unmated tangent neighbour's old UV curve, which can
  be the old corner. The shim predicts where vertices land (one shift every
  face around agrees on, or the point the moved planes and straight edges
  around it pin down) and moves the original brep's vertices there for the
  operation, restoring them after; patch `c282328` makes such a vertex a trim
  candidate.
- A closed intersection curve was trimmed to its far arc (patch `fdc9c9b`).
- Surfaces are extended by five times the "offset distance", which
  `LocalOperation` sets to half the bounding radius; the shim sizes it from the
  move. A cone slid along its axis keeps its old height band, so the executive
  would trim it on its NURBS extension; the shim widens the band (never past the
  apex). A face left on the NURBS extension of an analytic offset (a sphere's)
  gets the analytic surface back when all its edges lie on it and none is a
  seam. So does a cylinder the executive stretched along its axis (a fillet
  beside a pushed wall), which it makes a surface of extrusion: SMLib's
  `IsAnalytic` counts those, but the next offset wouldn't take one for a
  blend, and moving the top after pushing a wall under fillets left the solid
  open.
- Results are checked: manifold, the same faces and edges, every gap within the
  brep's tolerance, and the brep's own volume within 2% of its tessellation's.
  Moves that make a face vanish come back open, inside out (a box's walls
  shrunk to zero height: manifold, negative volume, a 32 mm gap) or with trims
  that disagree with the edges, and are refused.

Results (`rust/smlib/tests/cases.rs`, all against the part built directly at
its new size and, where known, its analytic volume, within 1e-6): plain,
chamfered, filleted and drafted tops; several faces of plain, chamfered and
filleted boxes, including every face at once; a wall pushed under fillets,
alone and with the top, and offsets one after another (a wall, then the top,
then another wall); a fully rounded box's top; a drafted wall; hole walls,
boss sides and tops, a frustum's side both ways, a countersink, a sphere-cut
box's sphere faces moved out, a torus. Every outcome was the same in eight
processes. Refused cleanly: fillets sliding along drafted walls (the result
tessellates right but its own volume is 20-30% off), a sphere-cut box with its
sphere faces moved in or its top moved up (the executive leaves it open), and
blends on a picked curved face (a rounded hole rim; their new shape is not a
translation).

Faces are picked by index, as everywhere else in `smlib`; `Brep::find_faces`
turns a `FaceQuery` into indices. An offset changes no topology, so every face
and edge keeps its label, attribute and all: the moved faces, the blends that
ride with them and the faces that stretch (`tests/labels.rs`). A drafted wall
is a B-spline surface whose control points share a plane, which
`SmSurface::IsPlanar` can't confirm below about 1e-5, so `outward_plane` (the
labels work's planar-face normal) also accepts what the offset's `plane_of`
does.
