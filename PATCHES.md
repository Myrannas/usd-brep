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
| `a9a7a28` | `SmOffsetExecutive.cpp` | Found offsetting faces: `DoExtendedOffset` mates edges in the order of `m_vEUToEU`, an `unordered_map` keyed by edgeuse addresses, which also picks each edge's target side; the same offset gave different last digits, vertices and trim curves from one process to the next. | For `SM_OO_LOCAL_OPERATION`, mate them in the original brep's order: its edges in turn, each edge's edgeuses in turn. |
| `4d13c5e` | `SmTess.h`, `SmTess.cpp` | Found keeping tessellated faces across rebuilds: `Phase1SetupBrep` builds every face's surface cache and splits every edge, so a small edit still pays for the whole solid's edge preparation. | Add `SmTess::Phase1SetupFaces`, which sets up only the faces given (with caches for the faces sharing their edges, which split only those edges) exactly as `Phase1SetupBrep` would, and reports refinement or anything `Phase1SetupBrep` would recover from; and a `SplitFaceEdgesToUVSize` overload that splits only edges shared with given faces. Existing calls are unchanged. |
| `7efc9ef` | `SmTess.h`, `SmTess.cpp` | Found keeping tessellated faces across rebuilds: a face's surface cache, often the largest part of setting it up, is built again every tessellation though it depends only on the surface and the parameters. | `Phase1SetupFaces` takes an optional `SmTessSrfCacheStore`: faces take its caches for identical surfaces, and faces set up only to split shared edges leave theirs there. `SmTessSrfCache::SetFace` moves a cache to another face and tessellation. |
| `76a3dbc` | `SmMemBlockMgr.h`, `SmTess.h`, `SmTess.cpp` | Found keeping tessellated faces across rebuilds: Phase 2 subdivides and labels a face's surface cache, so a kept cache serves one face once, and faces whose surfaces stay while their boundaries move (beside a dragged fillet) build theirs in every tessellation. | Add `SmTessSrfCache::Copy`, a deep copy of a cache as `BuildTree` leaves it: the tree's node and node-data blocks byte for byte with their pointers mapped (`SmMemBlockMgr::CopyFrom`), B-spline subdivision surfaces with their exact knots; NULL for a cache holding anything else. `Phase1SetupFaces` gives the faces it tessellates copies and leaves originals in the store, which also gets copies of caches it builds when `SmTessSrfCacheStore::WantsCopy` says to. `SetFace` repoints nodes from the old face's surface to the new one's. |
| `a9b66ee` | `SmBrep.cpp` | Found imprinting: `RemoveTopologicalEdgesAndVertices` tests each edge with `SmEdge::IsTopological`, which can expand one face's surface to cover the other's as a side effect, before it looks for the edge in the keep list. An edge kept between two pieces of a curved face (a panel imprinted on a crate's round, merged back off the round at the other end) leaves a piece on the stretched surface: its box grows by 12 mm, its area integral fails and a later union over it fails with `SM_ERR`. | Look for each edge in the keep list before testing it. Without a keep list the merge runs as before. |
| `cd41f50` | `SmTorus.cpp` | Found revolving: `SmTorus::IsNurbSurfaceTorus` snaps the generator circle's start and end angles to 0 and 360 only within `SM_EFF_ZERO` (1e-12 degrees). A circle revolved a whole turn on a tilted plane came out with both 2.7e-12 degrees off zero, and the torus got an empty (u, v) domain. | Snap within `SM_EFF_ZERO_DEG` and turn an end no further than the start into a whole turn, as `RebuildSTEPFromNURBParameters` in the same file does. Nothing in the kernel calls it (its call in `CopyAndAddAnalytics` is commented out); the shim's revolve does. |

`source/SM_API` is unmodified. Each offset patch changes only
`SM_OO_LOCAL_OPERATION` (`SmBrep::LocalOperation` and the shim's face
offset); solid offsets and shells take the code paths they always did. The
tessellation patches only add methods and a class (`SmMemBlockMgr` gains a
block copy and block accessors): `SmApiTessellate` and `Phase1SetupBrep`
take the code paths they always did. The merge patch changes only merges
given a keep list: `RemoveTopologicalEdgesAndVertices()` without one, as
every boolean's post-process and `Brep::merge_faces` call it, is as it was.
The torus patch changes only `SmTorus::IsNurbSurfaceTorus`, which only the
shim calls.

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
  been mated yet, follow the executive's mating order, so results varied with
  distances. Patch `a81ed1b` fixes the side. The order itself was the
  edgeuses' addresses, so results varied between runs too (ASLR); patch
  `a9a7a28` puts it in the brep's order.
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
  trim curve, which is off in the last digits, and which face followed the
  mating order. So is a seam (a hole's, under a pocket floor that moves): its
  two copies are both uses of one edge of the copied face, and a curve on a
  closed surface's seam drops onto it as a pair of UV curves, one on each
  side, so each copy takes the one on its own original's side.
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
  brep's tolerance, and a volume that can be trusted (see "Notes from the
  rebuild trims" below). Moves that make a face vanish come back open, inside
  out (a box's walls shrunk to zero height: manifold, negative volume, a 32 mm
  gap) or with faces whose integrals disagree with their edges, and are
  refused.

Results (`rust/smlib/tests/cases.rs`, all against the part built directly at
its new size and, where known, its analytic volume, within 1e-6): plain,
chamfered, filleted and drafted tops; several faces of plain, chamfered and
filleted boxes, including every face at once; a wall pushed under fillets,
alone and with the top, and offsets one after another (a wall, then the top,
then another wall); a fully rounded box's top; a drafted wall; hole walls, boss
sides and tops, a frustum's side both ways, a countersink, a sphere-cut box's
sphere faces moved out, a torus. A hole's seam under a pocket floor moved
either way comes out bit for bit as in the part built at its new depth. Every
outcome came out bit-identical in six processes once patch `a9a7a28` fixed the
mating order; before it, the sphere-cut box and a top-rounded box moved at its
top and walls varied in the last digits. Refused cleanly: fillets sliding along
drafted walls (the result tessellates right but its own volume is 20-30% off),
a sphere-cut box with its sphere faces moved in or its top moved up (the
executive leaves it open), and blends on a picked curved face (a rounded hole
rim; their new shape is not a translation).

Faces are picked by index, as everywhere else in `smlib`; `Brep::find_faces`
turns a `FaceQuery` into indices. An offset changes no topology, so every face
and edge keeps its label, attribute and all: the moved faces, the blends that
ride with them and the faces that stretch (`tests/labels.rs`). A drafted wall
is a B-spline surface whose control points share a plane, which
`SmSurface::IsPlanar` can't confirm below about 1e-5, so `outward_plane` (the
labels work's planar-face normal) also accepts what the offset's `plane_of`
does.

## Notes from the rebuild trims (no patch needed)

Profiling a part rebuilt as one of its features was dragged (a pocket under a
fillet, a face offset and another fillet) put about 40% of each rebuild into
checks and reads that the results didn't need in full. Results are unchanged:
meshes, volumes and names come out bit for bit as before, apart from results
that already varied from run to run.

- The face offset's volume check used to compare the brep's own volume
  (`SmApiBrepComputeVolume`) with what its tessellation encloses, to 2%. The
  tessellation cost more than the offset's own surfaces. SMLib integrates
  volume as the sum over faces of the integral of (z - z0) n_z dA, so raising
  the origin by h changes each face's share by h times its projected area,
  and for a face that integrates right that is also what its edges enclose
  (Stokes). The check now integrates the volume about two origins a box
  height apart and holds each face's change to its edges. Only when a face
  disagrees by more than 1e-3 of the shares, or the volume isn't positive,
  does the tessellation decide, as before.
  - The failure the old check was built for, fillets sliding along drafted
    walls, is a face whose integral over its own trims disagrees with its own
    edges: the fillets' projected areas came out 7.6 times too small, though
    their trims lie on their edges and their UV loops close. The new test
    scores it 6e-2; every valid result in the tests scores 2.6e-4 or less
    (a torus; 3.5e-5 for everything else), and the old test still refuses it.
  - Each face's volume is integrated alone (`SM_PPF_VOLUME`), not with every
    moment. The integrator's step control follows the area whatever is asked
    for, so the volume is bit for bit what `SmApiBrepComputeVolume` gives.
- `Brep::edges` no longer integrates each edge's exact length
  (`SmApiEdgeComputeLength`, an ODE per edge), which nothing read and which
  was 9% of the rebuild. `Brep::edge_length` integrates one edge's on demand.

## Notes from the sketches work (no patch needed)

`Brep::planar_faces` builds faces from lines, arcs and circles in a plane, for
a sketch's extrude to sweep with `Brep::sweep_face`. The shim's
`sm_brep_create_planar_faces` makes the curves in the plane z = 0 with
`SmApiCreateLineSegment`, `SmApiCreateArc` and `SmApiCreateCircle` (SM_API
makes its arcs and circles there), the faces with `SmApiCreatePlanarFaces`,
and then moves the sheet into the caller's frame with `SmApiTransform`.

- Each new edge holds the input curve it was made from, so the shim labels
  edges by curve pointer, with the nearest curve as a fallback.
  `sweep_face`'s sides find their source edges by where their start edges
  lie, so a prism's side traces back to its curve through the edge's label.
- `SmSurface::CreatePlaneWith3DCurves` picks each face's normal from the
  curves, not from the loops' direction: in the tests, both of two nested
  faces came out facing -z. The shim turns such faces over with
  `SmFace::SwapUV` (which `SmPlane` follows by its swap-UV bit), so every
  face faces the frame's z and a positive sweep runs along it.
- `SmApiCreateArc` runs counter-clockwise between angles in [-360, 360]; the
  shim normalises clockwise arcs, since the faces don't mind which way a
  curve runs.
- `SmApiCreatePlanarFaces` doesn't refuse every bad arrangement: two
  crossing squares, or two touching at a corner, still make faces. Callers
  hand over loops that neither cross nor touch.
- `CreateCurveSweepFromFaces` sweeps each loop's curves as NURBS, so every
  side of a swept prism was a NURBS surface, planar for lines (as for
  `sweep_face` on a solid's face). See "Notes from the analytic sweeps" below.
- Two faces of one feature can merge too, such as the ends of an extrude's
  adjacent regions. `LabelAttribute::Merge` compared only feature ids, so the
  merged face kept whichever label the survivor had, and two flush boxes of
  one feature kept a mix of both boxes' indices in one operand order. It now
  keeps the lower (feature, index), whichever face survives.

## Notes from pushing faces through (no patch needed)

A face pushed in can run into features standing beside it: in game-cad's
`bugs.gcad` a crate's front wall pushed in 30 or 46 mm passes under two bosses
on the lid near the front edge. The lid trims back with the wall, and its
inner loops (the bosses' rims) would then cross its new boundary or fall
outside it: new topology. The local operation doesn't make it. It moves the
lid's edge anyway and returns a solid that passed every check: manifold, the
same faces and edges, no gaps, and a volume within 2% because the bosses are
small. At 46 mm `SmApiTessellate` then fails with `SM_ERR`, and at 30 mm it
tessellates nonsense.

- The checks now include `loops_stay_apart`: on every planar face, no two
  loops cross, and every inner loop lies inside the outer one. Loops are
  polylines in the face's plane through 24 points per edge. Curved faces are
  skipped, since their loops can wrap around the seam in (u, v). The broken
  results above fail it with `SHIM_ERR_CHANGED_TOPOLOGY`.
- `find_obstacles` looks, before the operation, for inner loops of the faces
  beside the moving ones (from `plan_offset`'s motions) that reach into the
  space a pushed plane sweeps. The obstacle is the set of faces reachable
  from such a loop without crossing its face. It's a boss when it stands
  entirely on that face's outer side, a hole when it's entirely inside.
- `push_through` sets each one aside with a box on its face's plane, just
  covering it, and checks that nothing else reaches into the box. A boss
  comes out as `part ∩ box` and the part keeps `part − box`; a hole is filled
  with `part ∪ box` and kept as `box − part`. Once `RemoveTopologicalEdgesAndVertices`
  has merged the box's coplanar cuts back into the face, the picked faces are
  found again by a `SM_AB_COPY` attribute and offset as usual. Each boss is
  then cut by a box behind the pushed plane and joined back, unless it lies
  wholly beyond the plane. Each hole is cut back out. The result is checked
  like an offset, except for its face and edge counts, and must have as many
  shells as before. A boss on a side wall near the pushed face's corner
  shows why: the round between them moves in, the wall pulls back from under
  what's left of the boss, and that would come back as a separate lump. It's
  refused instead.
- A boss cut flush merges its new face into the pushed face, and a cut
  cylinder's side comes back in two pieces, split where its seam was. Volumes
  match the part without the features plus what's left of them to 1e-11
  (`tests/cases.rs`).
- With nothing in the way, the offset takes the same path as before, bit for
  bit. `offset_faces_unchecked` keeps the bare local operation, for
  diagnosing the kernel and for testing callers with a broken solid.

## Notes from the analytic sweeps (no patch needed)

`Brep::sweep_face` gives each side swept from a line a plane and each side
swept from an arc or circle a cylinder, where the sweep made NURBS surfaces,
so callers see a sketch extrude's walls, and an extruded face's, as what they
are. The shim's `make_sides_analytic` does it after the sweep:

- `CreateCurveSweep` joins each loop's curves into one composite NURBS,
  sweeps that along the path and splits the surface into faces wherever it
  is no better than tangent continuous: one face per curve, a circle staying
  whole. The path direction comes out
  degree 1 with two poles: a line's side is a bilinear patch (degree 2 with
  three poles when a loop mixes lines and arcs), an arc's or circle's a
  rational degree (2, 1) patch.
- `SmBSplineSurface::CopyAndAddAnalytics`, which `SmPrimitiveCreation` uses
  for its planes, recognises both: `SmPlane::IsNurbSurfacePlane` checks the
  plane against the NURBS at ten points and reparameterises it onto the
  NURBS's domain, and `SmCone::IsNurbSurfaceCone` takes rational degree
  (2, 1) patches with a circle's knots, making an `SmCone` with a zero
  semi-angle. SMLib's analytic surfaces are B-splines too, so a side keeps its
  (u, v): the shim swaps the surface only when the copy's domain is the
  NURBS's, and puts back the edges' UV trim curves, which `SmFace::SetSurface`
  drops. Any result other than a plane or a cylinder (a surface of extrusion)
  is discarded and the side keeps its NURBS.
- Prisms keep their volumes to the last digit and their largest vertex, edge
  and face gaps (5e-15 for a rectangle with a round hole); booleans and face
  merges with the new sides work as before (`tests/planar_faces.rs`,
  `tests/labels.rs`).

## Notes from the imprint work (no patch needed)

`Brep::imprint` and `Brep::imprint_classified` (`sm_brep_imprint`),
`Brep::merge_faces_keeping` (`sm_brep_merge_faces_keeping`) and Offset Face's
steps are shim and Rust only. What the shim relies on:

- `SmMerge::ManifoldBoolean` with `SM_BO_IMPRINT` (or `bImprinting`) returns
  as soon as the intersection is inserted and the loops are checked, before
  any face is deleted and before the post-process that merges faces on one
  surface, so `bBooleanPostProcess` makes no difference: the imprint
  survives with it on or off. `bImprintAndClassifyFaces` returns before the
  post-process too, after sorting the faces into keep and delete lists,
  which it marks with `SM_AI_BOOLEAN_SAVE` and `SM_AI_BOOLEAN_DELETE`
  (reference attributes) instead of deleting any. The shim reads and removes
  the marks.
- No public call fits. `SmApiBooleanWithOptions` keeps or deletes the tool
  but doesn't say which tool face each new edge lies on, and
  `SmApiBooleanWithCurves` returns the intersection edges with the default
  options, under which `SM_BO_IMPRINT` deletes the tool. So the shim builds
  the `SmMerge` itself with `bKeepOtherBrep`, reads each common edge's mate on
  the tool (`SmTopologyIntersector::GetCommonEdges`) and the labels of the
  mate's faces, and then deletes the tool.
- Split faces keep their label on every piece, and the pieces' UV domains
  are recomputed to fit them, so they and their neighbours key differently
  for tessellation. Where the tool only touches a face, as a cylinder tangent
  to a plane does, the imprint still splits the face along the touching
  line.
- `RemoveTopologicalEdgesAndVertices` without a keep list removes every
  imprinted edge between two pieces of one surface. With one, it used to
  test every edge before skipping the kept ones, which stretched a kept
  curved piece's surface (patch `a9b66ee`). That includes
  `Brep::merge_faces` and the shim's own merges in Offset Face's
  push-through.
- `SmBrep::LocalOperation` can't move a face away from another piece of its
  own surface, as a panel an imprint leaves on a lid: the offset executive
  only stretches and trims neighbours, and a neighbour on the face's own
  surface never meets its new place. It reports a broken solid, and the
  push-through took the wrong region (a 3 mm sink of an imprinted panel on
  the sample crate took 357165 mm³ for 36000). So the shim's Offset Face
  recognises such moves first and builds them from the faces swept into
  prisms, joined or cut, with a wall where the face leaves its own surface
  (see `sm_brep_offset_faces`): planar faces only, every other neighbour
  standing along the move, and the prisms adding or taking exactly their
  own volume, or it refuses (`SM_SHIM_ERR_UNSUPPORTED_STEP`,
  `SM_SHIM_ERR_OBSTRUCTED`). Steps never reach the push-through.
- The push-through's merges (`RemoveTopologicalEdgesAndVertices` after the
  boxes that set obstacles aside, and in the plugs) merged every pair of
  faces on one surface, so an imprint anywhere on the part was gone after
  a face was pushed past a boss. Before pushing through, the shim now marks
  every edge that joins two faces on one surface (`same_surface`) with a
  copy attribute, which split edges pass to both pieces, and each of those
  merges keeps the marked edges (`merge_unmarked`): only the faces its own
  cuts split merge back. With no such edge it is the plain merge, as before.
- Each step wall takes the label of the edge it stands on, if that has one,
  as the moved face takes its own: a caller that labelled its imprinted
  edges can tell the walls apart. The sides beside every other neighbour
  stay unlabelled, since they merge into it and a label of their own could
  outrank its.
- A label's role holds up to 255 bytes (`SM_LABEL_ROLE_CAPACITY` 256, was
  64), so a role can carry a sketch region's name: an imprinted panel's
  piece is named by every shape it lies inside and outside of.

## Notes from the solids work

Revolve, Draft, Pipe and Hollow, for game-cad's features of those names.
Shim and Rust only, apart from the torus patch (`cd41f50`).

### Revolve

`Brep::revolve_face` (`sm_brep_revolve_face`) turns a planar face about an
axis in its plane, part of a turn or a whole one, into a closed solid, and
says what each face came from as `sweep_face` does (`Sweep`).
`SmApiCreateRotationalSweep` takes curves, not a face, and doesn't say which
side came from which curve, so the shim runs `SmTopologySweep` with
`SmRotationalSweepGeometry` itself, as `SmPrimitiveCreation`'s rotational
sweep does, and reads its maps: `GetToHigher` gives each edge's side and
`GetToSame` the face's turned copy.

- The face is copied alone into a new brep (`SmBrep::CopyFaces` without
  merging); deleting a solid's other faces from a copy of it left a face
  the sweep failed on (`SM_ERR`). The copy's edges are new curves, traced
  back to the face's edges by where their midpoints lie.
- The face is moved into a frame with the axis along z and itself in the xz
  plane, swept there, and the solid moved back, so the sweep and the
  recognisers below see the same arithmetic whatever plane the face is in.
- The axis must lie in the face's plane (direction to 1e-9, point to the
  brep's tolerance) and may not run through the face's inside, which
  `SmFace::CurveOnClassify` decides as `SmPrimitiveCreation` does (where it
  returns `SM_ERR_AXIS_INSIDE_FACE`). Along the boundary or through a vertex
  is fine. An edge on the axis sweeps nothing and has no side; a face
  touching the axis at a point turns into cones meeting at their apex, which
  SMLib counts as a manifold solid.
- Part of a turn sweeps the face: the face is the start, its copy the end.
- A whole turn can't: the face would meet its own copy, and the result isn't
  a manifold solid. Each loop's edges are swept into closed sides instead,
  the face deleted and the sides stitched (`SmBrep::StitchAndOrient`), as
  `SmPrimitiveCreation`'s whole turns do. `StitchAndOrient` makes every
  closed shell a solid region, so a hole's ring came out solid: the volume
  was the face's without its hole. Each loop is now turned in a brep of its
  own and the holes' rings cut from the outer loop's solid
  (`SmApiBooleanDifference`), which leaves closed voids as a box minus an
  enclosed ball does. The sides are traced through the booleans by a
  temporary label.
- The sweep makes planes, cylinders, cones and spheres analytic itself
  (`CopyAndAddAnalytics`), but a torus stays a surface of revolution:
  `CopyAndAddAnalytics` has its torus test commented out. The shim calls
  `SmTorus::IsNurbSurfaceTorus` and keeps the torus when it has the
  surface's domain and points (patch `cd41f50`).
- Results (`rust/smlib/tests/revolve.rs`, on a tilted plane, whole and part
  turns, both ways): tubes, discs, cones, frustums, balls, rings, a tube with
  a ring-shaped hole and a box's face turned about its own edge, all manifold,
  volumes within 1e-6 of Pappus' (1e-8 in practice), meshes winding outward
  round them. Each side is the edge's: plane, cylinder, cone, sphere or
  torus.
- Turns from 359.999 degrees to just under 360 fail (`SHIM_ERR_BAD_SWEEP`:
  the start and end nearly meet); 359.99 works. A 0.01 degree turn leaves
  its cylinders as surfaces of extrusion (0.001 and 0.1 degrees don't).

