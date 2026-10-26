# Rust crates for USD BRep's kernel

- `smlib-sys`: builds SMLib (`source/SMLib`) and SM_API (`source/SM_API`) with
  the `cc` crate, plus a C shim (`smlib-sys/shim`), and exposes the shim through
  committed bindgen bindings.
- `smlib`: the safe API. `Kernel`, `Brep` (primitives, placement, transforms,
  scales and mirrors, booleans, face sweeps, face offsets, rigid face moves,
  face merging, fillets, chamfers, volume, manifold check, face and edge
  queries, labels) and `Mesh`, a tessellation that keeps its links to
  the B-rep: per-triangle face ids, per-vertex surface (u, v), per-face
  surface parameters, per-edge polylines, and every face's and edge's label.

## Fillets and chamfers

Edges are picked by query, since indices change with every operation. An
`EdgeQuery` tests each edge's `EdgeGeometry` (sample points, the adjacent
faces and their labels, and their outward normals at its midpoint):
`MidpointOnPlane`, `InPlane`, `OnCircle`, `ParallelTo` (straight edges along a
direction), `BetweenFaces` (by normals), `Between` (by the faces' labels),
`Sharp` (by the angle between normals), `Indices` (positions a caller found in
`Brep::edges`, valid until the brep changes), and `And` / `Or` / `Not`.
`Brep::edge_length` integrates one edge's exact length, which costs more than
reading every edge's geometry.

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

## Placing solids on faces

- `Brep::faces` gives every face's surface, label, bounds and, for planar
  faces, outward normal. It's cheap; `Brep::face_centroid` integrates one
  face's area and centroid to 1e-7, which can take milliseconds on a face
  trimmed by many curves.
- `Brep::place` moves a solid rigidly into a `Frame` (`SmApiTransform`):
  build a primitive at the origin, then place it on a face's frame.
- `Brep::sweep_face` sweeps a planar face along its outward normal into a
  closed prism (`SmApiCreateCurveSweepFromFaces`), positive outward or
  negative inward. `Sweep::faces` says which prism face is the start cap,
  which the end cap, and which edge of the swept face each side came from.
  Union or subtract the prism to extrude or offset the face.
- `Brep::planar_faces` builds a sheet of planar faces from closed loops of
  lines, arcs and circles (`PlaneCurve`) in a `Frame`'s xy plane
  (`SmApiCreatePlanarFaces`): outer loops bound faces, loops inside them
  holes. Every face faces the frame's z, and each curve's edge carries the
  `BoundaryCurve`'s label, so the sides of a prism swept from a face trace
  back to their curves through `Sweep::faces`. This is a sketch's extrude.
- `Brep::revolve_face` turns a planar face about an axis in its plane, by
  any angle up to a whole turn either way, into a closed solid
  (`SmTopologySweep` with `SmRotationalSweepGeometry`). `Sweep::faces` says
  which face is the start, which the end (a whole turn has neither) and
  which edge each side was turned from; an edge on the axis has none. Sides
  are planes, cylinders, cones, spheres and tori. This is a sketch's
  revolve.
- `Brep::sweep_face_tapered` is `sweep_face` with draft: the sides lean in
  (positive) or out by an angle, corners sharp
  (`SmPrimitiveCreation::CreateTaperExtrude` per loop, holes cut out). Sides
  from lines are planes and from arcs cones; it refuses a draft whose sides
  meet before the end. This is a sketch's extrude with draft.
- `Brep::pipe` sweeps a circle, or a ring, along a path of lines and arcs
  in a `Frame`'s plane: a piece per curve (a disc swept or revolved, so
  walls are cylinders and tori), lines mitred where they meet at an angle,
  all joined. `Pipe::faces` says which face is the outer wall, the inner
  wall or a cap, and which path curve each wall was swept along.
- `Brep::shell` hollows a solid inward to a wall thickness, chosen faces
  left open (`SmApiShellBrepFull`). `Shell::faces` says which face stays
  outside, which is the rim of an opened face, which is the inner face
  offset from which face, and which is new. Faces that stay keep their
  labels. This is a part's Hollow.
- `Brep::merge_faces` merges faces lying on one surface and edges lying on
  one curve (`SmBrep::RemoveTopologicalEdgesAndVertices`). Booleans already
  do this along their intersection curves, so an offset face's walls come
  back as single faces either way.
- `Brep::combine` adds a copy of another brep's solids as lumps of their own
  (`SmBrep::MergeBrep`), with no boolean: for copies that neither touch nor
  overlap, such as a pattern's, it's their union in a fraction of the time,
  and one boolean then takes them all as one tool. Solids that touch or
  overlap need `union`.

## Moving, scaling and mirroring

- `Brep::transform(rotation, translation)` moves a brep rigidly
  (`SmApiTransform`): p goes to R p + t. `Brep::translate` and
  `Brep::rotate_about` (degrees about a line) are shorthands. R must be a
  rotation to 1e-9 (`SHIM_ERR_BAD_FRAME` otherwise). Every surface and curve
  keeps its type and (u, v), every face and edge its index and label, and the
  volume holds. A moved brep tessellates to the same triangles, its vertices
  where the motion takes the old ones to within `f32` rounding.
- `Brep::scale_about(point, factor)` scales uniformly about a point: planes,
  cylinders, cones, spheres and tori stay what they are, radii scale, angles
  don't, and the volume scales by the cube to 1e-8. Tolerances scale with the
  brep (see `PATCHES.md`).
- `Brep::scale_axes_about(point, factors)` scales each axis on its own; unequal
  factors turn every surface and curve into a NURBS first
  (`SmBrep::TurnToNURBS`), planes included.
- `Brep::scale_axes_keeping(point, factors, extrusions)` scales each axis on
  its own but keeps analytic whatever the scale leaves so: planes and lines
  always, a cylinder or cone along a scale axis whose other two factors are
  equal, spheres and tori under equal factors, circles in a plane the scale
  keeps round. Circles it stretches become ellipses. With `extrusions`, a
  cylinder stretched across its axis becomes an ellipse swept along a line
  (`SmSurfOfExtrusion`, reported as `Surface::Other`); otherwise it, like
  everything else the scale distorts, becomes a NURBS. A `ScaleReport` counts
  what became what. Faces keep their (u, v) domains, edges their intervals
  and UV trim curves, so it tessellates the same whatever ran first, which
  `scale_axes_about` doesn't.
- `Brep::mirror(point, normal)` reflects a brep across a plane. Every surface
  and curve keeps its kind and (u, v): planes, cylinders, cones, spheres,
  tori and surfaces of revolution (fillets on circular edges), lines,
  circles and ellipses. Faces' sides and the order of faces round each edge
  flip, so the result is a valid solid facing out, and every face and edge
  keeps its index and label. It costs about 0.5 ms on a 16-face part. A
  mirrored brep tessellates to the original's vertices reflected, every
  triangle wound the other way round; a planar face with many boundary
  points can split some quads along the other diagonal (a scale there and
  back does the same). `Brep::mirror_stock` is SMLib's own `SmBrep::Mirror`,
  for comparison: it turns every curve into a NURBS, reparameterises faces,
  rebuilds every UV trim curve and makes a torus a NURBS surface of
  revolution, in about 15 ms.
- `Brep::tolerance` reads the brep's tolerance; `Brep::write_file` and
  `Brep::read_file` write and read SMLib's binary format (labels don't come
  back).

`tests/transform.rs` and `tests/mirror.rs` check labels, surfaces, (u, v),
volumes, triangles and later fillets and booleans through each of them;
`tests/mirror.rs` also mirrors half a symmetric part and joins it back.

## Labels

`Brep::set_face_label` and `Brep::set_edge_label` attach a `Label` (a feature
id, a role string and an index) to a face or edge as an SMLib attribute
(`SmAttribute.h`). The kernel carries it from there:

- A face or edge an operation splits passes its label to both pieces;
  copies keep every label.
- When two labelled faces or edges merge, the result keeps the older label,
  the lower feature id (then the lower index, for two of one feature),
  whichever of the two the kernel keeps. The shim's
  attribute class overrides `SmAttribute::Merge` for this; with SMLib's
  default the result keeps the survivor's label, and which face survives
  depends on the order of the boolean's operands.
- Faces and edges an operation creates from nothing have no label: a
  fillet's or chamfer's new faces and their rails, and edges a boolean cuts.

`Brep::faces`, `Brep::edges` and `Mesh` report the labels. `tests/labels.rs`
checks them through booleans that split faces, fillets, chamfers, sweeps and
merges.

## Offset faces

`Brep::offset_faces(faces, distance)` moves the faces at the given indices
(from `Brep::faces`, or `Brep::find_faces` with a `FaceQuery`: `OnPlane`,
`Facing`, `Kind`, `Labelled`, `And` / `Or` / `Not`) along their outward
normals, the way Parasolid's face offset does: every other neighbour is
extended or trimmed along its own surface, so a chamfer grows or shrinks and
a drafted wall keeps its draft. Planes move flat, holes and bosses change
radius, cones slide along their axis, spheres and tori change radius. Fillets
and corner blends on a moved plane move with it, staying tangent to
everything they touch: they rise with a raised top, slide along a drafted
wall, follow two pushed faces at once.

It is `SmBrep::LocalOperation`, the kernel's own replace-a-face local
operation on the offset executive, reimplemented in the shim, with three
kernel patches (see `PATCHES.md`). With nothing in the way the topology
never changes, so every face and edge keeps its label and its index. The result is checked: a move that
would make a face vanish fails (`SHIM_ERR_CHANGED_TOPOLOGY`, or
`SHIM_ERR_BROKE_SOLID` when the kernel gets that far), and any result that
isn't manifold, has gaps over the brep's tolerance, or has a volume that
can't be trusted (faces whose integrals disagree with their edges, and a
tessellation that doesn't bear the volume out) fails with
`SHIM_ERR_BROKE_SOLID`.
Blends on a moved curved face (a rounded hole rim) fail with
`SHIM_ERR_UNSUPPORTED_FACE`, and so, for now, do fillets sliding along drafted
walls and a few moves of sphere faces that the kernel gets wrong. Every planar
face's loops are checked too: none may cross, and every inner loop must stay
inside the outer one.

A plane pushed in can run into features standing beside it, such as a boss or
a hole on a lid near the edge above a pushed wall. It cuts through them, as
Parasolid does: each is set aside, the face is offset without it, and it comes
back cut by the face's new plane. A boss is joined back as far as the plane
leaves it, and a hole is cut back out as a notch. Those faces keep their
labels. A cut boss's new face is unlabelled where it doesn't merge into the
moved face. `Brep::offset_obstacles` lists what an offset runs into.
`SHIM_ERR_OBSTRUCTED` means one of them can't be set aside: it doesn't stand
entirely on one side of a planar face, or something else stands in its way.
`Brep::offset_faces_unchecked` returns the kernel's bare result, for diagnosis
only.

`Brep::tapered_prism`, `Brep::cone` and `Brep::torus` build drafted prisms,
frustums and rings.

```sh
cargo run --release -p smlib --example offset_face   # timings
```

## Moving faces

`Brep::move_faces(faces, motion)` moves the faces at the given indices by a
`RigidMotion`: an optional `Turn` about an axis through a point (degrees,
anticlockwise looking down the axis from its tip), then a translation. Every
other neighbour is extended or trimmed along its own surface to meet them,
on the same local operation as `offset_faces`, so a face can lean at any
angle (tilted panels, angled facets, a draft one wall at a time) and a
feature's faces moved together slide it across the face it stands on. A
planar face moved within its own plane stays where it is.

A translation carries the blends on the moved faces as an offset does,
each sliding the shortest way that keeps it tangent to its neighbours; a
turn can't, and a round on a turned face is refused. Faces and edges keep
their order and labels, and `Moved::faces` says what happened to each:
`Moved`, `Carried` (a blend slid with them), `Followed` (extended or
trimmed to meet them) or `Kept`. The result is checked as an offset's is.
When the faces round the moved ones can't reach them (one would vanish or
cross another, the solid wouldn't close, a round can't follow, or the
kernel finds no way to meet them) the move is refused with
`Error::InvalidArgument`, its reason starting "the faces round them can't
reach them". A turned cylinder comes back open and is refused; turned
planes are exact. Nothing in the way is pushed through.

## Face-by-face tessellation

`Brep::tessellate` runs `SmApiTessellate`, which drives SMLib's tessellator
(`SmTess`) in three phases on a copy of the brep: phase 1 splits every edge
for the whole solid and builds each face's surface cache and boundary
polygon; phase 2 triangulates one face inside its boundary; the polygon
output turns that into 3D polygons. Phase 2 and the output read only their
face's cache, boundary, surface, trim curves and vertices.

`Brep::begin_tessellation` runs phase 1 as `SmApiTessellate` does and
returns a `Tessellation`; `Tessellation::mesh` then triangulates only the
faces asked for, each into a poly brep of its own. A face's triangles and
vertices come out bit for bit as the whole tessellation gives them, whichever
other faces are triangulated, so a caller can keep earlier faces' triangles
and stay watertight: every shared edge is split the same on both sides.

`Tessellation::face_keys` says when that's safe. A face's `FaceKey` is a
128-bit digest (MurmurHash3) of everything phase 2 and the output read for it:
the parameters, the face's surface and every edge and trim curve as SMLib
writes them to a file (`WriteToDB`, binary, so every bit counts), its domain,
flags and tolerances, its vertices, and its boundary polygon after phase 1.
Equal keys mean identical triangles. A neighbour's change that splits a
shared edge differently changes the boundary, and so the key, of both faces.

`Brep::tessellate_some` goes further and runs phase 1 for the faces asked for
only, through `SmTess::Phase1SetupFaces` (a kernel patch, `PATCHES.md`): it
builds surface caches for those faces and the faces sharing their edges,
splits only their edges, and builds only their boundaries. Phase 1 acts on one
face or one edge at a time until a boundary needs refining (its polygon
crosses itself), so then the faces come out as a whole phase 1 leaves them. Its
keys, from `Brep::tessellation_keys`, are taken before tessellating: the
parameters, the face's surface, domain, tolerances and loops (curves,
intervals, trim curves and vertices), and for each edge the surfaces, domains
and order of the faces across it. A face keyed as one an earlier
`tessellate_some` triangulated doesn't refine either, and gets the same
triangles. `tessellate_some` checks the faces it sets up and returns `None`
when one refines, for the caller to tessellate the whole solid instead.

A face's surface cache (`SmTessSrfCache`), a subdivision of its surface to the
parameters, depends only on the surface and the parameters, and is often the
largest part of setting a face up. With a `SurfaceCaches` store,
`tessellate_some` keeps the caches of faces it set up only to split shared
edges, and faces whose surfaces are identical take them later rather than
building their own (`SmTessSrfCacheStore`, the second kernel patch). In a drag
the faces around the moving ones keep theirs from step to step.

Triangulating a face subdivides and labels its cache, so a wanted face takes a
deep copy (`SmTessSrfCache::Copy`, the third kernel patch): the tree's nodes
and their data are fixed-size blocks without virtual functions, copied byte for
byte with their pointers mapped, and a B-spline surface's subdivision surfaces
are copied with their exact knots. The original stays in the store. A wanted
face that builds its cache leaves a copy there when one of the store's last two
calls built a cache for the same surface, so in a drag the faces whose
boundaries move while their surfaces stay take copies from the third step on,
and faces whose surfaces change every step cost no copies.

`tests/tessellation.rs` checks every face against `Brep::tessellate`, partial
meshes against whole ones, and kept faces standing in across edits, both ways
and with surface caches kept (copied for NURBS faces too), and that
`tessellate_some` declines faces whose boundaries need refining.

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
