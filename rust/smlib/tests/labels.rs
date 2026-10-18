//! Labels, placement and face sweeps on the real kernel: the pieces a feature
//! tree needs to put solids on faces and to keep naming faces through every
//! operation.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{
    Brep, EdgeBlend, EdgeQuery, Error, FaceCentroid, FaceGeometry, FilletOptions, Frame, Kernel,
    Label, LabelMatch, Surface, SweptFace, TessellationParams, Tolerance,
};

type Job = Box<dyn FnOnce(&Kernel) + Send>;

fn on_kernel<R: Send + 'static>(f: impl FnOnce(&Kernel) -> R + Send + 'static) -> R {
    static JOBS: OnceLock<mpsc::Sender<Job>> = OnceLock::new();
    let jobs = JOBS.get_or_init(|| {
        let (tx, rx) = mpsc::channel::<Job>();
        std::thread::Builder::new()
            .name("smlib-kernel".into())
            .spawn(move || {
                let kernel = Kernel::init().expect("kernel init");
                for job in rx {
                    let _ = panic::catch_unwind(AssertUnwindSafe(|| job(&kernel)));
                }
            })
            .unwrap();
        tx
    });
    let (reply, result) = mpsc::sync_channel(1);
    jobs.send(Box::new(move |k| {
        let _ = reply.send(f(k));
    }))
    .unwrap();
    result.recv().expect("kernel job panicked")
}

const BLOCK: u32 = 1;

/// A block's role for each outward normal.
fn axis_role(n: [f64; 3]) -> &'static str {
    const ROLES: [(&str, [f64; 3]); 6] = [
        ("+X", [1.0, 0.0, 0.0]),
        ("-X", [-1.0, 0.0, 0.0]),
        ("+Y", [0.0, 1.0, 0.0]),
        ("-Y", [0.0, -1.0, 0.0]),
        ("+Z", [0.0, 0.0, 1.0]),
        ("-Z", [0.0, 0.0, -1.0]),
    ];
    ROLES
        .iter()
        .find(|(_, axis)| dot(*axis, n) > 0.999)
        .map(|(role, _)| *role)
        .unwrap_or_else(|| panic!("{n:?} isn't along an axis"))
}

/// A box centred on `centre`, its faces labelled by outward normal and its
/// edges by the two faces they join.
fn labelled_box(k: &Kernel, feature: u32, centre: [f64; 3], size: [f64; 3]) -> Brep {
    let min = [0, 1, 2].map(|i| centre[i] - size[i] / 2.0);
    let mut solid = Brep::cuboid(k, min, size).unwrap();
    for face in solid.faces().unwrap() {
        let role = axis_role(face.normal.expect("box faces are planar"));
        solid
            .set_face_label(face.index, &Label::new(feature, role))
            .unwrap();
    }
    for edge in solid.edges().unwrap() {
        let [a, b] = edge.face_labels.map(|label| label.unwrap().role);
        let role = format!("{a}{b}");
        solid
            .set_edge_label(edge.index, &Label::new(feature, role))
            .unwrap();
    }
    solid
}

/// A cylinder of radius `r` and height `h` standing on the origin along +z,
/// labelled `side`, `top` and `bottom`, then placed in `frame`.
fn labelled_cylinder(k: &Kernel, feature: u32, r: f64, h: f64, frame: &Frame) -> Brep {
    let mut solid = Brep::cylinder(k, [0.0; 3], r, h).unwrap();
    for face in solid.faces().unwrap() {
        let role = match face.normal {
            None => "side",
            Some(n) if n[2] > 0.0 => "top",
            Some(_) => "bottom",
        };
        solid
            .set_face_label(face.index, &Label::new(feature, role))
            .unwrap();
    }
    solid.place(frame).unwrap();
    solid
}

/// A face with its area and centroid.
#[derive(Clone, Debug)]
struct Face {
    geometry: FaceGeometry,
    area: f64,
    centroid: [f64; 3],
}

impl std::ops::Deref for Face {
    type Target = FaceGeometry;

    fn deref(&self) -> &FaceGeometry {
        &self.geometry
    }
}

fn measure(solid: &Brep) -> Vec<Face> {
    solid
        .faces()
        .unwrap()
        .into_iter()
        .map(|geometry| {
            let FaceCentroid { area, centroid } = solid.face_centroid(geometry.index).unwrap();
            Face {
                geometry,
                area,
                centroid,
            }
        })
        .collect()
}

fn labels(faces: &[Face]) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for face in faces {
        let key = face
            .label
            .as_ref()
            .map_or_else(|| "-".to_owned(), ToString::to_string);
        *counts.entry(key).or_default() += 1;
    }
    counts
}

fn find<'a>(faces: &'a [Face], feature: u32, role: &str) -> Vec<&'a Face> {
    faces
        .iter()
        .filter(|f| {
            f.label
                .as_ref()
                .is_some_and(|l| l.feature == feature && l.role == role)
        })
        .collect()
}

fn assert_near(what: &str, actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    let gap = dot(sub(actual, expected), sub(actual, expected)).sqrt();
    assert!(
        gap <= tolerance,
        "{what}: {actual:?} isn't within {tolerance} of {expected:?}"
    );
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[test]
fn faces_report_area_centroid_and_outward_normal() {
    let faces = on_kernel(|k| {
        let block = labelled_box(k, BLOCK, [0.0; 3], [100.0, 60.0, 40.0]);
        let hole = Brep::cylinder(k, [20.0, 0.0, -20.0], 10.0, 40.0).unwrap();
        measure(&block.difference(hole).unwrap())
    });
    assert_eq!(faces.len(), 7);
    let top = find(&faces, BLOCK, "+Z")[0];
    assert_eq!(top.normal, Some([0.0, 0.0, 1.0]));
    // A 100 × 60 rectangle less a circle of radius 10 at x = 20.
    let hole_area = std::f64::consts::PI * 100.0;
    let area = 6000.0 - hole_area;
    // Centroids are integrated to a relative accuracy of 1e-7.
    assert!((top.area - area).abs() < 1e-7 * area, "{}", top.area);
    assert_near(
        "top centroid",
        top.centroid,
        [-20.0 * hole_area / area, 0.0, 20.0],
        1e-5,
    );
    let bottom = find(&faces, BLOCK, "-Z")[0];
    assert_eq!(bottom.normal, Some([0.0, 0.0, -1.0]));
    let wall = faces.iter().find(|f| f.label.is_none()).unwrap();
    assert!(matches!(wall.surface, Surface::Cylinder { .. }));
    assert_eq!(wall.normal, None, "curved faces have no normal");
    assert_near("hole wall centroid", wall.centroid, [20.0, 0.0, 0.0], 1e-5);
}

#[test]
fn place_moves_a_solid_into_a_frame() {
    // A cylinder along +z, placed on the +X face of a block: its axis follows
    // the frame's z axis, x × y.
    let frame = Frame::new([50.0, 5.0, -3.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
    assert_eq!(frame.z_axis(), [1.0, 0.0, 0.0]);
    assert_eq!(frame.to_world([1.0, 2.0, 3.0]), [53.0, 6.0, -1.0]);
    let faces = on_kernel(move |k| measure(&labelled_cylinder(k, 2, 8.0, 10.0, &frame)));
    let top = find(&faces, 2, "top")[0];
    assert_eq!(top.normal.map(|n| n.map(f64::round)), Some([1.0, 0.0, 0.0]));
    assert_near("top", top.centroid, [60.0, 5.0, -3.0], 1e-5);
    let bottom = find(&faces, 2, "bottom")[0];
    assert_near("bottom", bottom.centroid, [50.0, 5.0, -3.0], 1e-5);
    let Surface::Cylinder { axis, radius, .. } = find(&faces, 2, "side")[0].surface else {
        panic!("the side is a cylinder");
    };
    assert_near("axis", axis.map(f64::abs), [1.0, 0.0, 0.0], 1e-12);
    assert!((radius - 8.0).abs() < 1e-12);

    // Frames must be orthonormal.
    let skewed = Frame {
        origin: [0.0; 3],
        x_axis: [1.0, 0.0, 0.0],
        y_axis: [0.6, 0.8, 0.0],
    };
    let error = on_kernel(move |k| {
        Brep::cuboid(k, [0.0; 3], [1.0; 3])
            .unwrap()
            .place(&skewed)
            .unwrap_err()
    });
    assert!(matches!(error, Error::Kernel { .. }), "{error}");
}

#[test]
fn labels_survive_booleans_and_splits() {
    let (faces, edges, mesh) = on_kernel(|k| {
        let block = labelled_box(k, BLOCK, [0.0; 3], [100.0, 60.0, 40.0]);
        // A slot right across the top splits the top face in two.
        let slot = labelled_box(k, 2, [0.0, 0.0, 20.0], [20.0, 80.0, 20.0]);
        let slotted = block.difference(slot).unwrap();
        // A boss on one half, placed by a frame on the top.
        let frame = Frame::new([-30.0, 0.0, 20.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
        let boss = labelled_cylinder(k, 3, 6.0, 5.0, &frame);
        let part = slotted.union(boss).unwrap();
        // A copy keeps everything.
        let copy = part.try_clone().unwrap();
        drop(part);
        let mesh = copy.tessellate(&TessellationParams::default()).unwrap();
        (measure(&copy), copy.edges().unwrap(), mesh)
    });
    let counts = labels(&faces);
    assert_eq!(counts.get("-"), None, "every face is labelled: {counts:?}");
    let tops = find(&faces, BLOCK, "+Z");
    assert_eq!(tops.len(), 2, "the slot splits the top: {counts:?}");
    let mut xs: Vec<f64> = tops.iter().map(|f| f.centroid[0].round()).collect();
    xs.sort_by(f64::total_cmp);
    assert!(xs[0] < -10.0 && xs[1] > 10.0, "{xs:?}");
    // The slot's floor and walls, and the boss's side and top. The boss's
    // bottom and the slot's top and ends are gone.
    for (feature, role) in [(2, "-Z"), (2, "+X"), (2, "-X"), (3, "side"), (3, "top")] {
        assert_eq!(find(&faces, feature, role).len(), 1, "{feature}/{role}");
    }
    assert!(find(&faces, 3, "bottom").is_empty());
    assert!(find(&faces, 2, "+Z").is_empty());
    // The mesh carries the same labels in the same order.
    let mesh_labels: Vec<_> = mesh.faces.iter().map(|f| f.label.clone()).collect();
    let brep_labels: Vec<_> = faces.iter().map(|f| f.label.clone()).collect();
    assert_eq!(mesh_labels, brep_labels);
    // The block's vertical edges keep their labels; edges the booleans made
    // have none.
    let vertical = edges
        .iter()
        .filter(|e| {
            e.label
                .as_ref()
                .is_some_and(|l| l.role.contains("X") && l.role.contains("Y"))
        })
        .count();
    assert_eq!(vertical, 4);
    assert_eq!(
        mesh.edges
            .iter()
            .map(|e| e.label.clone())
            .collect::<Vec<_>>(),
        edges.iter().map(|e| e.label.clone()).collect::<Vec<_>>()
    );
}

#[test]
fn fillets_and_chamfers_keep_labels_and_add_unlabelled_faces() {
    let (filleted, chamfered) = on_kernel(|k| {
        let block = labelled_box(k, BLOCK, [0.0; 3], [100.0, 60.0, 40.0]);
        let frame = Frame::new([20.0, 0.0, 20.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
        let hole = labelled_cylinder(
            k,
            2,
            10.0,
            40.0,
            &Frame {
                origin: sub(frame.origin, [0.0, 0.0, 40.0]),
                ..frame
            },
        );
        let part = block.difference(hole).unwrap();
        let vertical = EdgeQuery::ParallelTo {
            direction: [0.0, 0.0, 1.0],
        }
        .and(EdgeQuery::Sharp {
            min_angle_deg: 30.0,
        });
        let part = part
            .fillet(
                &[EdgeBlend::round(vertical, 5.0)],
                &FilletOptions::default(),
            )
            .unwrap();
        let filleted = measure(&part);
        let rim = EdgeQuery::Between {
            a: LabelMatch::role(2, "side"),
            b: LabelMatch::role(BLOCK, "+Z"),
        };
        let part = part.chamfer(&rim, 2.0, &Tolerance::default()).unwrap();
        (filleted, measure(&part))
    });
    // Six block faces and the hole's wall keep their labels; four new rounds.
    assert_eq!(filleted.len(), 11);
    assert_eq!(labels(&filleted).get("-"), Some(&4));
    for face in filleted.iter().filter(|f| f.label.is_none()) {
        assert!(matches!(face.surface, Surface::Cylinder { .. }));
    }
    assert_eq!(chamfered.len(), 12);
    assert_eq!(labels(&chamfered).get("-"), Some(&5), "one new bevel");
    assert_eq!(find(&chamfered, 2, "side").len(), 1);
    assert_eq!(find(&chamfered, BLOCK, "+Z").len(), 1);
}

#[test]
fn sweeps_report_their_caps_and_sides() {
    for distance in [15.0, -15.0] {
        let (sweep, source, faces, volume) = on_kernel(move |k| {
            let block = labelled_box(k, BLOCK, [0.0; 3], [100.0, 60.0, 40.0]);
            let top = find(&measure(&block), BLOCK, "+Z")[0].index;
            let source = block.edges().unwrap();
            let sweep = block.sweep_face(top, distance).unwrap();
            let roles = sweep.faces.clone();
            let prism_faces = measure(&sweep.prism);
            let prism_volume = sweep.prism.volume(1e-9).unwrap();
            let mut prism = sweep.prism;
            for (index, role) in roles.iter().enumerate() {
                let label = match role {
                    SweptFace::Start => Label::new(9, "start"),
                    SweptFace::End => Label::new(9, "end"),
                    SweptFace::Side { .. } => Label::new(9, "side").with_index(index as u32),
                };
                prism.set_face_label(index, &label).unwrap();
            }
            let part = if distance > 0.0 {
                block.union(prism)
            } else {
                block.difference(prism)
            };
            let mut part = part.unwrap();
            part.merge_faces().unwrap();
            assert!(part.is_manifold().unwrap());
            let volume = part.volume(1e-9).unwrap();
            (
                (roles, prism_faces, prism_volume),
                source,
                measure(&part),
                volume,
            )
        });
        let (roles, prism_faces, prism_volume) = sweep;
        assert!(
            (prism_volume - 6000.0 * 15.0).abs() < 1e-6,
            "{prism_volume}"
        );
        assert_eq!(roles.iter().filter(|r| **r == SweptFace::Start).count(), 1);
        assert_eq!(roles.iter().filter(|r| **r == SweptFace::End).count(), 1);
        let start = roles.iter().position(|r| *r == SweptFace::Start).unwrap();
        let end = roles.iter().position(|r| *r == SweptFace::End).unwrap();
        assert_near("start", prism_faces[start].centroid, [0.0, 0.0, 20.0], 1e-6);
        assert_near(
            "end",
            prism_faces[end].centroid,
            [0.0, 0.0, 20.0 + distance],
            1e-6,
        );
        // Each side was swept from one of the top's four edges, and lies
        // above or below it.
        let mut sources = Vec::new();
        for (index, role) in roles.iter().enumerate() {
            let SweptFace::Side { source_edge } = role else {
                continue;
            };
            let edge = &source[source_edge.expect("every side has a source")];
            assert!(edge.face_labels.iter().flatten().any(|l| l.role == "+Z"));
            let c = prism_faces[index].centroid;
            let mid = edge.midpoint();
            assert_near(
                "side over its edge",
                [c[0], c[1], 0.0],
                [mid[0], mid[1], 0.0],
                1e-6,
            );
            sources.push(edge.index);
        }
        sources.sort();
        sources.dedup();
        assert_eq!(sources.len(), 4);

        // The walls grow or shrink as single faces and keep the block's
        // labels: the older label wins wherever pieces merge. The end cap is
        // the new top.
        assert!(
            (volume - (240_000.0 + 6000.0 * distance)).abs() < 1e-3,
            "{volume}"
        );
        assert_eq!(faces.len(), 6, "{:?}", labels(&faces));
        for role in ["+X", "-X", "+Y", "-Y", "-Z"] {
            assert_eq!(find(&faces, BLOCK, role).len(), 1, "{role}");
        }
        assert_eq!(find(&faces, 9, "end").len(), 1);
        let wall = find(&faces, BLOCK, "+X")[0];
        assert!((wall.area - 60.0 * (40.0 + distance)).abs() < 1e-6);
    }
}

#[test]
fn offsetting_a_face_between_rounds_keeps_every_face_single() {
    let faces = on_kernel(|k| {
        let block = labelled_box(k, BLOCK, [0.0; 3], [400.0, 300.0, 250.0]);
        let vertical = EdgeQuery::ParallelTo {
            direction: [0.0, 0.0, 1.0],
        };
        let mut part = block
            .fillet(
                &[EdgeBlend::round(vertical, 20.0)],
                &FilletOptions::default(),
            )
            .unwrap();
        for face in part.faces().unwrap() {
            if face.label.is_none() {
                part.set_face_label(
                    face.index,
                    &Label::new(2, "round").with_index(face.index as u32),
                )
                .unwrap();
            }
        }
        let top = find(&measure(&part), BLOCK, "+Z")[0].index;
        let sweep = part.sweep_face(top, 15.0).unwrap();
        let mut prism = sweep.prism;
        assert_eq!(sweep.faces.len(), 10, "two caps and eight sides");
        for (index, role) in sweep.faces.iter().enumerate() {
            // Sides take a newer feature's label, so merging leaves the
            // neighbours' labels in place.
            let label = match role {
                SweptFace::End => Label::new(BLOCK, "+Z"),
                _ => Label::new(3, "side").with_index(index as u32),
            };
            prism.set_face_label(index, &label).unwrap();
        }
        let mut part = part.union(prism).unwrap();
        part.merge_faces().unwrap();
        measure(&part)
    });
    assert_eq!(faces.len(), 10, "{:?}", labels(&faces));
    assert!(find(&faces, 3, "side").is_empty(), "every side merged");
    assert_eq!(find(&faces, 2, "round").len(), 4);
    let top = find(&faces, BLOCK, "+Z");
    assert_eq!(top.len(), 1);
    assert!((top[0].centroid[2] - 140.0).abs() < 1e-6);
    let wall = find(&faces, BLOCK, "+Y")[0];
    assert!(
        (wall.area - (400.0 - 40.0) * 265.0).abs() < 1e-6,
        "{}",
        wall.area
    );
}

#[test]
fn merged_faces_keep_the_older_label_whichever_survives() {
    // Two boxes side by side with flush tops, in both orders: the merged top,
    // bottom and sides take the lower feature id. The ends were never merged.
    for older_first in [true, false] {
        let faces = on_kernel(move |k| {
            let a = labelled_box(k, 1, [-20.0, 0.0, 0.0], [40.0, 30.0, 20.0]);
            let b = labelled_box(k, 7, [20.0, 0.0, 0.0], [40.0, 30.0, 20.0]);
            let mut part = if older_first { a.union(b) } else { b.union(a) }.unwrap();
            part.merge_faces().unwrap();
            measure(&part)
        });
        let counts = labels(&faces);
        let expected = ["1/+Y.0", "1/+Z.0", "1/-X.0", "1/-Y.0", "1/-Z.0", "7/+X.0"];
        assert_eq!(
            counts.keys().map(String::as_str).collect::<Vec<_>>(),
            expected,
            "older_first {older_first}"
        );
    }
}

#[test]
fn merged_faces_of_one_feature_keep_the_lower_index() {
    // As above, but both boxes are one feature's, told apart by index: the
    // merged faces take the lower index in either order.
    for lower_first in [true, false] {
        let faces = on_kernel(move |k| {
            let label = |solid: &mut Brep, index: u32| {
                for face in solid.faces().unwrap() {
                    let role = axis_role(face.normal.unwrap());
                    let label = Label::new(BLOCK, role).with_index(index);
                    solid.set_face_label(face.index, &label).unwrap();
                }
            };
            let mut a = Brep::cuboid(k, [-40.0, -15.0, -10.0], [40.0, 30.0, 20.0]).unwrap();
            let mut b = Brep::cuboid(k, [0.0, -15.0, -10.0], [40.0, 30.0, 20.0]).unwrap();
            label(&mut a, 3);
            label(&mut b, 8);
            let mut part = if lower_first { a.union(b) } else { b.union(a) }.unwrap();
            part.merge_faces().unwrap();
            measure(&part)
        });
        let counts = labels(&faces);
        let expected = ["1/+X.8", "1/+Y.3", "1/+Z.3", "1/-X.3", "1/-Y.3", "1/-Z.3"];
        assert_eq!(
            counts.keys().map(String::as_str).collect::<Vec<_>>(),
            expected,
            "lower_first {lower_first}"
        );
    }
}

#[test]
fn edge_queries_pick_by_direction_label_and_index() {
    let (picked, edges) = on_kernel(|k| {
        let block = labelled_box(k, BLOCK, [0.0; 3], [100.0, 60.0, 40.0]);
        let frame = Frame::new([0.0, 0.0, 20.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
        let part = block
            .union(labelled_cylinder(k, 4, 5.0, 10.0, &frame))
            .unwrap();
        let tolerance = Tolerance::default();
        let count = |query: EdgeQuery| part.count_edges(&query, &tolerance).unwrap();
        let picked = [
            count(EdgeQuery::ParallelTo {
                direction: [0.0, 0.0, -2.0],
            }),
            count(EdgeQuery::ParallelTo {
                direction: [1.0, 0.0, 0.0],
            }),
            count(EdgeQuery::Between {
                a: LabelMatch::role(4, "side"),
                b: LabelMatch::role(BLOCK, "+Z"),
            }),
            count(EdgeQuery::Between {
                a: LabelMatch {
                    feature: BLOCK,
                    role: None,
                    index: None,
                },
                b: LabelMatch {
                    feature: BLOCK,
                    role: None,
                    index: None,
                },
            }),
            count(EdgeQuery::Indices(vec![0, 3])),
        ];
        (picked, part.edges().unwrap())
    });
    // Four vertical box edges (the boss's seam is vertical too, but its
    // faces are the same curved face), four along x, the boss's rim, the
    // box's twelve edges, two by index.
    let seams = edges.iter().filter(|e| e.faces[0] == e.faces[1]).count();
    assert_eq!(picked, [4 + seams, 4, 1, 12, 2]);
}

#[test]
fn offsets_keep_every_label_in_place() {
    // Offsets change no topology, so every face and edge keeps its label, on
    // the face or edge it was on: the moved faces, the rounds that ride on
    // them and the walls that stretch.
    type Edges = Vec<(Option<Label>, [Option<Label>; 2])>;
    let edges_of = |solid: &Brep| -> Edges {
        let mut edges: Edges = solid
            .edges()
            .unwrap()
            .into_iter()
            .map(|e| {
                let mut faces = e.face_labels.clone();
                faces.sort();
                (e.label.clone(), faces)
            })
            .collect();
        edges.sort();
        edges
    };
    let (before, pushed, moved, edges_before, edges_pushed, edges_moved) = on_kernel(move |k| {
        let block = labelled_box(k, BLOCK, [0.0; 3], [100.0, 60.0, 40.0]);
        let hole = labelled_cylinder(
            k,
            2,
            10.0,
            40.0,
            &Frame::new([0.0, 0.0, -20.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap(),
        );
        let part = block.difference(hole).unwrap();
        let top_edges = EdgeQuery::MidpointOnPlane {
            point: [0.0, 0.0, 20.0],
            normal: [0.0, 0.0, 1.0],
        }
        .and(!EdgeQuery::Between {
            a: LabelMatch::role(2, "side"),
            b: LabelMatch::role(BLOCK, "+Z"),
        });
        let mut part = part
            .fillet(
                &[EdgeBlend::round(top_edges, 8.0)],
                &FilletOptions::default(),
            )
            .unwrap();
        for face in part.faces().unwrap() {
            if face.label.is_none() {
                part.set_face_label(
                    face.index,
                    &Label::new(3, "round").with_index(face.index as u32),
                )
                .unwrap();
            }
        }
        for edge in part.edges().unwrap() {
            if edge.label.is_none() {
                part.set_edge_label(
                    edge.index,
                    &Label::new(3, "rail").with_index(edge.index as u32),
                )
                .unwrap();
            }
        }
        let before = measure(&part);
        let edges_before = edges_of(&part);
        let wall = find(&before, BLOCK, "+X")[0].index;
        let part = part.offset_faces(&[wall], 10.0).unwrap();
        let pushed = measure(&part);
        let edges_pushed = edges_of(&part);
        let top = find(&pushed, BLOCK, "+Z")[0].index;
        let hole_wall = find(&pushed, 2, "side")[0].index;
        let part = part.offset_faces(&[top, hole_wall], 2.0).unwrap();
        let moved = measure(&part);
        let edges_moved = edges_of(&part);
        (
            before,
            pushed,
            moved,
            edges_before,
            edges_pushed,
            edges_moved,
        )
    });
    assert_eq!(before.len(), 11, "{:?}", labels(&before));
    assert!(!labels(&before).contains_key("-"), "every face labelled");
    assert_eq!(labels(&pushed), labels(&before));
    assert_eq!(labels(&moved), labels(&before));
    assert_eq!(
        edges_pushed, edges_before,
        "edges keep their labels and faces"
    );
    assert_eq!(
        edges_moved, edges_before,
        "edges keep their labels and faces"
    );

    // The labelled faces are where the moves put them.
    let wall = find(&pushed, BLOCK, "+X")[0];
    assert!(
        (wall.centroid[0] - 60.0).abs() < 1e-4,
        "{:?}",
        wall.centroid
    );
    let far = find(&pushed, BLOCK, "+Y")[0];
    assert!(
        far.area > find(&before, BLOCK, "+Y")[0].area + 10.0 * 32.0 - 1e-6,
        "the +Y wall stretched"
    );
    let top = find(&moved, BLOCK, "+Z")[0];
    assert!((top.centroid[2] - 22.0).abs() < 1e-4, "{:?}", top.centroid);
    match find(&moved, 2, "side")[0].surface {
        Surface::Cylinder { radius, .. } => assert!((radius - 8.0).abs() < 1e-9, "{radius}"),
        ref other => panic!("the hole's wall is {other:?}"),
    }
    // Each round rode with what it rounds. The push carried the +X one out
    // 10 and stretched the two along it by 10; the top's move lifted them
    // all by 2. (Centroids come out good to 1e-6 here, areas to a relative
    // 2e-6.)
    let close = |a: [f64; 3], b: [f64; 3]| (0..3).all(|i| (a[i] - b[i]).abs() < 1e-4);
    let same_area = |a: f64, b: f64| (a - b).abs() < 1e-5 * b;
    let rounds = before
        .iter()
        .filter(|f| f.label.as_ref().is_some_and(|l| l.role == "round"));
    for old in rounds {
        let same = |faces: &[Face]| faces.iter().find(|f| f.label == old.label).unwrap().clone();
        let (out, up) = (same(&pushed), same(&moved));
        let c = old.centroid;
        match old.surface {
            Surface::Cylinder { axis, .. } if axis[0].abs() > 0.5 => {
                assert!(
                    same_area(out.area, old.area + 10.0 * 4.0 * PI),
                    "{:?} stretched",
                    old.label
                );
            }
            _ if c[0] > 0.0 => assert!(
                close(out.centroid, [c[0] + 10.0, c[1], c[2]]),
                "{:?}",
                old.label
            ),
            _ => assert!(close(out.centroid, c), "{:?}", old.label),
        }
        let c = out.centroid;
        assert!(
            close(up.centroid, [c[0], c[1], c[2] + 2.0]),
            "{:?}: {c:?} -> {:?}",
            old.label,
            up.centroid
        );
        assert!(
            same_area(up.area, out.area),
            "{:?}: {} -> {}",
            old.label,
            out.area,
            up.area
        );
    }
}
