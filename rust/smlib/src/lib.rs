//! Safe Rust API over NVIDIA's SMLib B-rep solid modelling kernel (USD BRep).
//!
//! SMLib keeps process-global state and is not thread-safe, so this crate
//! confines it to one thread: [`Kernel::init`] binds the kernel to the first
//! thread that calls it, and every [`Brep`] is created through that thread's
//! [`Kernel`] handle. Neither type is `Send` or `Sync`, so no kernel object can
//! leave that thread. Tessellation output ([`Mesh`]) is plain data and can.
//!
//! Fillets and chamfers ([`Brep::fillet`], [`Brep::chamfer`]) pick their
//! edges with an [`EdgeQuery`]. Solids built at the origin can be placed on a
//! face's [`Frame`] ([`Brep::place`]), any solid moved, turned, scaled and
//! mirrored ([`Brep::transform`], [`Brep::scale_about`], [`Brep::mirror`]),
//! and planar faces swept into prisms ([`Brep::sweep_face`], with draft
//! [`Brep::sweep_face_tapered`]) or revolved ([`Brep::revolve_face`]),
//! whether a solid's or ones built from lines, arcs and circles in a plane
//! ([`Brep::planar_faces`]), and circles swept along paths into pipes
//! ([`Brep::pipe`]), or any closed profile ([`Brep::sweep_profile`]).
//! [`Brep::shell`] hollows a solid, leaving chosen faces
//! open. [`Brep::offset_faces`]
//! moves faces along their normals, extending or trimming their neighbours,
//! as Parasolid's face offset does, and [`Brep::move_faces`] moves or turns
//! them rigidly the same way. Faces and edges can carry a [`Label`] that the kernel keeps
//! through every later operation.
//!
//! The kernel is unit-agnostic; lengths are in whatever unit the model uses.
//! The [`TessellationParams`] defaults assume millimetres.
//!
//! ```no_run
//! use smlib::{Brep, EdgeBlend, EdgeQuery, FilletOptions, Kernel, TessellationParams};
//!
//! let kernel = Kernel::init()?;
//! let block = Brep::cuboid(&kernel, [-50.0, -30.0, -20.0], [100.0, 60.0, 40.0])?;
//! let pin = Brep::cylinder(&kernel, [0.0, 0.0, -25.0], 10.0, 50.0)?;
//! let part = block.difference(pin)?;
//! assert!(part.is_manifold()?);
//! // Round every sharp edge of the part, including the hole's rims.
//! let sharp = EdgeQuery::Sharp { min_angle_deg: 30.0 };
//! let part = part.fillet(&[EdgeBlend::round(sharp, 2.0)], &FilletOptions::default())?;
//! let mesh = part.tessellate(&TessellationParams::default())?;
//! println!("{} triangles over {} faces", mesh.triangle_count(), mesh.faces.len());
//! # Ok::<(), smlib::Error>(())
//! ```

#![warn(missing_docs)]
#![deny(unsafe_op_in_unsafe_fn)]

mod brep;
mod edges;
mod error;
mod faces;
mod fillet;
mod imprint;
mod kernel;
mod label;
mod mesh;
mod motion;
mod pipe;
mod shell;
mod tessellation;

pub use brep::Brep;
pub use edges::{EDGE_SAMPLES, EdgeGeometry, EdgeQuery, Tolerance};
pub use error::{Error, Status};
pub use faces::{
    BoundaryCurve, FaceCentroid, FaceGeometry, FaceQuery, Frame, PlaneCurve, Relationship,
    ScaleReport, SurfaceKind, Sweep, SweptFace,
};
pub use fillet::{Continuity, CrossSection, EdgeBlend, FilletOptions};
pub use imprint::{Imprint, ImprintEdge, UvCurve};
pub use kernel::Kernel;
pub use label::{Label, LabelMatch, MAX_ROLE_LEN};
pub use mesh::{Edge, Face, Mesh, Surface, TessellationParams};
pub use motion::{Moved, MovedFace, RigidMotion, Turn};
pub use pipe::{PathEnd, Pipe, PipeFace, ProfileFace, ProfileSweep};
pub use shell::{Shell, ShellFace};
pub use tessellation::{FaceKey, SurfaceCaches, Tessellation};
