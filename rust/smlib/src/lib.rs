//! Safe Rust API over NVIDIA's SMLib B-rep solid modelling kernel (USD BRep).
//!
//! SMLib keeps process-global state and is not thread-safe, so this crate
//! confines it to one thread: [`Kernel::init`] binds the kernel to the first
//! thread that calls it, and every [`Brep`] is created through that thread's
//! [`Kernel`] handle. Neither type is `Send` or `Sync`, so no kernel object can
//! leave that thread. Tessellation output ([`Mesh`]) is plain data and can.
//!
//! The kernel is unit-agnostic; lengths are in whatever unit the model uses.
//! The [`TessellationParams`] defaults assume millimetres.
//!
//! ```no_run
//! use smlib::{Brep, Kernel, TessellationParams};
//!
//! let kernel = Kernel::init()?;
//! let block = Brep::cuboid(&kernel, [-50.0, -30.0, -20.0], [100.0, 60.0, 40.0])?;
//! let pin = Brep::cylinder(&kernel, [0.0, 0.0, -25.0], 10.0, 50.0)?;
//! let part = block.difference(pin)?;
//! assert!(part.is_manifold()?);
//! let mesh = part.tessellate(&TessellationParams::default())?;
//! println!("{} triangles over {} faces", mesh.triangle_count(), mesh.faces.len());
//! # Ok::<(), smlib::Error>(())
//! ```

#![warn(missing_docs)]
#![deny(unsafe_op_in_unsafe_fn)]

mod brep;
mod error;
mod kernel;
mod mesh;

pub use brep::Brep;
pub use error::{Error, Status};
pub use kernel::Kernel;
pub use mesh::{Edge, Face, Mesh, Surface, TessellationParams};
