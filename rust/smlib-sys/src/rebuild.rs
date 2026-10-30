//! Rebuilding faces (`shim/rebuild.h`), written out by hand rather than by
//! bindgen, so `shim.h` and `bindings.rs` stay as they are. Keep it in step
//! with `rebuild.h`.

use crate::{sm_brep, sm_label, sm_status, sm_surface};

pub const SM_SHIM_ERR_NOT_RIM: sm_status = 2101;
pub const SM_SHIM_ERR_OPEN_LOOP: sm_status = 2102;
pub const SM_SHIM_ERR_BAD_TOOL: sm_status = 2103;
pub const SM_SHIM_ERR_SURFACE_FAILED: sm_status = 2104;
pub const SM_SHIM_ERR_LEFT_OPEN: sm_status = 2105;
pub const SM_SHIM_ERR_NOT_SOLID: sm_status = 2106;

pub const SM_REBUILD_FILL: i32 = 0;
pub const SM_REBUILD_COONS: i32 = 1;
pub const SM_REBUILD_CORNER: i32 = 2;
pub const SM_REBUILD_EDGE_BLEND: i32 = 3;
pub const SM_REBUILD_SKIN: i32 = 4;

pub const SM_MEETS_POSITION: i32 = 0;
pub const SM_MEETS_TANGENT: i32 = 1;
pub const SM_MEETS_CURVATURE: i32 = 2;

pub const SM_EXTEND_LINEAR: i32 = 0;
pub const SM_EXTEND_SMOOTH: i32 = 1;
pub const SM_EXTEND_NATURAL: i32 = 2;
pub const SM_EXTEND_RULED: i32 = 3;

pub const SM_SIDE_UMIN: i32 = 0;
pub const SM_SIDE_VMIN: i32 = 1;
pub const SM_SIDE_UMAX: i32 = 2;
pub const SM_SIDE_VMAX: i32 = 3;
pub const SM_SIDE_ALL: i32 = 4;

pub const SM_EDGE_CURVATURE_SAMPLES: usize = 12;

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct sm_rebuild_side {
    pub edge: u32,
    pub side: u32,
    pub meets: i32,
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct sm_rebuild_surface {
    pub tool: i32,
    pub sides: *const sm_rebuild_side,
    pub side_count: usize,
    pub order: i32,
    pub tangent_tol_deg: f64,
    pub fit_tolerance: f64,
    pub label: sm_label,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct sm_rebuild_report {
    pub solid: i32,
    pub stitched: u32,
    pub open_edges: u32,
    pub vertex_gap: f64,
    pub edge_gap: f64,
    pub open_vertex_gap: f64,
    pub open_edge_gap: f64,
    pub faces_made: u32,
    pub failed_surface: i64,
    pub surface_status: sm_status,
    pub surface_seconds: f64,
    pub stitch_seconds: f64,
    pub open_gap: f64,
    pub open_at: i64,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct sm_edge_curvature {
    pub valid: i32,
    pub a: [f64; SM_EDGE_CURVATURE_SAMPLES],
    pub b: [f64; SM_EDGE_CURVATURE_SAMPLES],
}

unsafe extern "C" {
    pub fn sm_brep_rebuild_faces(
        brep: *const sm_brep,
        faces: *const u32,
        face_count: usize,
        surfaces: *const sm_rebuild_surface,
        surface_count: usize,
        stitch_tolerance: f64,
        out: *mut *mut sm_brep,
        report: *mut sm_rebuild_report,
    ) -> sm_status;

    pub fn sm_brep_edge_curvatures(
        brep: *const sm_brep,
        out: *mut sm_edge_curvature,
        capacity: usize,
        count: *mut usize,
    ) -> sm_status;

    pub fn sm_surface_extend(
        surface: *const sm_surface,
        side: i32,
        distance: f64,
        shape: i32,
        out: *mut *mut sm_surface,
    ) -> sm_status;
}
