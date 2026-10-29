// Rebuilding faces for the smlib-sys shim: faces taken out of a solid, the
// gaps rebuilt from the rim's edges and sewn back in one call. See rebuild.h.

#include "rebuild.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <new>
#include <utility>
#include <vector>

#include <SmApiGeneral.h>

#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmBrep.h>
#include <SmCurve.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmFace.h>
#include <SmFaceuse.h>
#include <SmObject.h>
#include <SmStitch.h>
#include <SmSurface.h>
#include <SmTangentField.h>

namespace {

SmBrep* brep_of(sm_brep* brep) { return reinterpret_cast<SmBrep*>(brep); }
const SmBrep* brep_of(const sm_brep* brep) { return reinterpret_cast<const SmBrep*>(brep); }
const SmSurface* surface_of(const sm_surface* surface) { return reinterpret_cast<const SmSurface*>(surface); }
sm_surface* surface_handle(SmSurface* surface) { return reinterpret_cast<sm_surface*>(surface); }
sm_curve* curve_handle(SmCurve* curve) { return reinterpret_cast<sm_curve*>(curve); }
SmSurface* surface_from(sm_surface* surface) { return reinterpret_cast<SmSurface*>(surface); }

template <typename F>
sm_status guarded(F&& f) {
    try {
        return f();
    } catch (const std::bad_alloc&) {
        return SM_SHIM_ERR_OUT_OF_MEMORY;
    } catch (...) {
        return SM_SHIM_ERR_EXCEPTION;
    }
}

const SmContext& context() { return *SmApiGetOrCreateContext(); }

sm_status status_of(SmStatus status) { return static_cast<sm_status>(status); }

double seconds_since(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

// Deletes what it holds unless released.
template <typename T>
struct Owned {
    std::vector<T*> items;
    ~Owned() {
        for (T* item : items) delete item;
    }
    T* add(T* item) {
        items.push_back(item);
        return item;
    }
    void release() { items.clear(); }
};

// A rim edge: the edge, and its edgeuse on the kept face, running with the
// edge's curve.
struct Rim {
    SmEdge* edge = nullptr;
    SmEdgeuse* use = nullptr;
    SmFace* face = nullptr;
    int32_t meets = SM_MEETS_POSITION;
};

sm_status rim_of(const SmTArray<SmEdge*>& edges, const std::vector<SmFace*>& removed, const sm_rebuild_side& side,
                 Rim& rim) {
    if (side.edge >= edges.GetSize()) return SM_SHIM_ERR_BAD_EDGE;
    if (side.meets < SM_MEETS_POSITION || side.meets > SM_MEETS_CURVATURE) return SM_SHIM_ERR_BAD_TOOL;
    SmEdge* edge = edges[side.edge];
    SmTArray<SmEdgeuse*> uses;
    edge->GetEdgeuses(uses);
    SmEdgeuse* kept = nullptr;
    bool borders_removed = false;
    for (ULONG i = 0; i < uses.GetSize(); ++i) {
        SmFace* face = uses[i]->GetFace();
        if (!face) continue;
        if (std::find(removed.begin(), removed.end(), face) != removed.end()) {
            borders_removed = true;
        } else if (!kept || (kept->GetOrientation() == SM_OT_OPPOSITE && uses[i]->GetOrientation() == SM_OT_SAME)) {
            kept = uses[i];
        }
    }
    if (!kept || !borders_removed) return SM_SHIM_ERR_NOT_RIM;
    // As SmBrep::CreateCornerBlend does: the edgeuse of the same face that
    // runs with the edge's curve.
    if (kept->GetOrientation() == SM_OT_OPPOSITE && kept->GetMate()) kept = kept->GetMate();
    rim.edge = edge;
    rim.use = kept;
    rim.face = kept->GetFace();
    rim.meets = side.meets;
    return SM_SUCCESS;
}

SmPoint3d point_at(const SmCurve& curve, double t) {
    SmPoint3d p;
    curve.EvaluatePoint(t, p);
    return p;
}

SmVector3d tangent_at(const SmCurve& curve, double t) {
    SmVector3d d[2];
    curve.Evaluate(t, 1, TRUE, d);
    return d[1];
}

// A side of a loop as a fill or Coons patch takes it: the rim edge's curve,
// non-rational and on [0, 1], and its trim curve on the kept face, both
// running in loop order.
struct Bound {
    Rim rim;
    SmBSplineCurve* curve = nullptr;
    SmBSplineCurve* uv = nullptr;
};

// The edge's curve trimmed to it, non-rational and on [0, 1], with its
// trim curve on the kept face to match: the edgeuse's own, or for a
// rational curve (an arc), approximated as sm_CreateNonRationalCurve does
// for CreateCornerBlend, one dropped onto the face's surface.
sm_status bound_of(const Rim& rim, double tolerance, Owned<SmBSplineCurve>& held, Bound& bound) {
    bound.rim = rim;
    SmBSplineCurve* nurbs = rim.edge->CreateTrimmedNURBSCurve(context());
    if (!nurbs) return SM_ERR;
    held.add(nurbs);
    SmBSplineCurve* curve = nurbs;
    SmBSplineCurve* uv_curve = nullptr;
    if (nurbs->IsRational()) {
        SmBSplineCurve* fitted = nullptr;
        SmStatus status = nurbs->ApproximateRationalCurve(context(), tolerance, fitted);
        if (status != SM_SUCCESS || !fitted) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        held.add(fitted);
        curve = fitted;
        const SmSurface& surface = *rim.face->GetSurface();
        double dropped = 0.0;
        double deviation = 0.0;
        SmTArray<SmBSplineCurve*> uvs;
        status = surface.DropCurve(context(), surface.GetNaturalUVDomain(), *fitted, fitted->GetNaturalInterval(),
                                   tolerance, dropped, deviation, uvs);
        for (ULONG i = 0; i < uvs.GetSize(); ++i) held.add(uvs[i]);
        if (status != SM_SUCCESS || uvs.GetSize() != 1) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        uv_curve = uvs[0];
    } else {
        SmBSplineCurve* uv = nullptr;
        rim.use->GetOrCreateUVTrimCurve(uv);
        if (!uv) return SM_ERR;
        SmCurve* uv_copy = nullptr;
        if (uv->Copy(context(), uv_copy) != SM_SUCCESS || !uv_copy) return SM_ERR;
        uv_curve = dynamic_cast<SmBSplineCurve*>(uv_copy);
        if (!uv_curve) {
            delete uv_copy;
            return SM_ERR;
        }
        held.add(uv_curve);
    }
    const SmExtent1d unit(0.0, 1.0);
    curve->EditParameterization(unit);
    uv_curve->EditParameterization(unit);
    bound.curve = curve;
    bound.uv = uv_curve;
    return SM_SUCCESS;
}

void reverse(Bound& bound) {
    SmExtent1d unit(0.0, 1.0);
    bound.curve->ReverseParameterization(unit, unit);
    bound.uv->ReverseParameterization(unit, unit);
}

// Turns the bounds to run end to end round their loop, in the order given:
// each one's end on the next one's start, the last's on the first's.
sm_status chain(std::vector<Bound>& bounds, double tolerance) {
    const size_t n = bounds.size();
    if (n < 2) return SM_SHIM_ERR_BAD_TOOL;
    {
        const SmPoint3d end = point_at(*bounds[0].curve, 1.0);
        const SmPoint3d start = point_at(*bounds[0].curve, 0.0);
        const SmPoint3d a = point_at(*bounds[1].curve, 0.0);
        const SmPoint3d b = point_at(*bounds[1].curve, 1.0);
        const double forward = std::min(end.DistanceBetween(a), end.DistanceBetween(b));
        const double backward = std::min(start.DistanceBetween(a), start.DistanceBetween(b));
        if (backward < forward) reverse(bounds[0]);
    }
    for (size_t i = 1; i < n; ++i) {
        const SmPoint3d last = point_at(*bounds[i - 1].curve, 1.0);
        const double as_is = last.DistanceBetween(point_at(*bounds[i].curve, 0.0));
        const double turned = last.DistanceBetween(point_at(*bounds[i].curve, 1.0));
        if (turned < as_is) reverse(bounds[i]);
        if (std::min(as_is, turned) > tolerance) return SM_SHIM_ERR_OPEN_LOOP;
    }
    if (point_at(*bounds[n - 1].curve, 1.0).DistanceBetween(point_at(*bounds[0].curve, 0.0)) > tolerance) {
        return SM_SHIM_ERR_OPEN_LOOP;
    }
    return SM_SUCCESS;
}

// One side of a run of bounds on the same kept face, meeting it the same
// way, as the ends of the corner part's rounds are each two edges: their
// curves and trim curves joined end to end, on [0, 1]. The bounds run end
// to end round their loop already (`chain`).
sm_status join_runs(std::vector<Bound>& bounds, Owned<SmBSplineCurve>& held) {
    auto joined = [&](SmBSplineCurve* a, SmBSplineCurve* b, SmBSplineCurve*& out) -> sm_status {
        SmTArray<SmBSplineCurve*> pair;
        pair.Add(a);
        pair.Add(b);
        out = nullptr;
        const SmStatus status = SmBSplineCurve::CreateByJoining(context(), pair, nullptr, out);
        if (status != SM_SUCCESS || !out) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        held.add(out);
        out->EditParameterization(SmExtent1d(0.0, 1.0));
        return SM_SUCCESS;
    };
    bool merged = true;
    while (merged && bounds.size() > 2) {
        merged = false;
        for (size_t i = 0; i < bounds.size() && bounds.size() > 2; ++i) {
            const size_t next = (i + 1) % bounds.size();
            if (bounds[i].rim.face != bounds[next].rim.face || bounds[i].rim.meets != bounds[next].rim.meets) continue;
            Bound bound = bounds[i];
            sm_status status = joined(bounds[i].curve, bounds[next].curve, bound.curve);
            if (status == SM_SUCCESS) status = joined(bounds[i].uv, bounds[next].uv, bound.uv);
            if (status != SM_SUCCESS) return status;
            bounds[i] = bound;
            bounds.erase(bounds.begin() + static_cast<std::ptrdiff_t>(next));
            if (next == 0) std::rotate(bounds.begin(), bounds.end() - 1, bounds.end());
            merged = true;
            break;
        }
    }
    return SM_SUCCESS;
}

// A cross-derivative curve along the bound, compatible with it, from the
// kept face's tangent plane: SmTangentField blends the take-off direction
// from `start` (the start curve's tangent at its parameter, turned round for
// SM_OT_OPPOSITE) to `end`, and it's fitted with matching parameters to
// tolerance.
sm_status field_of(const Bound& bound, const SmCurve& start, double start_t, SmOrientType start_orient,
                   const SmCurve& end, double end_t, SmOrientType end_orient, double tolerance, double tangent_tol,
                   Owned<SmBSplineCurve>& held, SmBSplineCurve*& out) {
    SmTangentField* field = nullptr;
    SmStatus status = SmTangentField::Create(context(), 2.0 * tolerance, *bound.curve, *bound.uv,
                                             *bound.rim.face->GetSurface(), start, start_t, start_orient, end, end_t,
                                             end_orient, field);
    if (status != SM_SUCCESS || !field) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
    std::unique_ptr<SmTangentField> owned(field);
    // The field's vectors are as long as the curves beside it are fast, so a
    // tolerance in radians is a length relative to them.
    const double scale =
        std::max(tangent_at(start, start_t).Length(), tangent_at(end, end_t).Length());
    SmTArray<double> knots;
    bound.curve->GetKnots(knots);
    double gap = 0.0;
    SmBSplineCurve* fitted = nullptr;
    status = field->ApproximateCurve(context(), knots, tangent_tol * std::max(scale, 1e-9), gap, fitted, FALSE, TRUE,
                                     FALSE);
    if (status != SM_SUCCESS || !fitted) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
    held.add(fitted);
    out = fitted;
    return SM_SUCCESS;
}

// A cubic cross-derivative curve from `from` to `to` on [0, 1], leaving
// each end at the rate given (the corner's twist), for a side with no face
// beyond it: with no twists, as SmBSplineSurface::CreateCornerBlend makes
// one.
sm_status blend_of(const SmVector3d& from, const SmVector3d& to, const SmVector3d& from_twist,
                   const SmVector3d& to_twist, Owned<SmBSplineCurve>& held, SmBSplineCurve*& out) {
    SmPoint3d data[4];
    SmTArray<SmPoint3d> poles(4, data);
    poles.Add(from);
    poles.Add(from + from_twist / 3.0);
    poles.Add(to - to_twist / 3.0);
    poles.Add(to);
    double knot_data[2] = {0.0, 1.0};
    SmTArray<double> knots(2, knot_data, 2);
    ULONG mult_data[2] = {4, 4};
    SmTArray<ULONG> mults(2, mult_data, 2);
    SmBSplineCurve* curve = nullptr;
    const SmStatus status = SmBSplineCurve::CreateCanonical(context(), 3, 3, poles, SM_CF_UNSPECIFIED, mults, knots,
                                                            SM_KT_UNSPECIFIED, nullptr, nullptr, curve);
    if (status != SM_SUCCESS || !curve) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
    held.add(curve);
    out = curve;
    return SM_SUCCESS;
}

// Fills the loop with SmBSplineSurface::CreateNSidedPatch, a Tangent side
// taking its cross-derivatives from the kept face.
sm_status make_fill(std::vector<Bound>& bounds, double tolerance, double tangent_tol_deg,
                    std::vector<SmSurface*>& made) {
    const size_t n = bounds.size();
    if (n < 3 || n > 8) return SM_SHIM_ERR_BAD_TOOL;
    Owned<SmBSplineCurve> held;
    std::vector<sm_curve*> boundary(n);
    std::vector<sm_curve*> derivatives(n, nullptr);
    bool any = false;
    for (size_t i = 0; i < n; ++i) {
        boundary[i] = curve_handle(bounds[i].curve);
        if (bounds[i].rim.meets == SM_MEETS_POSITION) continue;
        const Bound& before = bounds[(i + n - 1) % n];
        const Bound& after = bounds[(i + 1) % n];
        SmBSplineCurve* field = nullptr;
        // Into the hole at each end: along the side beside it, away from
        // the corner.
        const sm_status status = field_of(bounds[i], *before.curve, 1.0, SM_OT_OPPOSITE, *after.curve, 0.0,
                                          SM_OT_SAME, tolerance, tangent_tol_deg * SM_PI / 180.0, held, field);
        if (status != SM_SUCCESS) return status;
        // N_FillNSidedHole takes each cross-derivative compatible with its
        // boundary: one degree and one knot vector.
        SmTArray<SmBSplineCurve*> pair;
        pair.Add(bounds[i].curve);
        pair.Add(field);
        const SmStatus compatible = SmBSplineCurve::MakeCurvesCompatible(pair, 2.0e-8);
        if (compatible != SM_SUCCESS) return status_of(compatible);
        derivatives[i] = curve_handle(field);
        any = true;
    }
    for (size_t i = 0; i < n; ++i) {
    }
    std::vector<sm_surface*> out(n, nullptr);
    size_t count = 0;
    const sm_status status = sm_surface_fill(boundary.data(), any ? derivatives.data() : nullptr, n,
                                             tangent_tol_deg, out.data(), &count);
    if (status != SM_SUCCESS) return status;
    for (size_t i = 0; i < count; ++i) made.push_back(surface_from(out[i]));
    return SM_SUCCESS;
}

// A bicubic Coons patch through four sides in loop order, laid out as
// SmBSplineSurface::CreateCornerBlend lays its four out: left (u min),
// bottom (v min), right (u max) and top (v max), each running with its
// parameter, and each cross-derivative the patch's derivative across it.
sm_status make_coons(std::vector<Bound>& bounds, double tolerance, double tangent_tol_deg,
                     std::vector<SmSurface*>& made) {
    if (bounds.size() != 4) return SM_SHIM_ERR_BAD_TOOL;
    enum { kLeft = 0, kBottom = 1, kRight = 2, kTop = 3 };
    // In loop order the left and top run against their parameters.
    reverse(bounds[kLeft]);
    reverse(bounds[kTop]);
    // Opposite sides compatible, as CreateCornerBlend makes them: on [0, 1]
    // already, so their trim curves still match them.
    for (auto [a, b] : {std::pair<int, int>{kLeft, kRight}, std::pair<int, int>{kBottom, kTop}}) {
        SmTArray<SmBSplineCurve*> pair;
        pair.Add(bounds[a].curve);
        pair.Add(bounds[b].curve);
        const SmStatus compatible = SmBSplineCurve::MakeCurvesCompatible(pair, 2.0e-8);
        if (compatible != SM_SUCCESS) return status_of(compatible);
    }
    Owned<SmBSplineCurve> held;
    SmBSplineCurve* fields[4] = {nullptr, nullptr, nullptr, nullptr};
    // As CreateCornerBlend: left and right run between the bottom's and the
    // top's ends, bottom and top between the left's and right's.
    auto start_of = [](int i) { return i % 2 == 1 ? kLeft : kBottom; };
    auto end_of = [](int i) { return i % 2 == 1 ? kRight : kTop; };
    auto at_of = [](int i) { return i >= 2 ? 1.0 : 0.0; };
    // The Tangent sides first, so a Position side between them can take
    // the corners' twists from their fields: the patch's ∂²S/∂u∂v there,
    // the rate a field beside it changes at the corner.
    for (int i = 0; i < 4; ++i) {
        if (bounds[i].rim.meets == SM_MEETS_POSITION) continue;
        const sm_status status = field_of(bounds[i], *bounds[start_of(i)].curve, at_of(i), SM_OT_SAME,
                                          *bounds[end_of(i)].curve, at_of(i), SM_OT_SAME, tolerance,
                                          tangent_tol_deg * SM_PI / 180.0, held, fields[i]);
        if (status != SM_SUCCESS) return status;
    }
    for (int i = 0; i < 4; ++i) {
        if (bounds[i].rim.meets != SM_MEETS_POSITION) continue;
        const int start = start_of(i);
        const int end = end_of(i);
        const double at = at_of(i);
        const SmVector3d none(0.0, 0.0, 0.0);
        const SmVector3d from_twist = fields[start] ? tangent_at(*fields[start], at) : none;
        const SmVector3d to_twist = fields[end] ? tangent_at(*fields[end], at) : none;
        const sm_status status = blend_of(tangent_at(*bounds[start].curve, at), tangent_at(*bounds[end].curve, at),
                                          from_twist, to_twist, held, fields[i]);
        if (status != SM_SUCCESS) return status;
    }
    SmTArray<SmBSplineCurve*> curves;
    SmTArray<SmBSplineCurve*> tangents;
    for (int i : {kLeft, kRight, kBottom, kTop}) {
        curves.Add(bounds[i].curve);
        tangents.Add(fields[i]);
    }
    SmBSplineSurface* surface = nullptr;
    const SmStatus status = SmBSplineSurface::CreateCoonsPatch(context(), tolerance, FALSE, curves, tangents, surface);
    if (status != SM_SUCCESS || !surface) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
    made.push_back(surface);
    return SM_SUCCESS;
}

sm_status make_corner(const SmBrep& brep, const std::vector<Rim>& rims, double tolerance, double tangent_tol_deg,
                      std::vector<SmSurface*>& made) {
    if (rims.size() < 3 || rims.size() > 4) return SM_SHIM_ERR_BAD_TOOL;
    SmTArray<SmEdgeuse*> uses;
    SmTArray<SmBoolean> tangent;
    for (const Rim& rim : rims) {
        uses.Add(rim.use);
        tangent.Add(rim.meets != SM_MEETS_POSITION ? TRUE : FALSE);
    }
    SmTArray<SmSurface*> surfaces;
    const SmStatus status =
        brep.CreateCornerBlend(context(), uses, &tangent, tolerance, tangent_tol_deg * SM_PI / 180.0, surfaces);
    if (status != SM_SUCCESS) {
        for (ULONG i = 0; i < surfaces.GetSize(); ++i) delete surfaces[i];
        return status_of(status);
    }
    for (ULONG i = 0; i < surfaces.GetSize(); ++i) made.push_back(surfaces[i]);
    return surfaces.GetSize() > 0 ? SM_SUCCESS : SM_ERR;
}

sm_status make_edge_blend(const std::vector<Rim>& rims, int32_t order, std::vector<SmSurface*>& made) {
    if (rims.size() != 2 || order < 0 || order > 3) return SM_SHIM_ERR_BAD_TOOL;
    SmSurface* surface = nullptr;
    const SmStatus status = SmBrep::CreateEdgeEdgeBlend(context(), static_cast<ULONG>(order), rims[0].use,
                                                        rims[1].use, 0, nullptr, surface);
    if (status != SM_SUCCESS || !surface) {
        delete surface;
        return status != SM_SUCCESS ? status_of(status) : SM_ERR;
    }
    made.push_back(surface);
    return SM_SUCCESS;
}

// The kept face's surface normal at a point of the bound, as the surface
// itself has it (not the faceuse's).
SmVector3d surface_normal(const SmSurface& surface, const SmPoint2d& uv) {
    SmVector3d d[4];
    surface.Evaluate(uv, 1, 1, TRUE, TRUE, FALSE, d);
    SmVector3d normal = d[2] * d[1];  // Du x Dv
    normal.Unitize();
    return normal;
}

sm_status make_skin(const std::vector<Rim>& rims, double tolerance, std::vector<SmSurface*>& made) {
    if (rims.size() != 2) return SM_SHIM_ERR_BAD_TOOL;
    Owned<SmBSplineCurve> held;
    Bound bounds[2];
    for (int i = 0; i < 2; ++i) {
        const sm_status status = bound_of(rims[i], tolerance, held, bounds[i]);
        if (status != SM_SUCCESS) return status;
    }
    // Run the second section the way the first runs.
    {
        double same = 0.0;
        double opposite = 0.0;
        for (int k = 0; k <= 10; ++k) {
            const double t = k / 10.0;
            same += point_at(*bounds[0].curve, t).DistanceBetween(point_at(*bounds[1].curve, t));
            opposite += point_at(*bounds[0].curve, t).DistanceBetween(point_at(*bounds[1].curve, 1.0 - t));
        }
        if (opposite < same) reverse(bounds[1]);
    }
    SmBSplineSurface* derivs[2] = {nullptr, nullptr};
    Owned<SmBSplineSurface> held_derivs;
    for (int i = 0; i < 2; ++i) {
        if (rims[i].meets == SM_MEETS_POSITION) continue;
        // Which side of the section the gap lies on: away from the kept face.
        const SmExtent1d interval = rims[i].edge->GetInterval();
        SmPoint3d at;
        SmVector3d into;
        SmVector3d along;
        SmVector3d faceuse_normal;
        SmPoint2d uv;
        SmStatus status =
            rims[i].use->EvaluateBinormal(interval.Evaluate(0.5), TRUE, at, into, &along, &faceuse_normal, &uv);
        if (status != SM_SUCCESS) return status_of(status);
        // Right of the section as it runs, which needn't be as its edgeuse
        // runs: the bound's own direction where the edge's middle is.
        along = tangent_at(*bounds[i].curve, 0.5);
        const SmSurface& surface = *rims[i].face->GetSurface();
        const SmVector3d right = along * surface_normal(surface, uv);
        const bool gap_on_the_right = right.Dot(-into) > 0.0;
        // The skin leaves the first section and arrives at the second, as
        // CreateBlendSurface samples its two curves: leaving off to the right,
        // or arriving in from the left.
        const bool leaving = i == 0;
        const SmBoolean side = (leaving ? gap_on_the_right : !gap_on_the_right) ? TRUE : FALSE;
        SmDerivSurfDefinition definition;
        definition.m_eTangencyType = SM_TT_DIR_TO_REFERENCE_CURVE;
        definition.m_pTangentReference = bounds[1 - i].curve;
        definition.m_eCurvatureType = SM_CT_NO_CURVATURE;
        SmExtent1d unit(0.0, 1.0);
        SmBSplineSurface* deriv = nullptr;
        status = definition.CreateDerivativeSurface(context(), surface, *bounds[i].curve, unit, side,
                                                    leaving ? TRUE : FALSE, 0.0, deriv);
        if (status != SM_SUCCESS || !deriv) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        held_derivs.add(deriv);
        derivs[i] = deriv;
    }
    SmTArray<SmBSplineCurve*> sections;
    sections.Add(bounds[0].curve);
    sections.Add(bounds[1].curve);
    SmBSplineSurface* skin = nullptr;
    const SmStatus status = SmBSplineSurface::CreateSkinnedSurface(
        context(), sections, FALSE, SM_SP_U, 0.0, nullptr, nullptr, FALSE, nullptr, derivs, skin, 3);
    if (status != SM_SUCCESS || !skin) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
    made.push_back(skin);
    return SM_SUCCESS;
}

// Labels the faces through shim.cpp's attribute, by their index in the
// brep's face list.
sm_status label_faces(SmBrep& brep, const SmTArray<SmFace*>& made, const sm_label& label) {
    if (!label.present) return SM_SUCCESS;
    SmTArray<SmFace*> all;
    brep.GetFaces(all);
    for (ULONG i = 0; i < made.GetSize(); ++i) {
        for (ULONG j = 0; j < all.GetSize(); ++j) {
            if (all[j] != made[i]) continue;
            const sm_status status = sm_brep_set_face_label(reinterpret_cast<sm_brep*>(&brep), j, &label);
            if (status != SM_SUCCESS) return status;
        }
    }
    return SM_SUCCESS;
}

// The edge_creases points along an edge, as shares of its interval.
const double kShares[SM_EDGE_CURVATURE_SAMPLES] = {0.0,  1e-3, 0.12, 0.25, 0.34,        0.5,
                                                   0.55, 0.75, 0.78, 0.93, 1.0 - 1e-3, 1.0};

// The normal curvature across the edge on the use's face at t, against the
// face's outward normal.
bool curvature_across(const SmEdgeuse& use, double t, double& out) {
    SmPoint3d at;
    SmVector3d into;
    SmVector3d faceuse_normal;
    SmPoint2d uv;
    if (use.EvaluateBinormal(t, TRUE, at, into, nullptr, &faceuse_normal, &uv) != SM_SUCCESS) return false;
    const SmSurface* surface = use.GetFace() ? use.GetFace()->GetSurface() : nullptr;
    if (!surface) return false;
    SmVector3d direction;
    SmVector3d normal;
    double curvature = 0.0;
    if (surface->EvaluateNormalSection(uv, TRUE, TRUE, into, direction, curvature, normal) != SM_SUCCESS) return false;
    out = normal.Dot(faceuse_normal) < 0.0 ? -curvature : curvature;
    return std::isfinite(out);
}

}  // namespace

extern "C" {

sm_status sm_brep_rebuild_faces(const sm_brep* brep, const uint32_t* faces, size_t face_count,
                                const sm_rebuild_surface* surfaces, size_t surface_count, double stitch_tolerance,
                                sm_brep** out, sm_rebuild_report* report) {
    if (!brep || !out || !report || (face_count > 0 && !faces) || (surface_count > 0 && !surfaces)) {
        return SM_SHIM_ERR_NULL_ARGUMENT;
    }
    *report = sm_rebuild_report{};
    report->open_vertex_gap = -1.0;
    report->open_edge_gap = -1.0;
    report->failed_surface = -1;
    report->surface_status = SM_SUCCESS;
    if (face_count == 0 || surface_count == 0 || !std::isfinite(stitch_tolerance)) return SM_ERR_INVALID_INPUT;
    for (size_t i = 0; i < surface_count; ++i) {
        if (!surfaces[i].sides || surfaces[i].side_count == 0) return SM_SHIM_ERR_NULL_ARGUMENT;
    }
    return guarded([&]() -> sm_status {
        sm_brep* handle = nullptr;
        sm_status status = sm_brep_copy(brep, &handle);
        if (status != SM_SUCCESS) return status;
        struct Freed {
            sm_brep* handle;
            ~Freed() {
                if (handle) sm_brep_free(handle);
            }
        } freed{handle};
        SmBrep& copy = *brep_of(handle);

        // The copy's lists are in the original's order; check the picked
        // edges are where they were.
        SmTArray<SmEdge*> edges;
        copy.GetEdges(edges);
        SmTArray<SmEdge*> original_edges;
        brep_of(brep)->GetEdges(original_edges);
        if (edges.GetSize() != original_edges.GetSize()) return SM_ERR;
        SmTArray<SmFace*> all_faces;
        copy.GetFaces(all_faces);
        std::vector<SmFace*> removed;
        for (size_t i = 0; i < face_count; ++i) {
            if (faces[i] >= all_faces.GetSize()) return SM_SHIM_ERR_BAD_FACE;
            SmFace* face = all_faces[faces[i]];
            if (std::find(removed.begin(), removed.end(), face) != removed.end()) return SM_SHIM_ERR_BAD_FACE;
            removed.push_back(face);
        }

        SmExtent3d box;
        copy.CalculateBoundingBox(box, TRUE);
        const double size = box.GetSize().Length();
        // How far apart a loop's sides may end, and how near the kernel
        // fits the surfaces' curves.
        const double loop_gap = std::max(1e-6 * size, 1e-3);
        const double approx = std::max(copy.GetTolerance().val, 1e-7 * size);

        // Make every surface while the faces are all still there.
        const auto surfaces_start = std::chrono::steady_clock::now();
        std::vector<std::vector<SmSurface*>> made(surface_count);
        struct DropMade {
            std::vector<std::vector<SmSurface*>>& made;
            ~DropMade() {
                for (auto& list : made) {
                    for (SmSurface* surface : list) delete surface;
                }
            }
        } drop_made{made};
        for (size_t s = 0; s < surface_count; ++s) {
            const sm_rebuild_surface& request = surfaces[s];
            std::vector<Rim> rims(request.side_count);
            for (size_t i = 0; i < request.side_count; ++i) {
                const SmEdge* picked = original_edges[std::min<ULONG>(request.sides[i].edge, edges.GetSize() - 1)];
                if (request.sides[i].edge < edges.GetSize()) {
                    SmPoint3d was;
                    SmPoint3d is;
                    picked->GetCurve()->EvaluatePoint(picked->GetInterval().Evaluate(0.5), was);
                    SmEdge* now = edges[request.sides[i].edge];
                    now->GetCurve()->EvaluatePoint(now->GetInterval().Evaluate(0.5), is);
                    if (was.DistanceBetween(is) > loop_gap) return SM_ERR;
                }
                status = rim_of(edges, removed, request.sides[i], rims[i]);
                if (status != SM_SUCCESS) return status;
            }
            const double tangent_tol = request.tangent_tol_deg > 0.0 ? request.tangent_tol_deg : 0.01;
            const double fit = request.fit_tolerance > 0.0 ? request.fit_tolerance : approx;
            switch (request.tool) {
                case SM_REBUILD_FILL:
                case SM_REBUILD_COONS: {
                    Owned<SmBSplineCurve> held;
                    std::vector<Bound> bounds(rims.size());
                    for (size_t i = 0; i < rims.size() && status == SM_SUCCESS; ++i) {
                        status = bound_of(rims[i], fit, held, bounds[i]);
                    }
                    if (status == SM_SUCCESS) status = chain(bounds, loop_gap);
                    if (status == SM_SUCCESS) status = join_runs(bounds, held);
                    if (status == SM_SUCCESS) {
                        status = request.tool == SM_REBUILD_FILL ? make_fill(bounds, fit, tangent_tol, made[s])
                                                                 : make_coons(bounds, fit, tangent_tol, made[s]);
                    }
                    break;
                }
                case SM_REBUILD_CORNER:
                    status = make_corner(copy, rims, fit, tangent_tol, made[s]);
                    break;
                case SM_REBUILD_EDGE_BLEND:
                    status = make_edge_blend(rims, request.order, made[s]);
                    break;
                case SM_REBUILD_SKIN:
                    status = make_skin(rims, fit, made[s]);
                    break;
                default:
                    status = SM_SHIM_ERR_BAD_TOOL;
            }
            if (status != SM_SUCCESS) {
                report->failed_surface = static_cast<int64_t>(s);
                report->surface_status = status;
                report->surface_seconds = seconds_since(surfaces_start);
                return status == SM_SHIM_ERR_BAD_TOOL || status == SM_SHIM_ERR_OPEN_LOOP ? status
                                                                                         : SM_SHIM_ERR_SURFACE_FAILED;
            }
        }
        report->surface_seconds = seconds_since(surfaces_start);

        // Take the faces out, add the surfaces and sew them in.
        const auto stitch_start = std::chrono::steady_clock::now();
        SmTArray<SmFace*> taken;
        for (SmFace* face : removed) taken.Add(face);
        SmTemporaryChangeValue<SmBoolean> editing(copy.m_bEditingEnabled, TRUE);
        SmStatus kernel = copy.RemoveFaces(taken);
        if (kernel != SM_SUCCESS) return status_of(kernel);
        for (size_t s = 0; s < surface_count; ++s) {
            std::vector<SmSurface*>& list = made[s];
            while (!list.empty()) {
                SmSurface* surface = list.back();
                list.pop_back();
                const SmExtent2d domain = surface->GetNaturalUVDomain();
                SmTArray<SmFace*> new_faces;
                // The brep owns the surface from here, whether it succeeds or not.
                kernel = copy.CreateFacesFromSurface(surface, domain, SM_CT_C0, new_faces);
                if (kernel != SM_SUCCESS) return status_of(kernel);
                report->faces_made += new_faces.GetSize();
                status = label_faces(copy, new_faces, surfaces[s].label);
                if (status != SM_SUCCESS) return status;
            }
        }

        ULONG glued = 0;
        ULONG lamina = 0;
        double vertex_gap = 0.0;
        double edge_gap = 0.0;
        if (stitch_tolerance > 0.0) {
            SmStitchCallback callback;
            SmStitch stitch(callback, stitch_tolerance, TRUE, FALSE);
            stitch.m_bFastEdgeCompare = TRUE;
            stitch.m_bSplitEdgesWithVertices = TRUE;
            stitch.m_bMakingManifoldSolid = TRUE;
            double open_vertex = SM_BIG_DOUBLE;
            double open_edge = SM_BIG_DOUBLE;
            kernel = stitch.DoStitching(&copy, nullptr, nullptr, glued, lamina, vertex_gap, edge_gap, &open_vertex,
                                        &open_edge);
            if (kernel != SM_SUCCESS) return status_of(kernel);
            report->open_vertex_gap = open_vertex < SM_BIG_DOUBLE ? open_vertex : -1.0;
            report->open_edge_gap = open_edge < SM_BIG_DOUBLE ? open_edge : -1.0;
        } else {
            SmBoolean produced = FALSE;
            kernel = SmStitch::StitchIntoSolid(&copy, produced, glued, vertex_gap, edge_gap);
            if (kernel != SM_SUCCESS) return status_of(kernel);
            SmTArray<SmEdge*> after;
            copy.GetEdges(after);
            for (ULONG i = 0; i < after.GetSize(); ++i) {
                if (after[i]->IsLamina()) ++lamina;
            }
        }
        report->stitched = glued;
        report->open_edges = lamina;
        report->vertex_gap = vertex_gap;
        report->edge_gap = edge_gap;
        if (lamina > 0) {
            report->stitch_seconds = seconds_since(stitch_start);
            return SM_SHIM_ERR_LEFT_OPEN;
        }
        SmBoolean maybe_open = FALSE;
        kernel = copy.OrientTrimmedSurfaces(TRUE, maybe_open, FALSE, nullptr);
        if (kernel != SM_SUCCESS) return status_of(kernel);
        int32_t solid = 0;
        status = sm_brep_is_manifold_solid(handle, &solid);
        report->stitch_seconds = seconds_since(stitch_start);
        if (status != SM_SUCCESS) return status;
        report->solid = solid;
        if (!solid) return SM_SHIM_ERR_NOT_SOLID;
        *out = handle;
        freed.handle = nullptr;
        return SM_SUCCESS;
    });
}

sm_status sm_brep_edge_curvatures(const sm_brep* brep, sm_edge_curvature* out, size_t capacity, size_t* count) {
    if (!brep || !count || (capacity > 0 && !out)) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        SmTArray<SmEdge*> edges;
        brep_of(brep)->GetEdges(edges);
        const size_t written = std::min(capacity, static_cast<size_t>(edges.GetSize()));
        for (size_t i = 0; i < written; ++i) {
            sm_edge_curvature& entry = out[i];
            entry = sm_edge_curvature{};
            const SmEdge* edge = edges[static_cast<ULONG>(i)];
            SmTArray<SmEdgeuse*> uses;
            edge->GetEdgeuses(uses);
            if (uses.GetSize() < 1) continue;
            const SmEdgeuse* use = uses[0];
            const SmEdgeuse* other = use->GetRadial();
            if (!other || other == use || !use->GetFace() || !other->GetFace() || use->GetFace() == other->GetFace()) {
                continue;
            }
            const SmExtent1d interval = edge->GetInterval();
            bool valid = true;
            for (int k = 0; k < SM_EDGE_CURVATURE_SAMPLES && valid; ++k) {
                const double t = interval.Evaluate(kShares[k]);
                valid = curvature_across(*use, t, entry.a[k]) && curvature_across(*other, t, entry.b[k]);
            }
            entry.valid = valid ? 1 : 0;
        }
        *count = edges.GetSize();
        return SM_SUCCESS;
    });
}

sm_status sm_surface_extend(const sm_surface* surface, int32_t side, double distance, int32_t shape,
                            sm_surface** out) {
    if (!surface || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (!(distance > 0.0) || !std::isfinite(distance)) return SM_ERR_INVALID_INPUT;
    SmContinuityType continuity;
    switch (shape) {
        case SM_EXTEND_LINEAR:
            continuity = SM_CT_G1;
            break;
        case SM_EXTEND_SMOOTH:
            continuity = SM_CT_G1_G2;
            break;
        case SM_EXTEND_NATURAL:
            continuity = SM_CT_CINFINITY;
            break;
        case SM_EXTEND_RULED:
            continuity = SM_CT_G1R;
            break;
        default:
            return SM_ERR_INVALID_INPUT;
    }
    SmSurfParamType direction;
    switch (side) {
        case SM_SIDE_UMIN:
            direction = SM_SP_UMIN;
            break;
        case SM_SIDE_VMIN:
            direction = SM_SP_VMIN;
            break;
        case SM_SIDE_UMAX:
            direction = SM_SP_UMAX;
            break;
        case SM_SIDE_VMAX:
            direction = SM_SP_VMAX;
            break;
        case SM_SIDE_ALL:
            direction = SM_SP_BOTH;
            break;
        default:
            return SM_ERR_INVALID_INPUT;
    }
    return guarded([&]() -> sm_status {
        SmSurface* extended = nullptr;
        SmSurface* source = const_cast<SmSurface*>(surface_of(surface));
        const SmStatus status = side == SM_SIDE_ALL
                                    ? source->CreateExtendedSurface(context(), distance, continuity, extended)
                                    : source->CreateExtendedSurface(context(), direction, distance, continuity,
                                                                    extended);
        if (status != SM_SUCCESS || !extended) {
            delete extended;
            return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        }
        *out = surface_handle(extended);
        return SM_SUCCESS;
    });
}

}  // extern "C"
