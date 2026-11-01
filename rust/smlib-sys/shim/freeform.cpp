// Free-form curves and surfaces for the smlib-sys shim: NURBS curves and
// surfaces from poles or from other curves, made into faces and stitched into
// solids. See the matching section of shim.h.

#include "shim.h"

#include <cmath>
#include <cstring>
#include <memory>
#include <new>
#include <vector>

#include <SmApiCurves.h>
#include <SmApiGeneral.h>
#include <SmApiQueries.h>
#include <SmApiSurfaces.h>

#include <SmAttribute.h>
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmBrep.h>
#include <SmCurve.h>
#include <SmEdge.h>
#include <SmFace.h>
#include <SmStitch.h>
#include <SmSurface.h>
#include <SmVertex.h>

namespace {

SmBrep* brep_of(sm_brep* brep) { return reinterpret_cast<SmBrep*>(brep); }
const SmBrep* brep_of(const sm_brep* brep) { return reinterpret_cast<const SmBrep*>(brep); }
SmCurve* curve_of(sm_curve* curve) { return reinterpret_cast<SmCurve*>(curve); }
const SmCurve* curve_of(const sm_curve* curve) { return reinterpret_cast<const SmCurve*>(curve); }
sm_curve* curve_handle(SmCurve* curve) { return reinterpret_cast<sm_curve*>(curve); }
const SmSurface* surface_of(const sm_surface* surface) { return reinterpret_cast<const SmSurface*>(surface); }
sm_surface* surface_handle(SmSurface* surface) { return reinterpret_cast<sm_surface*>(surface); }

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

// Deletes kernel objects it holds unless released.
template <typename T>
struct Owned {
    std::vector<T*> items;
    ~Owned() {
        for (T* item : items) delete item;
    }
    void release() { items.clear(); }
};

// A copy of the curve as an SmBSplineCurve, or nullptr.
SmBSplineCurve* bspline_copy(const SmCurve* curve) {
    SmCurve* copy = nullptr;
    if (dynamic_cast<const SmBSplineCurve*>(curve)) {
        if (curve->Copy(context(), copy) != SM_SUCCESS) return nullptr;
    } else if (curve->CopyAnalyticAsNurb(context(), copy) != SM_SUCCESS) {
        return nullptr;
    }
    SmBSplineCurve* bspline = dynamic_cast<SmBSplineCurve*>(copy);
    if (!bspline) delete copy;
    return bspline;
}

// The unique knots of a whole knot vector and their multiplicities.
bool split_knots(const double* knots, size_t count, SmTArray<double>& unique, SmTArray<ULONG>& multiplicities) {
    for (size_t i = 0; i < count; ++i) {
        if (!std::isfinite(knots[i])) return false;
        if (i > 0 && knots[i] < knots[i - 1]) return false;
        if (i > 0 && knots[i] == knots[i - 1]) {
            multiplicities[multiplicities.GetSize() - 1] += 1;
        } else {
            unique.Add(knots[i]);
            multiplicities.Add(1);
        }
    }
    return unique.GetSize() >= 2;
}

void put(sm_vec3* out, const SmVector3d& v) {
    out->x = v.x;
    out->y = v.y;
    out->z = v.z;
}

// Labels the faces through shim.cpp's attribute, by their index in the
// brep's face list.
sm_status label_faces(SmBrep& brep, const SmTArray<SmFace*>& made, const sm_label* label) {
    if (!label || !label->present) return SM_SUCCESS;
    SmTArray<SmFace*> all;
    brep.GetFaces(all);
    for (ULONG i = 0; i < made.GetSize(); ++i) {
        for (ULONG j = 0; j < all.GetSize(); ++j) {
            if (all[j] != made[i]) continue;
            const sm_status status = sm_brep_set_face_label(reinterpret_cast<sm_brep*>(&brep), j, label);
            if (status != SM_SUCCESS) return status;
        }
    }
    return SM_SUCCESS;
}

}  // namespace

extern "C" {

sm_status sm_curve_create(int32_t degree, const sm_vec3* poles, const double* weights, size_t count,
                          const double* knots, size_t knot_count, sm_curve** out) {
    if (!poles || !knots || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (degree < 1 || count < static_cast<size_t>(degree) + 1 || knot_count != count + degree + 1) {
        return SM_ERR_INVALID_INPUT;
    }
    return guarded([&]() -> sm_status {
        SmTArray<SmPoint3d> points;
        for (size_t i = 0; i < count; ++i) points.Add(SmPoint3d(poles[i].x, poles[i].y, poles[i].z));
        SmTArray<double> unique;
        SmTArray<ULONG> multiplicities;
        if (!split_knots(knots, knot_count, unique, multiplicities)) return SM_ERR_INVALID_INPUT;
        SmTArray<double> weight_list;
        if (weights) {
            for (size_t i = 0; i < count; ++i) {
                if (!(weights[i] > 0.0) || !std::isfinite(weights[i])) return SM_ERR_INVALID_INPUT;
                weight_list.Add(weights[i]);
            }
        }
        SmBSplineCurve* curve = nullptr;
        const SmApiStatus status = SmApiCreateCanonicalCurve(points, unique, multiplicities,
                                                             weights ? &weight_list : nullptr,
                                                             static_cast<ULONG>(degree), SM_CF_UNSPECIFIED, curve);
        if (status != SM_SUCCESS || !curve) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        *out = curve_handle(curve);
        return SM_SUCCESS;
    });
}

sm_status sm_curve_interpolate(const sm_vec3* points, size_t count, int32_t degree, int32_t closed,
                               sm_curve** out) {
    if (!points || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (degree < 1 || count < 2) return SM_ERR_INVALID_INPUT;
    return guarded([&]() -> sm_status {
        SmTArray<SmPoint3d> list;
        for (size_t i = 0; i < count; ++i) list.Add(SmPoint3d(points[i].x, points[i].y, points[i].z));
        SmBSplineCurve* curve = nullptr;
        const SmStatus status = SmBSplineCurve::InterpolatePoints(
            context(), list, nullptr, static_cast<ULONG>(degree), nullptr, nullptr, closed ? TRUE : FALSE,
            SM_IT_CHORDLENGTH, curve);
        if (status != SM_SUCCESS || !curve) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        *out = curve_handle(curve);
        return SM_SUCCESS;
    });
}

sm_status sm_curve_of_edge(const sm_brep* brep, size_t edge, sm_curve** out) {
    if (!brep || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        SmTArray<SmEdge*> edges;
        const_cast<SmBrep*>(brep_of(brep))->GetEdges(edges);
        if (edge >= edges.GetSize()) return SM_SHIM_ERR_BAD_EDGE;
        const SmEdge* picked = edges[static_cast<ULONG>(edge)];
        if (!picked->GetCurve()) return SM_ERR;
        std::unique_ptr<SmBSplineCurve> curve(bspline_copy(picked->GetCurve()));
        if (!curve) return SM_ERR;
        SmExtent1d interval = picked->GetInterval();
        const SmStatus trimmed = curve->Trim(interval);
        if (trimmed != SM_SUCCESS) return status_of(trimmed);
        *out = curve_handle(curve.release());
        return SM_SUCCESS;
    });
}

sm_status sm_curve_constant_like(const sm_curve* like, sm_vec3 vector, sm_curve** out) {
    if (!like || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        const SmBSplineCurve* source = dynamic_cast<const SmBSplineCurve*>(curve_of(like));
        if (!source || source->IsRational()) return SM_ERR_INVALID_INPUT;
        std::unique_ptr<SmBSplineCurve> curve(bspline_copy(source));
        if (!curve) return SM_ERR;
        const SmPoint3d pole(vector.x, vector.y, vector.z);
        for (ULONG i = 0; i < curve->GetNumberControlPoints(); ++i) {
            const SmStatus set = curve->SetControlPoint(SM_CP_NON_RATIONAL, i, pole, 1.0);
            if (set != SM_SUCCESS) return status_of(set);
        }
        *out = curve_handle(curve.release());
        return SM_SUCCESS;
    });
}

sm_status sm_curve_copy(const sm_curve* curve, sm_curve** out) {
    if (!curve || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        SmCurve* copy = nullptr;
        const SmStatus status = curve_of(curve)->Copy(context(), copy);
        if (status != SM_SUCCESS || !copy) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        *out = curve_handle(copy);
        return SM_SUCCESS;
    });
}

sm_status sm_curve_reverse(sm_curve* curve) {
    if (!curve) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        SmCurve* own = curve_of(curve);
        SmExtent1d turned;
        return status_of(own->ReverseParameterization(own->GetNaturalInterval(), turned));
    });
}

sm_status sm_curve_interval(const sm_curve* curve, double* min, double* max) {
    if (!curve || !min || !max) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        const SmExtent1d interval = curve_of(curve)->GetNaturalInterval();
        *min = interval.GetMin();
        *max = interval.GetMax();
        return SM_SUCCESS;
    });
}

sm_status sm_curve_evaluate(const sm_curve* curve, double t, sm_vec3* point, sm_vec3* derivative) {
    if (!curve || !point) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        SmVector3d values[2];
        const SmStatus status = curve_of(curve)->Evaluate(t, 1, TRUE, values, FALSE);
        if (status != SM_SUCCESS) return status_of(status);
        put(point, values[0]);
        if (derivative) put(derivative, values[1]);
        return SM_SUCCESS;
    });
}

sm_status sm_curve_shape(const sm_curve* curve, int32_t* degree, size_t* poles, int32_t* rational) {
    if (!curve || !degree || !poles || !rational) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        const SmBSplineCurve* bspline = dynamic_cast<const SmBSplineCurve*>(curve_of(curve));
        if (!bspline) return SM_ERR_INVALID_INPUT;
        *degree = static_cast<int32_t>(bspline->GetDegree());
        *poles = bspline->GetNumberControlPoints();
        *rational = bspline->IsRational() ? 1 : 0;
        return SM_SUCCESS;
    });
}

sm_status sm_curve_net(const sm_curve* curve, sm_vec3* poles, double* weights, double* knots) {
    if (!curve || !poles || !knots) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        const SmBSplineCurve* bspline = dynamic_cast<const SmBSplineCurve*>(curve_of(curve));
        if (!bspline) return SM_ERR_INVALID_INPUT;
        const ULONG count = bspline->GetNumberControlPoints();
        for (ULONG i = 0; i < count; ++i) {
            SmPoint3d pole;
            double weight = 1.0;
            const SmStatus got = bspline->GetControlPoint(SM_CP_EUCLIDIAN_RATIONAL, i, pole, weight);
            if (got != SM_SUCCESS) return status_of(got);
            put(&poles[i], pole);
            if (weights) weights[i] = bspline->IsRational() ? weight : 1.0;
        }
        SmTArray<double> unique;
        SmTArray<ULONG> multiplicities;
        const SmStatus got = bspline->GetKnots(unique, &multiplicities);
        if (got != SM_SUCCESS) return status_of(got);
        const size_t wanted = count + bspline->GetDegree() + 1;
        size_t at = 0;
        for (ULONG i = 0; i < unique.GetSize(); ++i) {
            for (ULONG m = 0; m < multiplicities[i]; ++m) {
                if (at == wanted) return SM_ERR;
                knots[at++] = unique[i];
            }
        }
        return at == wanted ? SM_SUCCESS : SM_ERR;
    });
}

void sm_curve_free(sm_curve* curve) { delete curve_of(curve); }

sm_status sm_surface_create(int32_t degree_u, int32_t degree_v, const sm_vec3* poles, size_t count_u, size_t count_v,
                            const double* weights, const double* knots_u, size_t knot_count_u, const double* knots_v,
                            size_t knot_count_v, sm_surface** out) {
    if (!poles || !knots_u || !knots_v || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (degree_u < 1 || degree_v < 1 || knot_count_u != count_u + degree_u + 1 ||
        knot_count_v != count_v + degree_v + 1) {
        return SM_ERR_INVALID_INPUT;
    }
    return guarded([&]() -> sm_status {
        SmTArray<SmPoint3d> points;
        for (size_t i = 0; i < count_u * count_v; ++i) points.Add(SmPoint3d(poles[i].x, poles[i].y, poles[i].z));
        SmTArray<double> unique_u, unique_v;
        SmTArray<ULONG> multiplicities_u, multiplicities_v;
        if (!split_knots(knots_u, knot_count_u, unique_u, multiplicities_u) ||
            !split_knots(knots_v, knot_count_v, unique_v, multiplicities_v)) {
            return SM_ERR_INVALID_INPUT;
        }
        SmTArray<double> weight_list;
        if (weights) {
            for (size_t i = 0; i < count_u * count_v; ++i) weight_list.Add(weights[i]);
        }
        SmBSplineSurface* surface = nullptr;
        const SmStatus status = SmBSplineSurface::CreateCanonical(
            context(), static_cast<ULONG>(degree_u), static_cast<ULONG>(degree_v), points, SM_SF_UNSPECIFIED,
            multiplicities_u, multiplicities_v, unique_u, unique_v, SM_KT_UNSPECIFIED, weights ? &weight_list : nullptr,
            nullptr, surface);
        if (status != SM_SUCCESS || !surface) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        *out = surface_handle(surface);
        return SM_SUCCESS;
    });
}

namespace {

// Copies of the curves as SmBSplineCurves, owned by the holder.
sm_status bsplines(const sm_curve* const* curves, size_t count, Owned<SmBSplineCurve>& held,
                   SmTArray<SmBSplineCurve*>& out, bool allow_null = false) {
    for (size_t i = 0; i < count; ++i) {
        if (!curves[i]) {
            if (!allow_null) return SM_SHIM_ERR_NULL_ARGUMENT;
            out.Add(nullptr);
            continue;
        }
        SmBSplineCurve* copy = bspline_copy(curve_of(curves[i]));
        if (!copy) return SM_ERR;
        held.items.push_back(copy);
        out.Add(copy);
    }
    return SM_SUCCESS;
}

}  // namespace

sm_status sm_surface_loft(const sm_curve* const* sections, size_t count, int32_t degree, sm_surface** out) {
    if (!sections || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (count < 2 || degree < 1) return SM_ERR_INVALID_INPUT;
    return guarded([&]() -> sm_status {
        Owned<SmBSplineCurve> held;
        SmTArray<SmBSplineCurve*> curves;
        const sm_status copied = bsplines(sections, count, held, curves);
        if (copied != SM_SUCCESS) return copied;
        SmBSplineSurface* derivatives[2] = {nullptr, nullptr};
        SmBSplineSurface* surface = nullptr;
        const SmStatus status = SmBSplineSurface::CreateSkinnedSurface(
            context(), curves, FALSE, SM_SP_U, 0.0, nullptr, nullptr, FALSE, nullptr, derivatives, surface,
            static_cast<ULONG>(degree));
        if (status != SM_SUCCESS || !surface) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        *out = surface_handle(surface);
        return SM_SUCCESS;
    });
}

sm_status sm_surface_gordon(const sm_curve* const* u, size_t u_count, const sm_curve* const* v, size_t v_count,
                            sm_surface** out) {
    if (!u || !v || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (u_count < 2 || v_count < 2) return SM_ERR_INVALID_INPUT;
    return guarded([&]() -> sm_status {
        Owned<SmBSplineCurve> held;
        SmTArray<SmBSplineCurve*> u_curves, v_curves;
        sm_status copied = bsplines(u, u_count, held, u_curves);
        if (copied != SM_SUCCESS) return copied;
        copied = bsplines(v, v_count, held, v_curves);
        if (copied != SM_SUCCESS) return copied;
        SmBSplineSurface* surface = nullptr;
        const SmStatus status = SmBSplineSurface::CreateGordonSurface(context(), u_curves, v_curves, surface);
        if (status != SM_SUCCESS || !surface) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        *out = surface_handle(surface);
        return SM_SUCCESS;
    });
}

sm_status sm_surface_gordon_at(const sm_curve* const* u, const double* u_at, size_t u_count, const sm_curve* const* v,
                               const double* v_at, size_t v_count, int32_t degree_u, int32_t degree_v,
                               sm_surface** out) {
    if (!u || !u_at || !v || !v_at || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (u_count < 2 || v_count < 2 || degree_u < 1 || degree_v < 1) return SM_ERR_INVALID_INPUT;
    return guarded([&]() -> sm_status {
        Owned<SmBSplineCurve> held;
        SmTArray<SmBSplineCurve*> u_curves, v_curves;
        sm_status copied = bsplines(u, u_count, held, u_curves);
        if (copied != SM_SUCCESS) return copied;
        copied = bsplines(v, v_count, held, v_curves);
        if (copied != SM_SUCCESS) return copied;
        SmTArray<double> u_params, v_params;
        for (size_t i = 0; i < u_count; ++i) u_params.Add(u_at[i]);
        for (size_t i = 0; i < v_count; ++i) v_params.Add(v_at[i]);
        SmBSplineSurface* surface = nullptr;
        const SmStatus status = SmBSplineSurface::CreateGordonSurfaceAt(
            context(), u_curves, u_params, v_curves, v_params, static_cast<ULONG>(degree_u),
            static_cast<ULONG>(degree_v), surface);
        if (status != SM_SUCCESS || !surface) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        *out = surface_handle(surface);
        return SM_SUCCESS;
    });
}

sm_status sm_surface_coons(const sm_curve* const* boundary, const sm_curve* const* derivatives, double tolerance,
                           sm_surface** out) {
    if (!boundary || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        Owned<SmBSplineCurve> held;
        SmTArray<SmBSplineCurve*> curves, derivs;
        sm_status copied = bsplines(boundary, 4, held, curves);
        if (copied != SM_SUCCESS) return copied;
        if (derivatives) {
            copied = bsplines(derivatives, 4, held, derivs);
            if (copied != SM_SUCCESS) return copied;
        }
        SmBSplineSurface* surface = nullptr;
        const SmStatus status = SmBSplineSurface::CreateCoonsPatch(
            context(), tolerance > 0.0 ? tolerance : 1e-4, derivatives ? FALSE : TRUE, curves, derivs, surface);
        if (status != SM_SUCCESS || !surface) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        *out = surface_handle(surface);
        return SM_SUCCESS;
    });
}

sm_status sm_surface_fill(const sm_curve* const* boundary, const sm_curve* const* derivatives, size_t count,
                          double angle_tol_deg, sm_surface** out, size_t* out_count) {
    if (!boundary || !out || !out_count) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (count < 3) return SM_ERR_INVALID_INPUT;
    return guarded([&]() -> sm_status {
        Owned<SmBSplineCurve> held;
        SmTArray<SmBSplineCurve*> curves, derivs;
        sm_status copied = bsplines(boundary, count, held, curves);
        if (copied != SM_SUCCESS) return copied;
        for (ULONG i = 0; i < curves.GetSize(); ++i) {
            if (curves[i]->IsRational()) return SM_ERR_INVALID_INPUT;
        }
        if (derivatives) {
            copied = bsplines(derivatives, count, held, derivs, true);
            if (copied != SM_SUCCESS) return copied;
        } else {
            for (size_t i = 0; i < count; ++i) derivs.Add(nullptr);
        }
        SmTArray<SmSurface*> surfaces;
        const SmStatus status =
            SmBSplineSurface::CreateNSidedPatch(context(), curves, derivs, angle_tol_deg, surfaces);
        if (status != SM_SUCCESS) {
            for (ULONG i = 0; i < surfaces.GetSize(); ++i) delete surfaces[i];
            return status_of(status);
        }
        if (surfaces.GetSize() > count) {
            for (ULONG i = 0; i < surfaces.GetSize(); ++i) delete surfaces[i];
            return SM_ERR;
        }
        for (ULONG i = 0; i < surfaces.GetSize(); ++i) out[i] = surface_handle(surfaces[i]);
        *out_count = surfaces.GetSize();
        return SM_SUCCESS;
    });
}

sm_status sm_surface_domain(const sm_surface* surface, double min[2], double max[2]) {
    if (!surface || !min || !max) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        const SmExtent2d domain = surface_of(surface)->GetNaturalUVDomain();
        min[0] = domain.GetMin().x;
        min[1] = domain.GetMin().y;
        max[0] = domain.GetMax().x;
        max[1] = domain.GetMax().y;
        return SM_SUCCESS;
    });
}

sm_status sm_surface_evaluate(const sm_surface* surface, double u, double v, sm_vec3* point, sm_vec3* normal) {
    if (!surface || !point) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        SmVector3d at;
        SmApiStatus status = SmApiEvaluateSurfacePoint(surface_of(surface), SmVector2d(u, v), at);
        if (status != SM_SUCCESS) return status_of(status);
        put(point, at);
        if (normal) {
            SmVector3d n;
            status = SmApiEvaluateSurfaceNormal(surface_of(surface), SmVector2d(u, v), n);
            if (status != SM_SUCCESS) return status_of(status);
            put(normal, n);
        }
        return SM_SUCCESS;
    });
}

namespace {

// The surface as an SmBSplineSurface: itself, or an analytic one converted
// into held.
const SmBSplineSurface* as_bspline(const SmSurface* surface, std::unique_ptr<SmSurface>& held) {
    if (const SmBSplineSurface* bspline = dynamic_cast<const SmBSplineSurface*>(surface)) return bspline;
    SmSurface* copy = nullptr;
    if (surface->CopyAnalyticAsNurb(context(), copy) != SM_SUCCESS || !copy) return nullptr;
    held.reset(copy);
    return dynamic_cast<const SmBSplineSurface*>(copy);
}

// The whole knot vector along one parameter, poles + degree + 1 knots.
bool whole_knots(const SmBSplineSurface& surface, SmSurfParamType along, double* out) {
    SmTArray<double> unique;
    SmTArray<ULONG> multiplicities;
    if (surface.GetKnots(along, unique, &multiplicities) != SM_SUCCESS) return false;
    const ULONG wanted = surface.GetNumberControlPoints(along) + surface.GetDegree(along) + 1;
    ULONG at = 0;
    for (ULONG i = 0; i < unique.GetSize(); ++i) {
        for (ULONG m = 0; m < multiplicities[i]; ++m) {
            if (at == wanted) return false;
            out[at++] = unique[i];
        }
    }
    return at == wanted;
}

}  // namespace

sm_status sm_surface_shape(const sm_surface* surface, int32_t* degree_u, int32_t* degree_v, size_t* poles_u,
                           size_t* poles_v, int32_t* rational) {
    if (!surface || !degree_u || !degree_v || !poles_u || !poles_v || !rational) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        std::unique_ptr<SmSurface> held;
        const SmBSplineSurface* bspline = as_bspline(surface_of(surface), held);
        if (!bspline) return SM_ERR_INVALID_INPUT;
        *degree_u = static_cast<int32_t>(bspline->GetDegree(SM_SP_U));
        *degree_v = static_cast<int32_t>(bspline->GetDegree(SM_SP_V));
        *poles_u = bspline->GetNumberControlPoints(SM_SP_U);
        *poles_v = bspline->GetNumberControlPoints(SM_SP_V);
        *rational = bspline->IsRational() ? 1 : 0;
        return SM_SUCCESS;
    });
}

sm_status sm_surface_net(const sm_surface* surface, sm_vec3* poles, double* weights, double* knots_u,
                         double* knots_v) {
    if (!surface || !poles || !knots_u || !knots_v) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        std::unique_ptr<SmSurface> held;
        const SmBSplineSurface* bspline = as_bspline(surface_of(surface), held);
        if (!bspline) return SM_ERR_INVALID_INPUT;
        const ULONG count_u = bspline->GetNumberControlPoints(SM_SP_U);
        const ULONG count_v = bspline->GetNumberControlPoints(SM_SP_V);
        for (ULONG i = 0; i < count_u; ++i) {
            for (ULONG j = 0; j < count_v; ++j) {
                SmPoint3d pole;
                double weight = 1.0;
                const SmStatus got = bspline->GetControlPoint(SM_CP_EUCLIDIAN_RATIONAL, i, j, pole, weight);
                if (got != SM_SUCCESS) return status_of(got);
                put(&poles[i * count_v + j], pole);
                if (weights) weights[i * count_v + j] = bspline->IsRational() ? weight : 1.0;
            }
        }
        if (!whole_knots(*bspline, SM_SP_U, knots_u) || !whole_knots(*bspline, SM_SP_V, knots_v)) return SM_ERR;
        return SM_SUCCESS;
    });
}

sm_status sm_surface_closest_point(const sm_surface* surface, sm_vec3 point, sm_vec3* closest, double uv[2],
                                   double* distance) {
    if (!surface || !closest || !uv || !distance) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        SmPoint3d nearest;
        SmTArray<SmPoint2d> uvs;
        double gap = 0.0;
        const SmApiStatus status =
            SmApiSurfaceClosestPoint(surface_of(surface), SmPoint3d(point.x, point.y, point.z), nearest, uvs, gap);
        if (status != SM_SUCCESS) return status_of(status);
        if (uvs.GetSize() == 0) return SM_ERR;
        put(closest, nearest);
        uv[0] = uvs[0].x;
        uv[1] = uvs[0].y;
        *distance = gap;
        return SM_SUCCESS;
    });
}

sm_status sm_surface_of_face(const sm_brep* brep, size_t face, sm_surface** out) {
    if (!brep || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        SmTArray<SmFace*> faces;
        const_cast<SmBrep*>(brep_of(brep))->GetFaces(faces);
        if (face >= faces.GetSize()) return SM_SHIM_ERR_BAD_FACE;
        const SmSurface* surface = faces[static_cast<ULONG>(face)]->GetSurface();
        if (!surface) return SM_ERR;
        SmSurface* copy = nullptr;
        const SmStatus status = surface->Copy(context(), copy);
        if (status != SM_SUCCESS || !copy) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        *out = surface_handle(copy);
        return SM_SUCCESS;
    });
}

void sm_surface_free(sm_surface* surface) { delete reinterpret_cast<SmSurface*>(surface); }

sm_status sm_brep_create_empty(sm_brep** out) {
    if (!out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        SmBrep* brep = new (context()) SmBrep();
        *out = reinterpret_cast<sm_brep*>(brep);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_add_surface(sm_brep* brep, const sm_surface* surface, const sm_label* label, size_t* faces) {
    if (!brep || !surface || !faces) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        SmSurface* copy = nullptr;
        SmStatus status = surface_of(surface)->Copy(context(), copy);
        if (status != SM_SUCCESS || !copy) return status != SM_SUCCESS ? status_of(status) : SM_ERR;
        const SmExtent2d domain = copy->GetNaturalUVDomain();
        SmTArray<SmFace*> made;
        status = brep_of(brep)->CreateFacesFromSurface(copy, domain, SM_CT_C1, made);
        if (status != SM_SUCCESS) return status_of(status);
        *faces = made.GetSize();
        return label_faces(*brep_of(brep), made, label);
    });
}

sm_status sm_brep_add_planar_face(sm_brep* brep, const sm_curve* const* loop, size_t count, const sm_label* label) {
    if (!brep || !loop) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (count == 0) return SM_ERR_INVALID_INPUT;
    return guarded([&]() -> sm_status {
        SmTArray<SmCurve*> curves;
        for (size_t i = 0; i < count; ++i) {
            if (!loop[i]) {
                for (ULONG j = 0; j < curves.GetSize(); ++j) delete curves[j];
                return SM_SHIM_ERR_NULL_ARGUMENT;
            }
            SmCurve* copy = nullptr;
            if (curve_of(loop[i])->Copy(context(), copy) != SM_SUCCESS || !copy) {
                for (ULONG j = 0; j < curves.GetSize(); ++j) delete curves[j];
                return SM_ERR;
            }
            curves.Add(copy);
        }
        const double tolerance = brep_of(brep)->GetTolerance().val;
        SmTArray<SmFace*> made;
        const SmStatus status = brep_of(brep)->CreatePlanarFacesWith3DCurves(
            nullptr, curves, tolerance > 0.0 ? tolerance : 1e-6, made);
        if (status != SM_SUCCESS) return status_of(status);
        if (made.GetSize() == 0) return SM_ERR;
        return label_faces(*brep_of(brep), made, label);
    });
}

sm_status sm_brep_stitch(sm_brep* brep, double tolerance, int32_t* solid) {
    if (!brep || !solid) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        const SmStatus status = brep_of(brep)->StitchAndOrient(tolerance > 0.0 ? tolerance : -0.1);
        if (status != SM_SUCCESS) return status_of(status);
        return sm_brep_is_manifold_solid(brep, solid);
    });
}

sm_status sm_brep_stitch_into_solid(sm_brep* brep, int32_t* solid, uint32_t* stitched, double* vertex_gap,
                                    double* edge_gap) {
    if (!brep || !solid || !stitched || !vertex_gap || !edge_gap) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&]() -> sm_status {
        SmBoolean produced = FALSE;
        ULONG glued = 0;
        const SmStatus status =
            SmStitch::StitchIntoSolid(brep_of(brep), produced, glued, *vertex_gap, *edge_gap);
        if (status != SM_SUCCESS) return status_of(status);
        *solid = produced ? 1 : 0;
        *stitched = static_cast<uint32_t>(glued);
        return SM_SUCCESS;
    });
}

}  // extern "C"
