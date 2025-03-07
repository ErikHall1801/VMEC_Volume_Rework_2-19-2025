#pragma once
#include "vmec_defines.h"
#include "vmec_types.h"



namespace vmec
{
    void m_sandbox_trivial(Vector3 &L, const Vector3 &p0, const Vector3 &p1, const realNumber &dx);

    void m_sandbox_single(Vector3 &L, const Vector3 &p0, const Vector3 &p1, const Vector3 &w0);

    void m_sandbox_single_scatter_ray(Vector3 &L_s, const Vector3 &x_s, const Vector3 &w_s);

    void m_sandbox_trace_null_collision_ray(Vector3 &L, const Vector3 &x0, const Vector3 &w0, const Vector3 &AABB_min, const Vector3 &AABB_max, realNumber &running_majorant);

    bool m_AABB_In_Out(const Vector3 &p0, const Vector3 &AABB_min, const Vector3 &AABB_max);

    Vector3 m_norm_rand_dir();

    Vector3 m_norm_HG_dir(const realNumber &epsilon, const realNumber &g, const Vector3 &w0);

    void m_sandbox_volume(const Vector3 &x_k, Vector3 &sigma_a, Vector3 &sigma_s, Vector3 &sigma_t, Vector3 &L_e, Vector3 &floating_majorant);

    int m_sandbox_mandelbulb(Vector3 &out, const Vector3 &p, const int &imax);

    void m_point_light(Vector3 &pos, Vector3 &col);

    void m_directional_light(const Vector3 &x_s, Vector3 &x_e, Vector3 &col);

    void m_Isotropic_Phase_BSDF(realNumber &pdf);

    void m_Anisotropic_Phase_BSDF(realNumber &pdf, const realNumber &g, const realNumber &costheta);

    void m_Drains_Phase_BSDF(realNumber &pdf, const realNumber &alpha, const realNumber& g, const realNumber &costheta);

    void m_Mie_Phase_BSDF(realNumber &pdf, const realNumber &d, const realNumber &costheta);
}

namespace Magik
{
    void trivial(Vector4 &RGBA_out, const Vector3 &rayOrig, const Vector3 &rayDir);
    void single(Vector4 &RGBA_out, const Vector3 &rayOrig, const Vector3 &rayDir);
    void multiple(Vector4 &RGBA_out, const Vector3 &rayOrig, const Vector3 &rayDir);
    void trivial_relativistic();
    void multiple_relativistic();
}