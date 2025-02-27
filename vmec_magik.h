#pragma once
#include "vmec_defines.h"
#include "vmec_types.h"



namespace vmec
{
    void m_sandbox_trivial(Vector3 &L, const Vector3 &p0, const Vector3 &p1, const realNumber &dx);

    void m_sandbox_single(Vector3 &L, const Vector3 &p0, const Vector3 &p1, const Vector3 &w0, const realNumber &dx);

    void m_sandbox_single_scatter_ray(Vector3 &L_s, const Vector3 &x_s, const Vector3 &w_s, const realNumber &ds);

    void m_sandbox_single_MFP(Vector3 &L, const Vector3 &p0, const Vector3 &p1, const Vector3 &w0);

    void m_sandbox_single_scatter_ray_MFP(Vector3 &L_s, const Vector3 &x_s, const Vector3 &w_s);

    void m_sandbox_estimate_transmittance(Vector3 &T, const Vector3 &x_s, const Vector3 &x_e, const Vector3 &w0);

    void m_sandbox_volume(const Vector3 &x_k, Vector3 &sigma_a, Vector3 &sigma_s, Vector3 &sigma_t, Vector3 &L_e);

    int m_sandbox_mandelbulb(const Vector3 &p, const int &imax);

    void m_point_light(Vector3 &pos, Vector3 &col);

    void m_directional_light(const Vector3 &x_s, Vector3 &x_e, Vector3 &col);

    void m_Isotropic_Phase_BSDF(realNumber &pdf);

    void m_Anisotropic_Phase_BSDF(realNumber &pdf, const realNumber &g, const realNumber & costheta);

    void m_Mean_Free_Path_dx(realNumber &dx, const Vector3& sigma_t, const realNumber& min_dx, const realNumber& max_dx);
}

namespace Magik
{
    void trivial(Vector4 &RGBA_out, const Vector3 &rayOrig, const Vector3 &rayDir);
    void single(Vector4 &RGBA_out, const Vector3 &rayOrig, const Vector3 &rayDir);
    void single_MFP(Vector4 &RGBA_out, const Vector3 &rayOrig, const Vector3 &rayDir);
    void multiple();
    void trivial_relativistic();
    void single_relativistic();
    void multiple_relativistic();
}