#include "vmec_common.h"
#include "vmec_magik.h"
#include "vmec_magik.h"

namespace vmec
{
    /*

    Couple of thoughts before we begin.

    There are a couple of considerations regarding the translation from flat to curved spacetime. The Pixar paper implicitly assumes flat spatial geometry. 

    The most immediate matter is that of the geodesic itself. The geodesic is not straight in Cartesian coordinates, so would the phase function apply between steps ?
    I dont think it should, because the geodesic is a straight line in the physical sense. So Magik should threat the entire geodesic as a continues straight line.

    To make sure everything does its job, we will need a couple of debug tools specifically for Magik. 
        For instance, a tool that calculates the divergence between Geodesic and light path. That is, some way of determining how much the geodesic and path walked by Magik overlap
        Ideally it is 1:1, but if we have any fault in our logic a value other than 1 implies that. 
        One way to do this is to divide the final euclidian lengths. The geodesic and lightpath should have the exact same length. Well, not exact same length. But close to it. 

    Regarding the syntax; Its m_ for functions from this file

    Lets write it out

    L = direct_radiance + incoming_radiance;
            direct_radiance = integral_t_d( T(t)*sigma_a(x_t)*L_e(x_t,w)*dt )
            incoming_radiance = integral_t_d( T(t)*sigma_s(x_t)*L_s(x_t,w)*dt )

        where;
            T(t) = exp(-sigma_t(x_t)*dt)
            sigma_t(x_t) = sigma_a(x_t) + sigma_s(x_t)
            L_e(x_t,w) = Emission at position, invarient of direction
            L_s(x_t,w) = integral_S²( fp(x,w,w´)*L(x,w)*dw )
                L(x,w)
                    For Single Scattering;   L(x,w) = direct_radiance (where w is towards the light source(s))
                    For Multiple Scattering; L(x,w) = integral_t_MFPL(x_t)( direct_radiance ) (Where the integral runs for N scatter events and the direction changes according to the BSDF each time)
                        MFPL(x_t) = 1.0 / sigma_t(x);

    Lets consider what trivial, singel and multiple are supposed to do, which parts of the L equation above do they solve ? 

    trivial; L = direct_radiance // For the case where -> integral_t_d( T(t)*sigma_a(x_t)*L_e(x_t,w)*dt )
    single; L = direct_radiance + incoming_radiance // For the case of only one scatter event towards a descrete, point, light source
    multiple; L = direct_radiance + incoming_radiance //

    Note; T(), sigma_t(), sigma_s() and sigma_a() are not floats, they are colors
        So T(t) for instance is actually
        T(t) = (exp(-sigma_t.x(x_t)*dt) , exp(-sigma_t.y(x_t)*dt) , exp(-sigma_t.z(x_t)*dt))

    It is important to keep these functions organized. We use m_ for Magik, well what about mT_ for "Magik Trivial", mS_ "Magik Single", mM "Magik Multiple" and then mMR (Magik Multiple Relativistic)
        So; mT_T(t), mT_sigma_t(), mT_sigma_s(), mT_sigma_a()
        and then later
            mMR_T(t) . . . mMR_sigma_a()

    Another thing to consider is that the sigma_a and sigma_s functions will be identical for all sandbox cases. Just the mandelbulb plus some emissive sphere.
    So its probably best if this one function, void, just returns 3 values, scatter, absorption and emission.

    */

    //Sandbox Scene Specific functions
        //trivial sandbox - direct radiance only
        void m_sandbox_trivial(Vector3 &L, const Vector3 &p0, const Vector3 &p1, const realNumber &dx)
        {
            realNumber s = euclidianDistance(p0,p1), k = 0.0;
            Vector3 x0 = p0, w = (p1-p0) / VectorLength(p1-p0), x1 = p0, BL;
            int n_steps = int(s/dx);
            realNumber fractionaldx = (s - realNumber(n_steps*dx))*dx;

            Vector3 T(1.0,1.0,1.0), sigma_a(0.0,0.0,0.0), sigma_s(0.0,0.0,0.0), sigma_t(0.0,0.0,0.0), L_e(0.0,0.0,0.0);

            for(int i = 0; i <= n_steps; i++)
            {
                k = realNumber(i)*dx;
                if(i == n_steps)
                {
                    k += fractionaldx;
                }

                //Volume Sample
                m_sandbox_volume(x1,sigma_a,sigma_s,sigma_t,L_e);

                //Transmittance in Interval
                T.set(std::exp(-sigma_t.x*dx) * T.x , std::exp(-sigma_t.y*dx) * T.y , std::exp(-sigma_t.z*dx) * T.z);

                //Direct Light
                L.set(L.x + (T.x * sigma_a.x * L_e.x * dx) , L.y + (T.y * sigma_a.y * L_e.y * dx) , L.z + (T.z * sigma_a.z * L_e.z * dx));

                x1 = x0 + w*k;
            }
        }

        //Single scattering (direct + single indreict radiance) - Constant dx
            //single sandbox - direct + single indirect radiance
            void m_sandbox_single(Vector3 &L, const Vector3 &p0, const Vector3 &p1, const Vector3 &w0, const realNumber &dx)
            {
                realNumber s = euclidianDistance(p0,p1), k = 0.0;
                Vector3 x0 = p0, w = (p1-p0) / VectorLength(p1-p0), x1 = p0, BL;
                int n_steps = int(s/dx);
                realNumber fractionaldx = (s - realNumber(n_steps*dx))*dx;

                Vector3 T(1.0,1.0,1.0), sigma_a(0.0,0.0,0.0), sigma_s(0.0,0.0,0.0), sigma_t(0.0,0.0,0.0), L_e(0.0,0.0,0.0), L_s(0.0,0.0,0.0), direct_radiance(0.0,0.0,0.0), indirect_single_radiance(0.0,0.0,0.0);

                //L = direct_radiance + incoming_single_radiance;
                //incoming_single_radiance = integral_t_d( T(t)*sigma_s(x_t)*L_s(x_t,w)*dt )
                //L_s(x_t,w) = integral_S²( fp(x,w,w´)*L(x,w)*dw )

                for(int i = 0; i <= n_steps; i++)
                {
                    k = realNumber(i)*dx;
                    if(i == n_steps)
                    {
                        k += fractionaldx;
                    }

                    //Volume Sample
                    m_sandbox_volume(x1,sigma_a,sigma_s,sigma_t,L_e);

                    //Transmittance in Interval
                    T.set(std::exp(-sigma_t.x*dx) * T.x , std::exp(-sigma_t.y*dx) * T.y , std::exp(-sigma_t.z*dx) * T.z);

                    //Direct Radiance
                    direct_radiance.set( (T.x * sigma_a.x * L_e.x * dx) , (T.y * sigma_a.y * L_e.y * dx) , (T.z * sigma_a.z * L_e.z * dx) );

                    //Indirect Single Radiance
                    m_sandbox_single_scatter_ray(L_s,x1,w0,0.1);
                    indirect_single_radiance.set( (T.x * sigma_s.x * L_s.x * dx) , (T.y * sigma_s.y * L_s.y * dx) , (T.z * sigma_s.z * L_s.z * dx) );

                    //Total Radiance 
                    L.set(L.x + direct_radiance.x + indirect_single_radiance.x , L.y + direct_radiance.y + indirect_single_radiance.y , L.z + direct_radiance.z + indirect_single_radiance.z);

                    x1 = x0 + w*k;
                }
            }

            //Evaluate Single scatter event
            void m_sandbox_single_scatter_ray(Vector3 &L_s, const Vector3 &x_s, const Vector3 &w_s, const realNumber &ds)
            {
                //L_s(x_t,w) = integral_S( fp(x,w,w´)*L(x,w)*dw )
                
                int n_steps = 0;
                realNumber cosTheta = 0.0, s = 0.0, k = 0.0, pdf = 0.0;
                Vector3 x_t = x_s, x_e(0.0,0.0,0.0), w_e(0.0,0.0,0.0), light_col(0.0,0.0,0.0); //_e means "end", whereas _s means "start"
                //m_point_light(x_e,light_col);
                m_directional_light(x_s,x_e,light_col);

                //Compute Scatter direction (w_e) and Angle (cosTheta)
                w_e = x_e - x_s;
                w_e = w_e / VectorLength(w_e);
                cosTheta = dotProduct(w_e,w_s);

                //Compute travel length and steps
                s = euclidianDistance(x_e,x_s);
                n_steps = int(s / ds);

                //Declare scatter variables
                Vector3 T(1.0,1.0,1.0), sigma_a(0.0,0.0,0.0), sigma_s(0.0,0.0,0.0), sigma_t(0.0,0.0,0.0), L_e(0.0,0.0,0.0), direct_radiance(0.0,0.0,0.0);

                //Scatter
                for(int i = 0; i <= n_steps; i++)
                {
                    k = realNumber(i)*ds;

                    //Volume Sample
                    m_sandbox_volume(x_t,sigma_a,sigma_s,sigma_t,L_e);

                    //Transmittance along ray
                    T.set(std::exp(-sigma_t.x*ds) * T.x , std::exp(-sigma_t.y*ds) * T.y , std::exp(-sigma_t.z*ds) * T.z);

                    //Accumilate Direct Radiance
                    direct_radiance.set( direct_radiance.x + (T.x * sigma_a.x * L_e.x * ds) , direct_radiance.y + (T.y * sigma_a.y * L_e.y * ds) , direct_radiance.z + (T.z * sigma_a.z * L_e.z * ds) );

                    x_t = x_s + w_e*k;
                }

                //Add Directional light
                direct_radiance.set( direct_radiance.x + (light_col.x*T.x) , direct_radiance.y + (light_col.y*T.y), direct_radiance.z + (light_col.z*T.z) );

                //m_Isotropic_Phase_BSDF(pdf);
                m_Anisotropic_Phase_BSDF(pdf,0.7,cosTheta);
                L_s = direct_radiance*pdf;
            }

        //Single scattering (direct + single indirect radiance) - Mean Free Path dx
            //single scattering
            void m_sandbox_single_MFP(Vector3 &L, const Vector3 &p0, const Vector3 &p1, const Vector3 &w0)
            {
                bool track_back = false;
                realNumber s = euclidianDistance(p0,p1), k = 0.0, dx = 0.0, dx_base = 0.15, dx_scale = 0.0, T_max = 1.0;
                Vector3 x0 = p0, w = (p1-p0) / VectorLength(p1-p0), x1 = p0, BL;

                Vector3 T(1.0,1.0,1.0), T_estimate(1.0,1.0,1.0), simga_t_old(0.0,0.0,0.0), sigma_a(0.0,0.0,0.0), sigma_s(0.0,0.0,0.0), sigma_t(0.0,0.0,0.0), L_e(0.0,0.0,0.0), L_s(0.0,0.0,0.0), direct_radiance(0.0,0.0,0.0), indirect_single_radiance(0.0,0.0,0.0);

                while(s >= k)
                {
                    //Volume Sample
                    m_sandbox_volume(x1,sigma_a,sigma_s,sigma_t,L_e);
                    
                    //Determine dx
                    T_max = std::max(std::max(T.x,T.y),T.z);
                    dx_scale = (dx_base + ((1.0-T_max)*(1.0-T_max)*(1.0-T_max)));
                    dx = dx_scale + ((nrandom()-0.5) * dx_scale);
                    k += dx;

                    //Transmittance in Interval
                    T.set(std::exp(-sigma_t.x*dx) * T.x , std::exp(-sigma_t.y*dx) * T.y , std::exp(-sigma_t.z*dx) * T.z);

                    //Direct Radiance
                    direct_radiance.set( (T.x * sigma_a.x * L_e.x * dx) , (T.y * sigma_a.y * L_e.y * dx) , (T.z * sigma_a.z * L_e.z * dx) );

                    //Indirect Single Radiance
                    m_sandbox_single_scatter_ray_MFP(L_s,x1,w0);
                    indirect_single_radiance.set( (T.x * sigma_s.x * L_s.x * dx) , (T.y * sigma_s.y * L_s.y * dx) , (T.z * sigma_s.z * L_s.z * dx) );

                    //Total Radiance 
                    L.set(L.x + direct_radiance.x + indirect_single_radiance.x , L.y + direct_radiance.y + indirect_single_radiance.y , L.z + direct_radiance.z + indirect_single_radiance.z);

                    x1 = x0 + w*k;
                }
            }

            //evaluate single scatter event
            void m_sandbox_single_scatter_ray_MFP(Vector3 &L_s, const Vector3 &x_s, const Vector3 &w_s)
            {
                //L_s(x_t,w) = integral_S( fp(x,w,w´)*L(x,w)*dw )
                
                bool track_back = false;
                realNumber cosTheta = 0.0, s = 0.0, k = 0.0, pdf = 0.0, ds = 0.0, ds_base = 0.15, ds_scale = 0.0, T_max = 1.0;
                Vector3 x_t = x_s, x_e(0.0,0.0,0.0), w_e(0.0,0.0,0.0), light_col(0.0,0.0,0.0); //_e means "end", whereas _s means "start"
                //m_point_light(x_e,light_col);
                m_directional_light(x_s,x_e,light_col);

                //Compute Scatter direction (w_e) and Angle (cosTheta)
                w_e = x_e - x_s;
                w_e = w_e / VectorLength(w_e);
                cosTheta = dotProduct(w_e,w_s);

                //Compute travel length and steps
                s = euclidianDistance(x_e,x_s);

                //Declare scatter variables
                Vector3 T(1.0,1.0,1.0), T_estimate(1.0,1.0,1.0), simga_t_old(0.0,0.0,0.0), sigma_a(0.0,0.0,0.0), sigma_s(0.0,0.0,0.0), sigma_t(0.0,0.0,0.0), L_e(0.0,0.0,0.0), direct_radiance(0.0,0.0,0.0);

                //Scatter
                while(s >= k)
                {
                    //Volume Sample
                    m_sandbox_volume(x_t,sigma_a,sigma_s,sigma_t,L_e);

                    //Determine dx
                    T_max = std::max(std::max(T.x,T.y),T.z);
                    ds_scale = (ds_base + ((1.0-T_max)*(1.0-T_max)*(1.0-T_max)));
                    ds = ds_scale + ((nrandom()-0.5) * ds_scale);
                    k += ds;

                    //Transmittance along ray
                    T.set(std::exp(-sigma_t.x*ds) * T.x , std::exp(-sigma_t.y*ds) * T.y , std::exp(-sigma_t.z*ds) * T.z);

                    //Accumilate Direct Radiance
                    direct_radiance.set( direct_radiance.x + (T.x * sigma_a.x * L_e.x * ds) , direct_radiance.y + (T.y * sigma_a.y * L_e.y * ds) , direct_radiance.z + (T.z * sigma_a.z * L_e.z * ds) );

                    x_t = x_s + w_e*k;
                }

                //Add Directional light
                direct_radiance.set( direct_radiance.x + (light_col.x*T.x) , direct_radiance.y + (light_col.y*T.y), direct_radiance.z + (light_col.z*T.z) );

                //m_Isotropic_Phase_BSDF(pdf);
                m_Anisotropic_Phase_BSDF(pdf,0.7,cosTheta);
                L_s = direct_radiance*pdf;
            }

            //Transmittance Estimator
            void m_sandbox_estimate_transmittance(Vector3 &T, const Vector3 &x_s, const Vector3 &x_e, const Vector3 &w0)
            {
                realNumber s = euclidianDistance(x_s,x_e), dx_scale = 1.0, dx = 0.0, k = 0.0;
                Vector3 x0 = x_s, x1(0.0,0.0,0.0), sigma_a(0.0,0.0,0.0), sigma_s(0.0,0.0,0.0), sigma_t(0.0,0.0,0.0), L_e(0.0,0.0,0.0);

                while(s >= k)
                {
                    //Volume Sample
                    m_sandbox_volume(x1,sigma_a,sigma_s,sigma_t,L_e);

                    //Step size
                    dx = ((nrandom() - 0.5) * (dx_scale/2.0)) + dx_scale;
                    k += dx;

                    //Transmittance along ray
                    T.set(std::exp(-sigma_t.x*dx) * T.x , std::exp(-sigma_t.y*dx) * T.y , std::exp(-sigma_t.z*dx) * T.z);

                    x1 = x0 + w0*k;
                }
            }


        //Mandelbulb
        int m_sandbox_mandelbulb(const Vector3 &p, const int &imax)
        {
            realNumber x = p.x, y = p.y, z = p.z, xnew, ynew, znew, n = 4.0, r, theta, phi;
            int i;
            
            for(i = 0; i < imax; i++)
            {
                r = std::sqrt(x*x+y*y+z*z);
                theta = std::atan2(std::sqrt(x*x+y*y),z);
                phi = std::atan2(y,x);

                xnew = std::pow(r,n) * std::sin(theta*n) * std::cos(phi*n) + p.x;
                ynew = std::pow(r,n) * std::sin(theta*n) * std::sin(phi*n) + p.y;
                znew = std::pow(r,n) * cos(theta*n) + p.z;

                if(xnew*xnew + ynew*ynew + znew*znew > imax)
                {
                    return i;
                }

                x = xnew;
                y = ynew;
                z = znew;
            }

            return imax;
        }

        //Point Light
        void m_point_light(Vector3 &pos, Vector3 &col)
        {
            //Set Color and Intensity
            realNumber intensity = 24.0;
            col.set(1.0,1.0,1.0);
            col = col*intensity;

            //Set Pos
            pos.set(-9.0,9.0,8.0);
        }

        //Directional light
        void m_directional_light(const Vector3 &x_s, Vector3 &x_e, Vector3 &col)
        {
            //Light Direction
            Vector3 
                w_light(-0.8,-0.3,0.0),
                col_light(1.0,1.0,1.0);
            
            realNumber
                scale_light = 24.0;
            
            w_light = w_light / VectorLength(w_light);

            //Both Intersections
            Vector3
                AABBmax(10.0,10.0,10.0),
                AABBmin(-10.0,-10.0,-10.0),
                hit1(0.0,0.0,0.0),
                hit0(0.0,0.0,0.0),
                w0(0.0,0.0,0.0),
                w1(0.0,0.0,0.0);

            intersectAABB(x_s,w_light,AABBmax,AABBmin,hit1,hit0);

            //Directions
            w0 = hit1 - x_s;
            w0 = w0 / VectorLength(w0);
            w1 = hit0 - x_s;
            w1 = w1 / VectorLength(w1);

            if(dotProduct(w0,w_light) < 0.0)
            {
                x_e = hit1;
            }
            else if(dotProduct(w1,w_light) < 0.0)
            {
                x_e = hit0;
            }

            //Col
            col = col_light*scale_light;
        }

        //Sandbox Volume
        void m_sandbox_volume(const Vector3 &x_k, Vector3 &sigma_a, Vector3 &sigma_s, Vector3 &sigma_t, Vector3 &L_e)
        {
            int mandelbulb = m_sandbox_mandelbulb(x_k*0.15,8);
            realNumber r = VectorLength(x_k);
            Vector3 scatter,absorption;

            //Handel Emissive Media
            if(r < 18.0)
            {
                Vector3 emissionColor(1.0,1.0,1.0);
                L_e = emissionColor*std::exp(-r*0.25)*0.0;
            }

            //Handel Absorption and Scattering
            if(mandelbulb < 8)
            {
                //Outside Mandelbulb
                absorption.set(0.01,0.015,0.02);
                scatter.set(0.015,0.02,0.025);
                sigma_a = absorption;
                sigma_s = scatter;
                sigma_t = sigma_a + sigma_s;
            }
            else
            {
                //Inside Mandelbulb
                absorption.set(0.2,0.3,0.4);
                scatter.set(0.95,0.95,0.95);
                sigma_a = absorption;
                sigma_s = scatter;
                sigma_t = sigma_a + sigma_s;
            }

            /*realNumber j = 5.0;

            absorption.set(0.1,0.2,0.3);
            scatter.set(0.2,0.5,0.8);

            if(r < j)
            {
                sigma_a = absorption * 1.0;
                sigma_s = scatter * 9.0;
                sigma_t = sigma_a + sigma_s;
            }

            if(r >= j)
            {
                realNumber scale = std::exp(-(r-j)*3.0);
                sigma_a = absorption * scale;
                sigma_s = scatter * scale;
                sigma_t = sigma_a + sigma_s;
            }*/
        }
    
    
    //General Functions
        //Mean Free Path Length
        void m_Mean_Free_Path_dx(realNumber &dx, const Vector3& sigma_t, const realNumber& min_dx, const realNumber& max_dx)
        {
            realNumber MFP_dx = 1.0 / (VectorLength(sigma_t) + (2.0 * (nrandom()-0.5) * (max_dx/2.0)));

            if(MFP_dx > max_dx)
            {
                dx = max_dx;
            }
            else if(MFP_dx < min_dx)
            {
                dx = min_dx;
            }
            else 
            {
                dx = MFP_dx;
            }
        }
    
        //Isotropic Phase BSDF
        void m_Isotropic_Phase_BSDF(realNumber &pdf)
        {
            pdf = 0.25 / M_PI;
        }

        //Anisotropic Phase BSDF
        void m_Anisotropic_Phase_BSDF(realNumber &pdf, const realNumber &g, const realNumber & costheta)
        {
            pdf = (1.0/(4.0*M_PI)) * ((1.0-g*g)/(std::pow(1.0+g*g-2.0*g*(costheta),3.0/2.0)));
        }

        //Normalized random direction
        //static Vector3 m_norm_rand_dir()
        //{

        //}
}


namespace Magik
{
    //Main
        //Non Relativistic - All Sandbox
        void trivial(Vector4 &RGBA_out, const Vector3 &rayOrig, const Vector3 &rayDir)
        {
            //Do Trivial
            Vector3
            AABBmax(10.0,10.0,10.0),
            AABBmin(-10.0,-10.0,-10.0),
            hit1(0.0,0.0,0.0),
            hit0(0.0,0.0,0.0),
            backgroundCol(0.0,0.0,0.0);
    
            if(intersectAABB(rayOrig,rayDir,AABBmax,AABBmin,hit1,hit0))
            {        
                Vector3 L(0.0,0.0,0.0);
                m_sandbox_trivial(L,hit0,hit1,0.1);
                backgroundCol.set(backgroundCol.x + L.x , backgroundCol.y + L.y , backgroundCol.z + L.z);
            }
        
            RGBA_out.set(backgroundCol.x,backgroundCol.y,backgroundCol.z,0.0);
        }

        void single(Vector4 &RGBA_out, const Vector3 &rayOrig, const Vector3 &rayDir)
        {
            //Do Single
            Vector3
            AABBmax(10.0,10.0,10.0),
            AABBmin(-10.0,-10.0,-10.0),
            hit1(0.0,0.0,0.0),
            hit0(0.0,0.0,0.0),
            backgroundCol(0.0,0.0,0.0);
    
            if(intersectAABB(rayOrig,rayDir,AABBmax,AABBmin,hit1,hit0))
            {        
                Vector3 L(0.0,0.0,0.0);
                m_sandbox_single(L,hit0,hit1,rayDir,0.1);
                backgroundCol.set(backgroundCol.x + L.x , backgroundCol.y + L.y , backgroundCol.z + L.z);
            }
        
            RGBA_out.set(backgroundCol.x,backgroundCol.y,backgroundCol.z,0.0);
        }

        void single_MFP(Vector4 &RGBA_out, const Vector3 &rayOrig, const Vector3 &rayDir)
        {
            //Do Single - Mean Free Path
            Vector3
            AABBmax(10.0,10.0,10.0),
            AABBmin(-10.0,-10.0,-10.0),
            hit1(0.0,0.0,0.0),
            hit0(0.0,0.0,0.0),
            backgroundCol(0.0,0.0,0.0);
    
            if(intersectAABB(rayOrig,rayDir,AABBmax,AABBmin,hit1,hit0))
            {        
                Vector3 L(0.0,0.0,0.0);
                //m_sandbox_estimate_transmittance(L,hit0,hit1,rayDir);
                m_sandbox_single_MFP(L,hit0,hit1,rayDir);
                backgroundCol.set(backgroundCol.x + L.x , backgroundCol.y + L.y , backgroundCol.z + L.z);
            }
        
            RGBA_out.set(backgroundCol.x,backgroundCol.y,backgroundCol.z,0.0);
        }

        void multiple()
        {
            //Do Multiple
        }


        //Relativisitc
        void trivial_relativistic()
        {
            //Do Trivial but relativisitc
        }

        void single_relativistic()
        {
            //Do Single but relativisitc
        }

        void multiple_relativistic()
        {
            //Do Multiple but relativisitc
        }
}
