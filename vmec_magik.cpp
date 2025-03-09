#include "vmec_common.h"
#include "vmec_magik.h"
#include "vmec_magik.h"

namespace vmec
{
    //Sandbox Scene Specific functions
        //trivial sandbox - direct radiance only
        void m_sandbox_trivial(Vector3 &L, const Vector3 &p0, const Vector3 &p1, const realNumber &dx)
        {
            realNumber s = euclidianDistance(p0,p1), k = 0.0;
            Vector3 x0 = p0, w = (p1-p0) / VectorLength(p1-p0), x1 = p0, BL,floating_majorant;
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
                m_sandbox_volume(x1,sigma_a,sigma_s,sigma_t,L_e,floating_majorant);

                //Transmittance in Interval
                T.set(std::exp(-sigma_t.x*dx) * T.x , std::exp(-sigma_t.y*dx) * T.y , std::exp(-sigma_t.z*dx) * T.z);

                //Direct Light
                L.set(L.x + (T.x * sigma_a.x * L_e.x * dx) , L.y + (T.y * sigma_a.y * L_e.y * dx) , L.z + (T.z * sigma_a.z * L_e.z * dx));

                x1 = x0 + w*k;
            }
        }

        //Single scattering (direct + single indirect radiance) - Mean Free Path dx
            //single scattering
            void m_sandbox_single(Vector3 &L, const Vector3 &p0, const Vector3 &p1, const Vector3 &w0)
            {
                bool track_back = false;
                realNumber s = euclidianDistance(p0,p1), k = 0.0, dx = 0.15;
                Vector3 x0 = p0, w = (p1-p0) / VectorLength(p1-p0), x1 = p0, BL,floating_majorant;

                Vector3 T(1.0,1.0,1.0), T_estimate(1.0,1.0,1.0), simga_t_old(0.0,0.0,0.0), sigma_a(0.0,0.0,0.0), sigma_s(0.0,0.0,0.0), sigma_t(0.0,0.0,0.0), L_e(0.0,0.0,0.0), L_s(0.0,0.0,0.0), direct_radiance(0.0,0.0,0.0), indirect_single_radiance(0.0,0.0,0.0);

                while(s >= k)
                {
                    //Volume Sample
                    m_sandbox_volume(x1,sigma_a,sigma_s,sigma_t,L_e,floating_majorant);
                    
                    //Parameter along path
                    k += (dx + (nrandom()-0.5)*(dx/4.0));

                    //Transmittance in Interval
                    T.set(std::exp(-sigma_t.x*dx) * T.x , std::exp(-sigma_t.y*dx) * T.y , std::exp(-sigma_t.z*dx) * T.z);

                    //Direct Radiance
                    direct_radiance.set( (T.x * sigma_a.x * L_e.x * dx) , (T.y * sigma_a.y * L_e.y * dx) , (T.z * sigma_a.z * L_e.z * dx) );

                    //Indirect Single Radiance
                    m_sandbox_single_scatter_ray(L_s,x1,w0);
                    indirect_single_radiance.set( (T.x * sigma_s.x * L_s.x * dx) , (T.y * sigma_s.y * L_s.y * dx) , (T.z * sigma_s.z * L_s.z * dx) );

                    //Total Radiance 
                    L.set(L.x + direct_radiance.x + indirect_single_radiance.x , L.y + direct_radiance.y + indirect_single_radiance.y , L.z + direct_radiance.z + indirect_single_radiance.z);

                    x1 = x0 + w*k;
                }
            }

            //evaluate single scatter event
            void m_sandbox_single_scatter_ray(Vector3 &L_s, const Vector3 &x_s, const Vector3 &w_s)
            {
                //L_s(x_t,w) = integral_S( fp(x,w,w´)*L(x,w)*dw )
                
                bool track_back = false;
                realNumber cosTheta = 0.0, s = 0.0, k = 0.0, pdf = 0.0, ds = 0.15;
                Vector3 x_t = x_s, x_e(0.0,0.0,0.0), w_e(0.0,0.0,0.0), light_col(0.0,0.0,0.0),floating_majorant; //_e means "end", whereas _s means "start"
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
                    m_sandbox_volume(x_t,sigma_a,sigma_s,sigma_t,L_e,floating_majorant);

                    //Parameter along Path
                    k += (ds + (nrandom()-0.5)*(ds/4.0));

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
                //m_Mie_Phase_BSDF(pdf,5.0,cosTheta);
                L_s = direct_radiance*pdf;
            }
        
        //Multiple Scattering - Null Collision
            //Null Collision
            void m_sandbox_trace_null_collision_ray(Vector3 &L, const Vector3 &x0, const Vector3 &w0, const Vector3 &AABB_min, const Vector3 &AABB_max, realNumber &running_majorant)
            {
                for(int i = 0; i <= 2; i++)
                {
                    int depth = 0;
                    realNumber xi, t, multiplier = 1.0, P_e, P_s, P_n, majorant = 1.05, g = 0.7, pdf = 1.0, costheta;
                    Vector3 xt = x0, wt = w0, w_temp, sigma_a, sigma_s, sigma_t, sigma_n, L_e, L_s(0.0,0.0,0.0), floating_majorant;
    
                    while(settings.magik_sandbox_multiple_scatter_events > depth)
                    {
                        //Random Number
                        xi = nrandom();
    
                        //Parameter along Path from distribution
                        t = -std::log(1.0-xi)/majorant;
    
                        //Realization sample position
                        xt = xt + wt*t;
    
                        if(m_AABB_In_Out(xt,AABB_min,AABB_max))
                        {
                            break;
                        }
    
                        //Volume sample
                        m_sandbox_volume(xt,sigma_a,sigma_s,sigma_t,L_e,floating_majorant);
                        majorant = 1.05; //floating_majorant[i];
                        sigma_n[i] = majorant - sigma_t[i];
    
                        //Probabilities
                        P_e = sigma_a[i]/majorant;
                        P_s = sigma_s[i]/majorant;
                        P_n = sigma_n[i]/majorant;
    
                        //New random number
                        xi = nrandom();
    
                        //Determine next event
                        if(xi < P_e)
                        {
                            L_s[i] = L_s[i] + (L_e[i] * (sigma_a[i]/majorant) * multiplier);
                            break;
                        }
                        else if(xi < 1.0 - P_n)
                        {
                            depth++;
                            w_temp = wt;
                            wt = m_norm_HG_dir(nrandom(),g,wt);
                            costheta = dotProduct(w_temp,wt);
                            m_Anisotropic_Phase_BSDF(pdf,g,costheta);
                            multiplier *= sigma_s[i]/majorant * pdf;
                        }
                    }
    
                    L[i] = L[i] + (L_s[i]*(1.0/realNumber(settings.magik_sandbox_monte_carlo_samples)));
                }
            }

        //Mandelbulb
        int m_sandbox_mandelbulb(Vector3 &out, const Vector3 &p, const int &imax)
        {
            Matrix3 m = rotateRay(0.0,0.0,-90.0);
            Vector3 x_rot = p;
            x_rot = m*x_rot;
            realNumber x = x_rot.x, y = x_rot.y, z = x_rot.z, xnew, ynew, znew, n = 4.0, r, theta, phi;
            int i;
            
            for(i = 0; i < imax; i++)
            {
                r = std::sqrt(x*x+y*y+z*z);
                theta = std::atan2(std::sqrt(x*x+y*y),z);
                phi = std::atan2(y,x);

                xnew = std::pow(r,n) * std::sin(theta*n) * std::cos(phi*n) + x_rot.x;
                ynew = std::pow(r,n) * std::sin(theta*n) * std::sin(phi*n) + x_rot.y;
                znew = std::pow(r,n) * cos(theta*n) + x_rot.z;

                if(xnew*xnew + ynew*ynew + znew*znew > imax)
                {
                    return i;
                    out.set(xnew,ynew,znew);
                }

                x = xnew;
                y = ynew;
                z = znew;
            }

            return imax;
            out.set(xnew,ynew,znew);
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
                w_light(-1.0,-1.0,-0.1),
                col_light(1.0,1.0,1.0);
            
            realNumber
                scale_light = 512.0;
            
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
        void m_sandbox_volume(const Vector3 &x_k, Vector3 &sigma_a, Vector3 &sigma_s, Vector3 &sigma_t, Vector3 &L_e, Vector3 &floating_majorant)
        {
            //Vector3 m_out;
            //int mandelbulb = m_sandbox_mandelbulb(m_out,x_k*0.15,8);
            
            //Blender Shader Copy
            //Emissive spheres
            Vector3
                pos_blue(1.5,-2.2,-5.0),
                pos_red(-6.0,5.22,0.0),
                col_blue(0.0,0.0,1.0*1560000.0),
                col_red(1560000.0,1560000.0,1560000.0);

            realNumber
                r_blue = 2.0,
                r_red = 2.0;

            //Blue
            if(euclidianDistance(pos_blue,x_k) <= r_blue)
            {
                //L_e = col_blue;
            }

            //Red
            if(euclidianDistance(pos_red,x_k) <= r_red)
            {
                L_e = col_red;
            }

            //Density
            Vector3
                absorption(0.01,0.01,0.01),
                scatter(0.05,0.05,0.05);
                sigma_a = absorption;
                sigma_s = scatter;
                sigma_t = sigma_a + sigma_s;
                floating_majorant = scatter;

            realNumber
                r_sphere = 5.0;

            if(VectorLength(x_k) < r_sphere) //Mandelbulb !(mandelbulb < 8)
            {
                absorption.set(0.01,0.01,0.01);
                scatter.set(0.305,0.51,1.02);
                sigma_a = absorption;
                sigma_s = scatter;
                sigma_t = sigma_a + sigma_s;
                floating_majorant = scatter;
            }
        }

        //In-Out check
        bool m_AABB_In_Out(const Vector3 &p0, const Vector3 &AABB_min, const Vector3 &AABB_max)
        {
            if(((p0.x > AABB_max.x) || (p0.y > AABB_max.y) || (p0.z > AABB_max.z)) || ((p0.x < AABB_min.x) || (p0.y < AABB_min.y) || (p0.z < AABB_min.z)))
            {
                return true; //Inisde
            }

            return false; //Outside
        }

    //General Functions
        //Unified Volume Function
        void m_invoke_unified_volume(const Vector3 x_k, const BVObjects& BV, Vector3 &sigma_a, Vector3 &sigma_s, Vector3 &sigma_t, Vector3 &L_e)
        {
            
        }

        //Isotropic Phase BSDF
        void m_Isotropic_Phase_BSDF(realNumber &pdf)
        {
            pdf = 0.25 / M_PI;
        }

        //Anisotropic Phase BSDF
        void m_Anisotropic_Phase_BSDF(realNumber &pdf, const realNumber &g, const realNumber &costheta)
        {
            pdf = (1.0/(4.0*M_PI)) * ((1.0-(g*g))/(std::pow(1.0+(g*g)-(2.0*g)*(costheta),3.0/2.0)));
        }

        //Drains Phase BSDF
        void m_Drains_Phase_BSDF(realNumber &pdf, const realNumber &alpha, const realNumber& g, const realNumber &costheta)
        {
            realNumber pdf1;
            m_Anisotropic_Phase_BSDF(pdf1,g,costheta);
            pdf = (0.25/M_PI) * (pdf1) * ((1.0+alpha*(costheta*costheta))/(1+((alpha*(1.0+2.0*g*g))/3.0)));
        }

        //Mie scattering approximation
        void m_Mie_Phase_BSDF(realNumber &pdf, const realNumber &d, const realNumber &costheta)
        {
            realNumber g_HG, g_D, alpha_D, w_D, Drains1, Drains2, d_scale;
            d_scale = d/(1000000.0);
            g_HG = std::exp(-(0.0990567/(d_scale-1.67154)));
            g_D = std::exp(-(2.20679/(d_scale+3.91029))-0.428934);
            alpha_D = std::exp(3.62489 - (8.29288/(d_scale+5.52825)));
            w_D = std::exp(-(0.599085/(d_scale-0.641583))-0.665888);

            m_Drains_Phase_BSDF(Drains1,0.0,g_HG,costheta);
            m_Drains_Phase_BSDF(Drains2,alpha_D,g_D,costheta);
            pdf = (1.0-w_D)*Drains1 + w_D*Drains2;
        }

        //Normalized random direction
        Vector3 m_norm_rand_dir()
        {
            realNumber x,y,z,phi,proj;
            Vector3 out;
            
            z = (nrandom()-0.5)*2.0;
            phi = nrandom()*M_PI*2.0;
            proj = std::sqrt(1.0-z*z);
            x = proj*std::cos(phi);
            y = proj*std::sin(phi);

            out.set(x,y,z);

            return out;
        }

        //HG Importance Sampler
        Vector3 m_norm_HG_dir(const realNumber &epsilon, const realNumber &g, const Vector3 &w0)
        {
            realNumber cosTheta, theta, phi, sphi;
            Vector3 w1;
            
            cosTheta = (1.0/(2.0*g)) * ( 1+g*g - ( ((1-g*g)/(1-g+2*g*epsilon))*((1-g*g)/(1-g+2*g*epsilon)) ) );
            theta = std::acos(cosTheta);
            phi = nrandom()*2.0*M_PI;
            sphi = std::sqrt(1.0-w0.y*w0.y);

            w1.x = w0.x*cosTheta + ((std::sin(theta)*(w0.x*w0.y*std::cos(phi)-w0.z*std::sin(phi)))/sphi);
            w1.z = w0.z*cosTheta + ((std::sin(theta)*(w0.z*w0.y*std::cos(phi)+w0.x*std::sin(phi)))/sphi);
            w1.y = w0.y*cosTheta - sphi * std::cos(theta)*std::cos(phi);

            w1 = w1 / VectorLength(w1);

            return w1;
        }
}



namespace Magik
{
    //Fonts; Accolade-Medium for MAGIK and Heebo for VMEC in general
    
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
            backgroundCol(0.01,0.01,0.01);
    
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
                m_sandbox_single(L,hit0,hit1,rayDir);
                backgroundCol.set(backgroundCol.x + L.x , backgroundCol.y + L.y , backgroundCol.z + L.z);
            }
        
            RGBA_out.set(backgroundCol.x,backgroundCol.y,backgroundCol.z,0.0);
        }

        void multiple(Vector4 &RGBA_out, const Vector3 &rayOrig, const Vector3 &rayDir)
        {
            //Do Multiple
            Vector3
            AABBmax(10.0,10.0,10.0),
            AABBmin(-10.0,-10.0,-10.0),
            hit1(0.0,0.0,0.0),
            hit0(0.0,0.0,0.0),
            backgroundCol(0.0,0.0,0.0);

            realNumber
            initial_majorant = 0.1;
    
            if(intersectAABB(rayOrig,rayDir,AABBmax,AABBmin,hit1,hit0))
            {        
                Vector3 L(0.0,0.0,0.0);
                for(int i = 0; i < settings.magik_sandbox_monte_carlo_samples; i++)
                {
                    m_sandbox_trace_null_collision_ray(L,hit0,rayDir,AABBmin,AABBmax,initial_majorant);
                }
                backgroundCol.set(backgroundCol.x + L.x , backgroundCol.y + L.y , backgroundCol.z + L.z);
            }
        
            RGBA_out.set(backgroundCol.x,backgroundCol.y,backgroundCol.z,0.0);
        }


        //Relativisitc
        void trivial_relativistic()
        {
            //Do Trivial but relativisitc
        }

        void multiple_relativistic()
        {
            //Do Multiple but relativisitc
        }
}