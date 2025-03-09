#include "vmec_common.h"

namespace vmec
{

//Metric Tensor Components
    realNumber metric_delta(const realNumber r)
    {
        return (r*r) - (2.0*r) + (settings.a*settings.a);
    }


    realNumber metric_sigma(const realNumber r, const realNumber theta)
    {
        return (r*r) + (settings.a*settings.a) * (std::cos(theta) * std::cos(theta));
    }


    realNumber metric_lambda(const realNumber r, const realNumber theta)
    {
        return (((r*r) + (settings.a*settings.a))*((r*r) + (settings.a*settings.a))) - (settings.a*settings.a * metric_delta(r) * (std::sin(theta)*std::sin(theta)));
    }


    realNumber metric_gPhiPhi(const realNumber r, const realNumber theta)
    {
        return (metric_lambda(r,theta)/metric_sigma(r,theta)) * std::sin(theta)*std::sin(theta);
    }


    realNumber metric_gPhi(const realNumber r, const realNumber theta)
    {
        return -((2.0*settings.a*r) / metric_sigma(r,theta)) * std::sin(theta)*std::sin(theta);
    }


    realNumber metric_gtt(const realNumber r, const realNumber theta)
    {
        return -(1.0-((2.0*r) / metric_sigma(r,theta)));
    }


    realNumber metric_grr(const realNumber r, const realNumber theta)
    {
        return metric_sigma(r,theta) / metric_delta(r);
    }


    realNumber metric_gThetaTheta(const realNumber r, const realNumber theta)
    {
        return metric_sigma(r,theta);
    }


    realNumber metric_gtPhi(const realNumber r, const realNumber theta)
    {
        return -((2.0*settings.a*r) / metric_sigma(r,theta)) * (std::sin(theta)*std::sin(theta));
    }


//Ray-Disk Intersection
bool rayDiskIntersection(const Vector3 &n, const Vector3 &p0, const Vector3 &rayOrig, const Vector3 &rayDir, const realNumber &dx, Vector3 &hit)
{
    realNumber denom = dotProduct(n,rayDir);
    Vector3 p010 = p0 - rayOrig;
    realNumber t = dotProduct(p010,n) / denom;
    if(t >= 0.0 && t <= dx)
    {
        Vector3 p = rayOrig + rayDir * t; //Intersection point
        Vector3 v = p - p0;
        realNumber d2 = dotProduct(v,v);
        hit = p;
        return (d2 <= settings.checker_disk_major_radius*settings.checker_disk_major_radius && d2 >= settings.checker_disk_minor_radius*settings.checker_disk_minor_radius);
    }
    return false;
}


//Ray-Plane Intersection
bool rayPlaneIntersection(const Vector3 &n, const Vector3 &p0, const Vector3 &rayOrig, const Vector3 &rayDir, const realNumber &dx, Vector3 &hit)
{
    realNumber denom = dotProduct(n,rayDir);
    Vector3 p010 = p0 - rayOrig;
    realNumber t = dotProduct(p010,n) / denom;
    if(t >= 0.0 && t <= dx)
    {
        hit = rayOrig + rayDir * t; //Intersection point
        return true;
    }
    return false;
}


//Sub Grid
void axis_aligned_gird(const Vector3 &rayOrig, const Vector3 &hit, Vector3 &col, const realNumber &grid_spacing, bool &hit_grid)
{
    bool is_grid = false, is_X = false, is_Z = false;
    Vector3 modPos(std::abs(((std::fmod(std::abs(hit.x),grid_spacing)*(1.0/grid_spacing))-0.5)*2.0) , std::abs(((std::fmod(std::abs(hit.y),grid_spacing)*(1.0/grid_spacing))-0.5)*2.0) , std::abs(((std::fmod(std::abs(hit.z),grid_spacing)*(1.0/grid_spacing))-0.5)*2.0));
    realNumber fade = 1.0, r = VectorLength(hit), grid_width = 0.0015, scale_width = (1.0-((grid_width*euclidianDistance(rayOrig,hit))/grid_spacing));

    if(modPos.x > scale_width)
    {
        is_grid = true;
        if(std::abs(hit.x) < scale_width) { is_X = true; }
    }

    if(modPos.z > scale_width)
    {
        is_grid = true;
        if(std::abs(hit.z) < scale_width) { is_Z = true; }
    }

    if(is_grid)
    {
        hit_grid = true;

        if(r >= (settings.celestial_sphere_radius*0.75))
        {
            fade = std::exp(-(r-settings.celestial_sphere_radius*0.75)*0.01);
        }

        if(is_X)
        {
            col.set(0.95*fade,0.1*fade,0.1*fade);
        }
        else if(is_Z)
        {
            col.set(0.1*fade,0.95*fade,0.1*fade);
        }
        else
        {
            col.set(0.05*fade,0.05*fade,0.05*fade);
        }
    }
}


//XZ Gird
void XZ_world_gird(const Vector3 &rayOrig, const Vector3 &hit, Vector3 &col, bool &hit_grid)
{
    //Check main grid
    axis_aligned_gird(rayOrig,hit,col,settings.grid_spacing_major,hit_grid);
    
    //Check minor grid
    if(!hit_grid)
    {
        axis_aligned_gird(rayOrig,hit,col,settings.grid_spacing_minor,hit_grid);
        col = col * std::exp(-std::abs(rayOrig.y)*0.015)*0.5;
    }
}


//Ray-Sphere Intersection
bool raySphereIntersection(const Vector3 &rayOrig, const Vector3 &rayDir, const Vector3 &pos, const realNumber &radius, const realNumber &dx)
{
    Vector3 L = pos - rayOrig;
    realNumber tca = dotProduct(L,rayDir);
    realNumber d2 = dotProduct(L,L) - tca*tca;
    if(d2 > radius*radius) return false;
    realNumber thc = std::sqrt((radius*radius) - d2);
    realNumber t0 = 0.0, t1 = 0.0;
    t0 = tca - thc;
    t1 = tca + thc;
    if(t0 > t1) std::swap(t0,t1);
    if(t0 < 0.0)
    {
        t0 = t1;
        if(t0 < 0.0) return false;
    }

    if(t0 > dx) return false;

    //t0 is now the distance
    return true;
}


//ray-box intersection
bool intersectAABB(const Vector3 &rayOrig, const Vector3 &rayDir, const Vector3 &max, const Vector3 &min, Vector3 &hitIn, Vector3 &hitOut)
{
    realNumber tmin = (min.x - rayOrig.x) / rayDir.x;
    realNumber tmax = (max.x - rayOrig.x) / rayDir.x;

    if(tmin > tmax) std::swap(tmin,tmax);

    realNumber tymin = (min.y - rayOrig.y) / rayDir.y;
    realNumber tymax = (max.y - rayOrig.y) / rayDir.y;

    if(tymin > tymax) std::swap(tymin, tymax);

    if((tmin > tymax) || (tymin > tmax))
        return false;

    if(tymin > tmin) tmin = tymin;
    if(tymax < tmax) tmax = tymax;

    realNumber tzmin = (min.z - rayOrig.z) / rayDir.z;
    realNumber tzmax = (max.z - rayOrig.z) / rayDir.z;

    if(tzmin > tzmax) std::swap(tzmin,tzmax);

    if((tmin > tzmax) || (tzmin > tmax))
        return false;

    if(tzmin > tmin) tmin = tzmin;
    if(tzmax < tmax) tmax = tzmax;

    hitIn = rayOrig + rayDir*tmax;
    hitOut = rayOrig + rayDir*tmin;

    return true;
}


//Camera proper time to global
void camera_proper_to_globa_time(RayProperties &pointProperties, bool &rs_hit)
{
    realNumber hStepSum = 0.0;
    
    while(settings.frame_time_step*settings.time_ramp > hStepSum && !rs_hit)
    {
        if(pointProperties.hStep > 1e-5) {pointProperties.hStep = 1e-5;}
        pointProperties = rayRKF45(pointProperties);
    
        if(isinf(pointProperties.hStep))
        {
            pointProperties.hStep = 1e-6;
        }
    
        hStepSum += pointProperties.hStep;
    
        if(pointProperties.r < settings.rs*settings.rs_scale*2.0)
        {
            rs_hit = true;
        }
    }
}


//RNG
    static std::default_random_engine rng;
    realNumber nrandom()
    {
        std::uniform_real_distribution<double> dist(0.0,1.0);
        return dist(rng);
    }


//Hash function
    unsigned int hashFunction(unsigned int seed) {
        seed ^= seed >> 16;
        seed *= 0x85ebca6b;
        seed ^= seed >> 13;
        seed *= 0xc2b2ae35;
        seed ^= seed >> 16;
        return seed;
    }


//Deterministic RNG
    realNumber ndetermin(const int seed)
    {
        std::default_random_engine generator(hashFunction(seed));
        std::uniform_real_distribution<double> dist(0.0,1.0);
        return dist(generator);
    }


//Sign function
    realNumber sign(const realNumber x)
    {
        if(x == 0 || x > 0)
        {
            return 1.0;
        }
        else
        {
            return -1.0;
        }
    }


//Linear interpolation Function
    realNumber lerp(const realNumber c, const realNumber b, const realNumber t)
    {
        return c+t*(b-c);
    }


//Distance between two points
    realNumber euclidianDistance(const Vector3& v1, const Vector3& v2)
    {
        return std::sqrt( ((v1.x - v2.x)*(v1.x - v2.x)) + ((v1.y - v2.y)*(v1.y - v2.y)) + ((v1.z - v2.z)*(v1.z - v2.z)) );
    }


//Distance between two points squared
    realNumber euclidianDistance2(const Vector3& v1, const Vector3& v2)
    {
        return ((v1.x - v2.x)*(v1.x - v2.x)) + ((v1.y - v2.y)*(v1.y - v2.y)) + ((v1.z - v2.z)*(v1.z - v2.z));
    }


//Rotate 2D
    Vector3 rotate2D(const Vector3& rayOrig, const realNumber angle)
    {
        Matrix2 m;
        Vector3 pRot;
        Vector2 A2;
        realNumber s, c;

        s = std::sin(angle);
        c = std::cos(angle);
        A2.x = rayOrig.x;
        A2.y = rayOrig.z;
        m.M00 = c;
        m.M01 = -s;
        m.M10 = s;
        m.M11 = c;
        A2 = m*A2;
        pRot.x = A2.x;
        pRot.y = rayOrig.y;
        pRot.z = A2.y;

        return pRot;
    }


//Rotate ray
//Note: input angles in degrees
//rotation is anticlockwise (right-hand rule), see usual definition on e.g. Wikipedia
    Matrix3 rotateRay(const realNumber rX, const realNumber rY, const realNumber rZ)
    {
        Matrix3 mX, mY, mZ, mRot;
        realNumber thetaX = rX * (M_PI/180.), thetaY = rY * (M_PI/180.), thetaZ = rZ * (M_PI/180.);

        mX.M00 = 1.0;
        mX.M01 = 0.0;
        mX.M02 = 0.0;
        mX.M10 = 0.0;
        mX.M11 = std::cos(thetaX);
        mX.M12 = -std::sin(thetaX);
        mX.M20 = 0.0;
        mX.M21 = std::sin(thetaX);
        mX.M22 = std::cos(thetaX);

        mY.M00 = std::cos(thetaY);
        mY.M01 = 0.0;
        mY.M02 = std::sin(thetaY);
        mY.M10 = 0.0;
        mY.M11 = 1.0;
        mY.M12 = 0.0;
        mY.M20 = -std::sin(thetaY);
        mY.M21 = 0.0;
        mY.M22 = std::cos(thetaY);

        mZ.M00 = std::cos(thetaZ);
        mZ.M01 = -std::sin(thetaZ);
        mZ.M02 = 0.0;
        mZ.M10 = std::sin(thetaZ);
        mZ.M11 = std::cos(thetaZ);
        mZ.M12 = 0.0;
        mZ.M20 = 0.0;
        mZ.M21 = 0.0;
        mZ.M22 = 1.0;
        mRot = mX*mY*mZ;

        return mRot;
    }


//Convert Spherical to Cartesian function
    Vector3 sphericalToCartesian(const realNumber r, const realNumber theta, const realNumber phi)
    {
        Vector3 out;
        out.x = r * std::sin(theta) * std::cos(phi); //X
        out.y = r * std::cos(theta); //Y
        out.z = r * std::sin(theta) * std::sin(phi); //Z
        return out;
    }


//Cartesian to Boyer-Lindquist Coordinates
    Vector3 cartesianToBL(const Vector3& rayOrig)
    {
        Vector3 out;
        realNumber w, r, theta, phi;

        w = (rayOrig.x*rayOrig.x + rayOrig.y*rayOrig.y + rayOrig.z*rayOrig.z) - settings.a*settings.a;
        r = std::sqrt(0.5 * (w + std::sqrt(w*w + (4.0f*settings.a*settings.a*rayOrig.y*rayOrig.y))));
        theta = std::acos(rayOrig.y / r);
        phi = std::atan2(rayOrig.z,rayOrig.x);

        out.set(r,theta,phi);
        return out;
    }


//BL coords to Cartesian
    Vector3 BoyerLindquistToCartesian(const realNumber r, const realNumber theta, const realNumber phi)
    {
        Vector3 out;
        out.x = std::sqrt((r*r) + (settings.a*settings.a)) * std::sin(theta) * std::cos(phi); //X
        out.y = r*std::cos(theta); //Y
        out.z = std::sqrt((r*r) + (settings.a*settings.a)) * std::sin(theta) * std::sin(phi); //Z
        return out;
    }


//Quadradic Bezier Function
    realNumber quadraticBezier(const realNumber x, const realNumber P0, const realNumber P1)
    {
        if(0 <= x && x <= 1.0)
        {
            realNumber term1 = 1.0-x;
            return term1*term1*P0 + 2.0*term1*x*P1 + x*x*P1;
        }
        else
        {
            return 0.0;
        }
    }


//Bezier Fade
    realNumber bezierFade(const realNumber x, const realNumber u1, const realNumber h, const realNumber N, const realNumber k)
    {
        if(x <= u1*k)
        {
            //Case 1; Constant function
            return h;
        }
        else
        {
            realNumber Args = 0.0, s = 0.0;
            s = (u1-u1*k)/2.0;

            //Case 2; Bezier part
            if(x <= (u1 + u1*k) / 2.0)
            {
                //Case 2a; Ease-In Bezier
                Args = (-x+s+u1*k)/s;
                return std::pow(quadraticBezier(Args,0.5,1.0),N)*h;
            }
            else
            {
                //Case 2b; Ease-Out Bezier
                Args = (x-u1+s)/s;
                return std::pow(quadraticBezier(Args,0.5,0.0),N)*h;
            }
        }
    }


//Falloff function
    realNumber gx(const realNumber N0, const realNumber N1, const realNumber s, const realNumber r, const realNumber x)
    {
        realNumber h = std::tan((s)*(M_PI/180.0))*r;
        return  std::pow(-std::pow(x,N0)*(1.0f/std::pow(h,N0))+1.0,N1);
    }


//Camera Initialization
    localProperties initializeCamera(const Vector3& camPos, const Vector3& camMomentum)
    {
        //Variables
        localProperties out;
        Vector3 blCoordsCam, cameraMomentum;
        realNumber tCam, rCam, thetaCam, phiCam, deltaCam, sigmaCam, u0Camera;

        //Convert Camera cartesian origin into Boyer-Lindqust coordinates and assign variables
        blCoordsCam = cartesianToBL(camPos);
        tCam = 0.0;
        rCam = blCoordsCam.x;
        thetaCam = blCoordsCam.y;
        phiCam = blCoordsCam.z;
        deltaCam = metric_delta(rCam), sigmaCam = metric_sigma(rCam,thetaCam);

        //2.3 Derive and assing camera´s 4-velocity u0 component
        u0Camera = deriveu0(camMomentum.x,camMomentum.y,camMomentum.z,rCam,thetaCam,deltaCam,sigmaCam,-1.0);

        //Load Output
        out.t = tCam;
        out.r = rCam;
        out.theta = thetaCam;
        out.phi = phiCam;
        out.u0 = u0Camera;
        out.u1 = camMomentum.x;
        out.u2 = camMomentum.y;
        out.u3 = camMomentum.z;
        return out;
    }


//Ray Initalization
RayProperties initializeRay(const Vector3& rayOrig, const Vector3& rayDir, const localProperties& cameraProperties)
{
    //2.1 Declare variables
    RayProperties out;
    Vector4 ZAMODir, photonMomentum;
    Vector3 blCoordsRay;
    realNumber tRay, rRay, thetaRay, phiRay, deltaRay, sigmaRay; //hStep and dx need to be initalized as a non 0 value

    //2.5 Convert cartesian ray origin to Boyer-Lindquist coordinates and assign variables
    blCoordsRay = cartesianToBL(rayOrig);
    tRay = 0.0;
    rRay = blCoordsRay.x;
    thetaRay = blCoordsRay.y;
    phiRay = blCoordsRay.z;
    deltaRay = metric_delta(rRay), sigmaRay = metric_sigma(rRay,thetaRay);

    //2.6 ray direction (and position), rest frame at infinity, to local Zero Angular Momentum Observer (ZAMO)
    ZAMODir = restframeToZAMO(rayDir * 0.05,cameraProperties); //0.05

    //2.7 ZAMO to global reference frame 4-Velocity vector and assign to initial
    photonMomentum = ZAMOToGlobal(ZAMODir,rRay,thetaRay,phiRay,deltaRay,sigmaRay);
    
    out.t = tRay;
    out.r = rRay;
    out.theta = thetaRay;
    out.phi = phiRay;
    out.u0 = photonMomentum.x;
    out.u1 = photonMomentum.y;
    out.u2 = photonMomentum.z;
    out.u3 = photonMomentum.t;
    out.hStep = 0.0;
    return out;
}


//In bound check
void pointInBounds(BVObjects& BoundingVolumes, const realNumber& r, const realNumber& theta)
{
    //Accretion disk Case
    if((settings.bv_disk_toggle || settings.vis_disk_BV) && std::sqrt(r*r + settings.a*settings.a)*std::sin(theta) < settings.bv_disk_major_radius && std::sqrt(r*r + settings.a*settings.a)*std::sin(theta) > settings.bv_disk_minor_radius && std::abs(r*std::cos(theta)) < std::tan(settings.bv_disk_mean_slope * (M_PI/180.0))*(r+settings.bv_disk_radius_offset))
    {
        BoundingVolumes.disk = true;
    }
    else
    {
        BoundingVolumes.disk = false;
    }

    //Astrophysical Jet case
    if((settings.bv_jet_toggle || settings.vis_jet_BV) && std::abs(r*std::cos(theta)) < settings.bv_jet_height && sqrt(r*r + settings.a*settings.a)*std::sin(theta)-settings.bv_jet_minor_radius < std::abs(r*std::cos(theta))*std::tan(settings.bv_jet_deflection*(M_PI/180.0)))
    {
        BoundingVolumes.jet = true;
    }
    else
    {
        BoundingVolumes.jet = false;
    }

    //Ambient Medium case
    if((settings.bv_ambient_medium_toggle || settings.vis_ambient_BV) && r < settings.bv_ambient_medium_radius)
    {
        BoundingVolumes.ambient = true;
    }
    else
    {
        BoundingVolumes.ambient = false;
    }
}


//Equations of motion for the Kerr metric
static RayProperties rayDerivatives(const RayProperties& incomingValues)
{
    RayProperties outgoingFunction;
    realNumber u0Dot = 0.0, u1Dot = 0.0, u2Dot = 0.0, u3Dot = 0.0;
    realNumber t = incomingValues.t;
    realNumber r = incomingValues.r;
    realNumber theta = incomingValues.theta;
    realNumber phi = incomingValues.phi;
    realNumber u0 = incomingValues.u0;
    realNumber u1 = incomingValues.u1;
    realNumber u2 = incomingValues.u2;
    realNumber u3 = incomingValues.u3;
    realNumber r2 = r*r, r3 = r2*r, r4 = r2*r2, a2 = settings.a*settings.a, a3 = a2*settings.a, a4 = a2*a2, cosTheta = std::cos(theta), sinTheta = std::sin(theta), cosTheta2 = cosTheta*cosTheta, sinTheta2 = sinTheta*sinTheta, sinTheta3 = sinTheta2*sinTheta, s = (a2*cosTheta2+r2), s2 = s*s;

    u0Dot = -2.0*u0*u1*(-settings.a*r*(4.0*settings.a*r2*sinTheta2/s2 - 2.0*settings.a*sinTheta2/s)/(2.0*a2*r*sinTheta2 - 2.0*a2*r - 2.0*r3 + a4*cosTheta2 + a2*r2*cosTheta2 + a2*r2 + r4) + 0.5*(-4.0*r2/s2 + 2.0/s)*(-2.0*a2*r*sinTheta2 - a4*cosTheta2 - a2*r2*cosTheta2 - a2*r2 - r4)/(2.0*a2*r*sinTheta2 - 2.0*a2*r - 2.0*r3 + a4*cosTheta2 + a2*r2*cosTheta2 + a2*r2 + r4)) - 2.0*u0*u2*(2.0*a2*r*(-2.0*a2*r*sinTheta2 - a4*cosTheta2 - a2*r2*cosTheta2 - a2*r2 - r4)*sinTheta*cosTheta/(s2*(2.0*a2*r*sinTheta2 - 2.0*a2*r - 2.0*r3 + a4*cosTheta2 + a2*r2*cosTheta2 + a2*r2 + r4)) - settings.a*r*(-4.0*a3*r*sinTheta3*cosTheta/s2 - 4.0*settings.a*r*sinTheta*cosTheta/s)/(2.0*a2*r*sinTheta2 - 2.0*a2*r - 2.0*r3 + a4*cosTheta2 + a2*r2*cosTheta2 + a2*r2 + r4)) - 2.0*u1*u3*(-settings.a*r*(-4.0*a2*r2*sinTheta2/s2 + 2.0*a2*sinTheta2/s + 2.0*r)*sinTheta2/(2.0*a2*r*sinTheta2 - 2.0*a2*r - 2.0*r3 + a4*cosTheta2 + a2*r2*cosTheta2 + a2*r2 + r4) + 0.5*(4.0*settings.a*r2*sinTheta2/s2 - 2.0*settings.a*sinTheta2/s)*(-2.0*a2*r*sinTheta2 - a4*cosTheta2 - a2*r2*cosTheta2 - a2*r2 - r4)/(2.0*a2*r*sinTheta2 - 2.0*a2*r - 2.0*r3 + a4*cosTheta2 + a2*r2*cosTheta2 + a2*r2 + r4)) - 2.0*u2*u3*(-settings.a*r*((4.0*a4*r*sinTheta3*cosTheta/s2 + 4.0*a2*r*sinTheta*cosTheta/s)*sinTheta2 + 2.0*(2.0*a2*r*sinTheta2/s + a2 + r2)*sinTheta*cosTheta)/(2.0*a2*r*sinTheta2 - 2.0*a2*r - 2.0*r3 + a4*cosTheta2 + a2*r2*cosTheta2 + a2*r2 + r4) + 0.5*(-4.0*a3*r*sinTheta3*cosTheta/s2 -4.0*settings.a*r*sinTheta*cosTheta/s)*(-2.0*a2*r*sinTheta2 - a4*cosTheta2 - a2*r2*cosTheta2 - a2*r2 - r4)/(2.0*a2*r*sinTheta2 - 2.0*a2*r - 2.0*r3 + a4*cosTheta2 + a2*r2*cosTheta2 + a2*r2 + r4));
    u1Dot = 2.0*a2*u1*u2*sinTheta*cosTheta/s + r*u2*u2*(-2.0*r + a2 + r2)/s - 0.5*u0*u0*(4.0*r2/s2 - 2.0/s)*(-2.0*r + a2 + r2)/s -u0*u3*(-4.0*settings.a*r2*sinTheta2/s2 + 2.0*settings.a*sinTheta2/s)*(-2.0*r + a2 + r2)/s - 0.5*u1*u1*(2.0*r/(-2.0*r + a2 + r2) + (2.0 - 2.0*r)*s/((-2.0*r + a2 + r2)*(-2.0*r + a2 + r2)))*(-2.0*r + a2 + r2)/s + 0.5*u3*u3*(-2.0*r + a2 + r2)*(-4.0*a2*r2*sinTheta2/s2 + 2.0*a2*sinTheta2/s + 2.0*r)*sinTheta2/s;
    u2Dot = 2.0*a2*r*u0*u0*sinTheta*cosTheta/(s2*s) - a2*u1*u1*sinTheta*cosTheta/((s)*(-2.0*r + a2 + r2)) + a2*u2*u2*sinTheta*cosTheta/(s) - 2.0*r*u1*u2/(s) -u0*u3*(4.0*a3*r*sinTheta3*cosTheta/s2 + 4.0*settings.a*r*sinTheta*cosTheta/(s))/(s) - 0.5*u3*u3*(-(4.0*a4*r*sinTheta3*cosTheta/s2 + 4.0*a2*r*sinTheta*cosTheta/(s))*sinTheta2 - 2.0*(2.0*a2*r*sinTheta2/(s) + a2 + r2)*sinTheta*cosTheta)/(s);
    u3Dot = -2.0*u0*u1*(-settings.a*r*(-4.0*r2/s2 + 2.0/(s))/(2.0*a2*r*sinTheta2 - 2.0*a2*r - 2.0*r3 + a4*cosTheta2 + a2*r2*cosTheta2 + a2*r2 + r4) + 0.5*(4.0*settings.a*r2*sinTheta2/s2 - 2.0*settings.a*sinTheta2/(s))*(-2.0*r + s)/(2.0*a2*r*sinTheta2*sinTheta2 - 2.0*a2*r*sinTheta2 - 2.0*r3*sinTheta2 + a4*sinTheta2*cosTheta2 + a2*r2*sinTheta2*cosTheta2 + a2*r2*sinTheta2 + r4*sinTheta2)) - 2.0*u0*u2*(-4.0*a3*r2*sinTheta*cosTheta/(s2*(2.0*a2*r*sinTheta2 - 2.0*a2*r - 2.0*r3 + a4*cosTheta2 + a2*r2*cosTheta2 + a2*r2 + r4)) + 0.5*(-4.0*a3*r*sinTheta3*cosTheta/s2 - 4.0*settings.a*r*sinTheta*cosTheta/(s))*(-2.0*r + s)/(2.0*a2*r*sinTheta2*sinTheta2 - 2.0*a2*r*sinTheta2 - 2.0*r3*sinTheta2 + a4*sinTheta2*cosTheta2 + a2*r2*sinTheta2*cosTheta2 + a2*r2*sinTheta2 + r4*sinTheta2)) - 2.0*u1*u3*(-settings.a*r*(4.0*settings.a*r2*sinTheta2/s2 - 2.0*settings.a*sinTheta2/(s))/(2.0*a2*r*sinTheta2 - 2.0*a2*r - 2.0*r3 + a4*cosTheta2 + a2*r2*cosTheta2 + a2*r2 + r4) + 0.5*(-2.0*r + s)*(-4.0*a2*r2*sinTheta2/s2 + 2.0*a2*sinTheta2/(s) + 2.0*r)*sinTheta2/(2.0*a2*r*sinTheta2*sinTheta2 - 2.0*a2*r*sinTheta2 - 2.0*r3*sinTheta2 + a4*sinTheta2*cosTheta2 + a2*r2*sinTheta2*cosTheta2 + a2*r2*sinTheta2 + r4*sinTheta2)) - 2.0*u2*u3*(-settings.a*r*(-4.0*a3*r*sinTheta3*cosTheta/s2 - 4.0*settings.a*r*sinTheta*cosTheta/(s))/(2.0*a2*r*sinTheta2 - 2.0*a2*r - 2.0*r3 + a4*cosTheta2 + a2*r2*cosTheta2 + a2*r2 + r4) + 0.5*((4.0*a4*r*sinTheta3*cosTheta/s2 + 4.0*a2*r*sinTheta*cosTheta/(s))*sinTheta2 + 2.0*(2.0*a2*r*sinTheta2/(s) + a2 + r2)*sinTheta*cosTheta)*(-2.0*r + s)/(2.0*a2*r*sinTheta2*sinTheta2 - 2.0*a2*r*sinTheta2 - 2.0*r3*sinTheta2 + a4*sinTheta2*cosTheta2 + a2*r2*sinTheta2*cosTheta2 + a2*r2*sinTheta2 + r4*sinTheta2));

    outgoingFunction.t = u0;
    outgoingFunction.r = u1;
    outgoingFunction.theta = u2;
    outgoingFunction.phi = u3;
    outgoingFunction.u0 = u0Dot;
    outgoingFunction.u1 = u1Dot;
    outgoingFunction.u2 = u2Dot;
    outgoingFunction.u3 = u3Dot;
    outgoingFunction.hStep = 0.0;
    return outgoingFunction;
}


//Restframe to ZAMO
Vector4 restframeToZAMO(const Vector3& rayDir, const localProperties& cameraProperties)
{
    Matrix4 mCamera;
    Vector4 ray4, out4;
    realNumber rCamera = cameraProperties.r;
    realNumber thetaCamera = cameraProperties.theta;
    realNumber phiCamera = cameraProperties.phi;
    realNumber deltaCamera = metric_delta(rCamera);
    realNumber sigmaCamera = metric_sigma(rCamera,thetaCamera);
    realNumber u0Camera = cameraProperties.u0;
    realNumber u1Camera = cameraProperties.u1;
    realNumber u2Camera = cameraProperties.u2;
    realNumber u3Camera = cameraProperties.u3;
    realNumber Lambda, ut, ur, uTheta, uPhi, vr, vTheta, vPhi, vX, vY, vZ, v2, gamma;

    Lambda = metric_lambda(rCamera,thetaCamera);
    ut = std::sqrt((deltaCamera*sigmaCamera) / Lambda) * u0Camera;
    ur = std::sqrt(sigmaCamera / deltaCamera) * u1Camera;
    uTheta = std::sqrt(sigmaCamera) * u2Camera;
    uPhi = (std::sin(thetaCamera) * std::sqrt(Lambda / sigmaCamera) * u3Camera) - (((2.0*settings.a*rCamera*std::sin(thetaCamera)) / std::sqrt(Lambda*sigmaCamera)) * u0Camera);
    vr = ur / ut;
    vTheta = uTheta / ut;
    vPhi = uPhi / ut;
    vX = (1.0 / std::sqrt((rCamera*rCamera) + (settings.a*settings.a) * (std::cos(thetaCamera)*std::cos(thetaCamera)))) * ( ((rCamera*std::sin(thetaCamera)*std::cos(phiCamera)*vr) + (std::sqrt((rCamera*rCamera) + (settings.a*settings.a))*std::cos(thetaCamera)*std::cos(phiCamera)*vTheta)) ) - (std::sin(phiCamera)*vPhi);
    vY = (1.0 / std::sqrt((rCamera*rCamera) + (settings.a*settings.a) * (std::cos(thetaCamera)*std::cos(thetaCamera)))) * ( (std::sqrt((rCamera*rCamera) + (settings.a*settings.a)) * std::cos(thetaCamera) * vr) - (rCamera*std::sin(thetaCamera)*vTheta) );
    vZ = (1.0 / std::sqrt((rCamera*rCamera) + (settings.a*settings.a) * (std::cos(thetaCamera)*std::cos(thetaCamera)))) * ( (rCamera*std::sin(thetaCamera)*std::sin(phiCamera)*vr) + (std::sqrt((rCamera*rCamera) + (settings.a*settings.a)) * std::cos(thetaCamera)*std::sin(phiCamera)*vTheta) ) + (std::cos(phiCamera)*vPhi);
    v2 = vX*vX + vY*vY + vZ*vZ;

    gamma = 1.0 / std::sqrt(1.0-v2);
    mCamera.M00 = 1.0 + (gamma - 1.0) * ((vX * vX) / (v2));
    mCamera.M01 = (gamma - 1.0) * ((vX * vY) / v2);
    mCamera.M02 = (gamma - 1.0) * ((vX * vZ) / v2);
    mCamera.M03 = -gamma * vX;
    mCamera.M10 = (gamma - 1.0) * ((vX * vY) / v2);
    mCamera.M11 = 1.0 + (gamma - 1.0) * ((vY * vY) / (v2));
    mCamera.M12 = (gamma - 1.0) * ((vY * vZ) / v2);
    mCamera.M13 = -gamma * vY;
    mCamera.M20 = (gamma - 1.0) * ((vX * vZ) / v2);
    mCamera.M21 = (gamma - 1.0) * ((vY * vZ) / v2);
    mCamera.M22 = 1.0 + (gamma - 1.0) * ((vZ * vZ) / (v2));
    mCamera.M23 = -gamma * vZ;
    mCamera.M30 = -gamma * vX;
    mCamera.M31 = -gamma * vY;
    mCamera.M32 = -gamma * vZ;
    mCamera.M33 = gamma;

    ray4.x = rayDir.x;
    ray4.y = rayDir.y;
    ray4.z = rayDir.z;
    ray4.t = std::sqrt(rayDir.x * rayDir.x + rayDir.y * rayDir.y + rayDir.z * rayDir.z);
    ray4 = mCamera*ray4;
    out4.x = ray4.x;
    out4.y = ray4.y;
    out4.z = ray4.z;
    out4.t = ray4.t;
    return out4;
}


//Local Velocity relative to ZAMO
realNumber localZAMOVelocity(const localProperties& positionProperties)
{
    realNumber rCamera = positionProperties.r;
    realNumber thetaCamera = positionProperties.theta;
    realNumber phiCamera = positionProperties.phi;
    realNumber deltaCamera = metric_delta(rCamera);
    realNumber sigmaCamera = metric_sigma(rCamera,thetaCamera);
    realNumber u0Camera = positionProperties.u0;
    realNumber u1Camera = positionProperties.u1;
    realNumber u2Camera = positionProperties.u2;
    realNumber u3Camera = positionProperties.u3;
    realNumber Lambda, ut, ur, uTheta, uPhi, vr, vTheta, vPhi, vX, vY, vZ, v2, gamma;

    Lambda = metric_lambda(rCamera,thetaCamera);
    ut = std::sqrt((deltaCamera*sigmaCamera) / Lambda) * u0Camera;
    ur = std::sqrt(sigmaCamera / deltaCamera) * u1Camera;
    uTheta = std::sqrt(sigmaCamera) * u2Camera;
    uPhi = (std::sin(thetaCamera) * std::sqrt(Lambda / sigmaCamera) * u3Camera) - (((2.0*settings.a*rCamera*std::sin(thetaCamera)) / std::sqrt(Lambda*sigmaCamera)) * u0Camera);
    vr = ur / ut;
    vTheta = uTheta / ut;
    vPhi = uPhi / ut;
    vX = (1.0 / std::sqrt((rCamera*rCamera) + (settings.a*settings.a) * (std::cos(thetaCamera)*std::cos(thetaCamera)))) * ( ((rCamera*std::sin(thetaCamera)*std::cos(phiCamera)*vr) + (std::sqrt((rCamera*rCamera) + (settings.a*settings.a))*std::cos(thetaCamera)*std::cos(phiCamera)*vTheta)) ) - (std::sin(phiCamera)*vPhi);
    vY = (1.0 / std::sqrt((rCamera*rCamera) + (settings.a*settings.a) * (std::cos(thetaCamera)*std::cos(thetaCamera)))) * ( (std::sqrt((rCamera*rCamera) + (settings.a*settings.a)) * std::cos(thetaCamera) * vr) - (rCamera*std::sin(thetaCamera)*vTheta) );
    vZ = (1.0 / std::sqrt((rCamera*rCamera) + (settings.a*settings.a) * (std::cos(thetaCamera)*std::cos(thetaCamera)))) * ( (rCamera*std::sin(thetaCamera)*std::sin(phiCamera)*vr) + (std::sqrt((rCamera*rCamera) + (settings.a*settings.a)) * std::cos(thetaCamera)*std::sin(phiCamera)*vTheta) ) + (std::cos(phiCamera)*vPhi);
    v2 = vX*vX + vY*vY + vZ*vZ;
    return (v2);
}


//ZAMO to Global reference frame
Vector4 ZAMOToGlobal(const Vector4& rayDir4, realNumber r, realNumber theta, realNumber phi, realNumber delta, realNumber sigma)
{
    Vector4 uOut;
    realNumber rayr, raytheta, rayphi, rayt, Lambda, u0, u1, u2, u3;

    rayr = (1.0 / std::sqrt((r*r) + (settings.a*settings.a) * (std::cos(theta)*std::cos(theta)))) * ((r*std::sin(theta)*std::cos(phi)*rayDir4.x) + (r*std::sin(theta)*std::sin(phi)*rayDir4.z) + (std::sqrt((r*r) + (settings.a*settings.a)) * std::cos(theta) * rayDir4.y));
    raytheta = (1.0 / std::sqrt((r*r) + (settings.a*settings.a) * (std::cos(theta)*std::cos(theta)))) * ( (std::sqrt((r*r) + (settings.a*settings.a)) * std::cos(theta) * std::cos(phi) * rayDir4.x) + (std::sqrt((r*r) + (settings.a*settings.a)) * std::cos(theta) * std::sin(phi) * rayDir4.z) - (r*std::sin(theta)*rayDir4.y) );
    rayphi = (-std::sin(phi) * rayDir4.x) + (std::cos(phi) * rayDir4.z);
    rayt = rayDir4.t;
    Lambda = metric_lambda(r,theta);
    u0 = std::sqrt( Lambda / (delta*sigma) ) * rayt;
    u1 = std::sqrt( delta / sigma ) * rayr;
    u2 = (1.0 / std::sqrt(sigma)) * raytheta;
    u3 = 2.0*settings.a*r/std::sqrt(Lambda*delta*sigma)*rayt + std::sqrt(sigma/Lambda)/std::sin(theta)*rayphi;

    uOut.x = u0;
    uOut.y = u1;
    uOut.z = u2;
    uOut.t = u3;
    return uOut;
}


//Compute u0
realNumber deriveu0(realNumber u1, realNumber u2, realNumber u3, realNumber r, realNumber theta, realNumber delta, realNumber sigma, realNumber My)
{
    realNumber gtPhi, gtt, grr, gThetaTheta, Lambda, gPhiPhi, rootTerm;

    gtPhi = metric_gtPhi(r,theta);
    gtt = metric_gtt(r,theta);
    grr = metric_grr(r,theta);
    gThetaTheta = sigma;
    Lambda = metric_lambda(r,theta);
    gPhiPhi = metric_gPhiPhi(r,theta);
    rootTerm = (((gtPhi/gtt) * (gtPhi/gtt))*(u3*u3)) - ((1.0 / gtt) * ( grr*(u1*u1) + gThetaTheta*(u2*u2) + gPhiPhi*(u3*u3) - My ));

    if(rootTerm < 0)
    {
        rootTerm = 0;
    }
    else
    {
        rootTerm = std::sqrt(rootTerm);
    }

    return ((-gtPhi / gtt)*u3) + rootTerm;
}


//Runge Kutta-Fehlberg Scheme
RayProperties rayRKF45(const RayProperties& incomingValues)
{
    realNumber rateOfChange, rateOfChangeC, minStep, hU, absError = 1e10, tol = 1e-5, newh = 0.0;
    RayProperties newy, error, k1, k2, k3, k4, k5, k6;
    hU = incomingValues.hStep;

    while(absError > tol)
    {
        k1 = rayDerivatives(incomingValues);
        k2 = rayDerivatives(incomingValues + (1.0/4.0)*k1*hU);
        k3 = rayDerivatives(incomingValues + (3.0/32.0)*k1*hU + (9.0/32.0)*k2*hU);
        k4 = rayDerivatives(incomingValues + (1932.0/2197.0)*k1*hU - (7200.0/2197.0)*k2*hU + (7296.0/2197.0)*k3*hU);
        k5 = rayDerivatives(incomingValues + (439.0/216.0)*k1*hU - 8.0*k2*hU + (3680.0/513.0)*k3*hU - (845.0/4104.0)*k4*hU);
        k6 = rayDerivatives(incomingValues - (8.0/27.0)*k1*hU + 2.0*k2*hU - (3544.0/2565.0)*k3*hU + (1859.0/4104.0)*k4*hU - (11.0/40.0)*k5*hU);


        newy = incomingValues + (16.0/135.0)*k1*hU + (6656.0/12825.0)*k3*hU + (28561.0/56430.0)*k4*hU - (9.0/50.0)*k5*hU + (2.0/55.0)*k6*hU;
        error = (-1.0/360.0)*k1*hU + (128.0/4275.0)*k3*hU + (2197.0/75240.0)*k4*hU - (1.0/50.0)*k5*hU - (2.0/55.0)*k6*hU;
        absError = std::sqrt((error.t*error.t) + (error.r*error.r) + (error.theta*error.theta) + (error.phi*error.phi) + (error.u0*error.u0) + (error.u1*error.u1) + (error.u2*error.u2) + (error.u3*error.u3) + (error.hStep*error.hStep));

        newh = 0.9*hU*std::pow((tol/absError),(1.0/5.0));

        if(absError > tol)
        {
            hU = newh;
        }
    }

    newy.hStep = newh;
    return newy;
}


//Redshift factor function
realNumber generalizedRedshiftFactor(const RayProperties& rayProperties, const localProperties& cameraProperties, const Vector4& initialU, const Vector3& Vel)
{
    realNumber delta, sigma, gtPhi, gtt, grr, gThetaTheta, lambda, gPhiPhi, vel0;

    delta = metric_delta(rayProperties.r);
    sigma = metric_sigma(rayProperties.r,rayProperties.theta);

    //Metric Tensor Components
    gtPhi = metric_gtPhi(rayProperties.r,rayProperties.theta);
    gtt = metric_gtt(rayProperties.r,rayProperties.theta);
    grr = sigma / delta;
    gThetaTheta = sigma;
    lambda = metric_lambda(rayProperties.r,rayProperties.theta);
    gPhiPhi = metric_gPhiPhi(rayProperties.r,rayProperties.theta);

    //Ergosphere case
    if(gtt < 0)
    {
        vel0 = -(gtPhi / gtt) * Vel.z + std::sqrt( ((gtPhi / gtt)*(gtPhi / gtt)) * (Vel.z * Vel.z) - (1.0/gtt) * (grr*(Vel.x * Vel.x)+gThetaTheta*(Vel.y * Vel.y)+gPhiPhi*(Vel.z * Vel.z) + 1.0));
    }
    else
    {
        vel0 = -(gtPhi / gtt) * Vel.z - std::sqrt( ((gtPhi / gtt)*(gtPhi / gtt)) * (Vel.z * Vel.z) - (1.0/gtt) * (grr*(Vel.x * Vel.x)+gThetaTheta*(Vel.y * Vel.y)+gPhiPhi*(Vel.z * Vel.z) + 1.0));
    }

    //Emitter frequency
    realNumber frequencyEmitter = (gtt*rayProperties.u0*vel0) + gtPhi * ((rayProperties.u0*Vel.z + rayProperties.u3*vel0)) + (grr*rayProperties.u1*Vel.x) + (gThetaTheta*rayProperties.u2*Vel.y) + (gPhiPhi*rayProperties.u3*Vel.z);

    //Camera Frequency
    realNumber deltaCam, sigmaCam, gtPhiCam, gttCam, grrCam, gThetaThetaCam, lambdaCam, gPhiPhiCam;
    deltaCam = metric_delta(cameraProperties.r);
    sigmaCam = metric_sigma(cameraProperties.r,cameraProperties.theta);
    gtPhiCam = metric_gPhi(cameraProperties.r,cameraProperties.theta);
    gttCam = metric_gtt(cameraProperties.r,cameraProperties.theta);
    grrCam = sigmaCam / deltaCam;
    gThetaThetaCam = sigmaCam;
    lambdaCam = metric_lambda(cameraProperties.r,cameraProperties.theta);
    gPhiPhiCam = metric_gPhiPhi(cameraProperties.r,cameraProperties.theta);

    //Final
    Vector4 intU = initialU; //For intU read FirstU0, FirstU.x, FirstU.y, FirstU.z
    Vector3 CameraVelocity;
    CameraVelocity.set(-cameraProperties.u1,-cameraProperties.u2,-cameraProperties.u3);
    realNumber frequencyCam = gttCam*intU.x*cameraProperties.u0 + gtPhiCam * (intU.x*CameraVelocity.z + intU.t * cameraProperties.u0) + grrCam * intU.y * CameraVelocity.x + gThetaThetaCam*intU.z*CameraVelocity.y + gPhiPhiCam*intU.t*CameraVelocity.z;
    realNumber Z = frequencyEmitter/frequencyCam;

    if (!std::isfinite(Z) || Z < 0){
        Z=1000.0;
    }
    return Z;
}


//Generalized time dilation formular for 4-Velocities
realNumber generalizedTimeDilationFactor(realNumber r, realNumber theta, const Vector3& uVel, realNumber mu)
{
    //Initials and pre computation
    realNumber sigma, delta, gtPhi, gtt, gThetaTheta, Lambda, gPhiPhi, grr, factor;

    //Metric Tensor components
    sigma = metric_sigma(r,theta);
    delta = metric_delta(r);
    gtPhi = metric_gtPhi(r,theta);
    gtt = metric_gtt(r,theta);
    gThetaTheta = sigma;
    Lambda = metric_lambda(r,theta);
    gPhiPhi = metric_gPhiPhi(r,theta);
    grr = sigma/delta;

    //Distinguish between ray being inside vs outside the Ergosphere
    if(gtt < 0)
    {
        factor = (-(gtPhi/gtt)*uVel.z) + std::sqrt( ((gtPhi/gtt)*uVel.z * (gtPhi/gtt)*uVel.z) - (1.0/gtt) * ( grr*uVel.x*uVel.x + gThetaTheta*uVel.y*uVel.y + gPhiPhi*uVel.z*uVel.z - mu ) );
    }
    else
    {
        factor = (-(gtPhi/gtt)*uVel.z) - std::sqrt( ((gtPhi/gtt)*uVel.z * (gtPhi/gtt)*uVel.z) - (1.0/gtt) * ( grr*uVel.x*uVel.x + gThetaTheta*uVel.y*uVel.y + gPhiPhi*uVel.z*uVel.z - mu ) );
    }

    //Error handling because Computers suck
    if(std::isnan(factor))
    {
        factor = 100.0;
    }

    return factor;
}


//Disk Velocity function
realNumber uPhi(realNumber r)
{
    return (sign(settings.a)*std::sqrt(r) / (r * std::sqrt( (r*r) - (3.0*r) + (2.0*std::abs(settings.a)*std::sqrt(r))))) * settings.disk_lattice_velocity_scale;
}


//Astrophysical Jet Velocity function
Vector3 uATJ(realNumber vtau, realNumber verticalDirection, realNumber rtau, realNumber thetatau)
{
    Vector3 out;

    //Precompute and Equation of Motion
    realNumber ftau = std::sqrt(rtau*rtau + settings.a*settings.a)*std::sin(thetatau);
    realNumber rtauSqrtTerm1 = std::sqrt( ((settings.vy*settings.vy*vtau*vtau - settings.a*settings.a + ftau*ftau)*(settings.vy*settings.vy*vtau*vtau - settings.a*settings.a + ftau*ftau)) + 4.0*settings.vy*settings.vy*settings.a*settings.a*vtau*vtau );
    realNumber rtauSqrtTerm2 = std::sqrt( settings.vy*settings.vy*vtau*vtau - settings.a*settings.a + ftau*ftau +  rtauSqrtTerm1 );

    //Velocity
    realNumber dfdtau = settings.vr / (2.0*ftau);
    realNumber term1 = (1.0/std::sqrt(2.0));
    realNumber term2 = 2.0*settings.vy*settings.vy*vtau + 2.0*dfdtau*ftau;
    realNumber term3 = 2.0*(settings.vy*settings.vy*vtau*vtau - settings.a*settings.a + ftau*ftau) * (2.0*settings.vy*settings.vy*vtau+2.0*dfdtau*ftau);
    realNumber term4 = 8.0*settings.vy*settings.vy*settings.a*settings.a*vtau;
    realNumber term5 = 2.0*rtauSqrtTerm1;
    realNumber term6 = 2.0*rtauSqrtTerm2;
    realNumber drdtau = term1 * ( (term2 + ( (term3+term4) / (term5) )) / (term6) );
    realNumber dthetadtau = -(1.0/(std::sqrt(1.0-( ((settings.vy*vtau)/(rtau))*((settings.vy*vtau)/(rtau)) )))) * ((settings.vy/rtau) - ( ((settings.vy*vtau)/(rtau*rtau))*(drdtau))) * sign(verticalDirection);
    realNumber dphidtau = omega(rtau,thetatau,drdtau,dthetadtau);

    out.x = drdtau;
    out.y = dthetadtau;
    out.z = dphidtau;
    return out;
}


//Wind velocity function
Vector3 uWind(realNumber r0, realNumber vtau, realNumber verticalDirection, realNumber r, realNumber theta)
{
    realNumber windVelocity = 24.0;
    Vector3 u1, u2;

    //Velocity Contributions
    u1.x = windVelocity;
    u1.y = 0.0;
    u1.z = 0.0;
    u2 = uATJ(vtau,verticalDirection,r,theta);

    //Factor
    realNumber k = 1.0 / (1.0 + std::exp(-1.0*(-r+r0)));

    return u1-(u2*k);
}


//Varying Phi Velocity for jet
realNumber omega(realNumber r, realNumber theta, realNumber drdtau, realNumber dthetadtau)
{
    realNumber Sigma = metric_sigma(r,theta);
    realNumber Delta = metric_delta(r);
    realNumber Lambda = metric_lambda(r,theta);
    realNumber gPhiPhi = metric_gPhiPhi(r,theta);
    realNumber gtPhi = metric_gtPhi(r,theta);
    realNumber gtt = metric_gtt(r,theta);
    realNumber grr = Sigma / Delta;
    realNumber gThetaTheta = Sigma;
    realNumber W = sign(settings.a)*((settings.w) / (std::pow(std::sqrt(r*r + settings.a*settings.a)*std::sin(theta),settings.alpha)+1.0_real));

    realNumber phiZAMO = -(gtPhi/gPhiPhi) * std::sqrt( (1.0_real+grr*drdtau*drdtau + gThetaTheta*dthetadtau*dthetadtau + gPhiPhi*W*W)/((gtPhi*gtPhi)/(gPhiPhi)-gtt) );

    return (sign(settings.a)*((settings.w) / (std::pow(std::sqrt(r*r + settings.a*settings.a)*std::sin(theta),settings.alpha)+1.0_real)) + phiZAMO);
}


//Fake / Unrealistic Temperature distroy
realNumber fauxTemperatureDistribution(realNumber T1, realNumber T2, realNumber r, realNumber r1, realNumber r2, realNumber N)
{
    return std::pow(std::abs(r-r2),N)*((T1-T2) / std::pow(std::abs(r1-r2),N))+T2;
}


//Stores Color matching functions for Wavelength to RGB conversion
static Vector3 colorMatchingFunction(int wavelength_index)
{
    realNumber m[243] = { 
        0.0014, 0.0000, 0.0065,
        0.0022, 0.0001, 0.0105,
        0.0042, 0.0001, 0.0201,
        0.0076, 0.0002, 0.0362,
        0.0143, 0.0004, 0.0679,
        0.0232, 0.0006, 0.1102,
        0.0435, 0.0012, 0.2074,
        0.0776, 0.0022, 0.3713,
        0.1344, 0.0040, 0.6456,
        0.2148, 0.0073, 1.0391,
        0.2839, 0.0116, 1.3856,
        0.3285, 0.0168, 1.6230,
        0.3483, 0.0230, 1.7471,
        0.3481, 0.0298, 1.7826,
        0.3362, 0.0380, 1.7721,
        0.3187, 0.0480, 1.7441,
        0.2908, 0.0600, 1.6692,
        0.2511, 0.0739, 1.5281,
        0.1954, 0.0910, 1.2876,
        0.1421, 0.1126, 1.0419,
        0.0956, 0.1390, 0.8130,
        0.0580, 0.1693, 0.6162,
        0.0320, 0.2080, 0.4652,
        0.0147, 0.2586, 0.3533,
        0.0049, 0.3230, 0.2720,
        0.0024, 0.4073, 0.2123,
        0.0093, 0.5030, 0.1582,
        0.0291, 0.6082, 0.1117,
        0.0633, 0.7100, 0.0782,
        0.1096, 0.7932, 0.0573,
        0.1655, 0.8620, 0.0422,
        0.2257, 0.9149, 0.0298,
        0.2904, 0.9540, 0.0203,
        0.3597, 0.9803, 0.0134,
        0.4334, 0.9950, 0.0087,
        0.5121, 1.0000, 0.0057,
        0.5945, 0.9950, 0.0039,
        0.6784, 0.9786, 0.0027,
        0.7621, 0.9520, 0.0021,
        0.8425, 0.9154, 0.0018,
        0.9163, 0.8700, 0.0017,
        0.9786, 0.8163, 0.0014,
        1.0263, 0.7570, 0.0011,
        1.0567, 0.6949, 0.0010,
        1.0622, 0.6310, 0.0008,
        1.0456, 0.5668, 0.0006,
        1.0026, 0.5030, 0.0003,
        0.9384, 0.4412, 0.0002,
        0.8544, 0.3810, 0.0002,
        0.7514, 0.3210, 0.0001,
        0.6424, 0.2650, 0.0000,
        0.5419, 0.2170, 0.0000,
        0.4479, 0.1750, 0.0000,
        0.3608, 0.1382, 0.0000,
        0.2835, 0.1070, 0.0000,
        0.2187, 0.0816, 0.0000,
        0.1649, 0.0610, 0.0000,
        0.1212, 0.0446, 0.0000,
        0.0874, 0.0320, 0.0000,
        0.0636, 0.0232, 0.0000,
        0.0468, 0.0170, 0.0000,
        0.0329, 0.0119, 0.0000,
        0.0227, 0.0082, 0.0000,
        0.0158, 0.0057, 0.0000,
        0.0114, 0.0041, 0.0000,
        0.0081, 0.0029, 0.0000,
        0.0058, 0.0021, 0.0000,
        0.0041, 0.0015, 0.0000,
        0.0029, 0.0010, 0.0000,
        0.0020, 0.0007, 0.0000,
        0.0014, 0.0005, 0.0000,
        0.0010, 0.0004, 0.0000,
        0.0007, 0.0002, 0.0000,
        0.0005, 0.0002, 0.0000,
        0.0003, 0.0001, 0.0000,
        0.0002, 0.0001, 0.0000,
        0.0002, 0.0001, 0.0000,
        0.0001, 0.0000, 0.0000,
        0.0001, 0.0000, 0.0000,
        0.0001, 0.0000, 0.0000,
        0.0000, 0.0000, 0.0000 };

    int index = wavelength_index * 3;
    Vector3 out;
    out.x = m[index];
    out.y = m[index + 1];
    out.z = m[index + 2];
    return out;
}


//Temperature, in Kelvin, to RGB
Vector4 plancksLaw(realNumber T, realNumber RedshiftFactor)
{
    realNumber intensity, lambda, cX = 0.0, cY = 0.0, cZ = 0.0, lightSpeed = 2.99792458e8, kR = 1.3806504e-23, hR = 6.62607015e-34, Temperature = T / RedshiftFactor, X,Z,R = 0.0,G = 0.0,B = 0.0, M = 1.0;
    Vector3 RGB_Intensity;
    Vector4 out;

    if(!std::isfinite(Temperature) || Temperature < 100.0){
        Temperature = 100.0;
    }

    for(int i = 0; i < 80; i++)
    {
        lambda = (i*5.0+380.0)*1E-9;
        intensity = ((2.0*hR*std::pow(lightSpeed,2)) / std::pow(lambda,5)) / (std::exp((hR*lightSpeed) / (lambda*kR*Temperature)) - 1.0);
        RGB_Intensity = colorMatchingFunction(i);
        cX += intensity * RGB_Intensity.x*5.0;
        cY += intensity * RGB_Intensity.y*5.0;
        cZ += intensity * RGB_Intensity.z*5.0;
    }

    X = cX/cY;
    Z = cZ/cY;
    R = 3.24096994*X - 1.53738318 - 0.49861076*Z;

    G = -0.96924364*X + 1.8759675 + 0.04155506*Z;
    B = 0.05563008*X - 0.20397696 + 1.05697151*Z;
    M = std::max(R,G);
    M = std::max(M,B);

    cY = std::max(cY,0.0);

    if(M > 1.0)
    {
        R = R / M;
        G = G / M;
        B = B / M;
    }

    R = std::max(R,0.0);
    G = std::max(G,0.0);
    B = std::max(B,0.0);

    out.x = R;
    out.y = G;
    out.z = B;
    out.t = cY;
    return out;
}


//Wind density function
Vector2 sampleWindDensity(realNumber r, const Vector3& WindBounds)
{
    realNumber rho0 = 5.0;
    realNumber rhoScalar = (1.0/(r*r*r))*rho0;
    realNumber rFade = 1.0/(1.0+std::exp(-10*(r-rho0)));

    realNumber Emission = rhoScalar*rFade, Absorption = 0.0;
    Vector2 out;
    out.x = Emission;
    out.y = Absorption;
    return out;
}


//Disk Density
Vector2 sampleAccretionDiskDensity(const RayProperties& rayProperties, const Vector3& rayOrig, noise_point* noise_points, realNumber noiseValue, realNumber constValue, realNumber derivativeValue)
{
    //Define Initials
    Vector2 out;
    Vector3 UpVector, DirVector;
    realNumber sampleAngle, initialSlope, angleFallOff, variableSlope, Emission = 0.0, Absorption = 0.0;

    //Define Up, Y Axis, direction relative to rayOrig
    UpVector.x = 0.0;
    UpVector.y = sign(rayOrig.y);
    UpVector.z = 0.0;
    DirVector = rayOrig / VectorLength(rayOrig);

    //Angle between Up vector and current rayOrig
    sampleAngle = dotProduct(UpVector,DirVector) * (180.0/M_PI);
    initialSlope = 10.0; //Represents the average slope of the disk

    //Disk slope falls off towards 0 near Event Horizon
    angleFallOff = (-exp(-1.0*rayProperties.r) + 1.0) * (-exp(-1.0*rayProperties.r) + 1.0);

    //Change the slope according to the constant vortexNoise output, that is the noise constant along Y
    variableSlope = (initialSlope + constValue*7.0f) * angleFallOff;

    //In Exact bounds test
    if(sampleAngle < variableSlope)
    {
        //Shader Initials
        realNumber inverseNoise, OuterEmissionFalloff, OuterAbsorptionFalloff, InnerAbsorptionFalloff, verticalFade;

        //To avoid harsh cutoffs at the Bounding region edges we fade the disk in and out
        OuterEmissionFalloff = 1.0 / (1.0 + std::exp(0.15 * (rayProperties.r - 44.0)));
        OuterAbsorptionFalloff = 1.0 / (1.0 + std::exp(0.1 * (rayProperties.r - 64.0)));
        InnerAbsorptionFalloff = 1.0 / (1.0 + std::exp(-0.1 * (rayProperties.r - 32.0)));
        verticalFade = gx(5.0,2.0,variableSlope,rayProperties.r,std::abs(rayOrig.y));

        //Basically the true "Shader" code, Emission and Absorption define the visual output, and its just kinda made up
        inverseNoise = abs(-noiseValue+1.0);
        Emission = ((noiseValue*noiseValue) + noiseValue + (derivativeValue*derivativeValue*96.0) + 0.1) * 0.5 * OuterEmissionFalloff * verticalFade;
        Absorption = ((inverseNoise*inverseNoise*inverseNoise)*2.0 + inverseNoise + 0.1) * OuterAbsorptionFalloff * verticalFade * 0.1;
    }

    out.x = Emission;
    out.y = Absorption;
    return out;
}


//Astrophysical Jet
Vector2 sampleAstrophysicalJetDensity(const RayProperties& rayProperties, const Vector3& AccretionDiskBounds, const Vector3& mJet)
{
    realNumber Emission = 0.0, Absorption = 0.0;
    realNumber term1 = std::sqrt(rayProperties.r * rayProperties.r + settings.a*settings.a);
    realNumber x = term1 * std::sin(rayProperties.theta);
    realNumber h = std::abs(rayProperties.r * std::cos(rayProperties.theta));
    realNumber maelstromNoiseField = std::abs(mJet.x);
    realNumber maelstormNoiseDerivative = std::abs(mJet.y);
    realNumber dY = x;
    realNumber inverseNoise = std::abs(-maelstromNoiseField + 1.0);

    //Fades
    //Vertical Fade
    realNumber fade1 = std::exp(-h * 0.05);

    //Bounding Region Edge Fade
    realNumber u1 = AccretionDiskBounds.z + h*std::tan(AccretionDiskBounds.y*(M_PI/180.0));
    realNumber fade2 = bezierFade(x,u1,1.0,3.0,0.3);

    //Y-Axis Cooridor Fade
    realNumber fade3 = 1.0;
    realNumber CorridorDepression = 1.0;
    if(h > (x*x)-CorridorDepression)
    {
        fade3 = 0.0;
    }
    else
    {
        realNumber dtoEdge = x - std::sqrt(h+CorridorDepression);
        realNumber CorridorFadeRadius = 4.0;
        fade3 = bezierFade(-dtoEdge+CorridorFadeRadius,CorridorFadeRadius,1.0,2.0,0.1);
    }

    realNumber mfade = fade1*fade2*fade3;
    realNumber aFade = 1.0/(1.0 + std::exp(-0.5*(rayProperties.r -24.0)));

    //Noise Fields
    realNumber EmissionScale = 16.0; //0.003
    realNumber AbsorptionScale = 2.0;

    Emission = (std::pow(maelstromNoiseField,9))*EmissionScale*mfade; //((9.0*noiseField*noiseField + 5.0*noiseField) + 0.01 + (verticalFade*noiseDerivative*0.0))*EmissionScale*mfade;
    Absorption = (aFade*maelstormNoiseDerivative*96.0)*AbsorptionScale*mfade; //(std::pow(std::pow(inverseNoise,5)*5.0,3))*AbsorptionScale*mfade

    Vector2 out;
    out.x = Emission;
    out.y = Absorption;
    return out;
}


//Radial Checkerboard
realNumber radialChecker(const realNumber r, const realNumber phi, const realNumber scale, const realNumber w, const realNumber x)
{
    //Radial strpes
    realNumber x1 = std::pow(2.0,std::round(x) * -1.0);
    realNumber p1 = std::fmod(((phi / M_PI) + 1.0) * 0.5,x1) * (1.0/x1);
    if(p1 > 0.5)
    {
        p1 = 1.0;
    }
    else
    {
        p1 = 0.0;
    }
    realNumber ip1 = std::abs(p1 - 1.0);

    //Radial bands
    realNumber r1 = r / scale;
    realNumber r2 = r1 + (w/2.0);
    realNumber b1 = std::fmod(r1,w) * (1.0/w);
    realNumber b2 = std::fmod(r2,w) * (1.0/w);

    if(b1 > 0.5)
    {
        b1 = 1.0;
    }
    else
    {
        b1 = 0.0;
    }

    if(b2 > 0.5)
    {
        b2 = 1.0;
    }
    else
    {
        b2 = 0.0;
    }

    b1 = b1 * p1;
    b2 = b2 * ip1;

    return b1+b2;
}

}
