#include "vmec_common.h"

#include <cstdlib>
#include <cstring>

namespace vmec
{
    /*

    Couple of thoughts before we begin.

    We will follow the paper to the letter, which means it might just be a good idea to set up a test enviornment within VMEC. One without any relativistic nonsense. A simple box raytracer where all we intersect is one cube.
    Which is then the boundin box for the test volume. The paper is not dealing with curved rays, so we should just do that. So we just have another render mode "toggleMagikSandboxScene"
    Of course at some point we will have to port the flat spacetime stuff to VMEC proper. I think we can do ourselfs a favor by not using the main geodesic for volume sampling. 
    Here is the idea, the Geodesic has a step size determined by the RKF45 method. This step size does not reflect what is needed for proper volume rendering. So why not trace a 2nd line
    Say RKF45 wants to take a step of 10, but Magik´s RNG wants 5. Ok so we walk along the geodesic and take the steps Magik wants. So say its 5,3,7. Well 5+3+7 > 10 so before the 3rd step we would continue the Geodesic with RKF45 until its long enough 
    for Magik to take over. I guess we do "Parallel Integration". So Magik is not dependent on RKF45 for its quality. Which also means we can take the restrictions out of the integrator. 

    */



    //Sandbox Scene Specific
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
}