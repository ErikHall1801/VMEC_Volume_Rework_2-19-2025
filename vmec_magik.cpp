#include "vmec_common.h"
#include "vmec_magik.h"

#include <cstdlib>
#include <cstring>

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

    */



    //Sandbox Scene Specific functions
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



    //Main
        //Non Relativistic - All Sandbox
        void Magik::trivial()
        {
            //Do Trivial
        }

        void Magik::single()
        {
            //Do Single
        }

        void Magik::multiple()
        {
            //Do Multiple
        }



        //Relativisitc
        void Magik::trivial_relativistic()
        {
            //Do Trivial but relativisitc
        }

        void Magik::single_relativistic()
        {
            //Do Single but relativisitc
        }

        void Magik::multiple_relativistic()
        {
            //Do Multiple but relativisitc
        }
}