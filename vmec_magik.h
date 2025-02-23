#pragma once
#include "vmec_defines.h"
#include "vmec_types.h"

#include <cstdlib>
#include <cstring>

namespace vmec
{
    bool intersectAABB(const Vector3 &rayOrig, const Vector3 &rayDir, const Vector3 &max, const Vector3 &min, Vector3 &hitIn, Vector3 &hitOut);

    class Magik
    {
        public:
        void trivial();
        void single();
        void multiple();

        void trivial_relativistic();
        void single_relativistic();
        void multiple_relativistic();
    };
}