#pragma once
#include "vmec_defines.h"
#include "vmec_types.h"

#include <random>
#include <time.h>
#include <cstdlib>
#include <cstring>

namespace vmec
{

//Generic mathematical functions
realNumber nrandom();
unsigned int hashFunction(unsigned int seed);
realNumber ndetermin(const int seed);
realNumber sign(const realNumber x);
realNumber lerp(const realNumber c, const realNumber b, const realNumber t);
realNumber euclidianDistance(const Vector3& v1, const Vector3& v2);
realNumber euclidianDistance2(const Vector3& v1, const Vector3& v2);
Vector3 rotate2D(const Vector3& rayOrig, const realNumber angle);
Matrix3 rotateRay(const realNumber rX, const realNumber rY, const realNumber rZ);

//Metric Tensor Components
realNumber metric_delta(const realNumber r);
realNumber metric_sigma(const realNumber r, const realNumber theta);
realNumber metric_lambda(const realNumber r, const realNumber theta);
realNumber metric_gPhiPhi(const realNumber r, const realNumber theta);
realNumber metric_gPhi(const realNumber r, const realNumber theta);
realNumber metric_gtt(const realNumber r, const realNumber theta);
realNumber metric_grr(const realNumber r, const realNumber theta);
realNumber metric_gThetaTheta(const realNumber r, const realNumber theta);
realNumber metric_gtPhi(const realNumber r, const realNumber theta);

//Conversions
Vector3 sphericalToCartesian(const realNumber r, const realNumber theta, const realNumber phi);
Vector3 cartesianToBL(const Vector3& rayOrig);
Vector3 BoyerLindquistToCartesian(const realNumber r, const realNumber theta, const realNumber phi);

//Specific mathematics
realNumber quadraticBezier(const realNumber x, const realNumber P0, const realNumber P1);
realNumber bezierFade(const realNumber x, const realNumber u1, const realNumber h, const realNumber N, const realNumber k);
realNumber gx(const realNumber N0, const realNumber N1, const realNumber s, const realNumber r, const realNumber x);
realNumber radialChecker(const realNumber r, const realNumber phi, const realNumber scale, const realNumber w, const realNumber x);
bool rayDiskIntersection(const Vector3 &n, const Vector3 &p0, const Vector3 &rayOrig, const Vector3 &rayDir, const realNumber &dx);

//Local initializations
localProperties initializeCamera(const Vector3& camPos, const Vector3& camMomentum);
RayProperties initializeRay(const Vector3& rayOrig, const Vector3& rayDir, const localProperties& cameraProperties);

//Ray operations
void pointInBounds(BVObjects& BoundingVolumes, const realNumber& r, const realNumber& theta);
RayProperties rayRKF45(const RayProperties& incomingValues, const BVObjects BV, realNumber Medium);

//Ray Initialization
Vector4 restframeToZAMO(const Vector3& rayDir, const localProperties& cameraProperties);
Vector4 ZAMOToGlobal(const Vector4& rayDir4, realNumber r, realNumber theta, realNumber phi, realNumber delta, realNumber sigma);
realNumber deriveu0(realNumber u1, realNumber u2, realNumber u3, realNumber r, realNumber theta, realNumber delta, realNumber sigma, realNumber My);
realNumber localZAMOVelocity(const localProperties& positionProperties);

//Relativistic physics
realNumber generalizedRedshiftFactor(const RayProperties& rayProperties, const localProperties& cameraProperties, const Vector4& initialU, const Vector3& Vel);
realNumber generalizedTimeDilationFactor(realNumber r, realNumber theta, const Vector3& uVel, realNumber mu);
realNumber uPhi(realNumber r);
Vector3 uATJ(realNumber vtau, realNumber verticalDirection, realNumber rtau, realNumber thetatau);
Vector3 uWind(realNumber r0, realNumber vtau, realNumber verticalDirection, realNumber r, realNumber theta);
realNumber omega(realNumber r, realNumber theta, realNumber drdtau, realNumber dthetadtau);

//Physical
realNumber fauxTemperatureDistribution(realNumber T1, realNumber T2, realNumber r, realNumber r1, realNumber r2, realNumber N);
Vector4 plancksLaw(realNumber T, realNumber RedshiftFactor);

//Volume
Vector2 sampleWindDensity(realNumber r, const Vector3& WindBounds);
Vector2 sampleAccretionDiskDensity(const RayProperties& rayProperties, const Vector3& rayOrig, noise_point* sorted_noise_points, realNumber noiseValue, realNumber constValue, realNumber derivativeValue);
Vector2 sampleAstrophysicalJetDensity(const RayProperties& rayProperties, const Vector3& AccretionDiskBounds, const Vector3& mJet);

    using usize = uint64_t;
    /// Allocate physical memory, not virtual.
    void*
    raw_allocate( usize size );

    /// Allocate physical memory, not virtual.
    template <typename T>
    T*
    allocate( usize size )
    {
        usize byte_size = (size * sizeof(T));
        void* result = malloc( byte_size );
        std::memset( result, 0, byte_size );
        return (T*)result;
    }

#ifndef CLING_INTERPRETER_H
#define print( ... ) log(  __FILE__, __VA_ARGS__ )
#endif

    template<typename... t_streamable>
    void
    FUNCTION log( const char* category, t_streamable... messages )
    {
        std::cout << "[" << category << "]      ";
        ((std::cout << messages  << " "), ...) << "\n";
    }

    #define log_flush() std::cout.flush()
    #define vmec_log( ... ) log( "VMEC", __VA_ARGS__)
}

using namespace vmec;
