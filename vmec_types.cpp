#include "vmec_defines.h"
#include "vmec_types.h"


RayProperties operator+(const RayProperties& rp1, const RayProperties& rp2) {
    RayProperties result;
    result.t = rp1.t + rp2.t;
    result.r = rp1.r + rp2.r;
    result.theta = rp1.theta + rp2.theta;
    result.phi = rp1.phi + rp2.phi;
    result.u0 = rp1.u0 + rp2.u0;
    result.u1 = rp1.u1 + rp2.u1;
    result.u2 = rp1.u2 + rp2.u2;
    result.u3 = rp1.u3 + rp2.u3;
    result.hStep = rp1.hStep + rp2.hStep;
    return result;
}

RayProperties operator-(const RayProperties& rp1, const RayProperties& rp2) {
    RayProperties result;
    result.t = rp1.t - rp2.t;
    result.r = rp1.r - rp2.r;
    result.theta = rp1.theta - rp2.theta;
    result.phi = rp1.phi - rp2.phi;
    result.u0 = rp1.u0 - rp2.u0;
    result.u1 = rp1.u1 - rp2.u1;
    result.u2 = rp1.u2 - rp2.u2;
    result.u3 = rp1.u3 - rp2.u3;
    result.hStep = rp1.hStep - rp2.hStep;
    return result;
}

RayProperties operator+(const RayProperties& rp1, const realNumber& c) {
    RayProperties result;
    result.t = rp1.t + c;
    result.r = rp1.r + c;
    result.theta = rp1.theta + c;
    result.phi = rp1.phi + c;
    result.u0 = rp1.u0 + c;
    result.u1 = rp1.u1 + c;
    result.u2 = rp1.u2 + c;
    result.u3 = rp1.u3 + c;
    result.hStep = rp1.hStep + c;
    return result;
}

RayProperties operator+(const realNumber& c, const RayProperties& rp1) {
    RayProperties result;
    result.t = rp1.t + c;
    result.r = rp1.r + c;
    result.theta = rp1.theta + c;
    result.phi = rp1.phi + c;
    result.u0 = rp1.u0 + c;
    result.u1 = rp1.u1 + c;
    result.u2 = rp1.u2 + c;
    result.u3 = rp1.u3 + c;
    result.hStep = rp1.hStep + c;
    return result;
}

RayProperties operator-(const RayProperties& rp1, const realNumber& c) {
    RayProperties result;
    result.t = rp1.t - c;
    result.r = rp1.r - c;
    result.theta = rp1.theta - c;
    result.phi = rp1.phi - c;
    result.u0 = rp1.u0 - c;
    result.u1 = rp1.u1 - c;
    result.u2 = rp1.u2 - c;
    result.u3 = rp1.u3 - c;
    result.hStep = rp1.hStep - c;
    return result;
}

RayProperties operator-(const realNumber& c, const RayProperties& rp1) {
    RayProperties result;
    result.t = c - rp1.t;
    result.r = c - rp1.r;
    result.theta = c - rp1.theta;
    result.phi = c - rp1.phi;
    result.u0 = c - rp1.u0;
    result.u1 = c - rp1.u1;
    result.u2 = c - rp1.u2;
    result.u3 = c - rp1.u3;
    result.hStep = c - rp1.hStep;
    return result;
}

RayProperties operator*(const RayProperties& rp1, const realNumber& c) {
    RayProperties result;
    result.t = rp1.t * c;
    result.r = rp1.r * c;
    result.theta = rp1.theta * c;
    result.phi = rp1.phi * c;
    result.u0 = rp1.u0 * c;
    result.u1 = rp1.u1 * c;
    result.u2 = rp1.u2 * c;
    result.u3 = rp1.u3 * c;
    result.hStep = rp1.hStep * c;
    return result;
}

RayProperties operator*(const realNumber& c, const RayProperties& rp1) {
    RayProperties result;
    result.t = rp1.t * c;
    result.r = rp1.r * c;
    result.theta = rp1.theta * c;
    result.phi = rp1.phi * c;
    result.u0 = rp1.u0 * c;
    result.u1 = rp1.u1 * c;
    result.u2 = rp1.u2 * c;
    result.u3 = rp1.u3 * c;
    result.hStep = rp1.hStep * c;
    return result;
}

RayProperties operator/(const RayProperties& rp1, const realNumber& c) {
    RayProperties result;
    result.t = rp1.t / c;
    result.r = rp1.r / c;
    result.theta = rp1.theta / c;
    result.phi = rp1.phi / c;
    result.u0 = rp1.u0 / c;
    result.u1 = rp1.u1 / c;
    result.u2 = rp1.u2 / c;
    result.u3 = rp1.u3 / c;
    result.hStep = rp1.hStep / c;
    return result;
}

Vector2 operator+(const Vector2& v1, const Vector2& v2) {
    Vector2 result;
    result.x = v1.x + v2.x;
    result.y = v1.y + v2.y;
    return result;
}

Vector3 operator+(const Vector3& v1, const Vector3& v2) {
    Vector3 result;
    result.x = v1.x + v2.x;
    result.y = v1.y + v2.y;
    result.z = v1.z + v2.z;
    return result;
}

Vector4 operator+(const Vector4& v1, const Vector4& v2) {
    Vector4 result;
    result.x = v1.x + v2.x;
    result.y = v1.y + v2.y;
    result.z = v1.z + v2.z;
    result.t = v1.t + v2.t;
    return result;
}

Vector2 operator+(const Vector2& v1, const realNumber& c) {
    Vector2 result;
    result.x = v1.x + c;
    result.y = v1.y + c;
    return result;
}

Vector3 operator+(const Vector3& v1, const realNumber& c) {
    Vector3 result;
    result.x = v1.x + c;
    result.y = v1.y + c;
    result.z = v1.z + c;
    return result;
}

Vector4 operator+(const Vector4& v1, const realNumber& c) {
    Vector4 result;
    result.x = v1.x + c;
    result.y = v1.y + c;
    result.z = v1.z + c;
    result.t = v1.t + c;
    return result;
}

Vector2 operator+(const realNumber& c, const Vector2& v1) {
    Vector2 result;
    result.x = v1.x + c;
    result.y = v1.y + c;
    return result;
}

Vector3 operator+(const realNumber& c, const Vector3& v1) {
    Vector3 result;
    result.x = v1.x + c;
    result.y = v1.y + c;
    result.z = v1.z + c;
    return result;
}

Vector4 operator+(const realNumber& c, const Vector4& v1) {
    Vector4 result;
    result.x = v1.x + c;
    result.y = v1.y + c;
    result.z = v1.z + c;
    result.t = v1.t + c;
    return result;
}

Vector2 operator-(const Vector2& v1, const Vector2& v2) {
    Vector2 result;
    result.x = v1.x - v2.x;
    result.y = v1.y - v2.y;
    return result;
}

Vector3 operator-(const Vector3& v1, const Vector3& v2) {
    Vector3 result;
    result.x = v1.x - v2.x;
    result.y = v1.y - v2.y;
    result.z = v1.z - v2.z;
    return result;
}

Vector4 operator-(const Vector4& v1, const Vector4& v2) {
    Vector4 result;
    result.x = v1.x - v2.x;
    result.y = v1.y - v2.y;
    result.z = v1.z - v2.z;
    result.t = v1.t - v2.t;
    return result;
}

Vector2 operator-(const Vector2& v1, const realNumber& c) {
    Vector2 result;
    result.x = v1.x - c;
    result.y = v1.y - c;
    return result;
}

Vector3 operator-(const Vector3& v1, const realNumber& c) {
    Vector3 result;
    result.x = v1.x - c;
    result.y = v1.y - c;
    result.z = v1.z - c;
    return result;
}

Vector4 operator-(const Vector4& v1, const realNumber& c) {
    Vector4 result;
    result.x = v1.x - c;
    result.y = v1.y - c;
    result.z = v1.z - c;
    result.t = v1.t - c;
    return result;
}

Vector2 operator*(const Vector2& v1, const realNumber& c) {
    Vector2 result;
    result.x = v1.x * c;
    result.y = v1.y * c;
    return result;
}

Vector3 operator*(const Vector3& v1, const realNumber& c) {
    Vector3 result;
    result.x = v1.x * c;
    result.y = v1.y * c;
    result.z = v1.z * c;
    return result;
}

Vector4 operator*(const Vector4& v1, const realNumber& c) {
    Vector4 result;
    result.x = v1.x * c;
    result.y = v1.y * c;
    result.z = v1.z * c;
    result.t = v1.t * c;
    return result;
}

Vector2 operator*(const realNumber& c, const Vector2& v1) {
    Vector2 result;
    result.x = v1.x * c;
    result.y = v1.y * c;
    return result;
}

Vector3 operator*(const realNumber& c, const Vector3& v1) {
    Vector3 result;
    result.x = v1.x * c;
    result.y = v1.y * c;
    result.z = v1.z * c;
    return result;
}

Vector4 operator*(const realNumber& c, const Vector4& v1) {
    Vector4 result;
    result.x = v1.x * c;
    result.y = v1.y * c;
    result.z = v1.z * c;
    result.t = v1.t * c;
    return result;
}

Vector2 operator/(const Vector2& v1, const realNumber& c) {
    Vector2 result;
    result.x = v1.x / c;
    result.y = v1.y / c;
    return result;
}

Vector3 operator/(const Vector3& v1, const realNumber& c) {
    Vector3 result;
    result.x = v1.x / c;
    result.y = v1.y / c;
    result.z = v1.z / c;
    return result;
}

Vector4 operator/(const Vector4& v1, const realNumber& c) {
    Vector4 result;
    result.x = v1.x / c;
    result.y = v1.y / c;
    result.z = v1.z / c;
    result.t = v1.t / c;
    return result;
}

realNumber dotProduct(const Vector2& v1, const Vector2& v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

realNumber dotProduct(const Vector3& v1, const Vector3& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

realNumber dotProduct(const Vector4& v1, const Vector4& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.t * v2.t;
}

realNumber VectorLength(const Vector2& v1) {
    return std::sqrt((v1.x * v1.x) + (v1.y * v1.y));
}

realNumber VectorLength(const Vector3& v1) {
    return std::sqrt((v1.x * v1.x) + (v1.y * v1.y) + (v1.z * v1.z));
}

realNumber VectorLength(const Vector4& v1) {
    return std::sqrt((v1.x * v1.x) + (v1.y * v1.y) + (v1.z * v1.z) + (v1.t * v1.t));
}

Vector2 operator*(const Matrix2& M1, const Vector2& v1) {
    Vector2 result;
    result.x = M1.M00 * v1.x + M1.M01 * v1.y;
    result.y = M1.M10 * v1.x + M1.M11 * v1.y;
    return result;
}

Vector3 operator*(const Matrix3& M1, const Vector3& v1) {
    Vector3 result;
    result.x = M1.M00 * v1.x + M1.M01 * v1.y + M1.M02 * v1.z;
    result.y = M1.M10 * v1.x + M1.M11 * v1.y + M1.M12 * v1.z;
    result.z = M1.M20 * v1.x + M1.M21 * v1.y + M1.M22 * v1.z;
    return result;
}

Vector4 operator*(const Matrix4& M1, const Vector4& v1) {
    Vector4 result;
    result.x = M1.M00 * v1.x + M1.M01 * v1.y + M1.M02 * v1.z + M1.M03 * v1.t;
    result.y = M1.M10 * v1.x + M1.M11 * v1.y + M1.M12 * v1.z + M1.M13 * v1.t;
    result.z = M1.M20 * v1.x + M1.M21 * v1.y + M1.M22 * v1.z + M1.M23 * v1.t;
    result.t = M1.M30 * v1.x + M1.M31 * v1.y + M1.M32 * v1.z + M1.M33 * v1.t;
    return result;
}

Matrix2 operator+(const Matrix2& M1, const realNumber& c) {
    Matrix2 result;
    result.M00 = M1.M00 + c;
    result.M01 = M1.M01 + c;
    result.M10 = M1.M10 + c;
    result.M11 = M1.M11 + c;
    return result;
}

Matrix3 operator+(const Matrix3& M1, const realNumber& c) {
    Matrix3 result;
    result.M00 = M1.M00 + c;
    result.M01 = M1.M01 + c;
    result.M02 = M1.M02 + c;
    result.M10 = M1.M10 + c;
    result.M11 = M1.M11 + c;
    result.M12 = M1.M12 + c;
    result.M20 = M1.M20 + c;
    result.M21 = M1.M21 + c;
    result.M22 = M1.M22 + c;
    return result;
}

Matrix4 operator+(const Matrix4& M1, const realNumber& c) {
    Matrix4 result;
    result.M00 = M1.M00 + c;
    result.M01 = M1.M01 + c;
    result.M02 = M1.M02 + c;
    result.M03 = M1.M03 + c;
    result.M10 = M1.M10 + c;
    result.M11 = M1.M11 + c;
    result.M12 = M1.M12 + c;
    result.M13 = M1.M13 + c;
    result.M20 = M1.M20 + c;
    result.M21 = M1.M21 + c;
    result.M22 = M1.M22 + c;
    result.M23 = M1.M23 + c;
    result.M30 = M1.M30 + c;
    result.M31 = M1.M31 + c;
    result.M32 = M1.M32 + c;
    result.M33 = M1.M33 + c;
    return result;
}

Matrix2 operator+(const realNumber& c, const Matrix2& M1) {
    Matrix2 result;
    result.M00 = M1.M00 + c;
    result.M01 = M1.M01 + c;
    result.M10 = M1.M10 + c;
    result.M11 = M1.M11 + c;
    return result;
}

Matrix3 operator+(const realNumber& c, const Matrix3& M1) {
    Matrix3 result;
    result.M00 = M1.M00 + c;
    result.M01 = M1.M01 + c;
    result.M02 = M1.M02 + c;
    result.M10 = M1.M10 + c;
    result.M11 = M1.M11 + c;
    result.M12 = M1.M12 + c;
    result.M20 = M1.M20 + c;
    result.M21 = M1.M21 + c;
    result.M22 = M1.M22 + c;
    return result;
}

Matrix4 operator+(const realNumber& c, const Matrix4& M1) {
    Matrix4 result;
    result.M00 = M1.M00 + c;
    result.M01 = M1.M01 + c;
    result.M02 = M1.M02 + c;
    result.M03 = M1.M03 + c;
    result.M10 = M1.M10 + c;
    result.M11 = M1.M11 + c;
    result.M12 = M1.M12 + c;
    result.M13 = M1.M13 + c;
    result.M20 = M1.M20 + c;
    result.M21 = M1.M21 + c;
    result.M22 = M1.M22 + c;
    result.M23 = M1.M23 + c;
    result.M30 = M1.M30 + c;
    result.M31 = M1.M31 + c;
    result.M32 = M1.M32 + c;
    result.M33 = M1.M33 + c;
    return result;
}

Matrix2 operator*(const Matrix2& M1, const realNumber& c) {
    Matrix2 result;
    result.M00 = M1.M00 * c;
    result.M01 = M1.M01 * c;
    result.M10 = M1.M10 * c;
    result.M11 = M1.M11 * c;
    return result;
}

Matrix3 operator*(const Matrix3& M1, const realNumber& c) {
    Matrix3 result;
    result.M00 = M1.M00 * c;
    result.M01 = M1.M01 * c;
    result.M02 = M1.M02 * c;
    result.M10 = M1.M10 * c;
    result.M11 = M1.M11 * c;
    result.M12 = M1.M12 * c;
    result.M20 = M1.M20 * c;
    result.M21 = M1.M21 * c;
    result.M22 = M1.M22 * c;
    return result;
}

Matrix4 operator*(const Matrix4& M1, const realNumber& c) {
    Matrix4 result;
    result.M00 = M1.M00 * c;
    result.M01 = M1.M01 * c;
    result.M02 = M1.M02 * c;
    result.M03 = M1.M03 * c;
    result.M10 = M1.M10 * c;
    result.M11 = M1.M11 * c;
    result.M12 = M1.M12 * c;
    result.M13 = M1.M13 * c;
    result.M20 = M1.M20 * c;
    result.M21 = M1.M21 * c;
    result.M22 = M1.M22 * c;
    result.M23 = M1.M23 * c;
    result.M30 = M1.M30 * c;
    result.M31 = M1.M31 * c;
    result.M32 = M1.M32 * c;
    result.M33 = M1.M33 * c;
    return result;
}

Matrix2 operator*(const realNumber& c, const Matrix2& M1) {
    Matrix2 result;
    result.M00 = M1.M00 * c;
    result.M01 = M1.M01 * c;
    result.M10 = M1.M10 * c;
    result.M11 = M1.M11 * c;
    return result;
}

Matrix3 operator*(const realNumber& c, const Matrix3& M1) {
    Matrix3 result;
    result.M00 = M1.M00 * c;
    result.M01 = M1.M01 * c;
    result.M02 = M1.M02 * c;
    result.M10 = M1.M10 * c;
    result.M11 = M1.M11 * c;
    result.M12 = M1.M12 * c;
    result.M20 = M1.M20 * c;
    result.M21 = M1.M21 * c;
    result.M22 = M1.M22 * c;
    return result;
}

Matrix4 operator*(const realNumber& c, const Matrix4& M1) {
    Matrix4 result;
    result.M00 = M1.M00 * c;
    result.M01 = M1.M01 * c;
    result.M02 = M1.M02 * c;
    result.M03 = M1.M03 * c;
    result.M10 = M1.M10 * c;
    result.M11 = M1.M11 * c;
    result.M12 = M1.M12 * c;
    result.M13 = M1.M13 * c;
    result.M20 = M1.M20 * c;
    result.M21 = M1.M21 * c;
    result.M22 = M1.M22 * c;
    result.M23 = M1.M23 * c;
    result.M30 = M1.M30 * c;
    result.M31 = M1.M31 * c;
    result.M32 = M1.M32 * c;
    result.M33 = M1.M33 * c;
    return result;
}

Matrix2 operator+(const Matrix2& M1, const Matrix2& M2) {
    Matrix2 result;
    result.M00 = M1.M00 + M2.M00;
    result.M01 = M1.M01 + M2.M01;
    result.M10 = M1.M10 + M2.M10;
    result.M11 = M1.M11 + M2.M11;
    return result;
}

Matrix3 operator+(const Matrix3& M1, const Matrix3& M2) {
    Matrix3 result;
    result.M00 = M1.M00 + M2.M00;
    result.M01 = M1.M01 + M2.M01;
    result.M02 = M1.M02 + M2.M02;
    result.M10 = M1.M10 + M2.M10;
    result.M11 = M1.M11 + M2.M11;
    result.M12 = M1.M12 + M2.M12;
    result.M20 = M1.M20 + M2.M20;
    result.M21 = M1.M21 + M2.M21;
    result.M22 = M1.M22 + M2.M22;
    return result;
}

Matrix4 operator+(const Matrix4& M1, const Matrix4& M2) {
    Matrix4 result;
    result.M00 = M1.M00 + M2.M00;
    result.M01 = M1.M01 + M2.M01;
    result.M02 = M1.M02 + M2.M02;
    result.M03 = M1.M03 + M2.M03;
    result.M10 = M1.M10 + M2.M10;
    result.M11 = M1.M11 + M2.M11;
    result.M12 = M1.M12 + M2.M12;
    result.M13 = M1.M13 + M2.M13;
    result.M20 = M1.M20 + M2.M20;
    result.M21 = M1.M21 + M2.M21;
    result.M22 = M1.M22 + M2.M22;
    result.M23 = M1.M23 + M2.M23;
    result.M30 = M1.M30 + M2.M30;
    result.M31 = M1.M31 + M2.M31;
    result.M32 = M1.M32 + M2.M32;
    result.M33 = M1.M33 + M2.M33;
    return result;
}

Matrix2 operator*(const Matrix2& M1, const Matrix2& M2) {
    Matrix2 result;
    result.M00 = M1.M00 * M2.M00 + M1.M01 * M2.M10;
    result.M01 = M1.M00 * M2.M01 + M1.M01 * M2.M11;
    result.M10 = M1.M10 * M2.M00 + M1.M11 * M2.M10;
    result.M11 = M1.M10 * M2.M01 + M1.M11 * M2.M11;
    return result;
}

Matrix3 operator*(const Matrix3& M1, const Matrix3& M2) {
    Matrix3 result;
    result.M00 = M1.M00 * M2.M00 + M1.M01 * M2.M10 + M1.M02 * M2.M20;
    result.M01 = M1.M00 * M2.M01 + M1.M01 * M2.M11 + M1.M02 * M2.M21;
    result.M02 = M1.M00 * M2.M02 + M1.M01 * M2.M12 + M1.M02 * M2.M22;
    result.M10 = M1.M10 * M2.M00 + M1.M11 * M2.M10 + M1.M12 * M2.M20;
    result.M11 = M1.M10 * M2.M01 + M1.M11 * M2.M11 + M1.M12 * M2.M21;
    result.M12 = M1.M10 * M2.M02 + M1.M11 * M2.M12 + M1.M12 * M2.M22;
    result.M20 = M1.M20 * M2.M00 + M1.M21 * M2.M10 + M1.M22 * M2.M20;
    result.M21 = M1.M20 * M2.M01 + M1.M21 * M2.M11 + M1.M22 * M2.M21;
    result.M22 = M1.M20 * M2.M02 + M1.M21 * M2.M12 + M1.M22 * M2.M22;
    return result;
}

Matrix4 operator*(const Matrix4& M1, const Matrix4& M2) {
    Matrix4 result;
    result.M00 = M1.M00 * M2.M00 + M1.M01 * M2.M10 + M1.M02 * M2.M20 + M1.M03 * M2.M30;
    result.M01 = M1.M00 * M2.M01 + M1.M01 * M2.M11 + M1.M02 * M2.M21 + M1.M03 * M2.M31;
    result.M02 = M1.M00 * M2.M02 + M1.M01 * M2.M12 + M1.M02 * M2.M22 + M1.M03 * M2.M32;
    result.M03 = M1.M00 * M2.M03 + M1.M01 * M2.M13 + M1.M02 * M2.M23 + M1.M03 * M2.M33;
    result.M10 = M1.M10 * M2.M00 + M1.M11 * M2.M10 + M1.M12 * M2.M20 + M1.M13 * M2.M30;
    result.M11 = M1.M10 * M2.M01 + M1.M11 * M2.M11 + M1.M12 * M2.M21 + M1.M13 * M2.M31;
    result.M12 = M1.M10 * M2.M02 + M1.M11 * M2.M12 + M1.M12 * M2.M22 + M1.M13 * M2.M32;
    result.M13 = M1.M10 * M2.M03 + M1.M11 * M2.M13 + M1.M12 * M2.M23 + M1.M13 * M2.M33;
    result.M20 = M1.M20 * M2.M00 + M1.M21 * M2.M10 + M1.M22 * M2.M20 + M1.M23 * M2.M30;
    result.M21 = M1.M20 * M2.M01 + M1.M21 * M2.M11 + M1.M22 * M2.M21 + M1.M23 * M2.M31;
    result.M22 = M1.M20 * M2.M02 + M1.M21 * M2.M12 + M1.M22 * M2.M22 + M1.M23 * M2.M32;
    result.M23 = M1.M20 * M2.M03 + M1.M21 * M2.M13 + M1.M22 * M2.M23 + M1.M23 * M2.M33;
    result.M30 = M1.M30 * M2.M00 + M1.M31 * M2.M10 + M1.M32 * M2.M20 + M1.M33 * M2.M30;
    result.M31 = M1.M30 * M2.M01 + M1.M31 * M2.M11 + M1.M32 * M2.M21 + M1.M33 * M2.M31;
    result.M32 = M1.M30 * M2.M02 + M1.M31 * M2.M12 + M1.M32 * M2.M22 + M1.M33 * M2.M32;
    result.M33 = M1.M30 * M2.M03 + M1.M31 * M2.M13 + M1.M32 * M2.M23 + M1.M33 * M2.M33;
    return result;
}