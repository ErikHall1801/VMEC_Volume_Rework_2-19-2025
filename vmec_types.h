#pragma once

#include "pch_vmec_primary.h"

#include "vmec_defines.h"

#pragma pack(push, 1)

struct pix_RGB {
    unsigned char R;
    unsigned char G;
    unsigned char B;
};

struct pix_realNumber {
    realNumber R;
    realNumber G;
    realNumber B;
};

struct rgba_f32 {
    f32 R;
    f32 G;
    f32 B;
    f32 A;
};

struct pixel_f32 {
    f32 R;
    f32 G;
    f32 B;
};

struct noise_point {
    realNumber r;
    realNumber t;
    realNumber p;
    realNumber v;
};

struct lattice_noise_point {
    realNumber r;
    realNumber t;
    realNumber p;
    realNumber v;
    realNumber W;
    realNumber s;
};

struct dynamic_lattice_point {
    realNumber r;
    realNumber theta;
    realNumber phi;
    realNumber tau;
};

//Properties structs
struct localProperties {
    realNumber t;
    realNumber r;
    realNumber theta;
    realNumber phi;
    realNumber u0;
    realNumber u1;
    realNumber u2;
    realNumber u3;
};

//reads tRay, rRay, thetaRay, phiRay, u0, u1, u2, u3, hStep
struct RayProperties {
    realNumber t;
    realNumber r;
    realNumber theta;
    realNumber phi;
    realNumber u0;
    realNumber u1;
    realNumber u2;
    realNumber u3;
    realNumber hStep;
};

struct BVObjects {
    bool disk;
    bool jet;
    bool ambient;
};

RayProperties operator+(const RayProperties& rp1, const RayProperties& rp2);

RayProperties operator-(const RayProperties& rp1, const RayProperties& rp2);

RayProperties operator+(const RayProperties& rp1, const realNumber& c);

RayProperties operator+(const realNumber& c, const RayProperties& rp1);

RayProperties operator-(const RayProperties& rp1, const realNumber& c);

RayProperties operator-(const realNumber& c, const RayProperties& rp1);

RayProperties operator*(const RayProperties& rp1, const realNumber& c);

RayProperties operator*(const realNumber& c, const RayProperties& rp1);

RayProperties operator/(const RayProperties& rp1, const realNumber& c);


//Vectors and matrices
struct Vector2 { realNumber x; realNumber y; };

struct Vector3 
{ 
    realNumber x;
    realNumber y;
    realNumber z;
    Vector3() : x(0), y(0), z(0) {} 
    
    Vector3(realNumber x_val, realNumber y_val, realNumber z_val) : x(x_val), y(y_val), z(z_val) {} 
    
    void set(realNumber x_val, realNumber y_val, realNumber z_val) {x = x_val; y = y_val; z = z_val;} 
    
    realNumber operator[](int i) const {switch (i) { case 0: return x; case 1: return y; case 2: return z; default: throw std::out_of_range("Illegal Index"); } }

    realNumber& operator[](int i) {switch (i) { case 0: return x; case 1: return y; case 2: return z; default: throw std::out_of_range("Illegal Index"); } }
};

struct ray_scene_information {
    bool ray_convergence;
    int i_rs;
    int i_celestial_sphere;
    
    bool hit_disk;
    int i_hit_disk;
    Vector3 hit_disk_pos;

    bool hit_jet;
    int i_hit_jet;
    Vector3 hit_jet_pos;

    bool hit_ambient;
    int i_hit_ambient;
    Vector3 hit_ambient_pos;
};

struct Vector4 { realNumber x; realNumber y; realNumber z; realNumber t; Vector4() : x(0), y(0), z(0), t(0) {} Vector4(realNumber x_val, realNumber y_val, realNumber z_val, realNumber t_val) : x(x_val), y(y_val), z(z_val), t(t_val) {} void set(realNumber x_val, realNumber y_val, realNumber z_val, realNumber t_val) {x = x_val; y = y_val; z = z_val; t = t_val;} };
/*                 | M00  M01  M02 |
   Matrix layout:  | M10  M11  M12 |
                   | M20  M21  M22 |   */
struct Matrix2 {
    realNumber M00; realNumber M01;
    realNumber M10; realNumber M11;
};
struct Matrix3 {
    realNumber M00; realNumber M01; realNumber M02;
    realNumber M10; realNumber M11; realNumber M12;
    realNumber M20; realNumber M21; realNumber M22;
};
struct Matrix4 {
    realNumber M00; realNumber M01; realNumber M02; realNumber M03;
    realNumber M10; realNumber M11; realNumber M12; realNumber M13;
    realNumber M20; realNumber M21; realNumber M22; realNumber M23;
    realNumber M30; realNumber M31; realNumber M32; realNumber M33;
};

struct SensorProperties {
    int maxWidth;
    int maxHeight;
    realNumber FocalLength;
    realNumber SensorSize;
    realNumber AngleX;
    realNumber AngleY;
    realNumber AngleZ;
    Vector3 camPos;
    Vector3 camMomentum;
};

Vector2 operator+(const Vector2& v1, const Vector2& v2);

Vector3 operator+(const Vector3& v1, const Vector3& v2);

Vector4 operator+(const Vector4& v1, const Vector4& v2);

Vector2 operator+(const Vector2& v1, const realNumber& c);

Vector3 operator+(const Vector3& v1, const realNumber& c);

Vector4 operator+(const Vector4& v1, const realNumber& c);

Vector2 operator+(const realNumber& c, const Vector2& v1);

Vector3 operator+(const realNumber& c, const Vector3& v1);

Vector4 operator+(const realNumber& c, const Vector4& v1);

Vector2 operator-(const Vector2& v1, const Vector2& v2);

Vector3 operator-(const Vector3& v1, const Vector3& v2);

Vector4 operator-(const Vector4& v1, const Vector4& v2);

Vector2 operator-(const Vector2& v1, const realNumber& c);

Vector3 operator-(const Vector3& v1, const realNumber& c);

Vector4 operator-(const Vector4& v1, const realNumber& c);

Vector2 operator*(const Vector2& v1, const realNumber& c);

Vector3 operator*(const Vector3& v1, const realNumber& c);

Vector4 operator*(const Vector4& v1, const realNumber& c);

Vector2 operator*(const realNumber& c, const Vector2& v1);

Vector3 operator*(const realNumber& c, const Vector3& v1);

Vector4 operator*(const realNumber& c, const Vector4& v1);

Vector2 operator/(const Vector2& v1, const realNumber& c);

Vector3 operator/(const Vector3& v1, const realNumber& c);

Vector4 operator/(const Vector4& v1, const realNumber& c);

realNumber dotProduct(const Vector2& v1, const Vector2& v2);

realNumber dotProduct(const Vector3& v1, const Vector3& v2);

realNumber dotProduct(const Vector4& v1, const Vector4& v2);

realNumber VectorLength(const Vector2& v1);

realNumber VectorLength(const Vector3& v1);

realNumber VectorLength(const Vector4& v1);

Vector3 VectorCross(const Vector3& v1, const Vector3& v2);

Vector2 operator*(const Matrix2& M1, const Vector2& v1);

Vector3 operator*(const Matrix3& M1, const Vector3& v1);

Vector4 operator*(const Matrix4& M1, const Vector4& v1);

Matrix2 operator+(const Matrix2& M1, const realNumber& c);

Matrix3 operator+(const Matrix3& M1, const realNumber& c);

Matrix4 operator+(const Matrix4& M1, const realNumber& c);

Matrix2 operator+(const realNumber& c, const Matrix2& M1);

Matrix3 operator+(const realNumber& c, const Matrix3& M1);

Matrix4 operator+(const realNumber& c, const Matrix4& M1);

Matrix2 operator*(const Matrix2& M1, const realNumber& c);

Matrix3 operator*(const Matrix3& M1, const realNumber& c);

Matrix4 operator*(const Matrix4& M1, const realNumber& c);

Matrix2 operator*(const realNumber& c, const Matrix2& M1);

Matrix3 operator*(const realNumber& c, const Matrix3& M1);

Matrix4 operator*(const realNumber& c, const Matrix4& M1);

Matrix2 operator+(const Matrix2& M1, const Matrix2& M2);

Matrix3 operator+(const Matrix3& M1, const Matrix3& M2);

Matrix4 operator+(const Matrix4& M1, const Matrix4& M2);

Matrix2 operator*(const Matrix2& M1, const Matrix2& M2);

Matrix3 operator*(const Matrix3& M1, const Matrix3& M2);

Matrix4 operator*(const Matrix4& M1, const Matrix4& M2);

#pragma pack(pop)
