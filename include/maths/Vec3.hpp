#pragma once
#include <algorithm>
#include <cmath>

struct Vec3
{
    float x = 0.f, y = 0.f, z = 0.f;
    Vec3() = default;

    Vec3(const float X, const float Y, const float Z): x(X), y(Y), z(Z)
    {
    }

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(const float s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(const float s) const { return {x / s, y / s, z / s}; }

    Vec3& operator+=(const Vec3& o)
    {
        x += o.x;
        y += o.y;
        z += o.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& o)
    {
        x -= o.x;
        y -= o.y;
        z -= o.z;
        return *this;
    }

    Vec3& operator*=(const float s)
    {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }
};

static float dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static float len2(const Vec3& v) { return dot(v, v); }
static float len(const Vec3& v) { return std::sqrt(len2(v)); }

static Vec3 normalize(const Vec3& v)
{
    const float l = len(v);
    return l > 1e-6f ? v / l : Vec3{0.f, 0.f, 0.f};
}

static Vec3 clampVec(const Vec3& v, const float m)
{
    if (const float l = len(v); l > m && l > 1e-6f) return v * (m / l);
    return v;
}

static float clampf(const float x, const float a, const float b) { return std::max(a, std::min(b, x)); }
static float lerpf(const float a, const float b, const float t) { return a + (b - a) * t; }

// cross product helper (cleaner than repeating formulas)
static Vec3 cross(const Vec3& a, const Vec3& b)
{
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}
