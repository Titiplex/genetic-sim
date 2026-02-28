#pragma once
#include <algorithm>
#include <cmath>

namespace evo
{
    struct Vec3
    {
        float x = 0.f, y = 0.f, z = 0.f;

        constexpr Vec3() = default;

        constexpr Vec3(const float X, const float Y, const float Z): x(X), y(Y), z(Z)
        {
        }

        constexpr Vec3 operator+(const Vec3 &o) const
        {
            return {x + o.x, y + o.y, z + o.z};
        }

        constexpr Vec3 operator-(const Vec3 &o) const
        {
            return {x - o.x, y - o.y, z - o.z};
        }

        constexpr Vec3 operator*(const float s) const
        {
            return {x * s, y * s, z * s};
        }

        constexpr Vec3 operator/(const float s) const
        {
            return {x / s, y / s, z / s};
        }

        Vec3 &operator+=(const Vec3 &o)
        {
            x += o.x;
            y += o.y;
            z += o.z;
            return *this;
        }

        Vec3 &operator-=(const Vec3 &o)
        {
            x -= o.x;
            y -= o.y;
            z -= o.z;
            return *this;
        }

        Vec3 &operator*=(const float s)
        {
            x *= s;
            y *= s;
            z *= s;
            return *this;
        }
    };

    inline float dot(const Vec3 &a, const Vec3 &b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    inline Vec3 cross(const Vec3 &a, const Vec3 &b)
    {
        return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
    }

    inline float len2(const Vec3 &v)
    {
        return dot(v, v);
    }

    inline float len(const Vec3 &v)
    {
        return std::sqrt(len2(v));
    }

    inline Vec3 normalize(const Vec3 &v)
    {
        const float l = len(v);
        return l > 1e-6f ? v / l : Vec3{0.f, 0.f, 0.f};
    }

    inline float clampf(const float x, const float a, const float b)
    {
        return std::max(a, std::min(b, x));
    }

    inline float lerpf(const float a, const float b, const float t)
    {
        return a + (b - a) * t;
    }

    inline Vec3 clampVec(const Vec3 &v, const float maxLen)
    {
        if (const float l = len(v); l > maxLen && l > 1e-6f)
            return v * (maxLen / l);
        return v;
    }

    struct Mat4
    {
        float m[16]{};

        static Mat4 identity()
        {
            Mat4 r;
            r.m[0] = r.m[5] = r.m[10] = r.m[15] = 1.f;
            return r;
        }
    };

    inline Mat4 perspective(const float fovyRad, const float aspect, const float znear, const float zfar)
    {
        Mat4        r{};
        const float t = std::tan(fovyRad * 0.5f);
        r.m[0]        = 1.0f / (aspect * t);
        r.m[5]        = 1.0f / t;
        r.m[10]       = -(zfar + znear) / (zfar - znear);
        r.m[11]       = -1.0f;
        r.m[14]       = -(2.0f * zfar * znear) / (zfar - znear);
        return r;
    }

    inline Mat4 lookAt(const Vec3 &eye, const Vec3 &center, const Vec3 &up)
    {
        const Vec3 f = normalize(center - eye);
        const Vec3 s = normalize(cross(f, up));
        const Vec3 u = cross(s, f);

        Mat4 r  = Mat4::identity();
        r.m[0]  = s.x;
        r.m[4]  = s.y;
        r.m[8]  = s.z;
        r.m[1]  = u.x;
        r.m[5]  = u.y;
        r.m[9]  = u.z;
        r.m[2]  = -f.x;
        r.m[6]  = -f.y;
        r.m[10] = -f.z;
        r.m[12] = -dot(s, eye);
        r.m[13] = -dot(u, eye);
        r.m[14] = dot(f, eye);
        return r;
    }

    inline Mat4 mul(const Mat4 &a, const Mat4 &b)
    {
        Mat4 r{};
        for (int c = 0; c < 4; ++c)
        {
            for (int r0 = 0; r0 < 4; ++r0)
            {
                r.m[c * 4 + r0] = a.m[0 * 4 + r0] * b.m[c * 4 + 0] + a.m[1 * 4 + r0] * b.m[c * 4 + 1] +
                                  a.m[2 * 4 + r0] * b.m[c * 4 + 2] + a.m[3 * 4 + r0] * b.m[c * 4 + 3];
            }
        }
        return r;
    }
} // namespace evo