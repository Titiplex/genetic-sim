#pragma once

#include <random>
#include <string>

#include "maths/Vec3.hpp"

static uint64_t fnv1a64(const std::string& s)
{
    uint64_t h = 1469598103934665603ULL;
    for (const unsigned char c : s)
    {
        h ^= static_cast<uint64_t>(c);
        h *= 1099511628211ULL;
    }
    return h;
}

static uint64_t mix64(uint64_t x)
{
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

static float hashToUnit(uint64_t x)
{
    x = mix64(x);
    return static_cast<float>(static_cast<double>(x >> 11) * (1.0 / 9007199254740992.0));
}

static float hashToRange(const uint64_t x, const float a, const float b)
{
    return lerpf(a, b, hashToUnit(x));
}

static Vec3 hashToVec3(const uint64_t base)
{
    const Vec3 v{
        hashToRange(base + 11, -1.f, 1.f),
        hashToRange(base + 23, -1.f, 1.f),
        hashToRange(base + 37, -1.f, 1.f)
    };
    return normalize(v);
}

static std::mt19937_64 g_rng(1234567ULL);

static float rndf(const float a = 0.f, const float b = 1.f)
{
    std::uniform_real_distribution<float> d(a, b);
    return d(g_rng);
}

static int rndi(const int a, const int b)
{
    std::uniform_int_distribution<int> d(a, b);
    return d(g_rng);
}

// ======================= Spatial Hash ==========================
struct IVec3
{
    int x, y, z;
};

static bool operator==(const IVec3& a, const IVec3& b) { return a.x == b.x && a.y == b.y && a.z == b.z; }

struct IVec3Hash
{
    size_t operator()(const IVec3& k) const noexcept
    {
        const uint64_t h = static_cast<uint64_t>(k.x * 73856093)
            ^ static_cast<uint64_t>(k.y * 19349663)
            ^ static_cast<uint64_t>(k.z * 83492791);
        return h;
    }
};
