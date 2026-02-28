#pragma once
#include <cstdint>
#include <string>
#include "Math.hpp"

namespace evo
{
    inline uint64_t fnv1a64(const std::string &s)
    {
        uint64_t h = 1469598103934665603ULL;
        for (const unsigned char c : s)
        {
            h ^= static_cast<uint64_t>(c);
            h *= 1099511628211ULL;
        }
        return h;
    }

    inline uint64_t mix64(uint64_t x)
    {
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27;
        x *= 0x94d049bb133111ebULL;
        x ^= x >> 31;
        return x;
    }

    inline float hashToUnit(uint64_t x)
    {
        x = mix64(x);
        return static_cast<float>(static_cast<double>(x >> 11) * (1.0 / 9007199254740992.0));
    }

    inline float hashToRange(const uint64_t x, const float a, const float b)
    {
        return lerpf(a, b, hashToUnit(x));
    }

    inline Vec3 hashToVec3(const uint64_t base)
    {
        const Vec3 v{
            hashToRange(base + 11, -1.f, 1.f),
            hashToRange(base + 23, -1.f, 1.f),
            hashToRange(base + 37, -1.f, 1.f),
        };
        return normalize(v);
    }
} // namespace evo