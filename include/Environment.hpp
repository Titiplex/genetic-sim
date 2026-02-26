#pragma once

#include "maths/Vec3.hpp"
#include <vector>
#include <maths/Hash.hpp>

// ======================= Environment (3 factors + biomass) ==========================
// 1) Nutrients N(x,y,z)
// 2) Light L(x,y,z)
// 3) Toxin X(x,y,z)
// + Biomass B from dead organisms (stored as pulses)
struct BiomassPulse
{
    Vec3  pos;
    float amount = 0.f;
    float sigma  = 5.f;
    float age    = 0.f;
    float ttl    = 80.f;
};

struct Environment
{
    float worldRadius = 80.f;
    float time        = 0.f;
    float fluidLevel  = 6.f;

    std::vector<BiomassPulse> biomassPulses;

    [[nodiscard]] float nutrient(const Vec3& p) const
    {
        const Vec3  c1{20.f * std::sin(time * 0.07f), 8.f * std::cos(time * 0.05f), 20.f * std::cos(time * 0.09f)};
        const Vec3  c2{-25.f * std::cos(time * 0.04f), -10.f * std::sin(time * 0.06f), 15.f * std::sin(time * 0.03f)};
        const float d1 = len2(p - c1);
        const float d2 = len2(p - c2);
        const float n  = 1.8f * std::exp(-d1 / (2.f * 18.f * 18.f))
            + 1.4f * std::exp(-d2 / (2.f * 14.f * 14.f));
        return clampf(n, 0.f, 2.5f);
    }

    [[nodiscard]] float light(const Vec3& p) const
    {
        const float y   = (p.y + worldRadius) / (2.f * worldRadius);
        const float day = 0.65f + 0.35f * std::sin(time * 0.3f);
        return clampf((0.1f + 1.2f * y) * day, 0.f, 1.8f);
    }

    [[nodiscard]] float toxin(const Vec3& p) const
    {
        const float r     = len(p);
        const float shell = std::exp(-((r - 58.f) * (r - 58.f)) / (2.f * 5.f * 5.f));
        const Vec3  hot{30.f * std::sin(time * 0.08f), 0.f, -30.f * std::cos(time * 0.06f)};
        const float hotspot = 1.6f * std::exp(-len2(p - hot) / (2.f * 10.f * 10.f));
        return clampf(0.8f * shell + hotspot, 0.f, 2.2f);
    }

    [[nodiscard]] float biomass(const Vec3& p) const
    {
        float b = 0.f;
        for (const auto& pulse : biomassPulses)
        {
            if (pulse.amount <= 0.f) continue;
            const float s2 = pulse.sigma * pulse.sigma;
            b += pulse.amount * std::exp(-len2(p - pulse.pos) / (2.f * s2));
        }
        return clampf(b, 0.f, 4.f);
    }

    [[nodiscard]] float groundHeight(const Vec3& p) const
    {
        const float ridge = 6.0f * std::sin(0.055f * p.x + 0.2f * std::sin(time * 0.08f))
            + 4.5f * std::cos(0.045f * p.z - 0.15f * std::cos(time * 0.07f));
        const float basin = -22.f + 0.08f * len(Vec3{p.x, 0.f, p.z});
        return basin + ridge;
    }

    [[nodiscard]] Vec3 fluidVelocity(const Vec3& p) const
    {
        const float yRatio = clampf((p.y + worldRadius) / (2.f * worldRadius), 0.f, 1.f);

        const Vec3 swirlCenter{
            24.f * std::sin(time * 0.08f),
            fluidLevel - 6.f + 4.f * std::sin(time * 0.06f),
            22.f * std::cos(time * 0.05f)
        };

        const Vec3 radial = p - swirlCenter;
        Vec3       tangent{-radial.z, 0.f, radial.x};
        tangent = normalize(tangent);

        const float vortex = std::exp(-len2(radial) / (2.f * 34.f * 34.f));
        const float stream = 0.8f * std::sin(0.03f * p.z + time * 0.9f)
            + 0.65f * std::cos(0.035f * p.x - time * 0.6f);

        return tangent * (1.8f * vortex)
            + Vec3{stream * (0.25f + 0.75f * yRatio), 0.15f * std::sin(0.04f * p.x + time * 0.7f), 0.f};
    }

    [[nodiscard]] Vec3 turbulence(const Vec3& p, const uint64_t seed) const
    {
        const float sx = std::sin(0.11f * p.x + time * 1.7f + hashToRange(seed, -3.f, 3.f));
        const float sy = std::sin(0.08f * p.y - time * 1.4f + hashToRange(seed + 19, -3.f, 3.f));
        const float sz = std::sin(0.10f * p.z + time * 1.9f + hashToRange(seed + 41, -3.f, 3.f));
        return Vec3{sx, sy, sz} * 0.22f;
    }

    [[nodiscard]] float buoyancy(const Vec3& p, const float bodyScale) const
    {
        const float depth      = fluidLevel - p.y;
        const float submersion = clampf(depth / 22.f, 0.f, 1.f);
        const float lift       = (0.20f + 0.55f * submersion) * (1.15f - 0.18f * bodyScale);
        const float gravity    = 0.15f + 0.05f * bodyScale;
        return lift - gravity;
    }

    [[nodiscard]] Vec3 gradN(const Vec3& p) const
    {
        constexpr float e = 0.8f;
        return Vec3{
            nutrient({p.x + e, p.y, p.z}) - nutrient({p.x - e, p.y, p.z}),
            nutrient({p.x, p.y + e, p.z}) - nutrient({p.x, p.y - e, p.z}),
            nutrient({p.x, p.y, p.z + e}) - nutrient({p.x, p.y, p.z - e})
        } / (2.f * e);
    }

    [[nodiscard]] Vec3 gradL(const Vec3& p) const
    {
        constexpr float e = 0.8f;
        return Vec3{
            light({p.x + e, p.y, p.z}) - light({p.x - e, p.y, p.z}),
            light({p.x, p.y + e, p.z}) - light({p.x, p.y - e, p.z}),
            light({p.x, p.y, p.z + e}) - light({p.x, p.y, p.z - e})
        } / (2.f * e);
    }

    [[nodiscard]] Vec3 gradX(const Vec3& p) const
    {
        constexpr float e = 0.8f;
        return Vec3{
            toxin({p.x + e, p.y, p.z}) - toxin({p.x - e, p.y, p.z}),
            toxin({p.x, p.y + e, p.z}) - toxin({p.x, p.y - e, p.z}),
            toxin({p.x, p.y, p.z + e}) - toxin({p.x, p.y, p.z - e})
        } / (2.f * e);
    }

    [[nodiscard]] Vec3 gradB(const Vec3& p) const
    {
        constexpr float e = 0.8f;
        return Vec3{
            biomass({p.x + e, p.y, p.z}) - biomass({p.x - e, p.y, p.z}),
            biomass({p.x, p.y + e, p.z}) - biomass({p.x, p.y - e, p.z}),
            biomass({p.x, p.y, p.z + e}) - biomass({p.x, p.y, p.z - e})
        } / (2.f * e);
    }

    void step(const float dt)
    {
        time += dt;
        for (auto& p : biomassPulses)
        {
            p.age += dt;
            // decay amount over time
            p.amount *= std::exp(-0.018f * dt);
        }

        biomassPulses.erase(
            std::remove_if(biomassPulses.begin(), biomassPulses.end(),
                           [](const BiomassPulse& p)
                           {
                               return p.age > p.ttl || p.amount < 0.001f;
                           }),
            biomassPulses.end()
        );
    }

    void depositBiomass(const Vec3& pos, float amount, float sigma = 5.f)
    {
        if (amount <= 0.f) return;
        biomassPulses.push_back({pos, amount, sigma, 0.f, 80.f + rndf(0.f, 40.f)});
    }

    void consumeBiomassNear(const Vec3& pos, float requested)
    {
        if (requested <= 0.f || biomassPulses.empty()) return;

        // remove from nearby pulses first
        for (auto& pulse : biomassPulses)
        {
            const float d2        = len2(pos - pulse.pos);
            const float influence = std::exp(-d2 / (2.f * pulse.sigma * pulse.sigma));
            if (influence < 0.02f) continue;

            const float take = std::min(pulse.amount, requested * influence * 0.35f);
            pulse.amount -= take;
        }
    }
};
