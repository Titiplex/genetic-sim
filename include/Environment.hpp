#pragma once

#include "maths/Vec3.hpp"
#include <algorithm>
#include <cstdint>
#include <vector>
#include <maths/Hash.hpp>

enum class Scenario : uint8_t
{
    Baseline,
    ClimateShock,
    HabitatFragmentation,
    InvasiveSpecies
};

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
    float temperatureBaseline = 0.52f;
    float humidityBaseline = 0.55f;
    Scenario scenario = Scenario::Baseline;

    std::vector<BiomassPulse> biomassPulses;

    [[nodiscard]] float seasonal() const { return 0.5f + 0.5f * std::sin(time * 0.035f); }
    [[nodiscard]] float droughtCycle() const { return 0.5f + 0.5f * std::sin(time * 0.022f + 1.4f); }

    [[nodiscard]] float temperature(const Vec3& p) const
    {
        float lat = 0.5f + 0.5f * (p.z / worldRadius);
        float alt = clampf((p.y + worldRadius) / (2.f * worldRadius), 0.f, 1.f);
        float t = temperatureBaseline + 0.35f * seasonal() + 0.2f * lat - 0.28f * alt;
        if (scenario == Scenario::ClimateShock) t += 0.28f * std::sin(time * 0.12f);
        return clampf(t, 0.f, 1.2f);
    }

    [[nodiscard]] float humidity(const Vec3& p) const
    {
        float h = humidityBaseline + 0.35f * shorelineBlend(p) + 0.15f * std::cos(0.03f * p.x + time * 0.3f) - 0.35f * droughtCycle();
        if (scenario == Scenario::ClimateShock) h -= 0.1f;
        return clampf(h, 0.f, 1.2f);
    }

    [[nodiscard]] float pH(const Vec3& p) const
    {
        const float acidPatch = 0.5f + 0.5f * std::sin(0.04f * p.x - 0.03f * p.z + time * 0.4f);
        return clampf(0.45f + 0.4f * acidPatch + 0.15f * (1.f - humidity(p)), 0.f, 1.f);
    }

    [[nodiscard]] float predationPressure(const Vec3& p) const
    {
        const Vec3 hunter{28.f * std::sin(time * 0.06f), -4.f, 24.f * std::cos(time * 0.05f)};
        float val = std::exp(-len2(p - hunter) / (2.f * 14.f * 14.f));
        if (scenario == Scenario::InvasiveSpecies) val += 0.35f;
        return clampf(val, 0.f, 1.5f);
    }

    [[nodiscard]] float toxin(const Vec3& p) const
    {
        const float r = len(p);
        const float shell = std::exp(-((r - 58.f) * (r - 58.f)) / (2.f * 5.f * 5.f));
        const Vec3 hot{30.f * std::sin(time * 0.08f), 0.f, -30.f * std::cos(time * 0.06f)};
        const float hotspot = 1.6f * std::exp(-len2(p - hot) / (2.f * 10.f * 10.f));
        const float dryBoost = 1.25f - 0.45f * shorelineBlend(p);
        return clampf((0.8f * shell + hotspot + 0.35f * (1.f - humidity(p))) * dryBoost, 0.f, 2.4f);
    }

    [[nodiscard]] float waterSurfaceHeight(const Vec3& p) const
    {
        return fluidLevel + 1.8f * std::sin(0.045f * p.x + time * 0.85f) + 1.2f * std::cos(0.04f * p.z - time * 0.65f);
    }

    [[nodiscard]] float shorelineBlend(const Vec3& p) const
    {
        const float shoreBand = 4.2f;
        return clampf((waterSurfaceHeight(p) - groundHeight(p) + shoreBand) / (2.f * shoreBand), 0.f, 1.f);
    }

    [[nodiscard]] float nutrient(const Vec3& p) const
    {
        const Vec3 c1{20.f * std::sin(time * 0.07f), 8.f * std::cos(time * 0.05f), 20.f * std::cos(time * 0.09f)};
        const Vec3 c2{-25.f * std::cos(time * 0.04f), -10.f * std::sin(time * 0.06f), 15.f * std::sin(time * 0.03f)};
        const float d1 = len2(p - c1);
        const float d2 = len2(p - c2);
        const float wetBoost = 0.7f + 0.6f * shorelineBlend(p);
        const float bloom = 0.5f + 0.5f * std::sin(time * 0.09f + 0.04f * p.x);
        const float n = wetBoost * bloom * (1.8f * std::exp(-d1 / (2.f * 18.f * 18.f)) + 1.4f * std::exp(-d2 / (2.f * 14.f * 14.f)));
        return clampf(n, 0.f, 3.0f);
    }

    [[nodiscard]] float light(const Vec3& p) const
    {
        const float y = (p.y + worldRadius) / (2.f * worldRadius);
        const float day = 0.65f + 0.35f * std::sin(time * 0.3f);
        const float underwater = clampf((waterSurfaceHeight(p) - p.y) / 18.f, 0.f, 1.f);
        const float attenuation = 1.f - 0.55f * underwater;
        return clampf((0.1f + 1.2f * y) * day * attenuation, 0.f, 1.8f);
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

    [[nodiscard]] float obstacleField(const Vec3& p) const
    {
        float ridge = std::abs(std::sin(0.09f * p.x) * std::cos(0.08f * p.z));
        if (scenario == Scenario::HabitatFragmentation) ridge += 0.35f;
        return clampf(ridge, 0.f, 1.f);
    }

    [[nodiscard]] float corridorField(const Vec3& p) const
    {
        const float c = std::exp(-std::pow(std::sin(0.015f * p.x + 0.023f * p.z), 2.f) * 7.f);
        return clampf(c, 0.f, 1.f);
    }

    [[nodiscard]] float groundHeight(const Vec3& p) const
    {
        const float ridge = 6.0f * std::sin(0.055f * p.x + 0.2f * std::sin(time * 0.08f)) + 4.5f * std::cos(0.045f * p.z - 0.15f * std::cos(time * 0.07f));
        const float basin = -22.f + 0.08f * len(Vec3{p.x, 0.f, p.z});
        return basin + ridge + 2.2f * obstacleField(p);
    }

    [[nodiscard]] Vec3 fluidVelocity(const Vec3& p) const
    {
        const float yRatio = clampf((p.y + worldRadius) / (2.f * worldRadius), 0.f, 1.f);
        const Vec3 swirlCenter{24.f * std::sin(time * 0.08f), fluidLevel - 6.f + 4.f * std::sin(time * 0.06f), 22.f * std::cos(time * 0.05f)};
        const Vec3 radial = p - swirlCenter;
        Vec3 tangent{-radial.z, 0.f, radial.x};
        tangent = normalize(tangent);
        const float vortex = std::exp(-len2(radial) / (2.f * 34.f * 34.f));
        const float stream = 0.8f * std::sin(0.03f * p.z + time * 0.9f) + 0.65f * std::cos(0.035f * p.x - time * 0.6f);
        return tangent * (1.8f * vortex) + Vec3{stream * (0.25f + 0.75f * yRatio), 0.15f * std::sin(0.04f * p.x + time * 0.7f), 0.f};
    }

    [[nodiscard]] Vec3 turbulence(const Vec3& p, const uint64_t seed) const
    {
        const float sx = std::sin(0.11f * p.x + time * 1.7f + hashToRange(seed, -3.f, 3.f));
        const float sy = std::sin(0.08f * p.y - time * 1.4f + hashToRange(seed + 19, -3.f, 3.f));
        const float sz = std::sin(0.10f * p.z + time * 1.9f + hashToRange(seed + 41, -3.f, 3.f));
        return Vec3{sx, sy, sz} * (0.22f + 0.15f * obstacleField(p));
    }

    [[nodiscard]] float buoyancy(const Vec3& p, const float bodyScale) const
    {
        const float depth = waterSurfaceHeight(p) - p.y;
        const float submersion = clampf(depth / 22.f, 0.f, 1.f);
        const float tempEffect = 0.9f + 0.2f * temperature(p);
        const float lift = (0.20f + 0.55f * submersion) * (1.15f - 0.18f * bodyScale) * tempEffect;
        const float gravity = 0.15f + 0.05f * bodyScale;
        return lift - gravity;
    }

    [[nodiscard]] Vec3 signalGradient(const Vec3& p, const float phase) const
    {
        constexpr float e = 0.7f;
        auto signal = [&](const Vec3& q)
        {
            return 0.6f * nutrient(q) + 0.4f * biomass(q) - 0.35f * toxin(q) + 0.25f * std::sin(phase + 0.03f * q.x);
        };
        return Vec3{signal({p.x + e, p.y, p.z}) - signal({p.x - e, p.y, p.z}), signal({p.x, p.y + e, p.z}) - signal({p.x, p.y - e, p.z}), signal({p.x, p.y, p.z + e}) - signal({p.x, p.y, p.z - e})} / (2.f * e);
    }

    [[nodiscard]] Vec3 gradN(const Vec3& p) const { constexpr float e = 0.8f; return Vec3{nutrient({p.x + e, p.y, p.z}) - nutrient({p.x - e, p.y, p.z}), nutrient({p.x, p.y + e, p.z}) - nutrient({p.x, p.y - e, p.z}), nutrient({p.x, p.y, p.z + e}) - nutrient({p.x, p.y, p.z - e})} / (2.f * e); }
    [[nodiscard]] Vec3 gradL(const Vec3& p) const { constexpr float e = 0.8f; return Vec3{light({p.x + e, p.y, p.z}) - light({p.x - e, p.y, p.z}), light({p.x, p.y + e, p.z}) - light({p.x, p.y - e, p.z}), light({p.x, p.y, p.z + e}) - light({p.x, p.y, p.z - e})} / (2.f * e); }
    [[nodiscard]] Vec3 gradX(const Vec3& p) const { constexpr float e = 0.8f; return Vec3{toxin({p.x + e, p.y, p.z}) - toxin({p.x - e, p.y, p.z}), toxin({p.x, p.y + e, p.z}) - toxin({p.x, p.y - e, p.z}), toxin({p.x, p.y, p.z + e}) - toxin({p.x, p.y, p.z - e})} / (2.f * e); }
    [[nodiscard]] Vec3 gradB(const Vec3& p) const { constexpr float e = 0.8f; return Vec3{biomass({p.x + e, p.y, p.z}) - biomass({p.x - e, p.y, p.z}), biomass({p.x, p.y + e, p.z}) - biomass({p.x, p.y - e, p.z}), biomass({p.x, p.y, p.z + e}) - biomass({p.x, p.y, p.z - e})} / (2.f * e); }

    void step(const float dt)
    {
        time += dt;
        fluidLevel = 4.8f + 2.8f * std::sin(time * 0.04f);
        for (auto& p : biomassPulses)
        {
            p.age += dt;
            p.amount *= std::exp(-0.018f * dt);
        }
        biomassPulses.erase(std::remove_if(biomassPulses.begin(), biomassPulses.end(), [](const BiomassPulse& p) { return p.age > p.ttl || p.amount < 0.001f; }), biomassPulses.end());
    }

    void depositBiomass(const Vec3& pos, float amount, float sigma = 5.f)
    {
        if (amount <= 0.f) return;
        biomassPulses.push_back({pos, amount, sigma, 0.f, 80.f + rndf(0.f, 40.f)});
    }

    void consumeBiomassNear(const Vec3& pos, float requested)
    {
        if (requested <= 0.f || biomassPulses.empty()) return;
        for (auto& pulse : biomassPulses)
        {
            if (requested <= 0.f) break;
            const float d2 = len2(pos - pulse.pos);
            const float influence = std::exp(-d2 / (2.f * pulse.sigma * pulse.sigma));
            if (influence < 0.02f) continue;
            const float take = std::min(pulse.amount, requested * influence * 0.35f);
            pulse.amount -= take;
            requested -= take;
        }
    }
};
