#include "evo/Environment.hpp"
#include <algorithm>
#include <cmath>

namespace evo
{
    float Environment::nutrient(const Vec3 &p) const
    {
        const Vec3 c1{20.f * std::sin(time * 0.07f), 8.f * std::cos(time * 0.05f),
                      20.f * std::cos(time * 0.09f)};
        const Vec3 c2{-25.f * std::cos(time * 0.04f), -10.f * std::sin(time * 0.06f),
                      15.f * std::sin(time * 0.03f)};
        const float d1 = len2(p - c1);
        const float d2 = len2(p - c2);
        const float n  = 1.8f * std::exp(-d1 / (2.f * 18.f * 18.f)) + 1.4f * std::exp(
                            -d2 / (2.f * 14.f * 14.f));
        return clampf(n, 0.f, 2.5f);
    }

    float Environment::light(const Vec3 &p) const
    {
        const float y   = (p.y + worldRadius) / (2.f * worldRadius);
        const float day = 0.65f + 0.35f * std::sin(time * 0.3f);
        return clampf((0.1f + 1.2f * y) * day, 0.f, 1.8f);
    }

    float Environment::toxin(const Vec3 &p) const
    {
        const float r     = len(p);
        const float shell = std::exp(-((r - 58.f) * (r - 58.f)) / (2.f * 5.f * 5.f));
        const Vec3  hot{30.f * std::sin(time * 0.08f), 0.f, -30.f * std::cos(time * 0.06f)};
        const float hotspot = 1.6f * std::exp(-len2(p - hot) / (2.f * 10.f * 10.f));
        return clampf(0.8f * shell + hotspot, 0.f, 2.2f);
    }

    float Environment::biomass(const Vec3 &p) const
    {
        float b = 0.f;
        for (const auto &pulse : m_pulses)
        {
            if (pulse.amount <= 0.f)
                continue;
            const float s2 = pulse.sigma * pulse.sigma;
            b += pulse.amount * std::exp(-len2(p - pulse.pos) / (2.f * s2));
        }
        return clampf(b, 0.f, 4.f);
    }

    static Vec3 gradCentral(const auto &f, const Vec3 &p)
    {
        constexpr float e = 0.8f;
        return Vec3{
                   f(Vec3{p.x + e, p.y, p.z}) - f(Vec3{p.x - e, p.y, p.z}),
                   f(Vec3{p.x, p.y + e, p.z}) - f(Vec3{p.x, p.y - e, p.z}),
                   f(Vec3{p.x, p.y, p.z + e}) - f(Vec3{p.x, p.y, p.z - e}),
               } / (2.f * e);
    }

    Vec3 Environment::gradN(const Vec3 &p) const
    {
        return gradCentral([&](const Vec3 &q)
        {
            return nutrient(q);
        }, p);
    }

    Vec3 Environment::gradL(const Vec3 &p) const
    {
        return gradCentral([&](const Vec3 &q)
        {
            return light(q);
        }, p);
    }

    Vec3 Environment::gradX(const Vec3 &p) const
    {
        return gradCentral([&](const Vec3 &q)
        {
            return toxin(q);
        }, p);
    }

    Vec3 Environment::gradB(const Vec3 &p) const
    {
        return gradCentral([&](const Vec3 &q)
        {
            return biomass(q);
        }, p);
    }

    void Environment::step(const float dt)
    {
        time += dt;
        for (auto &p : m_pulses)
        {
            p.age += dt;
            p.amount *= std::exp(-0.018f * dt);
        }
        std::erase_if(m_pulses,
                      [](const BiomassPulse &p)
                      {
                          return p.age > p.ttl || p.amount < 0.001f;
                      });
    }

    void Environment::depositBiomass(const Vec3 &pos, const float amount, const float sigma)
    {
        if (amount <= 0.f)
            return;
        m_pulses.push_back({pos, amount, sigma, 0.f, 80.f});
    }

    void Environment::consumeBiomassNear(const Vec3 &pos, const float requested)
    {
        if (requested <= 0.f)
            return;
        for (auto &pulse : m_pulses)
        {
            const float influence = std::exp(-len2(pos - pulse.pos) / (2.f * pulse.sigma * pulse.sigma));
            if (influence < 0.02f)
                continue;
            const float take = std::min(pulse.amount, requested * influence * 0.35f);
            pulse.amount -= take;
        }
    }
} // namespace evo