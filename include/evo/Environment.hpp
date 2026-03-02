#pragma once
#include <vector>
#include "Math.hpp"

#include <cstdint>

namespace evo
{
    struct BiomassPulse
    {
        Vec3  pos{};
        float amount = 0.f;
        float sigma  = 5.f;
        float age    = 0.f;
        float ttl    = 80.f;
    };

    class Environment
    {
        public:
            float worldRadius = 80.f;
            float time        = 0.f;

            [[nodiscard]] float nutrient(const Vec3 &p) const;
            [[nodiscard]] float light(const Vec3 &p) const;
            [[nodiscard]] float toxin(const Vec3 &p) const;

            [[nodiscard]] float biomass(const Vec3 &p) const;
            [[nodiscard]] Vec3  gradN(const Vec3 &p) const;
            [[nodiscard]] Vec3  gradL(const Vec3 &p) const;
            [[nodiscard]] Vec3  gradX(const Vec3 &p) const;
            [[nodiscard]] Vec3  gradB(const Vec3 &p) const;

            void step(float dt);

            void depositBiomass(const Vec3 &pos, float amount, float sigma = 5.f);
            void consumeBiomassNear(const Vec3 &pos, float requested);

            [[nodiscard]] const std::vector<BiomassPulse> &pulses() const
            {
                return m_pulses;
            }

            void reset(uint64_t seed)
            {
                time = 0.f;
                // reset biomass grid si tu en as une
                (void)seed;
            }

        private:
            std::vector<BiomassPulse> m_pulses;
    };
} // namespace evo