#pragma once
#include <random>

namespace evo
{
    class Rng
    {
        public:
            explicit Rng(const uint64_t seed = 1234567ULL): m_rng(seed)
            {
            }

            float uniform(const float a = 0.f, const float b = 1.f) const
            {
                std::uniform_real_distribution d(a, b);
                return d(m_rng);
            }

            int uniformInt(const int a, const int b) const
            {
                std::uniform_int_distribution d(a, b);
                return d(m_rng);
            }

        private:
            mutable std::mt19937_64 m_rng;
    };
} // namespace evo