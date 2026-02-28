#pragma once
#include <array>
#include <string>
#include <vector>
#include "Math.hpp"

namespace evo
{
    struct Gene
    {
        std::string allele;
    };

    struct Genome
    {
        std::vector<Gene> genes;
    };

    struct Phenotype
    {
        float nutrientAffinity = 1.0f;
        float photoAffinity    = 1.0f;
        float toxinResistance  = 0.2f;

        float aggression   = 0.1f;
        float scavenging   = 0.5f;
        float adhesionBias = 0.7f;
        float growthBias   = 0.5f;

        float maxAge         = 120.f;
        float baseMetabolism = 0.04f;
        float splitThreshold = 18.f;
        float mutationRate   = 0.06f;
        float cellRadius     = 0.65f;

        std::array<float, 16> w{};
        std::array<Vec3, 8>   basis{};
        std::array<float, 8>  basisGain{};

        Vec3 pigment{0.7f, 0.7f, 1.0f};
    };

    Phenotype interpretGenome(const Genome &g);

    std::string mutateAllele(const std::string &s, float rate, class Rng &rng);
    Genome      mutateGenome(const Genome &parent, float mutationRate, Rng &rng);
} // namespace evo