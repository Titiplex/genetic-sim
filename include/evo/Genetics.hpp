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

    struct CellMode
    {
        // Primitive gains (0..~2). 1.0 = neutre.
        float gUptakeN = 1.0f;
        float gPhoto   = 1.0f;
        float gDigestB = 1.0f;
        float gResistX = 1.0f;

        float gAttack   = 0.0f; // 0..1.5
        float gAdhesion = 1.0f; // multiplie stiffness/rest
        float gMaint    = 1.0f; // multiplie maintenance

        // Gating weights
        static constexpr int F = 8;
        std::array<float, F> gateW{};
        float                gateBias = 0.0f;

        // Rendering hint
        Vec3 tint{1.f, 1.f, 1.f};
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

        static constexpr int    K = 8;
        std::array<CellMode, K> modes{};
        float                   gateSharpness = 1.2f; // température softmax
    };

    Phenotype interpretGenome(const Genome &g);

    std::string mutateAllele(const std::string &s, float rate, const class Rng &rng);
    Genome      mutateGenome(const Genome &parent, float mutationRate, const Rng &rng);
} // namespace evo