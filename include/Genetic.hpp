#pragma once
#include <string>
#include <vector>
#include <array>
#include <maths/Hash.hpp>

// ======================= Genome -> Dynamic Interpretation ==========================
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
    float aggression       = 0.1f;
    float scavenging       = 0.5f; // digestion affinity (biomass)
    float adhesionBias     = 0.7f;
    float growthBias       = 0.5f;
    float maxAge           = 120.f;
    float baseMetabolism   = 0.04f;
    float splitThreshold   = 18.f;
    float mutationRate     = 0.06f;
    float cellRadius       = 0.65f;

    // Dynamic growth interpreter
    std::array<float, 16> w{};
    std::array<Vec3, 8>   basis{};
    std::array<float, 8>  basisGain{};

    Vec3 pigment{0.7f, 0.7f, 1.0f};
};

static Phenotype interpretGenome(const Genome& g)
{
    Phenotype p;
    if (g.genes.empty()) return p;

    uint64_t master = 0x9E3779B97F4A7C15ULL;
    for (size_t i = 0; i < g.genes.size(); ++i)
    {
        const uint64_t h = fnv1a64(g.genes[i].allele);
        master ^= mix64(h + (uint64_t)i * 0x9E3779B97F4A7C15ULL);
        master = mix64(master);
    }

    auto H = [&](const uint64_t k) { return mix64(master + k * 0x9E3779B97F4A7C15ULL); };

    p.nutrientAffinity = hashToRange(H(1), 0.1f, 2.2f);
    p.photoAffinity    = hashToRange(H(2), 0.0f, 2.0f);
    p.toxinResistance  = hashToRange(H(3), 0.0f, 1.8f);
    p.aggression       = hashToRange(H(4), 0.0f, 1.5f);
    p.scavenging       = hashToRange(H(5), 0.0f, 2.0f);
    p.adhesionBias     = hashToRange(H(6), 0.1f, 1.6f);
    p.growthBias       = hashToRange(H(7), 0.05f, 1.7f);
    p.maxAge           = hashToRange(H(8), 60.f, 240.f);
    p.baseMetabolism   = hashToRange(H(9), 0.015f, 0.12f);
    p.splitThreshold   = hashToRange(H(10), 10.f, 40.f);
    p.mutationRate     = hashToRange(H(11), 0.01f, 0.18f);
    p.cellRadius       = hashToRange(H(12), 0.4f, 1.0f);

    for (int i = 0; i < 16; ++i) p.w[i] = hashToRange(H(100 + i), -2.0f, 2.0f);
    for (int i = 0; i < 8; ++i)
    {
        p.basis[i]     = hashToVec3(H(200 + i));
        p.basisGain[i] = hashToRange(H(300 + i), -1.0f, 1.0f);
    }

    p.pigment = {
        clampf(0.35f + 0.65f * (p.aggression / 1.5f), 0.f, 1.f),
        clampf(0.25f + 0.75f * (p.nutrientAffinity / 2.2f), 0.f, 1.f),
        clampf(0.20f + 0.80f * (p.photoAffinity / 2.0f), 0.f, 1.f)
    };
    return p;
}

// ======================= Mutation ==========================

static std::string mutateAllele(const std::string& s, const float rate)
{
    static const std::string alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
    std::string              out      = s.empty() ? "a" : s;

    if (rndf() < rate && !out.empty())
    {
        const int i = rndi(0, static_cast<int>(out.size()) - 1);
        out[i]      = alphabet[rndi(0, static_cast<int>(alphabet.size()) - 1)];
    }
    if (rndf() < rate * 0.6f && out.size() < 12)
    {
        const int i = rndi(0, static_cast<int>(out.size()));
        out.insert(out.begin() + i, alphabet[rndi(0, static_cast<int>(alphabet.size()) - 1)]);
    }
    if (rndf() < rate * 0.4f && out.size() > 1)
    {
        const int i = rndi(0, static_cast<int>(out.size()) - 1);
        out.erase(out.begin() + i);
    }
    return out;
}

static Genome mutateGenome(const Genome& parent, const float mutationRate)
{
    Genome g = parent;
    if (g.genes.empty()) g.genes.push_back({"a"});

    for (auto& gene : g.genes)
    {
        gene.allele = mutateAllele(gene.allele, mutationRate);
    }

    if (rndf() < mutationRate * 0.35f && g.genes.size() < 16)
    {
        const std::string base = g.genes[rndi(0, static_cast<int>(g.genes.size()) - 1)].allele;
        g.genes.push_back({mutateAllele(base, 0.4f)});
    }
    if (rndf() < mutationRate * 0.2f && g.genes.size() > 1)
    {
        g.genes.erase(g.genes.begin() + rndi(0, static_cast<int>(g.genes.size()) - 1));
    }
    return g;
}
