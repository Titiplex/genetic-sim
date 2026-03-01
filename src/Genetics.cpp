#include "evo/Genetics.hpp"
#include "evo/Hash.hpp"
#include "evo/Random.hpp"

namespace evo
{
    Phenotype interpretGenome(const Genome &g)
    {
        Phenotype p;
        if (g.genes.empty())
            return p;

        uint64_t master = 0x9E3779B97F4A7C15ULL;
        for (size_t i = 0; i < g.genes.size(); ++i)
        {
            const uint64_t h = fnv1a64(g.genes[i].allele);
            master ^= mix64(h + i * 0x9E3779B97F4A7C15ULL);
            master = mix64(master);
        }

        const auto H = [&](const uint64_t k)
        {
            return mix64(master + k * 0x9E3779B97F4A7C15ULL);
        };

        p.gateSharpness = hashToRange(H(15), 0.6f, 3.0f);

        for (int k = 0; k < Phenotype::K; ++k)
        {
            auto HK = [&](const uint64_t kk)
            {
                return mix64(master ^ static_cast<uint64_t>(k) * 0xD1B54A32D192ED03ULL ^ kk);
            };

            CellMode m;

            // Gains raisonnables (évite explosion)
            m.gUptakeN = hashToRange(HK(1), 0.2f, 2.0f);
            m.gPhoto   = hashToRange(HK(2), 0.0f, 2.2f);
            m.gDigestB = hashToRange(HK(3), 0.0f, 2.0f);
            m.gResistX = hashToRange(HK(4), 0.2f, 2.4f);

            m.gAttack   = hashToRange(HK(5), 0.0f, 1.6f);
            m.gAdhesion = hashToRange(HK(6), 0.6f, 1.7f);
            m.gMaint    = hashToRange(HK(7), 0.7f, 1.8f);

            // Gating weights
            for (int i                          = 0; i < CellMode::F; ++i)
                m.gateW[static_cast<size_t>(i)] =
                    hashToRange(HK(100 + static_cast<uint64_t>(i)), -2.0f, 2.0f);

            m.gateBias = hashToRange(HK(200), -1.0f, 1.0f);

            // Tint for debugging/visualization (stable)
            m.tint = {
                hashToRange(HK(300), 0.4f, 1.0f),
                hashToRange(HK(301), 0.4f, 1.0f),
                hashToRange(HK(302), 0.4f, 1.0f),
            };

            p.modes[static_cast<size_t>(k)] = m;
        }

        p.nutrientAffinity = hashToRange(H(1), 0.1f, 2.2f);
        p.photoAffinity    = hashToRange(H(2), 0.0f, 2.0f);
        p.toxinResistance  = hashToRange(H(3), 0.0f, 1.8f);

        p.aggression   = hashToRange(H(4), 0.0f, 1.5f);
        p.scavenging   = hashToRange(H(5), 0.0f, 2.0f);
        p.adhesionBias = hashToRange(H(6), 0.1f, 1.6f);
        p.growthBias   = hashToRange(H(7), 0.05f, 1.7f);

        p.maxAge         = hashToRange(H(8), 60.f, 240.f);
        p.baseMetabolism = hashToRange(H(9), 0.015f, 0.12f);
        p.splitThreshold = hashToRange(H(10), 10.f, 40.f);
        p.mutationRate   = hashToRange(H(11), 0.01f, 0.18f);
        p.cellRadius     = hashToRange(H(12), 0.4f, 1.0f);

        for (int i = 0; i < 16; ++i)
            p.w[i] = hashToRange(H(100 + static_cast<uint64_t>(i)), -2.0f, 2.0f);
        for (int i = 0; i < 8; ++i)
        {
            p.basis[i]     = hashToVec3(H(200 + static_cast<uint64_t>(i)));
            p.basisGain[i] = hashToRange(H(300 + static_cast<uint64_t>(i)), -1.0f, 1.0f);
        }

        p.pigment = {
            clampf(0.35f + 0.65f * (p.aggression / 1.5f), 0.f, 1.f),
            clampf(0.25f + 0.75f * (p.nutrientAffinity / 2.2f), 0.f, 1.f),
            clampf(0.20f + 0.80f * (p.photoAffinity / 2.0f), 0.f, 1.f),
        };
        return p;
    }

    std::string mutateAllele(const std::string &s, const float rate, const Rng &rng)
    {
        static const std::string alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
        std::string              out      = s.empty() ? "a" : s;

        if (rng.uniform() < rate && !out.empty())
        {
            const int i                 = rng.uniformInt(0, static_cast<int>(out.size()) - 1);
            out[static_cast<size_t>(i)] = alphabet[static_cast<size_t>(rng.uniformInt(
                0, static_cast<int>(alphabet.size()) - 1))];
        }
        if (rng.uniform() < rate * 0.6f && out.size() < 12)
        {
            const int i = rng.uniformInt(0, static_cast<int>(out.size()));
            out.insert(out.begin() + i,
                       alphabet[static_cast<size_t>(rng.
                           uniformInt(0, static_cast<int>(alphabet.size()) - 1))]);
        }
        if (rng.uniform() < rate * 0.4f && out.size() > 1)
        {
            const int i = rng.uniformInt(0, static_cast<int>(out.size()) - 1);
            out.erase(out.begin() + i);
        }
        return out;
    }

    Genome mutateGenome(const Genome &parent, const float mutationRate, const Rng &rng)
    {
        Genome g = parent;
        if (g.genes.empty())
            g.genes.push_back({"a"});

        for (auto &[allele] : g.genes)
            allele = mutateAllele(allele, mutationRate, rng);

        if (rng.uniform() < mutationRate * 0.35f && g.genes.size() < 16)
        {
            const std::string base = g.genes[static_cast<size_t>(rng.uniformInt(
                0, static_cast<int>(g.genes.size()) - 1))].allele;
            g.genes.push_back({mutateAllele(base, 0.4f, rng)});
        }
        if (rng.uniform() < mutationRate * 0.2f && g.genes.size() > 1)
        {
            g.genes.erase(g.genes.begin() + rng.uniformInt(0, static_cast<int>(g.genes.size()) - 1));
        }
        return g;
    }
} // namespace evo