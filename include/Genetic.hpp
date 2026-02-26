#pragma once
#include <array>
#include <string>
#include <vector>
#include <maths/Hash.hpp>

struct Gene
{
    std::string allele;
    float       dominance = 0.5f;
    float       expressionBias = 0.f;
};

struct GenomeModule
{
    std::string       name;
    std::vector<Gene> genes;
};

struct Genome
{
    std::vector<GenomeModule> modules;
};

enum class LifeStage : uint8_t
{
    Juvenile,
    Adult,
    Senescent
};

struct Phenotype
{
    float nutrientAffinity = 1.0f;
    float photoAffinity    = 1.0f;
    float toxinResistance  = 0.2f;
    float aggression       = 0.1f;
    float scavenging       = 0.5f;
    float adhesionBias     = 0.7f;
    float growthBias       = 0.5f;
    float maxAge           = 120.f;
    float baseMetabolism   = 0.04f;
    float splitThreshold   = 18.f;
    float mutationRate     = 0.06f;
    float cellRadius       = 0.65f;
    float bodyDensity      = 1.0f;
    float motorDrive       = 0.6f;
    float fineness         = 0.5f;
    float immunity         = 0.4f;
    float diseaseResistance = 0.4f;
    float learningRate      = 0.2f;
    float reproductiveFlex  = 0.5f;
    float thermalTolerance  = 0.5f;
    float hydricTolerance   = 0.5f;
    float socialAffinity    = 0.5f;
    float microbiomeEfficiency = 0.5f;
    float circadianPeriod   = 1.0f;
    float oxygenTolerance   = 0.5f;
    float osmoregulation    = 0.5f;
    float parentalInvestment = 0.5f;
    float dispersalDrive    = 0.5f;

    std::array<float, 16> w{};
    std::array<Vec3, 8>   basis{};
    std::array<float, 8>  basisGain{};

    std::array<float, 12> motorW{};
    std::array<Vec3, 6>   actuatorAxis{};
    std::array<float, 6>  actuatorBias{};
    std::array<float, 6>  actuatorFreq{};

    std::array<float, 16> grnNodes{};
    std::array<float, 16 * 16> grnWeights{};

    Vec3 pigment{0.7f, 0.7f, 1.0f};
};

static Genome makeBaseGenome()
{
    Genome g;
    g.modules = {
        {"development", {{"devA", 0.5f, 0.0f}, {"devB", 0.5f, 0.0f}}},
        {"metabolism", {{"metA", 0.5f, 0.0f}, {"metB", 0.5f, 0.0f}}},
        {"perception", {{"perA", 0.5f, 0.0f}, {"perB", 0.5f, 0.0f}}},
        {"locomotion", {{"locA", 0.5f, 0.0f}, {"locB", 0.5f, 0.0f}}},
        {"reproduction", {{"repA", 0.5f, 0.0f}, {"repB", 0.5f, 0.0f}}}
    };
    return g;
}

static float expressionValue(const Gene& gene, const uint64_t seed, const float context)
{
    const float base = hashToRange(seed, -1.f, 1.f);
    return std::tanh(base + gene.expressionBias + context * gene.dominance);
}

static Phenotype interpretGenome(const Genome& g, const float temperature = 0.5f, const float resources = 0.5f, const float stress = 0.2f)
{
    Phenotype p;
    if (g.modules.empty()) return p;

    uint64_t master = 0x9E3779B97F4A7C15ULL;
    int gid = 0;
    for (const auto& mod : g.modules)
    {
        master ^= fnv1a64(mod.name);
        for (const auto& gene : mod.genes)
        {
            const uint64_t h = fnv1a64(gene.allele);
            master ^= mix64(h + static_cast<uint64_t>(gid++) * 0x9E3779B97F4A7C15ULL);
            master = mix64(master);
        }
    }

    auto H = [&](const uint64_t k) { return mix64(master + k * 0x9E3779B97F4A7C15ULL); };
    const float context = 0.45f * temperature + 0.35f * resources - 0.6f * stress;

    for (int i = 0; i < 16; ++i)
    {
        p.grnNodes[i] = hashToRange(H(900 + i), -0.7f, 0.7f) + context * 0.2f;
        for (int j = 0; j < 16; ++j) p.grnWeights[i * 16 + j] = hashToRange(H(1000 + i * 17 + j), -1.2f, 1.2f);
    }

    for (int s = 0; s < 3; ++s)
    {
        std::array<float, 16> next{};
        for (int i = 0; i < 16; ++i)
        {
            float acc = p.grnNodes[i];
            for (int j = 0; j < 16; ++j) acc += p.grnNodes[j] * p.grnWeights[i * 16 + j] * 0.12f;
            next[i] = std::tanh(acc + context * 0.35f);
        }
        p.grnNodes = next;
    }

    p.nutrientAffinity = clampf(1.1f + p.grnNodes[0], 0.1f, 2.5f);
    p.photoAffinity    = clampf(1.0f + p.grnNodes[1], 0.f, 2.1f);
    p.toxinResistance  = clampf(0.8f + p.grnNodes[2], 0.f, 2.0f);
    p.aggression       = clampf(0.5f + p.grnNodes[3], 0.f, 1.8f);
    p.scavenging       = clampf(0.8f + p.grnNodes[4], 0.f, 2.1f);
    p.adhesionBias     = clampf(0.8f + p.grnNodes[5], 0.1f, 1.8f);
    p.growthBias       = clampf(0.7f + p.grnNodes[6], 0.05f, 1.8f);
    p.maxAge           = clampf(140.f + p.grnNodes[7] * 80.f, 55.f, 260.f);
    p.baseMetabolism   = clampf(0.055f + p.grnNodes[8] * 0.04f, 0.015f, 0.16f);
    p.splitThreshold   = clampf(20.f + p.grnNodes[9] * 16.f, 8.f, 42.f);
    p.mutationRate     = clampf(0.08f + p.grnNodes[10] * 0.08f, 0.01f, 0.2f);
    p.cellRadius       = clampf(0.7f + p.grnNodes[11] * 0.25f, 0.35f, 1.1f);
    p.bodyDensity      = clampf(1.1f + p.grnNodes[12] * 0.6f, 0.5f, 1.9f);
    p.motorDrive       = clampf(0.7f + p.grnNodes[13] * 0.7f, 0.05f, 1.6f);
    p.fineness         = clampf(0.5f + p.grnNodes[14] * 0.5f, 0.1f, 1.0f);
    p.immunity         = clampf(0.5f + p.grnNodes[15] * 0.5f, 0.05f, 1.6f);
    p.diseaseResistance = clampf(0.7f * p.immunity + 0.4f * p.toxinResistance, 0.1f, 1.8f);
    p.learningRate = clampf(0.2f + 0.35f * (1.f - stress) + 0.2f * resources, 0.05f, 0.95f);
    p.reproductiveFlex = clampf(0.4f + 0.3f * resources + 0.2f * temperature - 0.15f * stress, 0.05f, 1.0f);
    p.thermalTolerance = clampf(0.55f + 0.45f * p.grnNodes[2] - 0.15f * p.grnNodes[8], 0.05f, 1.4f);
    p.hydricTolerance = clampf(0.55f + 0.4f * p.grnNodes[6] + 0.2f * p.grnNodes[12], 0.05f, 1.5f);
    p.socialAffinity = clampf(0.45f + 0.6f * p.grnNodes[4] - 0.3f * p.grnNodes[3], 0.f, 1.8f);
    p.microbiomeEfficiency = clampf(0.55f + 0.45f * p.grnNodes[15] + 0.2f * p.grnNodes[1], 0.05f, 1.8f);
    p.circadianPeriod = clampf(0.7f + 0.45f * (p.grnNodes[10] + 1.f), 0.45f, 1.7f);
    p.oxygenTolerance = clampf(0.75f + 0.5f * p.grnNodes[2] + 0.18f * p.grnNodes[1], 0.05f, 1.8f);
    p.osmoregulation = clampf(0.7f + 0.45f * p.grnNodes[14] + 0.2f * p.grnNodes[6], 0.05f, 1.8f);
    p.parentalInvestment = clampf(0.65f + 0.5f * p.grnNodes[9] - 0.15f * p.grnNodes[3], 0.05f, 1.8f);
    p.dispersalDrive = clampf(0.7f + 0.5f * p.grnNodes[13] + 0.1f * p.grnNodes[10], 0.05f, 1.8f);

    for (int i = 0; i < 16; ++i) p.w[i] = hashToRange(H(100 + i), -2.0f, 2.0f) + p.grnNodes[i] * 0.25f;
    for (int i = 0; i < 8; ++i)
    {
        p.basis[i]     = hashToVec3(H(200 + i));
        p.basisGain[i] = hashToRange(H(300 + i), -1.0f, 1.0f) + p.grnNodes[i] * 0.15f;
    }

    for (int i = 0; i < 12; ++i) p.motorW[i] = hashToRange(H(400 + i), -2.2f, 2.2f) + p.grnNodes[i % 16] * 0.2f;
    for (int i = 0; i < 6; ++i)
    {
        p.actuatorAxis[i] = hashToVec3(H(500 + i));
        p.actuatorBias[i] = hashToRange(H(600 + i), -1.0f, 1.0f);
        p.actuatorFreq[i] = hashToRange(H(700 + i), 0.15f, 2.4f);
    }

    p.pigment = {
        clampf(0.25f + 0.75f * (p.aggression / 1.8f), 0.f, 1.f),
        clampf(0.2f + 0.8f * (p.nutrientAffinity / 2.5f), 0.f, 1.f),
        clampf(0.15f + 0.85f * (p.photoAffinity / 2.1f), 0.f, 1.f)
    };
    return p;
}

static uint64_t genomeSignature(const Genome& g)
{
    uint64_t s = 0x9E3779B97F4A7C15ULL;
    for (const auto& m : g.modules)
    {
        s ^= mix64(fnv1a64(m.name) + 0xA24BAED4963EE407ULL);
        for (const auto& gene : m.genes)
        {
            s ^= mix64(fnv1a64(gene.allele) ^ static_cast<uint64_t>(gene.dominance * 1000.f));
            s = mix64(s);
        }
    }
    return s;
}

static std::string mutateAllele(const std::string& s, const float rate)
{
    static const std::string alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
    std::string out = s.empty() ? "a" : s;
    if (rndf() < rate && !out.empty()) out[rndi(0, static_cast<int>(out.size()) - 1)] = alphabet[rndi(0, static_cast<int>(alphabet.size()) - 1)];
    if (rndf() < rate * 0.6f && out.size() < 12) out.insert(out.begin() + rndi(0, static_cast<int>(out.size())), alphabet[rndi(0, static_cast<int>(alphabet.size()) - 1)]);
    if (rndf() < rate * 0.4f && out.size() > 1) out.erase(out.begin() + rndi(0, static_cast<int>(out.size()) - 1));
    return out;
}

static void mutateModuleStructural(GenomeModule& m, const float mutationRate)
{
    if (m.genes.empty()) m.genes.push_back({"a", 0.5f, 0.f});
    for (auto& gene : m.genes)
    {
        gene.allele = mutateAllele(gene.allele, mutationRate);
        gene.dominance = clampf(gene.dominance + rndf(-mutationRate, mutationRate), 0.05f, 0.95f);
        gene.expressionBias = clampf(gene.expressionBias + rndf(-mutationRate, mutationRate), -1.2f, 1.2f);
    }

    if (rndf() < mutationRate * 0.35f && m.genes.size() < 20)
    {
        const Gene base = m.genes[rndi(0, static_cast<int>(m.genes.size()) - 1)]; // duplication
        m.genes.push_back(base);
    }
    if (rndf() < mutationRate * 0.25f && m.genes.size() > 2)
    {
        m.genes.erase(m.genes.begin() + rndi(0, static_cast<int>(m.genes.size()) - 1)); // deletion
    }
    if (rndf() < mutationRate * 0.20f && m.genes.size() > 3)
    {
        const int a = rndi(0, static_cast<int>(m.genes.size()) - 2);
        const int b = rndi(a + 1, static_cast<int>(m.genes.size()) - 1);
        std::reverse(m.genes.begin() + a, m.genes.begin() + b); // inversion
    }
    if (rndf() < mutationRate * 0.15f && m.genes.size() > 1)
    {
        const int from = rndi(0, static_cast<int>(m.genes.size()) - 1);
        const int to = rndi(0, static_cast<int>(m.genes.size()) - 1);
        Gene ins = m.genes[from];
        ins.allele = mutateAllele(ins.allele, mutationRate * 2.f);
        m.genes.insert(m.genes.begin() + to, ins); // insertion
    }
}

static Genome mutateGenome(const Genome& parent, const float mutationRate)
{
    Genome g = parent;
    if (g.modules.empty()) g = makeBaseGenome();
    for (auto& mod : g.modules) mutateModuleStructural(mod, mutationRate);
    if (rndf() < mutationRate * 0.08f && g.modules.size() > 2)
    {
        std::swap(g.modules[rndi(0, static_cast<int>(g.modules.size()) - 1)], g.modules[rndi(0, static_cast<int>(g.modules.size()) - 1)]);
    }
    return g;
}

static Genome recombineGenomes(const Genome& a, const Genome& b, const float crossingRate)
{
    Genome out;
    const size_t moduleCount = std::min(a.modules.size(), b.modules.size());
    for (size_t m = 0; m < moduleCount; ++m)
    {
        GenomeModule nm;
        nm.name = a.modules[m].name;
        const auto& ga = a.modules[m].genes;
        const auto& gb = b.modules[m].genes;
        const size_t maxGenes = std::max(ga.size(), gb.size());
        if (maxGenes == 0)
        {
            nm.genes.push_back({"fallback", 0.5f, 0.f});
            out.modules.push_back(std::move(nm));
            continue;
        }
        const size_t cut = static_cast<size_t>(rndi(0, static_cast<int>(maxGenes) - 1));
        for (size_t i = 0; i < maxGenes; ++i)
        {
            const bool crossoverFromA = (i <= cut);
            const bool randomMix = rndf() > crossingRate;
            const bool pickA = randomMix ? (rndf() < 0.5f) : crossoverFromA;
            if (pickA)
            {
                if (i < ga.size()) nm.genes.push_back(ga[i]);
                else if (i < gb.size()) nm.genes.push_back(gb[i]);
            }
            else
            {
                if (i < gb.size()) nm.genes.push_back(gb[i]);
                else if (i < ga.size()) nm.genes.push_back(ga[i]);
            }
        }
        if (nm.genes.empty()) nm.genes.push_back({"fallback", 0.5f, 0.f});
        out.modules.push_back(std::move(nm));
    }
    if (out.modules.empty()) out = makeBaseGenome();
    return out;
}
