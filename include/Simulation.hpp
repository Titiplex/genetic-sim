#pragma once
#include <Cell.hpp>
#include <Environment.hpp>
#include <fstream>
#include <unordered_map>
#if defined(EVO_HAS_OPENMP)
#include <omp.h>
#endif

static uint64_t g_nextId = 1;

struct SimMetrics
{
    float heterozygosity = 0.f;
    float meanFitness = 0.f;
    float fitnessStd = 0.f;
    float trophicStability = 0.f;
    float cooperationIndex = 0.f;
    float diseaseLoad = 0.f;
    float climateStress = 0.f;
    float oxygenStress = 0.f;
    float salinityStress = 0.f;
    float acidStress = 0.f;
    float dormancyRatio = 0.f;
    float reproductionRate = 0.f;
    float meanGenomeComplexity = 0.f;
    float shannonDiversity = 0.f;
    float nicheDiversity = 0.f;
    float meanCellCount = 0.f;
    int extinctions = 0;
    int speciations = 0;
};

struct Sim
{
    Environment env;
    std::vector<Organism> organisms;

    float dt = 0.035f;
    float cellSpacing = 1.25f;
    int targetPopulation = 320;
    float worldRadius = 130.f;

    std::unordered_map<IVec3, std::vector<int>, IVec3Hash> grid;
    float gridSize = 7.f;
    mutable std::vector<int> nearbyScratch;

    SimMetrics metrics;
    int steps = 0;
    int birthsThisStep = 0;
    float cachedSeasonality = 1.f;
    std::unordered_map<uint64_t, int> lineageBirths;

    static Vec3 organismCenterOfMass(const Organism& o)
    {
        if (o.cells.empty()) return o.pos;
        Vec3 s{0, 0, 0};
        for (const auto& c : o.cells) s += o.pos + c.localPos;
        return s / static_cast<float>(o.cells.size());
    }

    static Vec3 organismLocalCOM(const Organism& o)
    {
        if (o.cells.empty()) return Vec3{0, 0, 0};
        Vec3 s{0, 0, 0};
        for (const auto& c : o.cells) s += c.localPos;
        return s / static_cast<float>(o.cells.size());
    }

    static void recenterLocalBody(Organism& o)
    {
        if (o.cells.empty()) return;
        const Vec3 lc = organismLocalCOM(o);
        for (auto& c : o.cells) c.localPos -= lc;
        o.pos += lc;
    }

    IVec3 cellOf(const Vec3& p) const
    {
        return IVec3{static_cast<int>(std::floor(p.x / gridSize)), static_cast<int>(std::floor(p.y / gridSize)), static_cast<int>(std::floor(p.z / gridSize))};
    }

    void buildGrid()
    {
        grid.clear();
        grid.reserve(organisms.size() * 2);
        for (int i = 0; i < static_cast<int>(organisms.size()); ++i)
            if (organisms[i].alive) grid[cellOf(organisms[i].pos)].push_back(i);
    }

    void nearbyOrganisms(const Vec3& p, std::vector<int>& out) const
    {
        out.clear();
        const auto [x, y, z] = cellOf(p);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz)
                {
                    const IVec3 q{x + dx, y + dy, z + dz};
                    if (auto it = grid.find(q); it != grid.end()) out.insert(out.end(), it->second.begin(), it->second.end());
                }
    }

    static float genomeCompatibility(const Organism& a, const Organism& b)
    {
        const uint64_t d = a.geneticSignature ^ b.geneticSignature;
        const float pop = static_cast<float>(__builtin_popcountll(d));
        return clampf(1.f - pop / 64.f, 0.f, 1.f);
    }

    float localDensity(const Vec3& p, const int selfIndex) const
    {
        nearbyOrganisms(p, nearbyScratch);
        const int count = static_cast<int>(std::count_if(nearbyScratch.begin(), nearbyScratch.end(), [&](const int idx) { return idx != selfIndex; }));
        return clampf(static_cast<float>(count) / 28.f, 0.f, 2.5f);
    }

    [[nodiscard]] float reproductionSeasonality() const
    {
        const float annual = 0.35f + 0.65f * (0.5f + 0.5f * std::sin(env.time * 0.08f));
        const float circadian = 0.75f + 0.25f * (0.5f + 0.5f * std::sin(env.time * 0.28f));
        const float lunar = 0.82f + 0.18f * (0.5f + 0.5f * std::sin(env.time * 0.018f + 0.7f));
        return annual * circadian * lunar;
    }

    [[nodiscard]] bool shouldRefreshPhenotype(const Organism& o, const float stress) const
    {
        const float cadence = clampf(0.12f + 0.55f * (1.f - o.pheno.learningRate) + 0.22f * stress, 0.08f, 0.95f);
        return o.phenotypeUpdateTimer >= cadence || stress > 1.5f;
    }

    float ecologicalCompatibility(const Organism& a, const Organism& b) const
    {
        const float marineA = 1.f - a.terrestrialAffinity;
        const float marineB = 1.f - b.terrestrialAffinity;
        const float habitat = 1.f - std::abs(marineA - marineB);
        const float social = 1.f - std::abs(a.socialDrive - b.socialDrive);
        return clampf(0.6f * habitat + 0.4f * social, 0.f, 1.f);
    }

    Vec3 schoolingForce(const Organism& o, const int selfIndex) const
    {
        nearbyOrganisms(o.pos, nearbyScratch);
        Vec3 align{0,0,0}, cohesion{0,0,0}, separation{0,0,0};
        int neighbors = 0;
        for (const int idx : nearbyScratch)
        {
            if (idx == selfIndex) continue;
            const Organism& n = organisms[static_cast<size_t>(idx)];
            if (!n.alive) continue;
            const Vec3 d = n.pos - o.pos;
            const float d2 = len2(d);
            if (d2 < 0.25f || d2 > 12.f * 12.f) continue;
            neighbors++;
            align += n.vel;
            cohesion += n.pos;
            separation -= d / std::max(0.5f, std::sqrt(d2));
        }
        if (neighbors == 0) return {0,0,0};
        const float inv = 1.f / static_cast<float>(neighbors);
        const Vec3 desiredAlign = normalize(align * inv) * 0.35f;
        const Vec3 desiredCohesion = normalize((cohesion * inv) - o.pos) * 0.28f;
        const Vec3 desiredSeparation = normalize(separation * inv) * 0.42f;
        const float socialWeight = clampf(o.pheno.socialAffinity * o.socialDrive, 0.f, 1.4f);
        return (desiredAlign + desiredCohesion + desiredSeparation) * socialWeight;
    }

    void seedInitial(const int n)
    {
        organisms.clear();
        organisms.reserve(n * 2);
        for (int i = 0; i < n; i++)
        {
            Organism o;
            o.id = g_nextId++;
            o.genome = makeBaseGenome();
            if (i > 0)
            {
                o.genome = mutateGenome(o.genome, 0.45f);
                if (rndf() < 0.5f) o.genome = mutateGenome(o.genome, 0.2f);
            }
            o.pheno = interpretGenome(o.genome, 0.5f, 0.5f, 0.2f);
            o.geneticSignature = genomeSignature(o.genome);
            o.matePreference = rndf(0.f, 1.f);
            o.socialDrive = rndf(0.f, 1.f);
            o.pos = {rndf(-55.f, 55.f), rndf(-30.f, 45.f), rndf(-55.f, 55.f)};
            o.vel = {rndf(-0.2f, 0.2f), rndf(-0.2f, 0.2f), rndf(-0.2f, 0.2f)};
            o.energy = rndf(6.f, 12.f);
            o.storedEnergy = rndf(1.f, 3.f);
            o.deme = rndi(0, 4);
            o.terrestrialAffinity = clampf(0.2f + 0.6f * o.pheno.hydricTolerance + rndf(-0.15f, 0.15f), 0.f, 1.f);
            o.cells.push_back({Vec3{0, 0, 0}, Vec3{0, 0, 0}, 1.f, 0.f, 1.f});
            organisms.push_back(std::move(o));
        }
    }

    Vec3 chooseGrowthDirection(const Organism& o) const
    {
        const Vec3 p = o.pos;
        const Vec3 gN = env.gradN(p), gL = env.gradL(p), gX = env.gradX(p), gB = env.gradB(p);
        const Vec3 localCom = organismLocalCOM(o);
        Vec3 principal{0, 0, 0}, shellOut{0, 0, 0};
        for (const auto& c : o.cells)
        {
            const Vec3 d = c.localPos - localCom;
            const float l = len(d);
            if (l > 1e-4f)
            {
                principal += normalize(d) * (0.1f + l);
                if (l > 0.75f) shellOut += normalize(d);
            }
        }
        principal = normalize(principal);
        shellOut = normalize(shellOut);
        Vec3 repulse{0, 0, 0};
        for (const auto& c : o.cells)
        {
            const Vec3 d = localCom - c.localPos;
            const float l = len(d);
            if (l > 1e-4f) repulse += d / (l * l + 0.2f);
        }
        repulse = normalize(repulse);
        const Vec3 randomDrift = normalize(Vec3{rndf(-1, 1), rndf(-1, 1), rndf(-1, 1)});
        const std::array<Vec3, 16> feat = {normalize(gN), normalize(gL), normalize(gX) * -1.f, normalize(gB), normalize(o.vel), principal, repulse, shellOut, randomDrift, normalize(p * -1.f), normalize(gN + gL * 0.5f - gX * 0.7f), normalize(gN + gB * 0.4f), normalize(gL - gX), Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}};
        Vec3 dir{0, 0, 0};
        for (int i = 0; i < 16; i++) dir += feat[i] * o.pheno.w[i];
        for (int i = 0; i < 8; i++) dir += o.pheno.basis[i] * o.pheno.basisGain[i] * (0.15f + 0.85f * rndf());
        if (len2(dir) < 1e-6f) dir = randomDrift;
        return normalize(dir);
    }

    void applySprings(Organism& o) const
    {
        for (const auto& [a, b, restLen, stiffness, damping] : o.springs)
        {
            if (a < 0 || b < 0 || a >= static_cast<int>(o.cells.size()) || b >= static_cast<int>(o.cells.size())) continue;
            Cell& A = o.cells[a];
            Cell& B = o.cells[b];
            const Vec3 d = B.localPos - A.localPos;
            const float l = len(d);
            if (l < 1e-6f) continue;
            const Vec3 n = d / l;
            const float relv = dot(B.localVel - A.localVel, n);
            const Vec3 force = n * (stiffness * (l - restLen) + damping * relv);
            A.localVel += force * (dt / std::max(0.2f, A.mass));
            B.localVel -= force * (dt / std::max(0.2f, B.mass));
        }
    }

    void applySoftCellRepulsion(Organism& o) const
    {
        const float minDist = std::max(0.55f, 1.15f * o.pheno.cellRadius);
        for (int i = 0; i < static_cast<int>(o.cells.size()); ++i)
            for (int j = i + 1; j < static_cast<int>(o.cells.size()); ++j)
            {
                Vec3 d = o.cells[j].localPos - o.cells[i].localPos;
                const float l = len(d);
                if (l < 1e-5f)
                {
                    Vec3 push = normalize(Vec3{rndf(-1, 1), rndf(-1, 1), rndf(-1, 1)});
                    if (len2(push) < 1e-6f) push = Vec3{1, 0, 0};
                    push *= 0.02f;
                    o.cells[i].localVel -= push;
                    o.cells[j].localVel += push;
                    continue;
                }
                if (l < minDist)
                {
                    const Vec3 f = (d / l) * (8.0f * (minDist - l));
                    o.cells[i].localVel -= f * dt;
                    o.cells[j].localVel += f * dt;
                }
            }
    }

    void integrateInternalBody(Organism& o) const
    {
        if (o.cells.empty()) return;
        if (o.macroMode)
        {
            const int stride = std::max(1, static_cast<int>(o.cells.size() / 24));
            for (int i = 0; i < static_cast<int>(o.cells.size()); i += stride)
            {
                o.cells[i].localVel *= 0.9f;
                o.cells[i].localPos += o.cells[i].localVel * dt;
                o.cells[i].age += dt;
            }
            recenterLocalBody(o);
            return;
        }
        applySprings(o);
        applySoftCellRepulsion(o);
        for (auto& c : o.cells)
        {
            c.localVel *= 0.95f;
            c.localPos += c.localVel * dt;
            c.age += dt;
        }
        recenterLocalBody(o);
    }

    static void connectNewCell(Organism& o, const int newCellIndex)
    {
        if (newCellIndex <= 0) return;
        std::vector<std::pair<float, int>> ds;
        ds.reserve(newCellIndex);
        for (int i = 0; i < newCellIndex; ++i) ds.emplace_back(len2(o.cells[i].localPos - o.cells[newCellIndex].localPos), i);
        int linkCount = std::min(1 + (o.pheno.adhesionBias > 0.7f) + (o.pheno.adhesionBias > 1.15f), static_cast<int>(ds.size()));
        if (linkCount < static_cast<int>(ds.size())) std::nth_element(ds.begin(), ds.begin() + linkCount, ds.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        for (int k = 0; k < linkCount; ++k)
        {
            Spring s;
            s.a = newCellIndex;
            s.b = ds[k].second;
            s.restLen = clampf(std::sqrt(ds[k].first), 0.45f, 2.6f);
            s.stiffness = 5.f + 6.f * o.pheno.adhesionBias;
            o.springs.push_back(s);
        }
    }

    static void pruneSprings(Organism& o)
    {
        o.springs.erase(std::remove_if(o.springs.begin(), o.springs.end(), [&](const Spring& s)
        {
            return s.a < 0 || s.b < 0 || s.a >= static_cast<int>(o.cells.size()) || s.b >= static_cast<int>(o.cells.size()) || s.a == s.b;
        }), o.springs.end());
    }

    void growOrganism(Organism& o) const
    {
        if (o.cells.size() >= 220 || o.energy < 2.f) return;
        const float macroPenalty = o.macroMode ? 0.45f : 1.f;
        const float desire = macroPenalty * o.pheno.growthBias * (0.35f + 0.65f * clampf(o.energy / o.pheno.splitThreshold, 0.f, 2.f));
        if (rndf() > clampf(desire * dt * 1.7f, 0.f, 0.30f)) return;
        const Vec3 dir = chooseGrowthDirection(o);
        int anchorIdx = 0;
        float bestL2 = -1.f;
        for (int i = 0; i < static_cast<int>(o.cells.size()); ++i)
            if (const float d2 = len2(o.cells[i].localPos); d2 > bestL2) { bestL2 = d2; anchorIdx = i; }
        if (rndf() < 0.55f) anchorIdx = rndi(0, static_cast<int>(o.cells.size()) - 1);
        const Vec3 anchor = o.cells[anchorIdx].localPos;
        Vec3 bestPos = anchor + dir * (cellSpacing * o.pheno.cellRadius);
        float bestScore = -1e9f;
        for (int k = 0; k < 10; k++)
        {
            const Vec3 candDir = normalize(dir + hashToVec3(static_cast<uint64_t>(k) + o.id + static_cast<uint64_t>(o.age * 37.f)) * rndf(0.0f, 0.9f));
            const Vec3 cand = anchor + candDir * (cellSpacing * o.pheno.cellRadius * rndf(0.9f, 1.45f));
            float minD = 1e9f, rep = 0.f;
            for (const auto& c : o.cells)
            {
                const float d = len(c.localPos - cand);
                minD = std::min(minD, d);
                if (d < 1.7f * o.pheno.cellRadius) rep += (1.7f * o.pheno.cellRadius - d);
            }
            const float score = minD - rep * 0.9f + rndf(-0.15f, 0.15f);
            if (score > bestScore) { bestScore = score; bestPos = cand; }
        }
        if (bestScore > 0.12f)
        {
            o.cells.push_back({bestPos, dir * rndf(0.03f, 0.15f), 1.f, 0.f, 1.f});
            connectNewCell(o, static_cast<int>(o.cells.size()) - 1);
            o.energy -= 0.55f + 0.035f * static_cast<float>(o.cells.size());
        }
    }

    void huntInteractions()
    {
        for (int i = 0; i < static_cast<int>(organisms.size()); ++i)
        {
            Organism& a = organisms[i];
            if (!a.alive) continue;
            nearbyOrganisms(a.pos, nearbyScratch);
            for (int j : nearbyScratch)
            {
                if (j <= i) continue;
                Organism& b = organisms[j];
                if (!b.alive) continue;
                const float d2 = len2(b.pos - a.pos);
                if (d2 > 16.f) continue;

                const float contagion = (a.pathogenLoad + b.pathogenLoad) * 0.0016f * dt;
                if (contagion > 0.f)
                {
                    a.pathogenLoad += contagion * (1.f - 0.45f * a.pheno.diseaseResistance);
                    b.pathogenLoad += contagion * (1.f - 0.45f * b.pheno.diseaseResistance);
                }

                const bool kinLike = (a.deme == b.deme) || genomeCompatibility(a, b) > 0.75f;
                if (kinLike && a.energy > 3.5f && b.energy < 2.5f)
                {
                    const float share = 0.02f * a.pheno.socialAffinity * dt * 20.f;
                    const float transfer = std::min(share, a.energy - 3.2f);
                    if (transfer > 0.f) { a.energy -= transfer; b.energy += transfer * 0.9f; }
                }

                if (kinLike && b.energy > 3.5f && a.energy < 2.5f)
                {
                    const float share = 0.02f * b.pheno.socialAffinity * dt * 20.f;
                    const float transfer = std::min(share, b.energy - 3.2f);
                    if (transfer > 0.f) { b.energy -= transfer; a.energy += transfer * 0.9f; }
                }

                auto attack = [&](Organism& attacker, Organism& defender)
                {
                    if (attacker.pheno.aggression < 0.30f || attacker.energy < 2.f) return;
                    const float huntRange = 2.5f + 0.08f * static_cast<float>(attacker.cells.size());
                    if (d2 > huntRange * huntRange) return;
                    const float aPower = attacker.pheno.aggression + 0.018f * static_cast<float>(attacker.cells.size());
                    const float bDef = 0.5f * defender.pheno.toxinResistance + 0.010f * static_cast<float>(defender.cells.size());
                    const float damage = clampf((aPower - 0.4f * bDef + 0.2f * env.predationPressure(defender.pos)) * dt * 4.2f, 0.f, 0.7f);
                    if (damage <= 0.f) return;
                    defender.energy -= damage;
                    defender.pathogenLoad += damage * 0.05f;
                    attacker.energy += damage * (0.45f + 0.15f * attacker.pheno.microbiomeEfficiency);
                };

                attack(a, b);
                attack(b, a);
            }
        }
    }

    void digestBiomass(Organism& o)
    {
        const Vec3 com = organismCenterOfMass(o);
        const float B = env.biomass(com);
        if (B <= 0.0001f) return;
        const float request = o.pheno.scavenging * (0.02f + 0.005f * static_cast<float>(o.cells.size())) * dt * 20.f;
        const float requested = std::min(B, request);
        const float consumed = env.consumeBiomassNear(com, requested);
        const float gain = consumed * 0.9f;
        const float microbiomeGain = 0.75f + 0.5f * o.pheno.microbiomeEfficiency * o.microbiomeHealth;
        o.energy += gain * microbiomeGain;
        o.microbiomeHealth = clampf(o.microbiomeHealth + gain * 0.002f - 0.001f * o.pathogenLoad, 0.05f, 1.4f);
    }

    Vec3 computeActuation(Organism& o, const Vec3& com, const Vec3& gN, const Vec3& gL, const Vec3& gX, const Vec3& gB) const
    {
        const float waterDepth = env.waterSurfaceHeight(com) - com.y;
        const float submerged = clampf(waterDepth / 12.f, 0.f, 1.f);
        const float onGround = clampf((env.groundHeight(com) + o.pheno.cellRadius - com.y) / 2.5f, 0.f, 1.f);
        const float noise = rndf(-0.1f, 0.1f);
        const std::array<float, 12> feat = {clampf(len(o.vel) / 4.5f, 0.f, 1.f), clampf(o.energy / std::max(2.f, o.pheno.splitThreshold), 0.f, 2.f), clampf(o.age / std::max(20.f, o.pheno.maxAge), 0.f, 2.f), clampf(len(gN) * 2.2f + noise, 0.f, 1.f), clampf(len(gL) * 2.2f + noise, 0.f, 1.f), clampf(len(gX) * 2.2f + noise, 0.f, 1.f), clampf(len(gB) * 2.2f + noise, 0.f, 1.f), submerged, onGround, clampf(static_cast<float>(o.cells.size()) / 28.f, 0.f, 2.f), clampf(env.shorelineBlend(com), 0.f, 1.f), 1.f};
        float command = 0.f;
        for (int i = 0; i < 12; ++i) command += feat[i] * o.pheno.motorW[i];
        command = std::tanh(command + o.neuralMemory * 0.4f);
        o.neuralMemory = clampf(o.neuralMemory * (1.f - 0.15f * o.pheno.learningRate) + command * 0.08f + o.lastReward * 0.03f, -1.f, 1.f);
        o.gaitPhase += dt * (0.25f + 0.75f * o.pheno.fineness);
        Vec3 thrust{0, 0, 0};
        for (int i = 0; i < 6; ++i)
        {
            const Vec3 axis = normalize(o.pheno.actuatorAxis[i]);
            const float osc = std::sin(o.gaitPhase * (0.8f + o.pheno.actuatorFreq[i]) + o.pheno.actuatorBias[i] * 3.14159f);
            const float amp = (0.35f + 0.65f * submerged + 0.55f * onGround) * o.pheno.motorDrive;
            thrust += axis * (command * osc * amp);
        }
        return thrust;
    }

    void updateLifeStage(Organism& o)
    {
        const float relAge = o.age / std::max(1.f, o.pheno.maxAge);
        if (relAge < 0.23f) o.stage = LifeStage::Juvenile;
        else if (relAge < 0.78f) o.stage = LifeStage::Adult;
        else o.stage = LifeStage::Senescent;
    }

    void applyDiseaseAndImmunity(Organism& o, const float stress, const float localPredation)
    {
        const float chronicStress = clampf(stress + 0.45f * o.epigeneticStress, 0.f, 3.f);
        o.pathogenLoad += (0.0025f + 0.006f * chronicStress + 0.003f * localPredation) * dt;
        const float immuneResponse = (o.pheno.immunity + 0.5f * o.immuneMemory + 0.22f * o.microbiomeHealth) * dt * 0.04f;
        o.pathogenLoad = std::max(0.f, o.pathogenLoad - immuneResponse);
        o.immuneMemory = clampf(o.immuneMemory * 0.995f + o.pathogenLoad * 0.0009f, 0.f, 1.5f);
        o.epigeneticStress = clampf(o.epigeneticStress * 0.996f + chronicStress * 0.002f - 0.0015f * o.pheno.learningRate, 0.f, 2.5f);
    }

    bool attemptReproduction(Organism& o, std::vector<Organism>& newborns)
    {
        const float season = cachedSeasonality;
        const float popPressure = clampf(static_cast<float>(organisms.size()) / std::max(1.f, static_cast<float>(targetPopulation)), 0.2f, 2.6f);
        const float localMates = localDensity(o.pos, -1);
        const float allee = clampf(localMates / 0.35f, 0.35f, 1.2f);
        const float reproductionStress = clampf(1.f - 0.3f * o.epigeneticStress - 0.2f * o.pathogenLoad, 0.2f, 1.1f);
        if (!(o.energy > o.pheno.splitThreshold && o.age > 6.f && o.cells.size() >= 2 && o.reproductionCooldown <= 0.f && o.dormancy < 0.65f && rndf() < (season * allee * reproductionStress) / popPressure)) return false;
        Organism child;
        child.id = g_nextId++;
        const bool sexual = rndf() < clampf(0.25f + 0.55f * o.pheno.reproductiveFlex - 0.25f * o.energyDebt, 0.05f, 0.85f);
        if (sexual)
        {
            int mateIdx = -1;
            nearbyOrganisms(o.pos, nearbyScratch);
            float bestScore = -1.f;
            for (int idx : nearbyScratch)
            {
                if (organisms[idx].id == o.id || !organisms[idx].alive || organisms[idx].energy <= 5.f) continue;
                const float mateDist2 = len2(organisms[idx].pos - o.pos);
                if (mateDist2 > 12.f * 12.f) continue;
                const float comp = genomeCompatibility(o, organisms[idx]);
                const float ecoComp = ecologicalCompatibility(o, organisms[idx]);
                const float prefAlign = 1.f - std::abs(o.matePreference - comp);
                const float inbreedingPenalty = (comp > 0.92f) ? (comp - 0.92f) * 3.f : 0.f;
                const float outbreedingPenalty = (comp < 0.12f) ? (0.12f - comp) * 2.5f : 0.f;
                const float score = 0.42f * comp + 0.3f * ecoComp + 0.28f * prefAlign - inbreedingPenalty - outbreedingPenalty;
                if (score > bestScore) { bestScore = score; mateIdx = idx; }
            }
            if (mateIdx >= 0)
            {
                Organism& mate = organisms[mateIdx];
                child.genome = recombineGenomes(o.genome, mate.genome, 0.45f);
                child.genome = mutateGenome(child.genome, 0.5f * (o.pheno.mutationRate + mate.pheno.mutationRate));
                mate.energy *= 0.93f;
                child.deme = (rndf() < 0.5f) ? o.deme : mate.deme;
            }
            else
            {
                child.genome = mutateGenome(o.genome, o.pheno.mutationRate);
                child.deme = o.deme;
            }
        }
        else
        {
            child.genome = mutateGenome(o.genome, o.pheno.mutationRate);
            child.deme = o.deme;
        }
        if (env.climateShock > 0.7f && rndf() < 0.22f) child.genome = mutateGenome(child.genome, 0.12f + 0.2f * env.climateShock);
        child.pheno = interpretGenome(child.genome, env.temperature(o.pos), env.nutrient(o.pos), env.toxin(o.pos));
        child.geneticSignature = genomeSignature(child.genome);
        child.matePreference = clampf((o.matePreference + rndf(-0.08f, 0.08f)), 0.f, 1.f);
        child.socialDrive = clampf((o.socialDrive + rndf(-0.08f, 0.08f)), 0.f, 1.f);
        child.microbiomeHealth = clampf(o.microbiomeHealth * rndf(0.7f, 1.05f), 0.05f, 1.2f);
        child.epigeneticStress = clampf(o.epigeneticStress * rndf(0.45f, 0.85f), 0.f, 2.f);
        child.pathogenLoad = clampf(o.pathogenLoad * rndf(0.03f, 0.12f), 0.f, 1.5f);
        child.terrestrialAffinity = clampf(o.terrestrialAffinity + rndf(-0.08f, 0.08f), 0.f, 1.f);
        child.pos = o.pos + normalize(Vec3{rndf(-1, 1), rndf(-1, 1), rndf(-1, 1)}) * rndf(2.f, 6.f);
        child.vel = o.vel + Vec3{rndf(-0.4f, 0.4f), rndf(-0.4f, 0.4f), rndf(-0.4f, 0.4f)};
        const float investment = clampf(0.24f + 0.22f * o.pheno.parentalInvestment, 0.16f, 0.62f);
        child.energy = o.energy * investment + o.offspringInvestmentReserve;
        o.energy *= (1.f - investment * 0.9f);
        o.offspringInvestmentReserve = 0.f;
        const int take = std::min<int>(std::max<int>(1, static_cast<int>(o.cells.size()) / 3), 12);
        std::vector<std::pair<float, int>> order;
        order.reserve(o.cells.size());
        for (int i = 0; i < static_cast<int>(o.cells.size()); ++i) order.emplace_back(len2(o.cells[i].localPos), i);
        if (take < static_cast<int>(order.size())) std::nth_element(order.begin(), order.begin() + take, order.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
        std::vector<int> idxs;
        idxs.reserve(take);
        for (int k = 0; k < take && k < static_cast<int>(order.size()); ++k) idxs.push_back(order[k].second);
        std::sort(idxs.begin(), idxs.end(), std::greater<int>());
        for (int idx : idxs)
        {
            Cell c = o.cells[idx]; c.localPos *= 0.6f; c.localVel = Vec3{0, 0, 0}; c.age = 0.f;
            child.cells.push_back(c); o.cells.erase(o.cells.begin() + idx);
        }
        pruneSprings(o);
        if (child.cells.empty()) child.cells.push_back({Vec3{0, 0, 0}, Vec3{0, 0, 0}, 1.f, 0.f, 1.f});
        Vec3 cc{0, 0, 0}; for (const auto& c : child.cells) cc += c.localPos; cc = cc / static_cast<float>(child.cells.size()); for (auto& c : child.cells) c.localPos -= cc;
        child.springs.clear();
        for (int i = 1; i < static_cast<int>(child.cells.size()); ++i)
        {
            Spring s; s.a = i; s.b = rndi(0, i - 1); s.restLen = clampf(len(child.cells[i].localPos - child.cells[s.b].localPos), 0.45f, 2.6f); s.stiffness = 5.f + 6.f * child.pheno.adhesionBias; child.springs.push_back(s);
        }
        const float newbornViability = clampf(0.85f + 0.25f * child.pheno.reproductiveFlex - 0.18f * child.epigeneticStress - 0.2f * child.pathogenLoad, 0.05f, 1.1f);
        if (rndf() > newbornViability)
        {
            env.depositBiomass(child.pos, 0.6f + 0.05f * static_cast<float>(child.cells.size()), 2.4f);
            o.reproductionCooldown = clampf(5.f + 2.f * (1.f - o.pheno.reproductiveFlex), 2.f, 14.f);
            return false;
        }
        o.reproductionCooldown = clampf(8.f + 5.f * (1.f - o.pheno.reproductiveFlex) + 2.5f * o.pheno.parentalInvestment, 3.f, 24.f);
        lineageBirths[child.geneticSignature]++;
        newborns.push_back(std::move(child));
        birthsThisStep++;
        return true;
    }

    void updateMacroState(Organism& o)
    {
        const size_t c = o.cells.size();
        if (!o.macroMode && c >= 48) o.macroMode = true;
        if (o.macroMode && c <= 30) o.macroMode = false;
        o.macroScale = clampf(0.85f + std::pow(static_cast<float>(std::max<size_t>(1, c)) / 24.f, 0.55f), 0.85f, 4.5f);
    }

    void updateOne(const int selfIndex, Organism& o, std::vector<Organism>& newborns)
    {
        if (!o.alive) return;
        o.age += dt;
        o.phenotypeUpdateTimer += dt;
        updateLifeStage(o);
        integrateInternalBody(o);
        Vec3 com = organismCenterOfMass(o);
        updateMacroState(o);
        const float bodyScale = clampf(static_cast<float>(o.cells.size()) / 18.f, 0.35f, 4.8f);
        const float N = env.nutrient(com), L = env.light(com), X = env.toxin(com);
        const float O2 = env.oxygen(com), S = env.salinity(com);
        const float localPredation = env.predationPressure(com);
        const float co2 = env.dissolvedCO2(com);
        const float vent = env.hydrothermalVent(com);
        const float temp = env.temperature(com);
        const float humid = env.humidity(com);
        float thermalMismatch = std::abs(temp - 0.55f) * (1.5f - o.pheno.thermalTolerance);
        float hydricMismatch = std::max(0.f, 0.5f - humid) * (1.5f - o.pheno.hydricTolerance);
        const float hypoxia = std::max(0.f, 0.7f - O2) * (1.4f - o.pheno.oxygenTolerance);
        const float osmotic = std::abs(S - 0.75f) * (1.4f - o.pheno.osmoregulation);
        const float acidStress = std::max(0.f, co2 - (0.55f * o.pheno.toxinResistance + 0.45f * o.pheno.immunity));
        float stress = clampf(0.5f * X + 0.35f * hydricMismatch + 0.25f * thermalMismatch + 0.28f * hypoxia + 0.2f * osmotic + 0.18f * acidStress + 0.2f * localPredation + 0.25f * env.climateShock + 0.12f * vent, 0.f, 2.8f);
        if (shouldRefreshPhenotype(o, stress))
        {
            o.pheno = interpretGenome(o.genome, temp, N, stress);
            o.phenotypeUpdateTimer = 0.f;
        }
        thermalMismatch = std::abs(temp - 0.55f) * (1.5f - o.pheno.thermalTolerance);
        hydricMismatch = std::max(0.f, 0.5f - humid) * (1.5f - o.pheno.hydricTolerance);
        stress = clampf(0.5f * X + 0.35f * hydricMismatch + 0.25f * thermalMismatch + 0.28f * hypoxia + 0.2f * osmotic + 0.18f * acidStress + 0.2f * localPredation + 0.25f * env.climateShock + 0.12f * vent, 0.f, 2.8f);
        if (o.stage == LifeStage::Juvenile)
        {
            o.pheno.motorDrive *= 0.98f;
            o.pheno.growthBias *= 1.01f;
        }
        else if (o.stage == LifeStage::Senescent)
        {
            o.pheno.baseMetabolism *= 1.01f;
            o.pheno.growthBias *= 0.98f;
        }

        const Vec3 gN = env.gradN(com), gL = env.gradL(com), gX = env.gradX(com), gB = env.gradB(com), gO2 = env.gradO2(com), gS = env.signalGradient(com, o.gaitPhase);
        const float circadian = 0.65f + 0.35f * std::sin(env.time * 0.3f / std::max(0.4f, o.pheno.circadianPeriod));
        Vec3 steer = normalize(gN) * (0.31f * o.pheno.nutrientAffinity) + normalize(gL) * (0.16f * o.pheno.photoAffinity) + normalize(gB) * (0.21f * o.pheno.scavenging) + normalize(gO2) * (0.18f * o.pheno.oxygenTolerance) + normalize(gS) * 0.1f - normalize(gX) * (0.48f * (1.2f - clampf(o.pheno.toxinResistance, 0.f, 1.2f))) + hashToVec3(o.id + static_cast<uint64_t>(o.age * 10.0f)) * 0.12f;
        if (len2(steer) < 1e-6f) steer = Vec3{rndf(-1, 1), rndf(-1, 1), rndf(-1, 1)};
        steer = normalize(steer);

        const Vec3 flow = env.fluidVelocity(com), gust = env.turbulence(com, o.id + static_cast<uint64_t>(o.age * 100.f));
        const Vec3 wind = env.windVelocity(com);
        const float buoyancy = env.buoyancy(com, bodyScale);
        const Vec3 actuation = computeActuation(o, com, gN, gL, gX, gB);
        const float submersion = clampf((env.waterSurfaceHeight(com) - com.y) / 9.f, 0.f, 1.f);
        const float landFactor = 1.f - submersion;
        const Biome biome = env.biomeAt(com);
        const float marineBias = env.biomeMarineBias(biome);
        const float habitatAffinity = clampf((1.f - std::abs(o.terrestrialAffinity - (1.f - marineBias))), 0.25f, 1.25f);
        const Vec3 obstaclePush = normalize(Vec3{-env.obstacleField(com) * com.x, 0.f, -env.obstacleField(com) * com.z});
        const float corridorBoost = env.corridorField(com);
        const float dispersal = clampf(o.pheno.dispersalDrive * (0.35f + localDensity(com, selfIndex)), 0.f, 2.2f);
        const Vec3 school = schoolingForce(o, selfIndex);

        o.vel += steer * (0.42f * dt * circadian);
        o.vel += school * (0.22f * dt);
        o.vel += flow * (0.34f * dt * submersion);
        o.vel += wind * (0.26f * dt * landFactor);
        o.vel += gust * dt;
        o.vel += actuation * (0.95f * dt * (0.7f + 0.3f * habitatAffinity));
        o.vel += obstaclePush * (0.12f * dt);
        o.vel += normalize(com * -1.f) * (0.09f * dt * dispersal);
        o.vel *= (0.96f + 0.05f * corridorBoost);
        o.vel.y += buoyancy * dt * (1.0f + 0.2f * (1.3f - o.pheno.bodyDensity));
        if (landFactor > 0.45f) o.vel.y -= 0.13f * dt * (0.4f + o.pheno.bodyDensity);
        const Vec3 relFluid = o.vel - flow;
        o.vel -= relFluid * (0.10f + 0.06f * bodyScale + 0.03f * o.pheno.bodyDensity) * dt;

        const float macroDrag = o.macroMode ? 0.76f : 1.f;
        o.vel = clampVec(o.vel, macroDrag * 5.6f / (1.0f + 0.024f * static_cast<float>(o.cells.size())));
        o.pos += o.vel * dt * 6.0f;

        const float floorY = env.groundHeight(o.pos) + o.pheno.cellRadius * 0.9f;
        if (o.pos.y < floorY) { o.pos.y = floorY; if (o.vel.y < 0.f) o.vel.y *= -0.18f; o.vel.x *= 0.92f; o.vel.z *= 0.92f; }
        if (const float r = len(o.pos); r > worldRadius) { const Vec3 inward = normalize(o.pos) * -1.f; o.vel += inward * dt * 4.f; o.pos = normalize(o.pos) * worldRadius; o.vel *= 0.8f; }

        applyDiseaseAndImmunity(o, stress, localPredation);

        const float crowdPenalty = clampf(localDensity(com, selfIndex) * (1.f - o.pheno.socialAffinity * 0.35f), 0.f, 1.6f);
        const float fertility = env.biomeFertility(biome);
        const float marineFit = clampf(1.f - std::abs((1.f - o.terrestrialAffinity) - marineBias), 0.35f, 1.35f);
        const float landFit = clampf(1.f - std::abs(o.terrestrialAffinity - (1.f - marineBias)), 0.35f, 1.35f);
        const float chemo = vent * (0.012f + 0.004f * static_cast<float>(o.cells.size())) * (0.35f + 0.65f * submersion);
        const float uptake = o.pheno.nutrientAffinity * (N + chemo) * fertility * (0.05f + 0.010f * static_cast<float>(o.cells.size())) * (1.f - 0.18f * crowdPenalty) * (submersion * marineFit + landFactor * landFit);
        const float photo = o.pheno.photoAffinity * L * (0.03f + 0.008f * static_cast<float>(o.cells.size())) * (0.75f + 0.25f * landFactor);
        const float toxHit = std::max(0.f, X - 0.7f * o.pheno.toxinResistance) * (0.04f + 0.005f * static_cast<float>(o.cells.size()));
        const float moveCost = 0.015f * len(o.vel) * (1.0f + 0.03f * static_cast<float>(o.cells.size())) * (1.1f - 0.2f * habitatAffinity);
        const float actuationCost = 0.010f * len(actuation) * (0.6f + 0.8f * o.pheno.motorDrive);
        const float maint = o.pheno.baseMetabolism * (1.0f + 0.055f * static_cast<float>(o.cells.size())) * (o.macroMode ? 0.82f : 1.f);
        const float ageDrain = 0.0008f * o.age;
        const float thermalCost = 0.02f * thermalMismatch * (1.f + o.pheno.bodyDensity);
        const float hydrationCost = 0.018f * hydricMismatch * (1.f + 0.5f * o.pheno.bodyDensity);
        const float immuneCost = 0.02f * o.pathogenLoad * (0.5f + o.pheno.immunity);
        const float hypoxiaCost = 0.028f * hypoxia * (1.f + o.pheno.bodyDensity);
        const float osmoticCost = 0.024f * osmotic;
        const float acidCost = 0.021f * acidStress;
        const float stressOxCost = 0.018f * o.oxidativeStress;
        const float microbiomeTax = 0.016f * (1.1f - o.microbiomeHealth);
        const float outOfWaterStress = landFactor * std::max(0.f, 0.55f - o.terrestrialAffinity) * 0.08f;
        const float epigeneticTax = 0.012f * o.epigeneticStress;
        const float dormancyTax = o.dormancy * (0.45f * uptake + 0.4f * photo);
        const float delta = (uptake + photo - dormancyTax - toxHit - moveCost - actuationCost - maint - ageDrain - thermalCost - hydrationCost - immuneCost - stressOxCost - microbiomeTax - hypoxiaCost - osmoticCost - acidCost - outOfWaterStress - epigeneticTax) * dt * 20.f;
        o.energy += delta;
        if (o.energy < 0.f)
        {
            o.energyDebt += -o.energy * 0.3f;
            o.storedEnergy = std::max(0.f, o.storedEnergy + o.energy * 0.2f);
        }
        else
        {
            o.storedEnergy += delta * 0.1f;
            o.energyDebt = std::max(0.f, o.energyDebt - delta * 0.05f);
        }
        if (delta < -0.04f)
        {
            o.starvationTimer = clampf(o.starvationTimer + dt * (0.6f + 0.4f * std::max(0.f, -delta)), 0.f, 12.f);
        }
        else
        {
            o.starvationTimer = std::max(0.f, o.starvationTimer - dt * 0.8f);
        }
        if (o.starvationTimer > 2.5f && o.cells.size() > 2 && rndf() < dt * 0.25f)
        {
            const float recovered = 0.25f + 0.07f * o.pheno.scavenging;
            o.cells.pop_back();
            pruneSprings(o);
            o.energy += recovered;
            env.depositBiomass(com, recovered * 0.5f, 2.8f);
        }
        o.oxidativeStress = clampf(o.oxidativeStress * 0.995f + std::max(0.f, moveCost + toxHit + 0.4f * env.climateShock - 0.5f * o.pheno.immunity) * 0.01f, 0.f, 2.f);
        o.reproductionCooldown = std::max(0.f, o.reproductionCooldown - dt * (0.7f + 0.6f * o.pheno.reproductiveFlex));
        o.heatStressMemory = clampf(o.heatStressMemory * 0.99f + thermalMismatch * 0.015f, 0.f, 2.f);
        o.hydrationStressMemory = clampf(o.hydrationStressMemory * 0.99f + hydricMismatch * 0.015f, 0.f, 2.f);
        o.microbiomeHealth = clampf(o.microbiomeHealth - dt * (0.0012f * X + 0.0008f * o.pathogenLoad), 0.03f, 1.5f);
        o.oxygenStressMemory = clampf(o.oxygenStressMemory * 0.992f + hypoxia * 0.03f, 0.f, 2.2f);
        o.salinityStressMemory = clampf(o.salinityStressMemory * 0.992f + osmotic * 0.03f, 0.f, 2.2f);
        o.acidStressMemory = clampf(o.acidStressMemory * 0.992f + acidStress * 0.028f, 0.f, 2.2f);
        const float dormancyTarget = clampf(0.35f * stress + 0.2f * o.energyDebt - 0.22f * o.pheno.reproductiveFlex, 0.f, 1.f);
        o.dormancy = clampf(o.dormancy * 0.985f + dormancyTarget * 0.015f, 0.f, 1.f);
        if (toxHit > 0.08f) env.depositToxin(com, toxHit * 0.04f, 4.5f + 0.08f * static_cast<float>(o.cells.size()));
        o.lastReward = clampf(delta, -1.f, 1.f);
        const float targetTerrestrial = (landFactor > 0.5f) ? 1.f : 0.f;
        o.terrestrialAffinity = clampf(o.terrestrialAffinity * 0.995f + targetTerrestrial * 0.005f + clampf(delta, -0.2f, 0.2f) * 0.002f, 0.f, 1.f);

        if (delta > 0.f && rndf() < 0.09f * dt * 20.f) env.depositNutrient(com, 0.04f * delta, 3.5f + 0.05f * static_cast<float>(o.cells.size()));
        digestBiomass(o);
        growOrganism(o);
        attemptReproduction(o, newborns);

        if (o.energy > o.pheno.splitThreshold * 0.75f) o.offspringInvestmentReserve = clampf(o.offspringInvestmentReserve + 0.02f * dt * o.pheno.parentalInvestment, 0.f, 2.6f);
        if (rndf() < dt * (0.004f + 0.012f * o.pheno.dispersalDrive + 0.006f * localDensity(com, selfIndex))) o.deme = (o.deme + rndi(1, 4)) % 5;
        const float senescence = std::max(0.f, o.age / std::max(1.f, o.pheno.maxAge) - 0.65f);
        const float mortalityRisk = clampf(0.0003f + 0.018f * senescence * senescence + 0.01f * o.pathogenLoad + 0.006f * o.starvationTimer + 0.004f * o.epigeneticStress, 0.f, 0.35f);
        if (o.energy <= -3.f || o.age > o.pheno.maxAge || o.pathogenLoad > 3.2f || o.heatStressMemory > 1.8f || o.oxygenStressMemory > 1.9f || rndf() < mortalityRisk * dt * 20.f) o.alive = false;
    }

    void computeMetrics()
    {
        if (organisms.empty()) return;
        float fitSum = 0.f, fitSq = 0.f;
        float coop = 0.f, disease = 0.f, climate = 0.f, oxygen = 0.f, salinity = 0.f, acid = 0.f, dormancy = 0.f, genomeComplexity = 0.f, meanCells = 0.f;
        std::unordered_map<int, int> demeCount;
#if defined(EVO_HAS_OPENMP)
#pragma omp parallel for reduction(+:fitSum,fitSq,coop,disease,climate,oxygen,salinity,acid,dormancy,genomeComplexity,meanCells) schedule(static)
#endif
        for (int i = 0; i < static_cast<int>(organisms.size()); ++i)
        {
            const auto& o = organisms[static_cast<size_t>(i)];
            const float fitness = std::max(0.f, o.energy + 0.5f * o.storedEnergy - o.energyDebt - 0.4f * o.pathogenLoad);
            fitSum += fitness;
            fitSq += fitness * fitness;
            coop += clampf(o.pheno.socialAffinity * o.socialDrive, 0.f, 2.f);
            disease += o.pathogenLoad;
            climate += 0.5f * (o.heatStressMemory + o.hydrationStressMemory);
            oxygen += o.oxygenStressMemory;
            salinity += o.salinityStressMemory;
            acid += o.acidStressMemory;
            dormancy += o.dormancy;
            for (const auto& m : o.genome.modules) genomeComplexity += static_cast<float>(m.genes.size());
            meanCells += static_cast<float>(o.cells.size());
        }
        for (const auto& o : organisms) demeCount[o.deme]++;
        const float n = static_cast<float>(organisms.size());
        metrics.meanFitness = fitSum / n;
        metrics.fitnessStd = std::sqrt(std::max(0.f, fitSq / n - metrics.meanFitness * metrics.meanFitness));
        metrics.heterozygosity = 0.f;
        for (const auto& [_, c] : demeCount)
        {
            const float p = c / n;
            metrics.heterozygosity += p * (1.f - p);
        }
        metrics.trophicStability = clampf(1.f - metrics.fitnessStd / std::max(0.1f, metrics.meanFitness + 1.f), 0.f, 1.f);
        metrics.shannonDiversity = 0.f;
        for (const auto& [sig, c] : lineageBirths)
        {
            (void)sig;
            const float p = static_cast<float>(c) / std::max(1.f, static_cast<float>(std::max(1, birthsThisStep)));
            if (p > 1e-5f) metrics.shannonDiversity -= p * std::log(p);
        }
        metrics.nicheDiversity = 0.f;
        for (const auto& [_, c] : demeCount)
        {
            const float p = c / n;
            if (p > 1e-5f) metrics.nicheDiversity -= p * std::log(p);
        }
        metrics.cooperationIndex = coop / std::max(1.f, n);
        metrics.diseaseLoad = disease / std::max(1.f, n);
        metrics.climateStress = climate / std::max(1.f, n);
        metrics.oxygenStress = oxygen / std::max(1.f, n);
        metrics.salinityStress = salinity / std::max(1.f, n);
        metrics.acidStress = acid / std::max(1.f, n);
        metrics.dormancyRatio = dormancy / std::max(1.f, n);
        metrics.meanGenomeComplexity = genomeComplexity / std::max(1.f, n);
        metrics.meanCellCount = meanCells / std::max(1.f, n);
        metrics.reproductionRate = static_cast<float>(birthsThisStep) / std::max(1.f, n);
    }

    void exportMetricsIfNeeded()
    {
        if (steps % 200 != 0) return;
        std::ofstream f("sim_metrics.csv", std::ios::app);
        f << steps << "," << organisms.size() << "," << metrics.heterozygosity << "," << metrics.meanFitness << "," << metrics.fitnessStd << "," << metrics.trophicStability << "," << metrics.cooperationIndex << "," << metrics.diseaseLoad << "," << metrics.climateStress << "," << metrics.oxygenStress << "," << metrics.salinityStress << "," << metrics.acidStress << "," << metrics.dormancyRatio << "," << metrics.reproductionRate << "," << metrics.meanGenomeComplexity << "," << metrics.extinctions << "," << metrics.speciations << "," << metrics.shannonDiversity << "," << metrics.nicheDiversity << "," << metrics.meanCellCount << "\n";
    }

    void step()
    {
        ++steps;
        birthsThisStep = 0;
        lineageBirths.clear();
        env.step(dt);
        cachedSeasonality = reproductionSeasonality();
        buildGrid();
        huntInteractions();
        std::vector<Organism> newborns;
        newborns.reserve(organisms.size() / 4 + 8);
        for (int i = 0; i < static_cast<int>(organisms.size()); ++i) updateOne(i, organisms[i], newborns);

        for (const auto& o : organisms)
            if (!o.alive)
            {
                const Vec3 corpse = organismCenterOfMass(o);
                env.depositBiomass(corpse, std::max(0.f, o.energy) + 0.45f * static_cast<float>(o.cells.size()), 4.f + 0.08f * static_cast<float>(o.cells.size()));
                if (o.pathogenLoad > 0.6f) env.depositToxin(corpse, 0.08f * o.pathogenLoad * static_cast<float>(o.cells.size()), 3.2f + 0.05f * static_cast<float>(o.cells.size()));
            }

        const size_t before = organisms.size();
        organisms.erase(std::remove_if(organisms.begin(), organisms.end(), [](const Organism& o) { return !o.alive; }), organisms.end());
        metrics.extinctions += static_cast<int>(before - organisms.size());
        for (auto& n : newborns) organisms.push_back(std::move(n));

        targetPopulation = static_cast<int>(clampf(260.f + 180.f * (1.f - env.climateShock) + 40.f * env.seasonal(), 180.f, 520.f));

        if (static_cast<int>(organisms.size()) < targetPopulation / 3)
        {
            int toAdd = std::max(0, targetPopulation / 2 - static_cast<int>(organisms.size()));
            for (int i = 0; i < toAdd; i++)
            {
                Organism o;
                o.id = g_nextId++;
                o.genome = mutateGenome(makeBaseGenome(), 0.6f);
                o.pheno = interpretGenome(o.genome, env.temperature({0,0,0}), env.nutrient({0,0,0}), env.toxin({0,0,0}));
                o.geneticSignature = genomeSignature(o.genome);
                o.matePreference = rndf(0.f, 1.f);
                o.socialDrive = rndf(0.f, 1.f);
                o.pos = {rndf(-45, 45), rndf(-25, 45), rndf(-45, 45)};
                o.vel = {rndf(-.1f, .1f), rndf(-.1f, .1f), rndf(-.1f, .1f)};
                o.energy = rndf(6.f, 10.f);
                o.deme = rndi(0, 4);
                o.terrestrialAffinity = clampf(0.2f + 0.6f * o.pheno.hydricTolerance + rndf(-0.2f, 0.2f), 0.f, 1.f);
                o.cells.push_back({Vec3{0, 0, 0}, Vec3{0, 0, 0}, 1.f, 0.f, 1.f});
                o.epigeneticStress = rndf(0.f, 0.4f);
                organisms.push_back(std::move(o));
            }
        }

        computeMetrics();
        exportMetricsIfNeeded();
    }
};
