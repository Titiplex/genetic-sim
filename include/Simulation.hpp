#pragma once
#include <Cell.hpp>
#include <Environment.hpp>
#include <unordered_map>

static uint64_t g_nextId = 1;

// ======================= Simulation ==========================
struct Sim
{
    Environment           env;
    std::vector<Organism> organisms;

    float dt               = 0.035f;
    float cellSpacing      = 1.25f;
    int   targetPopulation = 180;
    float worldRadius      = 80.f;

    std::unordered_map<IVec3, std::vector<int>, IVec3Hash> grid;
    float                                                  gridSize = 7.f;

    // ---------- helpers ----------
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
        return IVec3{
            static_cast<int>(std::floor(p.x / gridSize)),
            static_cast<int>(std::floor(p.y / gridSize)),
            static_cast<int>(std::floor(p.z / gridSize))
        };
    }

    void buildGrid()
    {
        grid.clear();
        for (int i = 0; i < static_cast<int>(organisms.size()); ++i)
        {
            if (!organisms[i].alive) continue;
            grid[cellOf(organisms[i].pos)].push_back(i);
        }
    }

    std::vector<int> nearbyOrganisms(const Vec3& p)
    {
        std::vector<int> out;
        const auto       [x, y, z] = cellOf(p);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz)
                {
                    const IVec3 q{x + dx, y + dy, z + dz};
                    if (auto it = grid.find(q); it != grid.end())
                        out.insert(out.end(), it->second.begin(), it->second.end());
                }
        return out;
    }

    // ---------- seeding ----------
    void seedInitial(const int n)
    {
        organisms.clear();
        organisms.reserve(n * 2);

        for (int i = 0; i < n; i++)
        {
            Organism o;
            o.id           = g_nextId++;
            o.genome.genes = {{"a"}};

            if (i > 0)
            {
                o.genome = mutateGenome(o.genome, 0.45f);
                if (rndf() < 0.5f) o.genome = mutateGenome(o.genome, 0.2f);
            }

            o.pheno  = interpretGenome(o.genome);
            o.pos    = {rndf(-25.f, 25.f), rndf(-25.f, 25.f), rndf(-25.f, 25.f)};
            o.vel    = {rndf(-0.2f, 0.2f), rndf(-0.2f, 0.2f), rndf(-0.2f, 0.2f)};
            o.energy = rndf(6.f, 12.f);

            // single starting cell
            o.cells.push_back({Vec3{0, 0, 0}, Vec3{0, 0, 0}, 1.f, 0.f, 1.f});

            organisms.push_back(std::move(o));
        }
    }

    // ---------- dynamic growth direction ----------
    Vec3 chooseGrowthDirection(const Organism& o) const
    {
        const Vec3 p  = o.pos;
        const Vec3 gN = env.gradN(p);
        const Vec3 gL = env.gradL(p);
        const Vec3 gX = env.gradX(p);
        const Vec3 gB = env.gradB(p);

        const Vec3 localCom = organismLocalCOM(o);

        Vec3 principal{0, 0, 0};
        Vec3 shellOut{0, 0, 0};
        for (const auto& c : o.cells)
        {
            const Vec3 d = c.localPos - localCom;
            if (const float l = len(d); l > 1e-4f)
            {
                principal += normalize(d) * (0.1f + l);
                if (l > 0.75f) shellOut += normalize(d);
            }
        }
        principal = normalize(principal);
        shellOut  = normalize(shellOut);

        Vec3 repulse{0, 0, 0};
        for (const auto& c : o.cells)
        {
            const Vec3 d = localCom - c.localPos;
            if (const float l = len(d); l > 1e-4f) repulse += d / (l * l + 0.2f);
        }
        repulse = normalize(repulse);

        const Vec3 randomDrift = normalize(Vec3{rndf(-1, 1), rndf(-1, 1), rndf(-1, 1)});

        const std::array<Vec3, 16> feat = {
            normalize(gN),
            normalize(gL),
            normalize(gX) * -1.f,
            normalize(gB),
            normalize(o.vel),
            principal,
            repulse,
            shellOut,
            randomDrift,
            normalize(p * -1.f),
            normalize(gN + gL * 0.5f - gX * 0.7f),
            normalize(gN + gB * 0.4f),
            normalize(gL - gX),
            Vec3{1, 0, 0},
            Vec3{0, 1, 0},
            Vec3{0, 0, 1}
        };

        Vec3 dir{0, 0, 0};
        for (int i = 0; i < 16; i++) dir += feat[i] * o.pheno.w[i];
        for (int i = 0; i < 8; i++) dir += o.pheno.basis[i] * o.pheno.basisGain[i] * (0.15f + 0.85f * rndf());

        if (len2(dir) < 1e-6f) dir = randomDrift;
        return normalize(dir);
    }

    // ---------- body mechanics ----------
    void applySprings(Organism& o) const
    {
        for (const auto& [a, b, restLen, stiffness, damping] : o.springs)
        {
            if (a < 0 || b < 0 || a >= static_cast<int>(o.cells.size()) || b >= static_cast<int>(o.cells.size()))
                continue;

            Cell& A = o.cells[a];
            Cell& B = o.cells[b];

            const Vec3  d = B.localPos - A.localPos;
            const float l = len(d);
            if (l < 1e-6f) continue;

            const Vec3  n    = d / l;
            const float relv = dot(B.localVel - A.localVel, n);

            const float forceMag = stiffness * (l - restLen) + damping * relv;
            const Vec3  force    = n * forceMag;

            A.localVel += force * (dt / std::max(0.2f, A.mass));
            B.localVel -= force * (dt / std::max(0.2f, B.mass));
        }
    }

    void applySoftCellRepulsion(Organism& o) const
    {
        const float minDist = std::max(0.55f, 1.15f * o.pheno.cellRadius);

        for (int i = 0; i < static_cast<int>(o.cells.size()); ++i)
        {
            for (int j = i + 1; j < static_cast<int>(o.cells.size()); ++j)
            {
                Vec3        d = o.cells[j].localPos - o.cells[i].localPos;
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
                    constexpr float k   = 8.0f;
                    const Vec3      n   = d / l;
                    const float     pen = minDist - l;
                    const Vec3      f   = n * (k * pen);
                    o.cells[i].localVel -= f * dt;
                    o.cells[j].localVel += f * dt;
                }
            }
        }
    }

    void integrateInternalBody(Organism& o) const
    {
        if (o.cells.empty()) return;

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
        for (int i = 0; i < newCellIndex; ++i)
        {
            ds.emplace_back(len2(o.cells[i].localPos - o.cells[newCellIndex].localPos), i);
        }
        std::sort(ds.begin(), ds.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        int linkCount = 1;
        if (o.pheno.adhesionBias > 0.7f) linkCount++;
        if (o.pheno.adhesionBias > 1.15f) linkCount++;

        linkCount = std::min(linkCount, static_cast<int>(ds.size()));

        for (int k = 0; k < linkCount; ++k)
        {
            const int j = ds[k].second;
            Spring    s;
            s.a         = newCellIndex;
            s.b         = j;
            s.restLen   = clampf(std::sqrt(ds[k].first), 0.45f, 2.6f);
            s.stiffness = 5.f + 6.0f * o.pheno.adhesionBias;
            s.damping   = 0.35f;
            o.springs.push_back(s);
        }
    }

    // ---------- growth ----------
    void growOrganism(Organism& o) const
    {
        if (o.cells.size() >= 64) return;
        if (o.energy < 2.0f) return;

        if (const float desire = o.pheno.growthBias * (0.35f + 0.65f *
            clampf(o.energy / o.pheno.splitThreshold, 0.f, 2.f)); rndf() > clampf(desire * dt * 1.7f, 0.f, 0.30f))
            return;

        const Vec3 dir = chooseGrowthDirection(o);

        // choose anchor biased toward surface
        int   anchorIdx = 0;
        float bestL2    = -1.f;
        for (int i = 0; i < static_cast<int>(o.cells.size()); ++i)
        {
            if (const float d2 = len2(o.cells[i].localPos); d2 > bestL2)
            {
                bestL2    = d2;
                anchorIdx = i;
            }
        }
        if (rndf() < 0.55f) anchorIdx = rndi(0, static_cast<int>(o.cells.size()) - 1);

        const Vec3 anchor    = o.cells[anchorIdx].localPos;
        Vec3       bestPos   = anchor + dir * (cellSpacing * o.pheno.cellRadius);
        float      bestScore = -1e9f;

        for (int k = 0; k < 10; k++)
        {
            const Vec3 candDir = normalize(
                dir + hashToVec3(static_cast<uint64_t>(k) + o.id + static_cast<uint64_t>(o.age * 37.f)) * rndf(
                    0.0f, 0.9f));
            const Vec3 cand = anchor + candDir * (cellSpacing * o.pheno.cellRadius * rndf(0.9f, 1.45f));

            float minD = 1e9f;
            float rep  = 0.f;
            for (const auto& c : o.cells)
            {
                const float d = len(c.localPos - cand);
                minD          = std::min(minD, d);
                if (d < 1.7f * o.pheno.cellRadius) rep += (1.7f * o.pheno.cellRadius - d);
            }

            if (const float score = minD - rep * 0.9f + rndf(-0.15f, 0.15f); score > bestScore)
            {
                bestScore = score;
                bestPos   = cand;
            }
        }

        if (bestScore > 0.12f)
        {
            Cell c;
            c.localPos = bestPos;
            c.localVel = dir * rndf(0.03f, 0.15f);
            c.energy   = 1.f;
            c.age      = 0.f;
            c.mass     = 1.f;
            o.cells.push_back(c);
            connectNewCell(o, static_cast<int>(o.cells.size()) - 1);

            o.energy -= 0.55f + 0.035f * static_cast<float>(o.cells.size());
        }
    }

    // ---------- ecology ----------
    void huntInteractions()
    {
        buildGrid();

        for (int i = 0; i < static_cast<int>(organisms.size()); ++i)
        {
            Organism& a = organisms[i];
            if (!a.alive) continue;
            if (a.pheno.aggression < 0.30f) continue;
            if (a.energy < 2.0f) continue;

            auto nearby = nearbyOrganisms(a.pos);
            for (const int j : nearby)
            {
                if (i == j) continue;
                Organism& b = organisms[j];
                if (!b.alive) continue;

                const Vec3  d  = b.pos - a.pos;
                const float r2 = len2(d);
                if (const float huntRange = 2.5f + 0.08f * static_cast<float>(a.cells.size()); r2 > huntRange *
                    huntRange)
                    continue;

                const float aPower = a.pheno.aggression + 0.018f * static_cast<float>(a.cells.size());
                const float bDef   = 0.5f * b.pheno.toxinResistance + 0.010f * static_cast<float>(b.cells.size());

                const float damage = clampf((aPower - 0.4f * bDef) * dt * 4.2f, 0.f, 0.7f);
                if (damage <= 0.f) continue;

                b.energy -= damage;
                a.energy += damage * 0.50f;
                a.energy -= env.toxin(b.pos) * (1.0f - clampf(a.pheno.toxinResistance / 2.0f, 0.f, 0.9f)) * dt * 0.4f;
            }
        }
    }

    void digestBiomass(Organism& o)
    {
        const Vec3  com = organismCenterOfMass(o);
        const float B   = env.biomass(com);
        if (B <= 0.0001f) return;

        const float request = o.pheno.scavenging * (0.02f + 0.005f * static_cast<float>(o.cells.size())) * dt * 20.f;
        if (request <= 0.f) return;

        const float gain = std::min(B, request) * 0.9f;
        o.energy += gain * 0.8f; // imperfect conversion
        env.consumeBiomassNear(com, gain);
    }


    Vec3 computeActuation(Organism& o, const Vec3& com, const Vec3& gN, const Vec3& gL, const Vec3& gX, const Vec3& gB) const
    {
        const float waterDepth = env.waterSurfaceHeight(com) - com.y;
        const float submerged = clampf(waterDepth / 12.f, 0.f, 1.f);
        const float onGround = clampf((env.groundHeight(com) + o.pheno.cellRadius - com.y) / 2.5f, 0.f, 1.f);

        const std::array<float, 12> feat = {
            clampf(len(o.vel) / 4.5f, 0.f, 1.f),
            clampf(o.energy / std::max(2.f, o.pheno.splitThreshold), 0.f, 2.f),
            clampf(o.age / std::max(20.f, o.pheno.maxAge), 0.f, 2.f),
            clampf(len(gN) * 2.2f, 0.f, 1.f),
            clampf(len(gL) * 2.2f, 0.f, 1.f),
            clampf(len(gX) * 2.2f, 0.f, 1.f),
            clampf(len(gB) * 2.2f, 0.f, 1.f),
            submerged,
            onGround,
            clampf(static_cast<float>(o.cells.size()) / 28.f, 0.f, 2.f),
            clampf(env.shorelineBlend(com), 0.f, 1.f),
            1.f
        };

        float command = 0.f;
        for (int i = 0; i < 12; ++i) command += feat[i] * o.pheno.motorW[i];
        command = std::tanh(command);

        o.gaitPhase += dt * (0.25f + 0.75f * o.pheno.fineness);

        Vec3 thrust{0.f, 0.f, 0.f};
        for (int i = 0; i < 6; ++i)
        {
            const Vec3 axis = normalize(o.pheno.actuatorAxis[i]);
            const float osc = std::sin(o.gaitPhase * (0.8f + o.pheno.actuatorFreq[i]) + o.pheno.actuatorBias[i] * 3.14159f);
            const float amp = (0.35f + 0.65f * submerged + 0.55f * onGround) * o.pheno.motorDrive;
            thrust += axis * (command * osc * amp);
        }

        return thrust;
    }

    // ---------- update organism ----------
    void updateOne(Organism& o, std::vector<Organism>& newborns)
    {
        if (!o.alive) return;

        o.age += dt;

        // body mechanics first
        integrateInternalBody(o);

        Vec3 com = organismCenterOfMass(o);
        const float bodyScale = clampf(static_cast<float>(o.cells.size()) / 18.f, 0.35f, 3.0f);

        const float N = env.nutrient(com);
        const float L = env.light(com);
        const float X = env.toxin(com);

        const Vec3 gN = env.gradN(com);
        const Vec3 gL = env.gradL(com);
        const Vec3 gX = env.gradX(com);
        const Vec3 gB = env.gradB(com);

        Vec3 steer =
            normalize(gN) * (0.38f * o.pheno.nutrientAffinity)
            + normalize(gL) * (0.25f * o.pheno.photoAffinity)
            + normalize(gB) * (0.20f * o.pheno.scavenging)
            - normalize(gX) * (0.50f * (1.2f - clampf(o.pheno.toxinResistance, 0.f, 1.2f)))
            + hashToVec3(o.id + static_cast<uint64_t>(o.age * 10.0f)) * 0.12f;

        if (len2(steer) < 1e-6f) steer = Vec3{rndf(-1, 1), rndf(-1, 1), rndf(-1, 1)};
        steer = normalize(steer);

        const Vec3 flow      = env.fluidVelocity(com);
        const Vec3 gust      = env.turbulence(com, o.id + static_cast<uint64_t>(o.age * 100.f));
        const float buoyancy = env.buoyancy(com, bodyScale);
        const Vec3 actuation = computeActuation(o, com, gN, gL, gX, gB);

        o.vel += steer * (0.48f * dt);
        o.vel += flow * (0.38f * dt);
        o.vel += gust * dt;
        o.vel += actuation * (0.95f * dt);
        o.vel.y += buoyancy * dt * (1.0f + 0.2f * (1.3f - o.pheno.bodyDensity));

        const Vec3 relFluid = o.vel - flow;
        const float drag    = 0.10f + 0.06f * bodyScale + 0.03f * o.pheno.bodyDensity;
        o.vel -= relFluid * drag * dt;

        const float terrainDrag = clampf((env.groundHeight(o.pos) + 1.0f - o.pos.y) / 2.0f, 0.f, 1.f);
        o.vel.x *= (1.f - 0.08f * terrainDrag * dt * 20.f);
        o.vel.z *= (1.f - 0.08f * terrainDrag * dt * 20.f);

        o.vel = clampVec(o.vel, 5.6f / (1.0f + 0.024f * static_cast<float>(o.cells.size())));
        o.pos += o.vel * dt * 6.0f;

        const float floorY = env.groundHeight(o.pos) + o.pheno.cellRadius * 0.9f;
        if (o.pos.y < floorY)
        {
            o.pos.y = floorY;
            if (o.vel.y < 0.f) o.vel.y *= -0.18f;
            o.vel.x *= 0.92f;
            o.vel.z *= 0.92f;
        }

        if (const float r = len(o.pos); r > worldRadius)
        {
            const Vec3 inward = normalize(o.pos) * -1.f;
            o.vel += inward * dt * 4.f;
            o.pos = normalize(o.pos) * worldRadius;
            o.vel *= 0.8f;
        }

        // Energy model with 3 factors + biomass digestion
        const float uptake = o.pheno.nutrientAffinity * N * (0.05f + 0.010f * static_cast<float>(o.cells.size()));
        const float photo  = o.pheno.photoAffinity * L * (0.03f + 0.008f * static_cast<float>(o.cells.size()));
        const float toxHit = std::max(0.f, X - 0.7f * o.pheno.toxinResistance) * (0.04f + 0.005f * static_cast<float>(o.
                                                                                                                      cells
                                                                                                                      .size()));

        const float moveCost = 0.015f * len(o.vel) * (1.0f + 0.03f * static_cast<float>(o.cells.size()));
        const float actuationCost = 0.010f * len(actuation) * (0.6f + 0.8f * o.pheno.motorDrive);
        const float turbulenceCost = 0.006f * len(gust) * (1.0f + 0.02f * static_cast<float>(o.cells.size()));
        const float maint    = o.pheno.baseMetabolism * (1.0f + 0.055f * static_cast<float>(o.cells.size()));
        const float ageDrain = 0.0008f * o.age;

        o.energy += (uptake + photo - toxHit - moveCost - actuationCost - turbulenceCost - maint - ageDrain) * dt * 20.f;

        digestBiomass(o);

        growOrganism(o);

        // Reproduction (budding / split)
        if (o.energy > o.pheno.splitThreshold && o.age > 6.f && o.cells.size() >= 2)
        {
            Organism child;
            child.id     = g_nextId++;
            child.genome = mutateGenome(o.genome, o.pheno.mutationRate);
            child.pheno  = interpretGenome(child.genome);

            child.pos    = o.pos + normalize(Vec3{rndf(-1, 1), rndf(-1, 1), rndf(-1, 1)}) * rndf(2.f, 6.f);
            child.vel    = o.vel + Vec3{rndf(-0.4f, 0.4f), rndf(-0.4f, 0.4f), rndf(-0.4f, 0.4f)};
            child.energy = o.energy * 0.36f;
            o.energy *= 0.54f;

            const int take = std::min<int>(std::max<int>(1, static_cast<int>(o.cells.size()) / 3), 12);

            // take farthest cells first (surface budding)
            std::vector<std::pair<float, int>> order;
            order.reserve(o.cells.size());
            for (int i = 0; i < static_cast<int>(o.cells.size()); ++i)
                order.emplace_back(len2(o.cells[i].localPos), i);

            std::sort(order.begin(), order.end(),
                      [](const auto& a, const auto& b) { return a.first > b.first; });

            std::vector<int> idxs;
            idxs.reserve(take);
            for (int k = 0; k < take && k < static_cast<int>(order.size()); ++k) idxs.push_back(order[k].second);
            std::sort(idxs.begin(), idxs.end(), std::greater<int>());

            for (const int idx : idxs)
            {
                Cell c = o.cells[idx];
                c.localPos *= 0.6f;
                c.localVel = Vec3{0, 0, 0};
                c.age      = 0.f;
                child.cells.push_back(c);
                o.cells.erase(o.cells.begin() + idx);
            }

            if (child.cells.empty())
            {
                child.cells.push_back({Vec3{0, 0, 0}, Vec3{0, 0, 0}, 1.f, 0.f, 1.f});
            }

            // recenter child body
            {
                Vec3 cc{0, 0, 0};
                for (const auto& c : child.cells) cc += c.localPos;
                cc = cc / static_cast<float>(child.cells.size());
                for (auto& c : child.cells) c.localPos -= cc;
            }

            // rebuild simple child springs
            child.springs.clear();
            for (int i = 1; i < static_cast<int>(child.cells.size()); ++i)
            {
                int    j = rndi(0, i - 1);
                Spring s;
                s.a         = i;
                s.b         = j;
                s.restLen   = clampf(len(child.cells[i].localPos - child.cells[j].localPos), 0.45f, 2.6f);
                s.stiffness = 5.f + 6.f * child.pheno.adhesionBias;
                s.damping   = 0.35f;
                child.springs.push_back(s);
            }

            // parent may have too few springs left -> prune invalid + optionally reconnect lightly
            o.springs.erase(std::remove_if(o.springs.begin(), o.springs.end(),
                                           [&](const Spring& s)
                                           {
                                               return s.a < 0 || s.b < 0 || s.a >= static_cast<int>(o.cells.size()) || s
                                                   .b >= static_cast<int>(o.
                                                                          cells
                                                                          .size());
                                           }),
                            o.springs.end());

            newborns.push_back(std::move(child));
        }

        // Death
        if (o.energy <= 0.f || o.age > o.pheno.maxAge)
        {
            o.alive = false;
        }
    }

    // ---------- step ----------
    void step()
    {
        env.step(dt);

        // predation before metabolism
        huntInteractions();

        std::vector<Organism> newborns;
        newborns.reserve(organisms.size() / 4 + 8);

        for (auto& o : organisms)
        {
            updateOne(o, newborns);
        }

        // dead => biomass deposits
        for (const auto& o : organisms)
        {
            if (!o.alive)
            {
                const Vec3  com    = organismCenterOfMass(o);
                const float amount = std::max(0.f, o.energy) + 0.45f * static_cast<float>(o.cells.size());
                env.depositBiomass(com, amount, 4.f + 0.08f * static_cast<float>(o.cells.size()));
            }
        }

        organisms.erase(
            std::remove_if(organisms.begin(), organisms.end(),
                           [](const Organism& o) { return !o.alive; }),
            organisms.end()
        );

        for (auto& n : newborns) organisms.push_back(std::move(n));

        // Population floor for testing
        if (static_cast<int>(organisms.size()) < targetPopulation / 3)
        {
            int toAdd = targetPopulation / 2 - static_cast<int>(organisms.size());
            toAdd     = std::max(0, toAdd);

            for (int i = 0; i < toAdd; i++)
            {
                Organism o;
                o.id           = g_nextId++;
                o.genome.genes = {{"a"}};
                o.genome       = mutateGenome(o.genome, 0.6f);
                o.pheno        = interpretGenome(o.genome);
                o.pos          = {rndf(-20, 20), rndf(-20, 20), rndf(-20, 20)};
                o.vel          = {rndf(-.1f, .1f), rndf(-.1f, .1f), rndf(-.1f, .1f)};
                o.energy       = rndf(6.f, 10.f);
                o.cells.push_back({Vec3{0, 0, 0}, Vec3{0, 0, 0}, 1.f, 0.f, 1.f});
                organisms.push_back(std::move(o));
            }
        }
    }
};
