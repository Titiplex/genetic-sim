#include "evo/Sim.hpp"
#include "evo/Hash.hpp"
#include <algorithm>

namespace evo
{
    Sim::Sim(const uint64_t seed): m_rng(seed)
    {
    }

    void Sim::seedInitial(const int n)
    {
        m_orgs.clear();
        m_orgs.reserve(static_cast<size_t>(n) * 2);

        for (int i = 0; i < n; ++i)
        {
            Organism o;
            o.id           = m_nextId++;
            o.genome.genes = {{"a"}};
            if (i > 0)
            {
                o.genome = mutateGenome(o.genome, 0.45f, m_rng);
                if (m_rng.uniform() < 0.5f)
                    o.genome = mutateGenome(o.genome, 0.2f, m_rng);
            }

            o.pheno  = interpretGenome(o.genome);
            o.pos    = {m_rng.uniform(-25.f, 25.f), m_rng.uniform(-25.f, 25.f), m_rng.uniform(-25.f, 25.f)};
            o.vel    = {m_rng.uniform(-0.2f, 0.2f), m_rng.uniform(-0.2f, 0.2f), m_rng.uniform(-0.2f, 0.2f)};
            o.energy = m_rng.uniform(6.f, 12.f);

            o.cells.push_back({Vec3{0, 0, 0}, Vec3{0, 0, 0}, 0.f, 1.f});
            updateCellModes(o);
            m_orgs.push_back(std::move(o));
        }
    }

    IVec3 Sim::cellOf(const Vec3 &p) const
    {
        return IVec3{
            static_cast<int>(std::floor(p.x / m_gridSize)),
            static_cast<int>(std::floor(p.y / m_gridSize)),
            static_cast<int>(std::floor(p.z / m_gridSize)),
        };
    }

    void Sim::buildGrid()
    {
        m_grid.clear();
        for (int i = 0; i < static_cast<int>(m_orgs.size()); ++i)
        {
            if (!m_orgs[static_cast<size_t>(i)].alive)
                continue;
            m_grid[cellOf(m_orgs[static_cast<size_t>(i)].pos)].push_back(i);
        }
    }

    std::vector<int> Sim::nearby(const Vec3 &p)
    {
        std::vector<int> out;
        const auto       [x, y, z] = cellOf(p);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz)
                {
                    const IVec3 q{x + dx, y + dy, z + dz};
                    if (auto it = m_grid.find(q); it != m_grid.end())
                        out.insert(out.end(), it->second.begin(), it->second.end());
                }
        return out;
    }

    Vec3 Sim::centerOfMass(const Organism &o)
    {
        if (o.cells.empty())
            return o.pos;
        Vec3 s{0, 0, 0};
        for (const auto &c : o.cells)
            s += o.pos + c.localPos;
        return s / static_cast<float>(o.cells.size());
    }

    Vec3 Sim::localCOM(const Organism &o)
    {
        if (o.cells.empty())
            return Vec3{0, 0, 0};
        Vec3 s{0, 0, 0};
        for (const auto &c : o.cells)
            s += c.localPos;
        return s / static_cast<float>(o.cells.size());
    }

    void Sim::recenterLocalBody(Organism &o)
    {
        if (o.cells.empty())
            return;
        const Vec3 lc = localCOM(o);
        for (auto &c : o.cells)
            c.localPos -= lc;
        o.pos += lc;
    }

    Vec3 Sim::chooseGrowthDirection(const Organism &o) const
    {
        const Vec3 p  = o.pos;
        const Vec3 gN = m_env.gradN(p);
        const Vec3 gL = m_env.gradL(p);
        const Vec3 gX = m_env.gradX(p);
        const Vec3 gB = m_env.gradB(p);

        const Vec3 lc = localCOM(o);

        Vec3 principal{0, 0, 0};
        Vec3 shellOut{0, 0, 0};
        for (const auto &c : o.cells)
        {
            const Vec3 d = c.localPos - lc;
            if (const float l = len(d); l > 1e-4f)
            {
                principal += normalize(d) * (0.1f + l);
                if (l > 0.75f)
                    shellOut += normalize(d);
            }
        }
        principal = normalize(principal);
        shellOut  = normalize(shellOut);

        Vec3 repulse{0, 0, 0};
        for (const auto &c : o.cells)
        {
            const Vec3 d = lc - c.localPos;
            if (const float l = len(d); l > 1e-4f)
                repulse += d / (l * l + 0.2f);
        }
        repulse = normalize(repulse);

        const Vec3 randomDrift = normalize(Vec3{m_rng.uniform(-1, 1), m_rng.uniform(-1, 1),
                                                m_rng.uniform(-1, 1)});

        const std::array feat = {
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
            Vec3{0, 0, 1},
        };

        Vec3 dir{0, 0, 0};
        for (int i = 0; i < 16; ++i)
            dir += feat[static_cast<size_t>(i)] * o.pheno.w[static_cast<size_t>(i)];
        for (int i = 0; i < 8; ++i)
            dir += o.pheno.basis[static_cast<size_t>(i)] * o.pheno.basisGain[static_cast<size_t>(i)] * 0.6f;

        if (len2(dir) < 1e-6f)
            dir = randomDrift;
        return normalize(dir);
    }

    void Sim::applySprings(Organism &o) const
    {
        for (const auto &[a, b, restLen, stiffness, damping] : o.springs)
        {
            if (a < 0 || b < 0 || a >= static_cast<int>(o.cells.size()) || b >= static_cast<int>(o.cells.
                                                                                                   size()))
                continue;

            Cell &A = o.cells[static_cast<size_t>(a)];
            Cell &B = o.cells[static_cast<size_t>(b)];

            const Vec3  d = B.localPos - A.localPos;
            const float l = len(d);
            if (l < 1e-6f)
                continue;

            const Vec3  n    = d / l;
            const float relv = dot(B.localVel - A.localVel, n);

            const float forceMag = stiffness * (l - restLen) + damping * relv;
            const Vec3  force    = n * forceMag;

            A.localVel += force * (dt / std::max(0.2f, A.mass));
            B.localVel -= force * (dt / std::max(0.2f, B.mass));
        }
    }

    void Sim::softRepulsion(Organism &o) const
    {
        const float minDist = std::max(0.55f, 1.15f * o.pheno.cellRadius);

        for (int i = 0; i < static_cast<int>(o.cells.size()); ++i)
        {
            for (int j = i + 1; j < static_cast<int>(o.cells.size()); ++j)
            {
                Vec3 d = o.cells[static_cast<size_t>(j)].localPos - o.cells[static_cast<size_t>(i)].localPos;
                const float l = len(d);

                if (l < 1e-5f)
                    continue;
                if (l < minDist)
                {
                    constexpr float k   = 8.0f;
                    const Vec3      n   = d / l;
                    const float     pen = minDist - l;
                    const Vec3      f   = n * (k * pen);
                    o.cells[static_cast<size_t>(i)].localVel -= f * dt;
                    o.cells[static_cast<size_t>(j)].localVel += f * dt;
                }
            }
        }
    }

    void Sim::integrateBody(Organism &o) const
    {
        if (o.cells.empty())
            return;
        applySprings(o);
        softRepulsion(o);

        for (auto &c : o.cells)
        {
            c.localVel *= 0.95f;
            c.localPos += c.localVel * dt;
            c.age += dt;
        }
        recenterLocalBody(o);
    }

    void Sim::connectNewCell(Organism &o, const int newIdx)
    {
        if (newIdx <= 0)
            return;

        const auto &mNew = o.pheno.modes[static_cast<size_t>(o.cells[static_cast<size_t>(newIdx)].modeId)];

        std::vector<std::pair<float, int>> ds;
        ds.reserve(static_cast<size_t>(newIdx));
        for (int i = 0; i < newIdx; ++i)
            ds.emplace_back(
                len2(o.cells[static_cast<size_t>(i)].localPos - o.cells[static_cast<size_t>(newIdx)].
                     localPos), i);

        std::ranges::sort(ds, [](const auto &a, const auto &b)
        {
            return a.first < b.first;
        });

        int linkCount = 1;
        if (o.pheno.adhesionBias > 0.7f)
            linkCount++;
        if (o.pheno.adhesionBias > 1.15f)
            linkCount++;
        linkCount = std::min(linkCount, static_cast<int>(ds.size()));

        for (int k = 0; k < linkCount; ++k)
        {
            const int j = ds[static_cast<size_t>(k)].second;
            Spring    s;
            s.a       = newIdx;
            s.b       = j;
            s.restLen = clampf(std::sqrt(ds[static_cast<size_t>(k)].first), 0.45f, 2.6f);

            s.stiffness = (5.f + 6.0f * o.pheno.adhesionBias) * clampf(mNew.gAdhesion, 0.6f, 1.8f);

            s.damping = 0.35f;
            o.springs.push_back(s);
        }
    }

    void Sim::grow(Organism &o) const
    {
        if (o.cells.size() >= 64)
            return;
        if (o.energy < 2.0f)
            return;

        const float desire = o.pheno.growthBias * (0.35f + 0.65f * clampf(
                                                       o.energy / o.pheno.splitThreshold, 0.f, 2.f));
        if (m_rng.uniform() > clampf(desire * dt * 1.7f, 0.f, 0.30f))
            return;

        const Vec3 dir = chooseGrowthDirection(o);

        int   anchorIdx = 0;
        float bestL2    = -1.f;
        for (int i = 0; i < static_cast<int>(o.cells.size()); ++i)
        {
            if (const float d2 = len2(o.cells[static_cast<size_t>(i)].localPos); d2 > bestL2)
            {
                bestL2    = d2;
                anchorIdx = i;
            }
        }
        if (m_rng.uniform() < 0.55f)
            anchorIdx = m_rng.uniformInt(0, static_cast<int>(o.cells.size()) - 1);

        const Vec3  anchor  = o.cells[static_cast<size_t>(anchorIdx)].localPos;
        const float spacing = 1.25f * o.pheno.cellRadius;

        Vec3  bestPos   = anchor + dir * spacing;
        float bestScore = -1e9f;

        for (int k = 0; k < 10; ++k)
        {
            const Vec3 candDir = normalize(dir + hashToVec3(o.id + static_cast<uint64_t>(k)) * 0.7f);
            const Vec3 cand    = anchor + candDir * spacing * m_rng.uniform(0.9f, 1.45f);

            float minD = 1e9f;
            float rep  = 0.f;
            for (const auto &c : o.cells)
            {
                const float d = len(c.localPos - cand);
                minD          = std::min(minD, d);
                if (d < 1.7f * o.pheno.cellRadius)
                    rep += (1.7f * o.pheno.cellRadius - d);
            }

            if (const float score = minD - rep * 0.9f + m_rng.uniform(-0.15f, 0.15f); score > bestScore)
            {
                bestScore = score;
                bestPos   = cand;
            }
        }

        if (bestScore > 0.12f)
        {
            Cell nc;
            nc.localPos = bestPos;
            nc.localVel = dir * 0.08f;
            nc.age      = 0.f;
            nc.mass     = 1.f;
            nc.modeId   = 0;
            o.cells.push_back(nc);

            const int  newIdx = static_cast<int>(o.cells.size()) - 1;
            const Vec3 wpNew = o.pos + o.cells[static_cast<size_t>(newIdx)].localPos;
            o.cells[static_cast<size_t>(newIdx)].modeId = chooseCellMode(
                o, wpNew, o.cells[static_cast<size_t>(newIdx)].localPos, 0.f,
                newIdx);

            connectNewCell(o, static_cast<int>(o.cells.size()) - 1);
            o.energy -= 0.55f + 0.035f * static_cast<float>(o.cells.size());
        }
    }

    void Sim::huntInteractions()
    {
        buildGrid();

        for (int i = 0; i < static_cast<int>(m_orgs.size()); ++i)
        {
            Organism &a = m_orgs[static_cast<size_t>(i)];
            if (!a.alive)
                continue;
            if (a.pheno.aggression < 0.30f)
                continue;
            if (a.energy < 2.0f)
                continue;

            for (const auto near = nearby(a.pos); const int j : near)
            {
                if (i == j)
                    continue;
                Organism &b = m_orgs[static_cast<size_t>(j)];
                if (!b.alive)
                    continue;

                float gAttack = 0.f;
                for (const auto &c : a.cells)
                    gAttack += a.pheno.modes[static_cast<size_t>(c.modeId)].gAttack;
                gAttack /= std::max(1.f, static_cast<float>(a.cells.size()));
                gAttack = clampf(gAttack, 0.f, 1.8f);

                if (const float huntRange = 2.5f + 0.08f * static_cast<float>(a.cells.size());
                    len2(b.pos - a.pos) > huntRange * huntRange)
                    continue;

                const float aPower = (a.pheno.aggression + 0.6f * gAttack) + 0.018f * static_cast<float>(a.
                                                                                                         cells
                                                                                                         .size());
                const float bDef = 0.5f * b.pheno.toxinResistance + 0.010f * static_cast<float>(b.cells.
                                                                                                  size());

                const float damage = clampf((aPower - 0.4f * bDef) * dt * 4.2f, 0.f, 0.7f);
                if (damage <= 0.f)
                    continue;

                b.energy -= damage;
                a.energy += damage * 0.50f;
                a.energy -= m_env.toxin(b.pos) * (1.0f - clampf(a.pheno.toxinResistance / 2.0f, 0.f, 0.9f)) *
                    dt * 0.4f;
            }
        }
    }

    void Sim::digestBiomass(Organism &o, float gDigestB)
    {
        const Vec3  com = centerOfMass(o);
        const float B   = m_env.biomass(com);
        if (B <= 0.0001f)
            return;

        gDigestB = clampf(gDigestB, 0.f, 2.2f);

        const float request    = gDigestB * (0.02f + 0.005f * static_cast<float>(o.cells.size())) * dt * 20.f;
        const float eaten      = std::min(B, request);
        const float efficiency = 0.5f + 0.35f * (gDigestB / 2.2f);
        o.energy += eaten * efficiency;
        m_env.consumeBiomassNear(com, eaten);
    }

    void Sim::updateOne(Organism &o, std::vector<Organism> &newborns)
    {
        if (!o.alive)
            return;

        o.age += dt;
        integrateBody(o);
        updateCellModes(o);

        float gUptakeN = 0.f, gPhoto = 0.f, gDigestB = 0.f, gResistX = 0.f, gMaint = 0.f;
        for (const auto &c : o.cells)
        {
            const auto &m = o.pheno.modes[static_cast<size_t>(c.modeId)];
            gUptakeN += m.gUptakeN;
            gPhoto += m.gPhoto;
            gDigestB += m.gDigestB;
            gResistX += m.gResistX;
            gMaint += m.gMaint;
        }
        const float inv = 1.0f / std::max(1.f, static_cast<float>(o.cells.size()));
        gUptakeN *= inv;
        gPhoto *= inv;
        gDigestB *= inv;
        gResistX *= inv;
        gMaint *= inv;

        const Vec3  com = centerOfMass(o);
        const float N   = m_env.nutrient(com);
        const float L   = m_env.light(com);
        const float X   = m_env.toxin(com);

        const Vec3 gN = m_env.gradN(com);
        const Vec3 gL = m_env.gradL(com);
        const Vec3 gX = m_env.gradX(com);
        const Vec3 gB = m_env.gradB(com);

        Vec3 steer =
            normalize(gN) * (0.38f * o.pheno.nutrientAffinity)
            + normalize(gL) * (0.25f * o.pheno.photoAffinity)
            + normalize(gB) * (0.15f * clampf(gDigestB, 0.0f, 2.2f))
            - normalize(gX) * (0.50f * (1.2f - clampf(o.pheno.toxinResistance, 0.f, 1.2f)))
            + hashToVec3(o.id + static_cast<uint64_t>(o.age * 10.0f)) * 0.12f;

        if (len2(steer) < 1e-6f)
            steer = Vec3{m_rng.uniform(-1, 1), m_rng.uniform(-1, 1), m_rng.uniform(-1, 1)};
        steer = normalize(steer);

        o.vel += steer * (0.6f * dt);
        o.vel = clampVec(o.vel, 4.0f / (1.0f + 0.02f * static_cast<float>(o.cells.size())));
        o.pos += o.vel * dt * 6.0f;

        if (const float r = len(o.pos); r > worldRadius)
        {
            const Vec3 inward = normalize(o.pos) * -1.f;
            o.vel += inward * dt * 4.f;
            o.pos = normalize(o.pos) * worldRadius;
        }

        const auto sizeF = static_cast<float>(o.cells.size());

        // nutriments
        const float uptake = o.pheno.nutrientAffinity * N
                             * (0.035f + 0.012f * sizeF)
                             * clampf(gUptakeN, 0.2f, 2.2f);

        // photo
        const float photo = o.pheno.photoAffinity * L
                            * (0.020f + 0.010f * sizeF)
                            * clampf(gPhoto, 0.0f, 2.4f);

        // tox
        const float localResist = clampf(o.pheno.toxinResistance * clampf(gResistX, 0.2f, 2.6f), 0.f, 3.0f);
        const float toxHit      = std::max(0.f, X - 0.7f * localResist)
                                  * (0.040f + 0.005f * sizeF);

        // maint
        const float maint = o.pheno.baseMetabolism
                            * (1.0f + 0.055f * sizeF)
                            * clampf(gMaint, 0.6f, 2.0f);

        const float moveCost = 0.015f * len(o.vel) * (1.0f + 0.03f * static_cast<float>(o.cells.size()));
        const float ageDrain = 0.0008f * o.age;

        o.energy += (uptake + photo - toxHit - moveCost - maint - ageDrain) * dt * 20.f;

        digestBiomass(o, gDigestB);
        grow(o);

        // reproduction
        if (o.energy > o.pheno.splitThreshold && o.age > 6.f && o.cells.size() >= 2)
        {
            Organism child;
            child.id     = m_nextId++;
            child.genome = mutateGenome(o.genome, o.pheno.mutationRate, m_rng);
            child.pheno  = interpretGenome(child.genome);

            child.pos = o.pos + normalize(Vec3{m_rng.uniform(-1, 1), m_rng.uniform(-1, 1),
                                               m_rng.uniform(-1, 1)}) * 4.f;
            child.vel    = o.vel;
            child.energy = o.energy * 0.36f;
            o.energy *= 0.54f;

            // take 1/3 cells
            const int take = std::min<int>(std::max<int>(1, static_cast<int>(o.cells.size()) / 3), 12);
            for (int k = 0; k < take && !o.cells.empty(); ++k)
            {
                const int idx = m_rng.uniformInt(0, static_cast<int>(o.cells.size()) - 1);
                child.cells.push_back(o.cells[static_cast<size_t>(idx)]);
                o.cells.erase(o.cells.begin() + idx);
            }
            if (child.cells.empty())
                child.cells.push_back({Vec3{0, 0, 0}, Vec3{0, 0, 0}, 0.f, 1.f});

            updateCellModes(child);

            // recentre child
            {
                Vec3 cc{0, 0, 0};
                for (const auto &c : child.cells)
                    cc += c.localPos;
                cc = cc / static_cast<float>(child.cells.size());
                for (auto &c : child.cells)
                    c.localPos -= cc;
            }

            // rebuild child springs lightly
            for (int i = 1; i < static_cast<int>(child.cells.size()); ++i)
                connectNewCell(child, i);

            // prune parent invalid springs
            std::erase_if(o.springs,
                          [&](const Spring &s)
                          {
                              return s.a < 0 || s.b < 0 || s.a >= static_cast<int>(o.cells.size()) ||
                                     s.b >= static_cast<int>(o.cells.size());
                          });

            newborns.push_back(std::move(child));
        }

        if (o.energy <= 0.f || o.age > o.pheno.maxAge)
            o.alive = false;
    }

    void Sim::step()
    {
        m_env.step(dt);

        huntInteractions();

        std::vector<Organism> newborns;
        newborns.reserve(m_orgs.size() / 4 + 8);

        for (auto &o : m_orgs)
            updateOne(o, newborns);

        // dead -> biomass
        for (const auto &o : m_orgs)
        {
            if (!o.alive)
            {
                const Vec3  com    = centerOfMass(o);
                const float amount = std::max(0.f, o.energy) + 0.45f * static_cast<float>(o.cells.size());
                m_env.depositBiomass(com, amount, 4.f + 0.08f * static_cast<float>(o.cells.size()));
            }
        }

        std::erase_if(m_orgs, [](const Organism &o)
        {
            return !o.alive;
        });
        for (auto &n : newborns)
            m_orgs.push_back(std::move(n));

        if (static_cast<int>(m_orgs.size()) < targetPopulation / 3)
        {
            const int toAdd = std::max(0, targetPopulation / 2 - static_cast<int>(m_orgs.size()));
            for (int i = 0; i < toAdd; ++i)
            {
                Organism o;
                o.id = m_nextId++;
                o.genome.genes = {{"a"}};
                o.genome = mutateGenome(o.genome, 0.6f, m_rng);
                o.pheno = interpretGenome(o.genome);
                o.pos = {m_rng.uniform(-20, 20), m_rng.uniform(-20, 20), m_rng.uniform(-20, 20)};
                o.vel = {m_rng.uniform(-0.1f, 0.1f), m_rng.uniform(-0.1f, 0.1f), m_rng.uniform(-0.1f, 0.1f)};
                o.energy = m_rng.uniform(6.f, 10.f);
                o.cells.push_back({Vec3{0, 0, 0}, Vec3{0, 0, 0}, 0.f, 1.f});
                m_orgs.push_back(std::move(o));
            }
        }
    }

    float Sim::localDensity(const Organism &o, const int idx)
    {
        const Vec3  p     = o.cells[static_cast<size_t>(idx)].localPos;
        int         count = 0;
        const float r     = 1.8f * o.pheno.cellRadius;
        const float r2    = r * r;

        for (int j = 0; j < static_cast<int>(o.cells.size()); ++j)
        {
            if (j == idx)
                continue;
            if (len2(o.cells[static_cast<size_t>(j)].localPos - p) < r2)
                ++count;
        }
        // normalize roughly into 0..1
        return clampf(static_cast<float>(count) / 10.0f, 0.f, 1.f);
    }

    uint8_t Sim::chooseCellMode(const Organism &o,
                                const Vec3 &    worldPos,
                                const Vec3 &    localPos,
                                const float     cellAge,
                                const int       cellIndex) const
    {
        const float N = m_env.nutrient(worldPos);
        const float L = m_env.light(worldPos);
        const float X = m_env.toxin(worldPos);
        const float B = m_env.biomass(worldPos);

        const Vec3  comL     = localCOM(o);
        const float dist     = len(localPos - comL);
        const float distNorm = clampf(dist / (6.0f + 0.15f * static_cast<float>(o.cells.size())), 0.f, 1.f);

        const float dens    = localDensity(o, cellIndex);
        const float ageNorm = clampf(cellAge / 25.f, 0.f, 1.f);

        const float f[CellMode::F] = {
            clampf(N / 2.5f, 0.f, 1.f),
            clampf(L / 1.8f, 0.f, 1.f),
            clampf(X / 2.2f, 0.f, 1.f),
            clampf(B / 4.0f, 0.f, 1.f),
            distNorm,
            dens,
            ageNorm,
            1.0f,
        };

        // Stable per-cell randomness to break ties (but still heritable-ish)
        const uint64_t base = mix64(o.id ^ static_cast<uint64_t>(cellIndex) * 0x9E3779B97F4A7C15ULL);

        // Softmax over modes
        float scores[Phenotype::K];
        float maxS = -1e9f;

        for (int k = 0; k < Phenotype::K; ++k)
        {
            const auto &m = o.pheno.modes[static_cast<size_t>(k)];

            float s = m.gateBias;
            for (int i = 0; i < CellMode::F; ++i)
                s += m.gateW[static_cast<size_t>(i)] * f[i];

            // small deterministic noise per (cell,k)
            s += hashToRange(base + static_cast<uint64_t>(k) * 131, -0.08f, 0.08f);

            scores[k] = s;
            maxS      = std::max(maxS, s);
        }

        const float sharp  = clampf(o.pheno.gateSharpness, 0.5f, 4.0f);
        float       expSum = 0.f;
        float       probs[Phenotype::K];

        for (int k = 0; k < Phenotype::K; ++k)
        {
            const float e = std::exp((scores[k] - maxS) * sharp);
            probs[k]      = e;
            expSum += e;
        }

        // Deterministic sampling (not purely argmax)
        const float r   = hashToUnit(base + 777);
        float       acc = 0.f;
        for (int k = 0; k < Phenotype::K; ++k)
        {
            acc += probs[k] / expSum;
            if (r <= acc)
                return static_cast<uint8_t>(k);
        }
        return Phenotype::K - 1;
    }

    void Sim::updateCellModes(Organism &o) const
    {
        for (int i = 0; i < static_cast<int>(o.cells.size()); ++i)
        {
            auto &     c  = o.cells[static_cast<size_t>(i)];
            const Vec3 wp = o.pos + c.localPos;
            c.modeId      = chooseCellMode(o, wp, c.localPos, c.age, i);
        }
    }

    Sim::RuntimeStats Sim::computeStats() const
    {
        RuntimeStats s;
        s.organisms = static_cast<int>(m_orgs.size());
        int cells   = 0;
        for (const auto &o : m_orgs)
        {
            cells += static_cast<int>(o.cells.size());
            for (const auto &c : o.cells)
                s.modeCounts[static_cast<size_t>(c.modeId)]++;
        }
        s.totalCells = cells;
        s.avgCells   = s.organisms > 0 ? static_cast<float>(cells) / static_cast<float>(s.organisms) : 0.f;
        return s;
    }
} // namespace evo