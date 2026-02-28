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
            const Vec3  d = c.localPos - lc;
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
            const Vec3  d = lc - c.localPos;
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
            if (a < 0 || b < 0 || a >= static_cast<int>(o.cells.size()) || b >= static_cast<int>(o.cells.size()))
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
                Vec3  d = o.cells[static_cast<size_t>(j)].localPos - o.cells[static_cast<size_t>(i)].localPos;
                const float l = len(d);

                if (l < 1e-5f)
                    continue;
                if (l < minDist)
                {
                    constexpr float k   = 8.0f;
                    const Vec3  n   = d / l;
                    const float pen = minDist - l;
                    const Vec3  f   = n * (k * pen);
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

        std::vector<std::pair<float, int>> ds;
        ds.reserve(static_cast<size_t>(newIdx));
        for (int i = 0; i < newIdx; ++i)
            ds.emplace_back(len2(o.cells[static_cast<size_t>(i)].localPos - o.cells[static_cast<size_t>(newIdx)].localPos), i);

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
            s.a         = newIdx;
            s.b         = j;
            s.restLen   = clampf(std::sqrt(ds[static_cast<size_t>(k)].first), 0.45f, 2.6f);
            s.stiffness = 5.f + 6.0f * o.pheno.adhesionBias;
            s.damping   = 0.35f;
            o.springs.push_back(s);
        }
    }

    void Sim::grow(Organism &o)
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
            o.cells.push_back({bestPos, dir * 0.08f, 0.f, 1.f});
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

                if (const float huntRange = 2.5f + 0.08f * static_cast<float>(a.cells.size()); len2(b.pos - a.pos) > huntRange * huntRange)
                    continue;

                const float aPower = a.pheno.aggression + 0.018f * static_cast<float>(a.cells.size());
                const float bDef   = 0.5f * b.pheno.toxinResistance + 0.010f * static_cast<float>(b.cells.size());

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

    void Sim::digestBiomass(Organism &o)
    {
        const Vec3  com = centerOfMass(o);
        const float B   = m_env.biomass(com);
        if (B <= 0.0001f)
            return;

        const float request = o.pheno.scavenging * (0.02f + 0.005f * static_cast<float>(o.cells.size())) * dt * 20.f;
        const float gain    = std::min(B, request) * 0.9f;

        o.energy += gain * 0.8f;
        m_env.consumeBiomassNear(com, gain);
    }

    void Sim::updateOne(Organism &o, std::vector<Organism> &newborns)
    {
        if (!o.alive)
            return;

        o.age += dt;
        integrateBody(o);

        const Vec3  com = centerOfMass(o);
        const float N   = m_env.nutrient(com);
        const float L   = m_env.light(com);
        const float X   = m_env.toxin(com);

        const Vec3 gN = m_env.gradN(com);
        const Vec3 gL = m_env.gradL(com);
        const Vec3 gX = m_env.gradX(com);
        const Vec3 gB = m_env.gradB(com);

        Vec3 steer = normalize(gN) * (0.38f * o.pheno.nutrientAffinity) + normalize(gL) * (
                         0.25f * o.pheno.photoAffinity) +
                     normalize(gB) * (0.20f * o.pheno.scavenging) -
                     normalize(gX) * (0.50f * (1.2f - clampf(o.pheno.toxinResistance, 0.f, 1.2f))) +
                     hashToVec3(o.id + static_cast<uint64_t>(o.age * 10.0f)) * 0.12f;

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

        const float uptake = o.pheno.nutrientAffinity * N * (0.05f + 0.010f * static_cast<float>(o.cells.size()));
        const float photo  = o.pheno.photoAffinity * L * (0.03f + 0.008f * static_cast<float>(o.cells.size()));
        const float toxHit =
            std::max(0.f, X - 0.7f * o.pheno.toxinResistance) * (0.04f + 0.005f * static_cast<float>(o.cells.size()));

        const float moveCost = 0.015f * len(o.vel) * (1.0f + 0.03f * static_cast<float>(o.cells.size()));
        const float maint    = o.pheno.baseMetabolism * (1.0f + 0.055f * static_cast<float>(o.cells.size()));
        const float ageDrain = 0.0008f * o.age;

        o.energy += (uptake + photo - toxHit - moveCost - maint - ageDrain) * dt * 20.f;

        digestBiomass(o);
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
} // namespace evo