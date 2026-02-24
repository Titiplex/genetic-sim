#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

// ======================= Math ==========================
struct Vec3
{
    float x = 0, y = 0, z = 0;
    Vec3() = default;

    Vec3(const float X, const float Y, const float Z): x(X), y(Y), z(Z)
    {
    }

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(const float s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(const float s) const { return {x / s, y / s, z / s}; }

    Vec3& operator+=(const Vec3& o)
    {
        x += o.x;
        y += o.y;
        z += o.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& o)
    {
        x -= o.x;
        y -= o.y;
        z -= o.z;
        return *this;
    }

    Vec3& operator*=(const float s)
    {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }
};

static float dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static float len2(const Vec3& v) { return dot(v, v); }
static float len(const Vec3& v) { return std::sqrt(len2(v)); }

static Vec3 normalize(const Vec3& v)
{
    const float l = len(v);
    return l > 1e-6f ? v / l : Vec3{0, 0, 0};
}

static Vec3 clampVec(const Vec3& v, const float m)
{
    if (const float l = len(v); l > m && l > 1e-6f) return v * (m / l);
    return v;
}

static float clampf(const float x, const float a, const float b) { return std::max(a, std::min(b, x)); }
static float lerpf(const float a, const float b, const float t) { return a + (b - a) * t; }

struct Mat4
{
    float m[16]{};

    static Mat4 identity()
    {
        Mat4 r;
        r.m[0] = r.m[5] = r.m[10] = r.m[15] = 1.f;
        return r;
    }
};

static Mat4 perspective(const float fovy, const float aspect, const float znear, const float zfar)
{
    Mat4        r{};
    const float t = std::tan(fovy * 0.5f);
    r.m[0]        = 1.0f / (aspect * t);
    r.m[5]        = 1.0f / t;
    r.m[10]       = -(zfar + znear) / (zfar - znear);
    r.m[11]       = -1.0f;
    r.m[14]       = -(2.0f * zfar * znear) / (zfar - znear);
    return r;
}

static Mat4 lookAt(const Vec3& eye, const Vec3& center, const Vec3& up)
{
    const Vec3 f = normalize(center - eye);
    const Vec3 s = normalize(Vec3{
        f.y * up.z - f.z * up.y,
        f.z * up.x - f.x * up.z,
        f.x * up.y - f.y * up.x
    });
    const auto u = Vec3{
        s.y * f.z - s.z * f.y,
        s.z * f.x - s.x * f.z,
        s.x * f.y - s.y * f.x
    };

    Mat4 r  = Mat4::identity();
    r.m[0]  = s.x;
    r.m[4]  = s.y;
    r.m[8]  = s.z;
    r.m[1]  = u.x;
    r.m[5]  = u.y;
    r.m[9]  = u.z;
    r.m[2]  = -f.x;
    r.m[6]  = -f.y;
    r.m[10] = -f.z;
    r.m[12] = -dot(s, eye);
    r.m[13] = -dot(u, eye);
    r.m[14] = dot(f, eye);
    return r;
}

static Mat4 mul(const Mat4& a, const Mat4& b)
{
    Mat4 r{};
    for (int c = 0; c < 4; c++)
    {
        for (int r0 = 0; r0 < 4; r0++)
        {
            r.m[c * 4 + r0] =
                a.m[0 * 4 + r0] * b.m[c * 4 + 0] +
                a.m[1 * 4 + r0] * b.m[c * 4 + 1] +
                a.m[2 * 4 + r0] * b.m[c * 4 + 2] +
                a.m[3 * 4 + r0] * b.m[c * 4 + 3];
        }
    }
    return r;
}

// ======================= RNG / Hash ==========================
static uint64_t fnv1a64(const std::string& s)
{
    uint64_t h = 1469598103934665603ULL;
    for (const unsigned char c : s)
    {
        h ^= static_cast<uint64_t>(c);
        h *= 1099511628211ULL;
    }
    return h;
}

static uint64_t mix64(uint64_t x)
{
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

static float hashToUnit(uint64_t x)
{
    x = mix64(x);
    return static_cast<float>(static_cast<double>(x >> 11) * (1.0 / 9007199254740992.0)); // [0,1)
}

static float hashToRange(const uint64_t x, const float a, const float b)
{
    return lerpf(a, b, hashToUnit(x));
}

static Vec3 hashToVec3(const uint64_t base)
{
    const Vec3 v{
        hashToRange(base + 11, -1.f, 1.f),
        hashToRange(base + 23, -1.f, 1.f),
        hashToRange(base + 37, -1.f, 1.f)
    };
    return normalize(v);
}

static std::mt19937_64 g_rng(1234567);

static float rndf(const float a = 0.f, const float b = 1.f)
{
    std::uniform_real_distribution<float> d(a, b);
    return d(g_rng);
}

static int rndi(const int a, const int b)
{
    std::uniform_int_distribution<int> d(a, b);
    return d(g_rng);
}

// ======================= Environment (3 factors) ==========================
// 1) Nutrients N(x,y,z)
// 2) Light L(x,y,z)
// 3) Toxin X(x,y,z)
struct Environment
{
    float worldRadius = 80.f;
    float time        = 0.f;

    [[nodiscard]] float nutrient(const Vec3& p) const
    {
        // Two moving nutrient "clouds"
        const Vec3  c1{20.f * std::sin(time * 0.07f), 8.f * std::cos(time * 0.05f), 20.f * std::cos(time * 0.09f)};
        const Vec3  c2{-25.f * std::cos(time * 0.04f), -10.f * std::sin(time * 0.06f), 15.f * std::sin(time * 0.03f)};
        const float d1 = len2(p - c1);
        const float d2 = len2(p - c2);
        const float n  = 1.8f * std::exp(-d1 / (2.f * 18.f * 18.f))
            + 1.4f * std::exp(-d2 / (2.f * 14.f * 14.f));
        return clampf(n, 0.f, 2.5f);
    }

    [[nodiscard]] float light(const Vec3& p) const
    {
        // Light stronger toward +Y and slight day/night oscillation
        const float y   = (p.y + worldRadius) / (2.f * worldRadius); // 0..1
        const float day = 0.65f + 0.35f * std::sin(time * 0.3f);
        return clampf((0.1f + 1.2f * y) * day, 0.f, 1.8f);
    }

    [[nodiscard]] float toxin(const Vec3& p) const
    {
        // Toxic shell regions + a moving hotspot
        const float r     = len(p);
        const float shell = std::exp(-((r - 58.f) * (r - 58.f)) / (2.f * 5.f * 5.f));
        const Vec3  hot{30.f * std::sin(time * 0.08f), 0.f, -30.f * std::cos(time * 0.06f)};
        const float hotspot = 1.6f * std::exp(-len2(p - hot) / (2.f * 10.f * 10.f));
        return clampf(0.8f * shell + hotspot, 0.f, 2.2f);
    }

    [[nodiscard]] Vec3 gradN(const Vec3& p) const
    {
        constexpr float e = 0.8f;
        return Vec3{
                nutrient({p.x + e, p.y, p.z}) - nutrient({p.x - e, p.y, p.z}),
                nutrient({p.x, p.y + e, p.z}) - nutrient({p.x, p.y - e, p.z}),
                nutrient({p.x, p.y, p.z + e}) - nutrient({p.x, p.y, p.z - e})
            }
            /
            2.f * e;
    }

    [[nodiscard]] Vec3 gradL(const Vec3& p) const
    {
        constexpr float e = 0.8f;
        return Vec3{
                light({p.x + e, p.y, p.z}) - light({p.x - e, p.y, p.z}),
                light({p.x, p.y + e, p.z}) - light({p.x, p.y - e, p.z}),
                light({p.x, p.y, p.z + e}) - light({p.x, p.y, p.z - e})
            }
            /
            2.f * e;
    }

    [[nodiscard]] Vec3 gradX(const Vec3& p) const
    {
        constexpr float e = 0.8f;
        return Vec3{
                toxin({p.x + e, p.y, p.z}) - toxin({p.x - e, p.y, p.z}),
                toxin({p.x, p.y + e, p.z}) - toxin({p.x, p.y - e, p.z}),
                toxin({p.x, p.y, p.z + e}) - toxin({p.x, p.y, p.z - e})
            }
            /
            2.f * e;
    }
};

// ======================= Genome -> Dynamic Interpretation ==========================
struct Gene
{
    std::string allele; // e.g. "a", then mutates/grows
};

struct Genome
{
    std::vector<Gene> genes;
};

struct Phenotype
{
    // No direct "gene X => trait Y". Everything derived from hashes over alleles.
    float nutrientAffinity = 1.0f;
    float photoAffinity    = 1.0f;
    float toxinResistance  = 0.2f;
    float aggression       = 0.1f; // ability/interest to hunt
    float growthBias       = 0.5f; // tendency to build structure
    float maxAge           = 120.f;
    float baseMetabolism   = 0.04f;
    float splitThreshold   = 18.f;
    float mutationRate     = 0.06f;
    float cellRadius       = 0.65f;

    // Dynamic directional-growth interpreter weights
    // growth_dir = Σ wi * feature_i(position, local, env, internal) ; not axis-hardcoded
    std::array<float, 12> w{};
    std::array<Vec3, 8>   basis{};
    std::array<float, 8>  basisGain{};

    // Visual
    Vec3 pigment{0.7f, 0.7f, 1.0f};
};

static Phenotype interpretGenome(const Genome& g)
{
    Phenotype p;
    if (g.genes.empty()) return p;

    // Fold all alleles into a master hash stream
    uint64_t master = 0x9E3779B97F4A7C15ULL;
    for (size_t i = 0; i < g.genes.size(); ++i)
    {
        const uint64_t h = fnv1a64(g.genes[i].allele);
        master ^= mix64(h + (uint64_t)i * 0x9E3779B97F4A7C15ULL);
        master = mix64(master);
    }

    auto H = [&](const uint64_t k) { return mix64(master + k * 0x9E3779B97F4A7C15ULL); };

    // Emergent-ish phenotype params from aggregated hash
    p.nutrientAffinity = hashToRange(H(1), 0.1f, 2.2f);
    p.photoAffinity    = hashToRange(H(2), 0.0f, 2.0f);
    p.toxinResistance  = hashToRange(H(3), 0.0f, 1.8f);
    p.aggression       = hashToRange(H(4), 0.0f, 1.4f);
    p.growthBias       = hashToRange(H(5), 0.05f, 1.5f);
    p.maxAge           = hashToRange(H(6), 60.f, 220.f);
    p.baseMetabolism   = hashToRange(H(7), 0.015f, 0.12f);
    p.splitThreshold   = hashToRange(H(8), 10.f, 34.f);
    p.mutationRate     = hashToRange(H(9), 0.01f, 0.18f);
    p.cellRadius       = hashToRange(H(10), 0.4f, 1.0f);

    // Growth-direction interpreter weights and random bases
    for (int i = 0; i < 12; ++i) p.w[i] = hashToRange(H(100 + i), -2.0f, 2.0f);
    for (int i = 0; i < 8; ++i)
    {
        p.basis[i]     = hashToVec3(H(200 + i));
        p.basisGain[i] = hashToRange(H(300 + i), -1.0f, 1.0f);
    }

    p.pigment = {
        hashToRange(H(400), 0.15f, 1.0f),
        hashToRange(H(401), 0.15f, 1.0f),
        hashToRange(H(402), 0.15f, 1.0f)
    };
    // Slightly tie pigment to metabolism mode for readability (not phenotype hardcoding)
    p.pigment.x = clampf(0.35f + 0.65f * (p.aggression / 1.4f), 0.f, 1.f);
    p.pigment.y = clampf(0.25f + 0.75f * (p.nutrientAffinity / 2.2f), 0.f, 1.f);
    p.pigment.z = clampf(0.2f + 0.8f * (p.photoAffinity / 2.0f), 0.f, 1.f);

    return p;
}

// ======================= Organisms / Cells ==========================
struct Cell
{
    Vec3  localPos; // position relative to organism center
    float energy = 1.f;
    float age    = 0.f;
};

struct Organism
{
    uint64_t  id = 0;
    Genome    genome;
    Phenotype pheno;

    Vec3  pos;
    Vec3  vel;
    float age    = 0.f;
    float energy = 8.f;
    bool  alive  = true;

    std::vector<Cell> cells;
};

static uint64_t g_nextId = 1;

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

    // mutate existing alleles
    for (auto& [allele] : g.genes)
    {
        allele = mutateAllele(allele, mutationRate);
    }

    // add/remove gene sometimes (dynamic interpretation scales with gene count)
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

// ======================= Spatial Hash (for hunting & local interactions) ==========================
struct IVec3
{
    int x, y, z;
};

static bool operator==(const IVec3& a, const IVec3& b) { return a.x == b.x && a.y == b.y && a.z == b.z; }

struct IVec3Hash
{
    size_t operator()(const IVec3& k) const noexcept
    {
        const uint64_t h = static_cast<uint64_t>(k.x * 73856093) ^ static_cast<uint64_t>(k.y * 19349663) ^ static_cast<
            uint64_t>(k.z * 83492791);
        return (size_t)h;
    }
};

// ======================= Simulation ==========================
struct Sim
{
    Environment           env;
    std::vector<Organism> organisms;

    float dt               = 0.04f;
    float cellSpacing      = 1.2f;
    int   targetPopulation = 200; // increase later
    float worldRadius      = 80.f;

    std::unordered_map<IVec3, std::vector<int>, IVec3Hash> grid;
    float                                                  gridSize = 6.f;

    void seedInitial(const int n)
    {
        organisms.clear();
        organisms.reserve(n * 2);

        for (int i = 0; i < n; i++)
        {
            Organism o;
            o.id           = g_nextId++;
            o.genome.genes = {{"a"}};
            // add slight diversity immediately
            if (i > 0)
            {
                o.genome = mutateGenome(o.genome, 0.45f);
                if (rndf() < 0.5f) o.genome = mutateGenome(o.genome, 0.2f);
            }
            o.pheno = interpretGenome(o.genome);
            o.pos   = {
                rndf(-25.f, 25.f),
                rndf(-25.f, 25.f),
                rndf(-25.f, 25.f)
            };
            o.vel    = {rndf(-0.2f, 0.2f), rndf(-0.2f, 0.2f), rndf(-0.2f, 0.2f)};
            o.energy = rndf(6.f, 12.f);
            o.cells.push_back({Vec3{0, 0, 0}, 1.f, 0.f}); // starts as single cell
            organisms.push_back(std::move(o));
        }
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
                    IVec3 q{x + dx, y + dy, z + dz};
                    if (auto it = grid.find(q); it != grid.end())
                    {
                        out.insert(out.end(), it->second.begin(), it->second.end());
                    }
                }
        return out;
    }

    static Vec3 organismCenterOfMass(const Organism& o)
    {
        if (o.cells.empty()) return o.pos;
        Vec3 s{0, 0, 0};
        for (auto& c : o.cells) s += o.pos + c.localPos;
        return s / static_cast<float>(o.cells.size());
    }

    Vec3 chooseGrowthDirection(const Organism& o) const
    {
        // Dynamic interpretation from genome-derived weights + environment + local structure.
        const Vec3 p  = o.pos;
        const Vec3 gN = env.gradN(p);
        const Vec3 gL = env.gradL(p);
        const Vec3 gX = env.gradX(p);

        // Structural biases (computed from current geometry, no axis hardcoding)
        Vec3 com = {0, 0, 0};
        for (auto& c : o.cells) com += c.localPos;
        if (!o.cells.empty()) com = com / static_cast<float>(o.cells.size());

        // approximate radial spread / anisotropy direction
        Vec3 principal{0, 0, 0};
        for (auto& c : o.cells) principal += normalize(c.localPos - com) * (0.1f + len(c.localPos - com));
        principal = normalize(principal);

        Vec3 repulse{0, 0, 0};
        for (auto& c : o.cells)
        {
            Vec3 d = com - c.localPos;
            if (const float l = len(d); l > 1e-4f) repulse += d / (l * l + 0.2f);
        }
        repulse = normalize(repulse);

        const Vec3 randomDrift = normalize(Vec3{rndf(-1, 1), rndf(-1, 1), rndf(-1, 1)});

        // Features basis
        const std::array<Vec3, 12> feat = {
            normalize(gN),
            normalize(gL),
            normalize(gX) * -1.f, // often useful to avoid toxin
            normalize(o.vel),
            principal,
            repulse,
            randomDrift,
            normalize(p * -1.f), // keep inside world
            normalize(gN + gL * 0.5f - gX * 0.7f),
            normalize(gN - gL),
            normalize(gL - gX),
            Vec3{0, 0, 0}
        };

        Vec3 dir{0, 0, 0};
        for (int i = 0; i < 12; i++) dir += feat[i] * o.pheno.w[i];
        for (int i = 0; i < 8; i++) dir += o.pheno.basis[i] * o.pheno.basisGain[i] * (0.15f + 0.85f * rndf());

        if (len2(dir) < 1e-6f) dir = randomDrift;
        return normalize(dir);
    }

    void growOrganism(Organism& o) const
    {
        if (o.cells.size() >= 48) return; // MVP cap (raise later)
        if (o.energy < 2.0f) return;

        if (const float desire = o.pheno.growthBias * (0.4f + 0.6f *
            clampf(o.energy / o.pheno.splitThreshold, 0.f, 2.f)); rndf() > clampf(desire * dt * 1.5f, 0.f, 0.25f))
            return;

        const Vec3 dir = chooseGrowthDirection(o);

        // Place new cell near existing surface, trying several candidates
        const Vec3 anchor    = o.cells[rndi(0, static_cast<int>(o.cells.size()) - 1)].localPos;
        Vec3       bestPos   = anchor + dir * (cellSpacing * o.pheno.cellRadius);
        float      bestScore = -1e9f;

        for (int k = 0; k < 8; k++)
        {
            Vec3 candDir = normalize(dir + hashToVec3(static_cast<uint64_t>(k) + o.id) * rndf(0.0f, 0.9f));
            Vec3 cand    = anchor + candDir * (cellSpacing * o.pheno.cellRadius * rndf(0.9f, 1.4f));

            // avoid overlap with existing cells
            float minD = 1e9f;
            float rep  = 0.f;
            for (auto& c : o.cells)
            {
                float d = len(c.localPos - cand);
                minD    = std::min(minD, d);
                if (d < 2.0f) rep += 2.0f - d;
            }
            if (const float score = minD - rep * 0.8f + rndf(-0.2f, 0.2f); score > bestScore)
            {
                bestScore = score;
                bestPos   = cand;
            }
        }

        if (bestScore > 0.2f)
        {
            o.cells.push_back({bestPos, 1.f, 0.f});
            o.energy -= 0.6f + 0.03f * static_cast<float>(o.cells.size());
        }
    }

    void huntInteractions()
    {
        buildGrid();

        for (int i = 0; i < static_cast<int>(organisms.size()); ++i)
        {
            Organism& a = organisms[i];
            if (!a.alive) continue;
            if (a.pheno.aggression < 0.35f) continue;
            if (a.energy < 2.0f) continue;

            auto nearby = nearbyOrganisms(a.pos);
            for (const int j : nearby)
            {
                if (i == j) continue;
                Organism& b = organisms[j];
                if (!b.alive) continue;

                Vec3        d         = b.pos - a.pos;
                const float r2        = len2(d);
                const float huntRange = 2.5f + 0.1f * static_cast<float>(a.cells.size());
                if (r2 > huntRange * huntRange) continue;

                // Simple predation rule: aggression + size advantage helps
                const float aPower = a.pheno.aggression + 0.02f * static_cast<float>(a.cells.size());
                const float bDef   = b.pheno.toxinResistance + 0.01f * static_cast<float>(b.cells.size());

                const float damage = clampf((aPower - 0.4f * bDef) * dt * 4.0f, 0.f, 0.6f);
                if (damage <= 0.f) continue;

                b.energy -= damage;
                a.energy += damage * 0.55f; // imperfect conversion
                // toxin retaliation cost
                a.energy -= env.toxin(b.pos) * (1.0f - clampf(a.pheno.toxinResistance / 2.0f, 0.f, 0.9f)) * dt * 0.5f;
            }
        }
    }

    void updateOne(Organism& o, std::vector<Organism>& newborns) const
    {
        if (!o.alive) return;

        o.age += dt;

        // Environment sensing at center of mass
        Vec3  com = organismCenterOfMass(o);
        float N   = env.nutrient(com);
        float L   = env.light(com);
        float X   = env.toxin(com);

        Vec3 gN = env.gradN(com);
        Vec3 gL = env.gradL(com);
        Vec3 gX = env.gradX(com);

        // Very simple movement controller emergent from phenotype+env (no hardcoded "genes => speed")
        Vec3 steer =
            normalize(gN) * (0.4f * o.pheno.nutrientAffinity)
            + normalize(gL) * (0.25f * o.pheno.photoAffinity)
            - normalize(gX) * (0.5f * (1.2f - clampf(o.pheno.toxinResistance, 0.f, 1.2f)))
            + hashToVec3(o.id + static_cast<uint64_t>(o.age * 10.0f)) * 0.15f;

        o.vel += steer * (0.6f * dt);
        o.vel = clampVec(o.vel, 4.0f / (1.0f + 0.02f * static_cast<float>(o.cells.size())));
        o.pos += o.vel * dt * 6.0f;

        // Soft world confinement
        if (float r = len(o.pos); r > worldRadius)
        {
            Vec3 inward = normalize(o.pos) * -1.f;
            o.vel += inward * dt * 4.f;
            o.pos = normalize(o.pos) * worldRadius;
        }

        // Energy model with 3 environmental factors
        float uptake = o.pheno.nutrientAffinity * N * (0.05f + 0.01f * static_cast<float>(o.cells.size()));
        float photo  = o.pheno.photoAffinity * L * (0.03f + 0.008f * static_cast<float>(o.cells.size()));
        float toxHit = std::max(0.f, X - 0.7f * o.pheno.toxinResistance) * (0.04f + 0.005f * static_cast<float>(o.cells.
                                                                                                                  size()));

        float moveCost = 0.015f * len(o.vel) * (1.0f + 0.03f * static_cast<float>(o.cells.size()));
        float maint    = o.pheno.baseMetabolism * (1.0f + 0.05f * static_cast<float>(o.cells.size()));
        float ageDrain = 0.0008f * o.age;

        o.energy += (uptake + photo - toxHit - moveCost - maint - ageDrain) * dt * 20.f;

        // Cell aging (for visual / later structural death)
        for (auto& c : o.cells) c.age += dt;

        // Growth (shape evolution in 3D, dynamic direction)
        growOrganism(o);

        // Reproduction (asexual split / budding)
        if (o.energy > o.pheno.splitThreshold && o.age > 6.f && o.cells.size() >= 2)
        {
            Organism child;
            child.id     = g_nextId++;
            child.genome = mutateGenome(o.genome, o.pheno.mutationRate);
            child.pheno  = interpretGenome(child.genome);

            child.pos    = o.pos + normalize(Vec3{rndf(-1, 1), rndf(-1, 1), rndf(-1, 1)}) * rndf(2.f, 5.f);
            child.vel    = o.vel + Vec3{rndf(-0.4f, 0.4f), rndf(-0.4f, 0.4f), rndf(-0.4f, 0.4f)};
            child.energy = o.energy * 0.38f;
            o.energy *= 0.54f;

            // Split off some cells (if many), otherwise child starts single-cell
            if (int take = std::min<int>(static_cast<int>(o.cells.size()) / 3, 10); take <= 0)
            {
                child.cells.push_back({Vec3{0, 0, 0}, 1.f, 0.f});
            }
            else
            {
                for (int k = 0; k < take; k++)
                {
                    int  idx   = rndi(0, static_cast<int>(o.cells.size()) - 1);
                    Cell c     = o.cells[idx];
                    c.localPos = c.localPos * 0.6f;
                    c.age      = 0.f;
                    child.cells.push_back(c);
                    o.cells.erase(o.cells.begin() + idx);
                }
                // Recenter child local coordinates
                Vec3 cc{0, 0, 0};
                for (auto& c : child.cells) cc += c.localPos;
                cc = cc / static_cast<float>(child.cells.size());
                for (auto& c : child.cells) c.localPos -= cc;
            }

            newborns.push_back(std::move(child));
        }

        // Death: starvation / old age
        if (o.energy <= 0.f || o.age > o.pheno.maxAge)
        {
            o.alive = false;
        }
    }

    void step()
    {
        env.time += dt;

        // predation first (local interactions)
        huntInteractions();

        std::vector<Organism> newborns;
        newborns.reserve(organisms.size() / 4 + 8);

        for (auto& o : organisms)
        {
            updateOne(o, newborns);
        }

        // Remove dead, optionally turn dead into nutrient pulses (simple feedback)
        for (const auto& o : organisms)
        {
            if (!o.alive)
            {
                // no explicit field grid; dead simply disappear in MVP
                // later: deposit biomass field / debris entities here
            }
        }
        organisms.erase(
            std::remove_if(organisms.begin(), organisms.end(), [](const Organism& o) { return !o.alive; }),
            organisms.end()
        );

        // Spawn newborns
        for (auto& n : newborns) organisms.push_back(std::move(n));

        // Population floor (prevents extinction while testing)
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
                o.cells.push_back({Vec3{0, 0, 0}, 1.f, 0.f});
                organisms.push_back(std::move(o));
            }
        }
    }
};

// ======================= OpenGL helpers ==========================
static GLuint compileShader(const GLenum type, const char* src)
{
    const GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok = 0;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok)
    {
        char log[2048];
        glGetShaderInfoLog(s, sizeof(log), nullptr, log);
        std::cerr << "Shader compile error:\n" << log << "\n";
    }
    return s;
}

static GLuint makeProgram(const char* vs, const char* fs)
{
    const GLuint v = compileShader(GL_VERTEX_SHADER, vs);
    const GLuint f = compileShader(GL_FRAGMENT_SHADER, fs);
    const GLuint p = glCreateProgram();
    glAttachShader(p, v);
    glAttachShader(p, f);
    glLinkProgram(p);
    GLint ok = 0;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok)
    {
        char log[2048];
        glGetProgramInfoLog(p, sizeof(log), nullptr, log);
        std::cerr << "Program link error:\n" << log << "\n";
    }
    glDeleteShader(v);
    glDeleteShader(f);
    return p;
}

// ======================= Render Buffers ==========================
struct RenderPoint
{
    float x, y, z;
    float r, g, b;
    float size;
};

int main()
{
    if (!glfwInit())
    {
        std::cerr << "Failed to init GLFW\n";
        return 1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(1400, 900, "EvoLife3D - Dynamic Genome Interpretation", nullptr, nullptr);
    if (!window)
    {
        std::cerr << "Failed to create window\n";
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cerr << "Failed to init GLAD\n";
        return 1;
    }

    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_DEPTH_TEST);

    auto VS = R"(
        #version 330 core
        layout(location=0) in vec3 aPos;
        layout(location=1) in vec3 aColor;
        layout(location=2) in float aSize;

        uniform mat4 uVP;
        out vec3 vColor;

        void main(){
            gl_Position = uVP * vec4(aPos, 1.0);
            gl_PointSize = aSize;
            vColor = aColor;
        }
    )";

    auto FS = R"(
        #version 330 core
        in vec3 vColor;
        out vec4 FragColor;

        void main(){
            vec2 p = gl_PointCoord * 2.0 - 1.0;
            float r2 = dot(p,p);
            if(r2 > 1.0) discard;
            float edge = smoothstep(1.0, 0.5, r2);
            FragColor = vec4(vColor * edge, 1.0);
        }
    )";

    GLuint program = makeProgram(VS, FS);
    GLint  uVP     = glGetUniformLocation(program, "uVP");

    GLuint vao = 0, vbo = 0;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(RenderPoint) * 1000, nullptr, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(RenderPoint), static_cast<void*>(nullptr));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(RenderPoint), reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(RenderPoint), reinterpret_cast<void*>(6 * sizeof(float)));

    Sim sim;
    sim.seedInitial(260);

    double       lastTime = glfwGetTime();
    double       accum    = 0.0;
    const double fixedDt  = sim.dt;

    bool  paused = false;
    float camYaw = 0.8f, camPitch = 0.5f, camDist = 180.f;

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window, 1);
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) paused = true;
        if (glfwGetKey(window, GLFW_KEY_ENTER) == GLFW_PRESS) paused = false;

        // camera controls
        if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) camYaw -= 0.02f;
        if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) camYaw += 0.02f;
        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) camPitch = clampf(camPitch + 0.02f, -1.4f, 1.4f);
        if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) camPitch = clampf(camPitch - 0.02f, -1.4f, 1.4f);
        if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) camDist = std::max(40.f, camDist - 2.f);
        if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) camDist = std::min(400.f, camDist + 2.f);

        double now   = glfwGetTime();
        double frame = now - lastTime;
        lastTime     = now;
        accum += frame;

        if (!paused)
        {
            int maxSubSteps = 4;
            int steps       = 0;
            while (accum >= fixedDt && steps < maxSubSteps)
            {
                sim.step();
                accum -= fixedDt;
                steps++;
            }
        }

        // Build render point list (cells + tiny environment probes)
        std::vector<RenderPoint> pts;
        size_t                   estimated = 0;
        for (auto& o : sim.organisms) estimated += o.cells.size();
        pts.reserve(estimated + 512);

        // Organisms / cells
        for (const auto& o : sim.organisms)
        {
            Vec3  baseColor = o.pheno.pigment;
            float ageT      = clampf(o.age / o.pheno.maxAge, 0.f, 1.f);

            for (const auto& c : o.cells)
            {
                Vec3 wp = o.pos + c.localPos;

                // Slight color variation by local environment (visualization aid)
                float n = sim.env.nutrient(wp);
                float l = sim.env.light(wp);
                float x = sim.env.toxin(wp);

                Vec3 col = {
                    clampf(baseColor.x + 0.2f * n + 0.2f * o.pheno.aggression - 0.15f * x, 0.f, 1.f),
                    clampf(baseColor.y + 0.25f * n - 0.05f * ageT, 0.f, 1.f),
                    clampf(baseColor.z + 0.25f * l - 0.12f * x, 0.f, 1.f)
                };

                float size = 4.0f + 2.0f * o.pheno.cellRadius;
                pts.push_back({wp.x, wp.y, wp.z, col.x, col.y, col.z, size});
            }

            // Draw COM marker
            Vec3 com = Sim::organismCenterOfMass(o);
            pts.push_back({com.x, com.y, com.z, 1.f, 1.f, 1.f, 2.5f});
        }

        // Sparse environment visualization points (for 3 factors)
        for (int i = 0; i < 220; i++)
        {
            Vec3 p{
                rndf(-sim.worldRadius, sim.worldRadius),
                rndf(-sim.worldRadius, sim.worldRadius),
                rndf(-sim.worldRadius, sim.worldRadius)
            };
            if (len(p) > sim.worldRadius) continue;
            float n = sim.env.nutrient(p);
            float l = sim.env.light(p);
            float x = sim.env.toxin(p);
            Vec3  col{
                clampf(0.2f + 0.55f * x, 0.f, 1.f), // red = toxin
                clampf(0.15f + 0.6f * n, 0.f, 1.f), // green = nutrient
                clampf(0.15f + 0.6f * l, 0.f, 1.f)  // blue = light
            };
            float size = 2.0f;
            pts.push_back({p.x, p.y, p.z, col.x, col.y, col.z, size});
        }

        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(pts.size() * sizeof(RenderPoint)), pts.data(),
                     GL_DYNAMIC_DRAW);

        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.02f, 0.02f, 0.03f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        Vec3 eye{
            camDist * std::cos(camPitch) * std::cos(camYaw),
            camDist * std::sin(camPitch),
            camDist * std::cos(camPitch) * std::sin(camYaw)
        };
        Mat4 P = perspective(60.f * 3.1415926f / 180.f, static_cast<float>(w) / static_cast<float>(std::max(h, 1)),
                             0.1f, 1000.f);
        Mat4 V   = lookAt(eye, Vec3{0, 0, 0}, Vec3{0, 1, 0});
        auto [m] = mul(P, V);

        glUseProgram(program);
        glUniformMatrix4fv(uVP, 1, GL_FALSE, m);
        glBindVertexArray(vao);
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(pts.size()));

        glfwSwapBuffers(window);

        // Console stats every ~1 sec
        static double statTimer = 0.0;
        statTimer += frame;
        if (statTimer > 1.0)
        {
            statTimer         = 0.0;
            size_t totalCells = 0;
            for (auto& o : sim.organisms) totalCells += o.cells.size();
            std::cout << "Organisms: " << sim.organisms.size()
                << " | Cells: " << totalCells
                << " | Time: " << sim.env.time
                << (paused ? " [PAUSED]" : "")
                << "\n";
        }
    }

    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
    glDeleteProgram(program);

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
