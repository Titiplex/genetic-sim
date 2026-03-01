#pragma once
#include <unordered_map>
#include <vector>
#include "Environment.hpp"
#include "Genetics.hpp"
#include "Math.hpp"
#include "Random.hpp"

namespace evo
{
    struct Cell
    {
        Vec3  localPos{};
        Vec3  localVel{};
        float age  = 0.f;
        float mass = 1.f;

        uint8_t modeId = 0; // 0..K-1
    };

    struct Spring
    {
        int   a         = -1, b = -1;
        float restLen   = 1.f;
        float stiffness = 8.f;
        float damping   = 0.4f;
    };

    struct Organism
    {
        uint64_t  id = 0;
        Genome    genome;
        Phenotype pheno;

        Vec3  pos{};
        Vec3  vel{};
        float age    = 0.f;
        float energy = 8.f;
        bool  alive  = true;

        std::vector<Cell>   cells;
        std::vector<Spring> springs;
    };

    struct IVec3
    {
        int x = 0, y = 0, z = 0;

        friend bool operator==(const IVec3 &a, const IVec3 &b)
        {
            return a.x == b.x && a.y == b.y && a.z == b.z;
        }
    };

    struct IVec3Hash
    {
        size_t operator()(const IVec3 &k) const noexcept
        {
            const uint64_t h = static_cast<uint64_t>(k.x * 73856093) ^ static_cast<uint64_t>(k.y * 19349663) ^
                               static_cast<uint64_t>(k.z * 83492791);
            return h;
        }
    };

    class Sim
    {
        public:
            explicit Sim(uint64_t seed = 1234567ULL);

            void seedInitial(int n);
            void step();

            const Environment &env() const
            {
                return m_env;
            }

            const std::vector<Organism> &organisms() const
            {
                return m_orgs;
            }

            float dt               = 0.035f;
            float worldRadius      = 80.f;
            int   targetPopulation = 180;

        private:
            Environment           m_env;
            std::vector<Organism> m_orgs;
            Rng                   m_rng;

            uint64_t m_nextId = 1;

            std::unordered_map<IVec3, std::vector<int>, IVec3Hash> m_grid;
            float                                                  m_gridSize = 7.f;

            IVec3            cellOf(const Vec3 &p) const;
            void             buildGrid();
            std::vector<int> nearby(const Vec3 &p);

            static Vec3 centerOfMass(const Organism &o);
            static Vec3 localCOM(const Organism &o);
            static void recenterLocalBody(Organism &o);

            Vec3 chooseGrowthDirection(const Organism &o) const;

            void applySprings(Organism &o) const;
            void softRepulsion(Organism &o) const;
            void integrateBody(Organism &o) const;

            static void connectNewCell(Organism &o, int newIdx);
            void        grow(Organism &o) const;

            void huntInteractions();
            void digestBiomass(Organism &o, float gDigestB);

            void updateOne(Organism &o, std::vector<Organism> &newborns);

            uint8_t chooseCellMode(const Organism &o, const Vec3 &worldPos, const Vec3 &localPos,
                                   float           cellAge, int   cellIndex) const;
            void         updateCellModes(Organism &o) const;
            static float localDensity(const Organism &o, int idx);
    };
} // namespace evo