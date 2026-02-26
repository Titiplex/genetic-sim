#pragma once
#include <cstdint>
#include <Genetic.hpp>

#include "maths/Vec3.hpp"

// ======================= Organisms / Cells / Springs ==========================
struct Cell
{
    Vec3  localPos;
    Vec3  localVel;
    float energy = 1.f;
    float age    = 0.f;
    float mass   = 1.f;
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

    Vec3  pos;
    Vec3  vel;
    float age    = 0.f;
    float energy = 8.f;
    bool  alive  = true;

    std::vector<Cell>   cells;
    std::vector<Spring> springs;
};
