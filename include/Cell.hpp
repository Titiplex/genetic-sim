#pragma once
#include <cstdint>
#include <Genetic.hpp>

#include "maths/Vec3.hpp"

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
    float gaitPhase = 0.f;
    float age    = 0.f;
    float energy = 8.f;
    float storedEnergy = 0.f;
    float energyDebt = 0.f;
    float oxidativeStress = 0.f;
    float pathogenLoad = 0.f;
    float immuneMemory = 0.f;
    float neuralMemory = 0.f;
    float lastReward = 0.f;
    float reproductionCooldown = 0.f;
    float microbiomeHealth = 0.6f;
    float heatStressMemory = 0.f;
    float hydrationStressMemory = 0.f;
    float matePreference = 0.f;
    float socialDrive = 0.f;
    float phenotypeUpdateTimer = 0.f;
    float offspringInvestmentReserve = 0.f;
    float oxygenStressMemory = 0.f;
    float salinityStressMemory = 0.f;
    float terrestrialAffinity = 0.35f;
    float macroScale = 1.f;
    float dormancy = 0.f;
    float acidStressMemory = 0.f;
    float nicheFidelity = 0.f;
    float starvationTimer = 0.f;
    float epigeneticStress = 0.f;
    uint64_t geneticSignature = 0;
    int   deme = 0;
    LifeStage stage = LifeStage::Juvenile;
    bool  alive  = true;
    bool  macroMode = false;

    std::vector<Cell>   cells;
    std::vector<Spring> springs;
};
