#pragma once

#include "../particles/structs.h"

#include <random>

struct ParticleSource {
    double distance;
    double divergence;

    virtual ParticleState* genParticleState(ParticleInfo* particleType, ParticleState* existing = nullptr) = 0;
};

struct SquareSource : ParticleSource {
    double energy;

    long x_extent, y_extent;

    ParticleState* genParticleState(ParticleInfo* particleType, ParticleState* existing = nullptr) override;
};

struct HelixSource : ParticleSource {
    double energy;
    double dphi;

    long N;

    ParticleState* genParticleState(ParticleInfo* particleType, ParticleState* existing = nullptr) override;
};

struct ScatterSource : ParticleSource {
    std::uniform_real_distribution<double> zrand;
    std::uniform_real_distribution<double> phirand;
    std::default_random_engine re;

    double energy;

    long N;

    ParticleState* genParticleState(ParticleInfo* particleType, ParticleState* existing = nullptr) override;
};

void initPosPointSource(ParticleState* particles, double distance);
