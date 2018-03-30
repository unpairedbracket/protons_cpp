#pragma once

#include "../particles/structs.h"

#include <random>

struct ParticleSource {
    double distance;
    double divergence;
    double energy;

    virtual ParticleState* createParticleState(ParticleInfo particleType) = 0;
    virtual void setParticleState(ParticleState* state) = 0;
};

struct RectangleSource : ParticleSource {
    long x_points, y_points;
    double x_size, y_size;

    ParticleState* createParticleState(ParticleInfo particleType) override;
    void setParticleState(ParticleState* state) override;
};

struct HelixSource : ParticleSource {
    double dphi;

    long N;

    ParticleState* createParticleState(ParticleInfo particleType) override;
    void setParticleState(ParticleState* state) override;
};

struct ScatterSource : ParticleSource {
    std::uniform_real_distribution<double> zrand;
    std::uniform_real_distribution<double> phirand;
    std::default_random_engine re;

    long N;

    ParticleState* createParticleState(ParticleInfo particleType) override;
    void setParticleState(ParticleState* state) override;
};

void initPosPointSource(ParticleState* particles, double distance);
