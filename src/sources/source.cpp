#include "source.h"

#include <cassert>
#include <cstdio>

#include "../util/math.h"
#include "../util/physical_constants.h"

// Speed as a function of relativistic kinetic energy
inline double speed(double energy, double mass) {
    double gamma = (1+energy/(mass*c*c));
    double beta = sqrt(1 - 1.0/(gamma*gamma));
    return beta*c;
}

ParticleState* SquareSource::createParticleState(ParticleInfo particleType) {
    ParticleState* state = new ParticleState();
    initParticleState(state, particleType, this->x_extent * this->y_extent);
    return state;
}

void SquareSource::setParticleState(ParticleState* state) {
    initPosPointSource(state, this->distance);

    assert(this->divergence < pi() / 4);

    double v = speed(this->energy, state->particleInfo.mass);
    for(long i = 0; i < state->N; i++) {
        long j = i / this->y_extent;
        long k = i % this->y_extent;
        double angleX = (2*j / ((double)this->x_extent-1) - 1) * this->divergence;
        double angleY = (2*k / ((double)this->y_extent-1) - 1) * this->divergence;

        state->vel[i].x = v * sin(angleX);
        state->vel[i].y = v * sin(angleY);
        state->vel[i].z = v * sqrt(cos(angleX + angleY) * cos(angleX - angleY));
    }
}

ParticleState* HelixSource::createParticleState(ParticleInfo particleType) {
    ParticleState* state = new ParticleState();
    initParticleState(state, particleType, this->N);
    return state;
}

void HelixSource::setParticleState(ParticleState* state) {
    initPosPointSource(state, this->distance);

    double z_min = cos(this->divergence);
    double dz = (1 - z_min) / (this->N - 1);

    double v = speed(this->energy, state->particleInfo.mass);
    for(long i = 0; i < state->N; i++) {
        double z = z_min + i * dz;
        double phi = i * this->dphi;
        double rho = sqrt(1 - z*z);

        state->vel[i].x = v * rho * cos(phi);
        state->vel[i].y = v * rho * sin(phi);
        state->vel[i].z = v * z;
    }
}

ParticleState* ScatterSource::createParticleState(ParticleInfo particleType) {
    ParticleState* state = new ParticleState();
    initParticleState(state, particleType, this->N);
    return state;
}

void ScatterSource::setParticleState(ParticleState* state) {
    initPosPointSource(state, this->distance);

    double v = speed(this->energy, state->particleInfo.mass);
    for(long i = 0; i < state->N; i++) {
        double z = zrand(re);
        double phi = phirand(re);
        double rho = sqrt(1 - z*z);

        state->vel[i].x = v * rho * cos(phi);
        state->vel[i].y = v * rho * sin(phi);
        state->vel[i].z = v * z;
    }
}

void initPosPointSource(ParticleState* particles, double distance) {
    for(long i = 0; i < particles->N; i++) {
        particles->pos[i].x = 0;
        particles->pos[i].y = 0;
        particles->pos[i].z = -distance;
    }
}

