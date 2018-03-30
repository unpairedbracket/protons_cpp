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

ParticleState* RectangleSource::createParticleState(ParticleInfo particleType) {
    ParticleState* state = new ParticleState();
    initParticleState(state, particleType, this->x_points * this->y_points);
    return state;
}

void RectangleSource::setParticleState(ParticleState* state) {
    initPosPointSource(state, this->distance);

    double v = speed(this->energy, state->particleInfo.mass);
    double Z = this->distance;
    for(long i = 0; i < state->N; i++) {
        long j = i / this->y_points;
        long k = i % this->y_points;
        double X = (j / ((double)this->x_points-1) - 0.5) * this->x_size;
        double Y = (k / ((double)this->y_points-1) - 0.5) * this->y_size;

        double length = sqrt(X*X + Y*Y + Z*Z);

        state->vel[i].x = v * X / length;
        state->vel[i].y = v * Y / length;
        state->vel[i].z = v * Z / length;
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

