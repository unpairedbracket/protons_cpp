#include "source.h"

#include <cassert>
#include <cstdio>

#include "../util/math.h"

ParticleState* SquareSource::genParticleState(ParticleInfo* particleType, ParticleState* existing) {
    ParticleState* state;
    if(existing) {
        state = existing;
    } else {
        state = new ParticleState();
        initParticleState(state, particleType, this->x_extent * this->y_extent);
    }

    initPosPointSource(state, this->distance);

    assert(this->divergence < pi() / 4);

    double speed = sqrt(2 * this->energy / particleType->mass);
    for(long i = 0; i < state->N; i++) {
        long j = i / this->y_extent;
        long k = i % this->y_extent;
        double angleX = (2*j / ((double)this->x_extent-1) - 1) * this->divergence;
        double angleY = (2*k / ((double)this->y_extent-1) - 1) * this->divergence;

        state->vel[i].x = speed * sin(angleX);
        state->vel[i].y = speed * sin(angleY);
        state->vel[i].z = speed * sqrt(cos(angleX + angleY) * cos(angleX - angleY));
    }

    return state;
}

ParticleState* HelixSource::genParticleState(ParticleInfo* particleType, ParticleState* existing) {
    ParticleState* state;
    if(existing) {
        state = existing;
    } else {
        state = new ParticleState();
        initParticleState(state, particleType, this->N);
    }

    initPosPointSource(state, this->distance);

    double z_min = cos(this->divergence);
    double dz = (1 - z_min) / (this->N - 1);

    double speed = sqrt(2 * this->energy / particleType->mass);
    for(long i = 0; i < state->N; i++) {
        double z = z_min + i * dz;
        double phi = i * this->dphi;
        double rho = sqrt(1 - z*z);

        state->vel[i].x = speed * rho * cos(phi);
        state->vel[i].y = speed * rho * sin(phi);
        state->vel[i].z = speed * z;
    }

    return state;
}

ParticleState* ScatterSource::genParticleState(ParticleInfo* particleType, ParticleState* existing) {
    ParticleState* state;
    if(existing) {
        state = existing;
    } else {
        state = new ParticleState();
        initParticleState(state, particleType, this->N);
    }

    initPosPointSource(state, this->distance);

    double z_min = cos(this->divergence);

    double speed = sqrt(2 * this->energy / particleType->mass);
    for(long i = 0; i < state->N; i++) {
        double z = zrand(re);
        double phi = phirand(re);
        double rho = sqrt(1 - z*z);

        state->vel[i].x = speed * rho * cos(phi);
        state->vel[i].y = speed * rho * sin(phi);
        state->vel[i].z = speed * z;
    }

    return state;
}

void initPosPointSource(ParticleState* particles, double distance) {
    for(long i = 0; i < particles->N; i++) {
        particles->pos[i].x = 0;
        particles->pos[i].y = 0;
        particles->pos[i].z = -distance;
    }
}
