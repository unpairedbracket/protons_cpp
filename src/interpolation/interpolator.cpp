#include "interpolator.h"

void Interpolator::initState(ParticleInfo type) {
    if(!this->interpParticles) this->interpParticles = this->interpSource->createParticleState(type);
    this->interpSource->setParticleState(this->interpParticles);
}

