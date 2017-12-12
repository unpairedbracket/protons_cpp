#pragma once

#include "../particles/structs.h"

struct FieldStructure {
    long N;
    Vector3 *E;
    Vector3 *B;

    Vector3 xaxis, yaxis, zaxis;
    double theta, phi;
    double min_z;

    virtual void initFields() = 0;
    virtual void getFields(ParticleState* state) = 0;
    virtual void invalidatePositions(ParticleState* state) = 0;

    void initFieldArrays(long N);
    void orientBeam(ParticleState* state);
    void deorientBeam(ParticleState* state);
};

