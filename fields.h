#pragma once

#include "structs.h"

struct FieldStructure {
    long N;
    Vector3 *E;
    Vector3 *B;

    virtual void initFields() = 0;
    virtual void getFields(ParticleState* state) = 0;
};
