#pragma once

#include "structs.h"

struct FieldStructure {
    long N;
    double *Ex, *Ey, *Ez;
    double *Bx, *By, *Bz;

    virtual void initFields() = 0;
    virtual void getFields(ParticleState* state) = 0;
};
