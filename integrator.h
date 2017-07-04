#pragma once

#include "structs.h"
#include "fields.h"

struct integrator {
    static const double a2[];
    static const double a3[];
    static const double a4[];
    static const double a5[];
    static const double a6[];
    static const double a7[];

    static const double b4[];
    static const double b5[];

    double *accelx, *accely, *accelz;
    double *dt;
    bool started = false;
};

integrator* makeIntegrator(long N);

void step(
        integrator* method,
        ParticleState* state, FieldStructure* field,
        void (*accelFunc)(ParticleInfo*, FieldStructure*, double*, double*, double*, double*, double*, double*, double*, double*, double*, long N)
    );
