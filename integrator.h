#pragma once

#include "structs.h"
#include "fields.h"

struct Integrator {
    virtual void init(long N) = 0;
    virtual void deinit() = 0;
    virtual void setInitTimestep(double dt) = 0;
    virtual void step(ParticleState* state, FieldStructure* field) = 0;
};

struct EulerIntegrator : Integrator {
    void (*accelFunc)(ParticleState*, FieldStructure*);

    double dt;

    void init(long N) override;
    void deinit() override;
    void setInitTimestep(double dt) override;
    void step(ParticleState* state, FieldStructure* field) override;
};

struct RK4Integrator : Integrator {
    void (*accelFunc)(ParticleState*, FieldStructure*);

    double dt;

    ParticleState *state1, *state2, *state3;

    RK4Integrator(ParticleState* state);

    void init(long N) override;
    void deinit() override;
    void setInitTimestep(double dt) override;
    void step(ParticleState* state, FieldStructure* field) override;
};

struct RKDPIntegrator : Integrator {
    static const double a2[];
    static const double a3[];
    static const double a4[];
    static const double a5[];
    static const double a6[];
    static const double a7[];

    static const double b4[];
    static const double b5[];

    void (*accelFunc)(ParticleState*, FieldStructure*);

    double *dt;

    long N;

    void init(long N) override;
    void deinit() override;
    void setInitTimestep(double dt) override;
    void step(ParticleState* state, FieldStructure* field) override;
};
