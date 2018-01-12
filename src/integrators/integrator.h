#pragma once

#include <cmath>

#include "../particles/structs.h"
#include "../fields/fields.h"
#include "../util/accel.h"

struct Integrator {
    virtual void deinit() = 0;
    virtual void setInitTimestep(double dt) = 0;
    virtual void setRelativistic(bool isRelativistic) = 0;
    virtual void initFromState(ParticleState* state) = 0;
    virtual void step(ParticleState* state, FieldStructure* field) = 0;
};

struct EulerIntegrator : Integrator {
    void (*accelFunc)(ParticleState*, FieldStructure*);

    double dt;

    void deinit() override;
    void setInitTimestep(double dt) override;
    void setRelativistic(bool isRelativistic) override;
    void initFromState(ParticleState* state) override;
    void step(ParticleState* state, FieldStructure* field) override;
};

struct RK4Integrator : Integrator {
    void (*accelFunc)(ParticleState*, FieldStructure*);

    double dt;

    ParticleState *state1, *state2, *state3;

    void deinit() override;
    void setInitTimestep(double dt) override;
    void setRelativistic(bool isRelativistic) override;
    void initFromState(ParticleState* state) override;
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

    //// Adaptation params
    // How large a proportion of the velocity are we allowed to be wrong by?
    double rtol = 1e-9;
    // err ~ O(h^5) so we take our step size down by (err/tol)^(1/5) if we're wrong
    double error_scale_power = 1.0/5.0;
    // How much are we allowed to increase step size by per non-failed step?
    double maxLengthen = 1.1;
    // How much are we allowed to decrease step size by on first failure?
    double maxFirstShorten = 10;
    // How much do we decrease step size by on subsequent failures?
    double maxOtherShorten = 2;
    // How short is a step allowed to be
    double dt_min = 5E-18; // seconds
    // Do we want to be loud about what's happening
    bool verbose = false;

    void (*accelFunc)(ParticleState*, FieldStructure*);

    double dt_init;
    double *dt;

    long N;

    ParticleState *state1, *state2, *state3, *state4, *state5, *state6, *state7;

    bool* lastIterationSuccess;

    void deinit() override;
    void setInitTimestep(double dt) override;
    void setRelativistic(bool isRelativistic) override;
    void initFromState(ParticleState* state) override;
    void step(ParticleState* state, FieldStructure* field) override;
};