#include "integrator.h"

#include <cstdio>

const double RKDPIntegrator::a2[] = {1.0/5.0};
const double RKDPIntegrator::a3[] = {3.0/40.0, 9.0/40.0};
const double RKDPIntegrator::a4[] = {44.0/45.0, -56.0/15.0, 32.0/9.0};
const double RKDPIntegrator::a5[] = {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0};
const double RKDPIntegrator::a6[] = {9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0};
const double RKDPIntegrator::a7[] = {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0};

const double RKDPIntegrator::b4[] = {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0};
const double RKDPIntegrator::b5[] = {5179.0/57600.0, 0.0, 7571.0/16695.0, 393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 1.0/40.0};

void EulerIntegrator::initFromState(ParticleState* state) {
}

void EulerIntegrator::deinit() {
}

void EulerIntegrator::setInitTimestep(double dt) {
    this->dt = dt;
}

void EulerIntegrator::setRelativistic(bool relativistic) {
    this->accelFunc = relativistic ? accel_relativistic : accel;
}

void EulerIntegrator::step(ParticleState* state, FieldStructure* field) {

    accelFunc(state, field);
    
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state->pos[j].x += state->vel[j].x * dt;
            state->pos[j].y += state->vel[j].y * dt;
            state->pos[j].z += state->vel[j].z * dt;

            state->vel[j].x += state->acc[j].x * dt;
            state->vel[j].y += state->acc[j].y * dt;
            state->vel[j].z += state->acc[j].z * dt;
        }
    }
}

void RK4Integrator::initFromState(ParticleState* state) {
    state1 = new ParticleState();
    state2 = new ParticleState();
    state3 = new ParticleState();

    shadowParticleState(state1, state);
    shadowParticleState(state2, state);
    shadowParticleState(state3, state);
}

void RK4Integrator::deinit() {
}

void RK4Integrator::setInitTimestep(double dt) {
    this->dt = dt;
}

void RK4Integrator::setRelativistic(bool relativistic) {
    this->accelFunc = relativistic ? accel_relativistic : accel;
}

void RK4Integrator::step(ParticleState* state, FieldStructure* field) {

    accelFunc(state, field);
    
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state1->pos[j].x = state->pos[j].x + 1/2 * state->vel[j].x * dt;
            state1->pos[j].y = state->pos[j].y + 1/2 * state->vel[j].y * dt;
            state1->pos[j].z = state->pos[j].z + 1/2 * state->vel[j].z * dt;

            state1->vel[j].x = state->vel[j].x + 1/2 * state->acc[j].x * dt;
            state1->vel[j].y = state->vel[j].y + 1/2 * state->acc[j].y * dt;
            state1->vel[j].z = state->vel[j].z + 1/2 * state->acc[j].z * dt;
        }
    }

    accelFunc(state1, field);
    
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state2->pos[j].x = state->pos[j].x + 1/2 * state1->vel[j].x * dt;
            state2->pos[j].y = state->pos[j].y + 1/2 * state1->vel[j].y * dt;
            state2->pos[j].z = state->pos[j].z + 1/2 * state1->vel[j].z * dt;

            state2->vel[j].x = state->vel[j].x + 1/2 * state1->acc[j].x * dt;
            state2->vel[j].y = state->vel[j].y + 1/2 * state1->acc[j].y * dt;
            state2->vel[j].z = state->vel[j].z + 1/2 * state1->acc[j].z * dt;
        }
    }

    accelFunc(state2, field);
    
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state3->pos[j].x = state->pos[j].x + state2->vel[j].x * dt;
            state3->pos[j].y = state->pos[j].y + state2->vel[j].y * dt;
            state3->pos[j].z = state->pos[j].z + state2->vel[j].z * dt;

            state3->vel[j].x = state->vel[j].x + state2->acc[j].x * dt;
            state3->vel[j].y = state->vel[j].y + state2->acc[j].y * dt;
            state3->vel[j].z = state->vel[j].z + state2->acc[j].z * dt;
        }
    }

    accelFunc(state3, field);
    
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state->pos[j].x += (state->vel[j].x + 2*state1->vel[j].x + 2*state2->vel[j].x + state3->vel[j].x) * dt / 6;
            state->pos[j].y += (state->vel[j].y + 2*state1->vel[j].y + 2*state2->vel[j].y + state3->vel[j].y) * dt / 6;
            state->pos[j].z += (state->vel[j].z + 2*state1->vel[j].z + 2*state2->vel[j].z + state3->vel[j].z) * dt / 6;

            state->vel[j].x += (state->acc[j].x + 2*state1->acc[j].x + 2*state2->acc[j].x + state3->acc[j].x) * dt / 6;
            state->vel[j].y += (state->acc[j].y + 2*state1->acc[j].y + 2*state2->acc[j].y + state3->acc[j].y) * dt / 6;
            state->vel[j].z += (state->acc[j].z + 2*state1->acc[j].z + 2*state2->acc[j].z + state3->acc[j].z) * dt / 6;
        }
    }

}

void RKDPIntegrator::initFromState(ParticleState* state) {
    this->lastIterationSuccess = new bool[state->N];
    this->dt = new double[state->N];
    this->N = state->N;

    state1 = new ParticleState();
    state2 = new ParticleState();
    state3 = new ParticleState();
    state4 = new ParticleState();
    state5 = new ParticleState();
    state6 = new ParticleState();
    state7 = new ParticleState();
    
    shadowParticleState(state1, state);
    shadowParticleState(state2, state);
    shadowParticleState(state3, state);
    shadowParticleState(state4, state);
    shadowParticleState(state5, state);
    shadowParticleState(state6, state);
    shadowParticleState(state7, state);

}

void RKDPIntegrator::reset() {
    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        this->dt[j] = dt_init;
    }
}

void RKDPIntegrator::deinit() {
    delete[] this->lastIterationSuccess;
    delete[] this->dt;
    //TODO: delete everything contained within
}

void RKDPIntegrator::setInitTimestep(double dt_0) {
    dt_init = dt_0;
}

void RKDPIntegrator::setRelativistic(bool relativistic) {
    this->accelFunc = relativistic ? accel_relativistic : accel;
}

void RKDPIntegrator::step(ParticleState* state, FieldStructure* field) {

    accelFunc(state, field);
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state1->pos[j].x = state->pos[j].x;
            state1->pos[j].y = state->pos[j].y;
            state1->pos[j].z = state->pos[j].z;

            state1->vel[j].x = state->vel[j].x;
            state1->vel[j].y = state->vel[j].y;
            state1->vel[j].z = state->vel[j].z;

            state1->acc[j].x = state->acc[j].x;
            state1->acc[j].y = state->acc[j].y;
            state1->acc[j].z = state->acc[j].z;
        }
    }
    
    // All of our input protons are running to begin with
    bool* idxrunning = new bool[state->N];
    
    bool nofailures = true;
    
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state2->pos[j].x = state1->pos[j].x + ( a2[0] * state1->vel[j].x ) * dt[j];
            state2->pos[j].y = state1->pos[j].y + ( a2[0] * state1->vel[j].y ) * dt[j];
            state2->pos[j].z = state1->pos[j].z + ( a2[0] * state1->vel[j].z ) * dt[j];

            state2->vel[j].x = state1->vel[j].x + ( a2[0] * state1->acc[j].x ) * dt[j];
            state2->vel[j].y = state1->vel[j].y + ( a2[0] * state1->acc[j].y ) * dt[j];
            state2->vel[j].z = state1->vel[j].z + ( a2[0] * state1->acc[j].z ) * dt[j];
        }
    }

    accelFunc(state2, field); // kv2 = f(x2, v2)

    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state3->pos[j].x = state1->pos[j].x + ( a3[0] * state1->vel[j].x
                                                  + a3[1] * state2->vel[j].x
                                                  ) * dt[j];
            
            state3->pos[j].y = state1->pos[j].y + ( a3[0] * state1->vel[j].y
                                                  + a3[1] * state2->vel[j].y
                                                  ) * dt[j];
            
            state3->pos[j].z = state1->pos[j].z + ( a3[0] * state1->vel[j].z
                                                  + a3[1] * state2->vel[j].z
                                                  ) * dt[j];

            state3->vel[j].x = state1->vel[j].x + ( a3[0] * state1->acc[j].x
                                                  + a3[1] * state2->acc[j].x
                                                  ) * dt[j];
            
            state3->vel[j].y = state1->vel[j].y + ( a3[0] * state1->acc[j].y
                                                  + a3[1] * state2->acc[j].y
                                                  ) * dt[j];
            
            state3->vel[j].z = state1->vel[j].z + ( a3[0] * state1->acc[j].z
                                                  + a3[1] * state2->acc[j].z
                                                  ) * dt[j];
        }
    }

    accelFunc(state3, field); // kv3 = f(x3, v3)

    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state4->pos[j].x = state1->pos[j].x + ( a4[0] * state1->vel[j].x
                                                  + a4[1] * state2->vel[j].x 
                                                  + a4[2] * state3->vel[j].x
                                                  ) * dt[j];

            state4->pos[j].y = state1->pos[j].y + ( a4[0] * state1->vel[j].y
                                                  + a4[1] * state2->vel[j].y 
                                                  + a4[2] * state3->vel[j].y
                                                  ) * dt[j];

            state4->pos[j].z = state1->pos[j].z + ( a4[0] * state1->vel[j].z
                                                  + a4[1] * state2->vel[j].z 
                                                  + a4[2] * state3->vel[j].z
                                                  ) * dt[j];

            state4->vel[j].x = state1->vel[j].x + ( a4[0] * state1->acc[j].x
                                                  + a4[1] * state2->acc[j].x 
                                                  + a4[2] * state3->acc[j].x
                                                  ) * dt[j];

            state4->vel[j].y = state1->vel[j].y + ( a4[0] * state1->acc[j].y
                                                  + a4[1] * state2->acc[j].y 
                                                  + a4[2] * state3->acc[j].y
                                                  ) * dt[j];

            state4->vel[j].z = state1->vel[j].z + ( a4[0] * state1->acc[j].z
                                                  + a4[1] * state2->acc[j].z 
                                                  + a4[2] * state3->acc[j].z
                                                  ) * dt[j];
        }
    }

    accelFunc(state4, field); // kv4 = f(x4, v4)

    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state5->pos[j].x = state1->pos[j].x + ( a5[0] * state1->vel[j].x
                                                  + a5[1] * state2->vel[j].x 
                                                  + a5[2] * state3->vel[j].x
                                                  + a5[3] * state4->vel[j].x
                                                  ) * dt[j];

            state5->pos[j].y = state1->pos[j].y + ( a5[0] * state1->vel[j].y
                                                  + a5[1] * state2->vel[j].y 
                                                  + a5[2] * state3->vel[j].y
                                                  + a5[3] * state4->vel[j].y
                                                  ) * dt[j];

            state5->pos[j].z = state1->pos[j].z + ( a5[0] * state1->vel[j].z
                                                  + a5[1] * state2->vel[j].z 
                                                  + a5[2] * state3->vel[j].z
                                                  + a5[3] * state4->vel[j].z
                                                  ) * dt[j];

            state5->vel[j].x = state1->vel[j].x + ( a5[0] * state1->acc[j].x
                                                  + a5[1] * state2->acc[j].x 
                                                  + a5[2] * state3->acc[j].x
                                                  + a5[3] * state4->acc[j].x
                                                  ) * dt[j];

            state5->vel[j].y = state1->vel[j].y + ( a5[0] * state1->acc[j].y
                                                  + a5[1] * state2->acc[j].y 
                                                  + a5[2] * state3->acc[j].y
                                                  + a5[3] * state4->acc[j].y
                                                  ) * dt[j];

            state5->vel[j].z = state1->vel[j].z + ( a5[0] * state1->acc[j].z
                                                  + a5[1] * state2->acc[j].z 
                                                  + a5[2] * state3->acc[j].z
                                                  + a5[3] * state4->acc[j].z
                                                  ) * dt[j];
        }
    }

    accelFunc(state5, field); // kv5 = f(x5, v5)

    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state6->pos[j].x = state1->pos[j].x + ( a6[0] * state1->vel[j].x
                                                  + a6[1] * state2->vel[j].x 
                                                  + a6[2] * state3->vel[j].x
                                                  + a6[3] * state4->vel[j].x
                                                  + a6[4] * state5->vel[j].x
                                                  ) * dt[j];

            state6->pos[j].y = state1->pos[j].y + ( a6[0] * state1->vel[j].y
                                                  + a6[1] * state2->vel[j].y 
                                                  + a6[2] * state3->vel[j].y
                                                  + a6[3] * state4->vel[j].y
                                                  + a6[4] * state5->vel[j].y
                                                  ) * dt[j];

            state6->pos[j].z = state1->pos[j].z + ( a6[0] * state1->vel[j].z
                                                  + a6[1] * state2->vel[j].z 
                                                  + a6[2] * state3->vel[j].z
                                                  + a6[3] * state4->vel[j].z
                                                  + a6[4] * state5->vel[j].z
                                                  ) * dt[j];

            state6->vel[j].x = state1->vel[j].x + ( a6[0] * state1->acc[j].x
                                                  + a6[1] * state2->acc[j].x 
                                                  + a6[2] * state3->acc[j].x
                                                  + a6[3] * state4->acc[j].x
                                                  + a6[4] * state5->acc[j].x
                                                  ) * dt[j];

            state6->vel[j].y = state1->vel[j].y + ( a6[0] * state1->acc[j].y
                                                  + a6[1] * state2->acc[j].y 
                                                  + a6[2] * state3->acc[j].y
                                                  + a6[3] * state4->acc[j].y
                                                  + a6[4] * state5->acc[j].y
                                                  ) * dt[j];

            state6->vel[j].z = state1->vel[j].z + ( a6[0] * state1->acc[j].z
                                                  + a6[1] * state2->acc[j].z 
                                                  + a6[2] * state3->acc[j].z
                                                  + a6[3] * state4->acc[j].z
                                                  + a6[4] * state5->acc[j].z
                                                  ) * dt[j];
        }
    }

    accelFunc(state6, field); // kv6 = f(x6, v6)

    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state7->pos[j].x = state1->pos[j].x + ( a7[0] * state1->vel[j].x
                                                  + a7[1] * state2->vel[j].x 
                                                  + a7[2] * state3->vel[j].x
                                                  + a7[3] * state4->vel[j].x
                                                  + a7[4] * state5->vel[j].x
                                                  + a7[5] * state6->vel[j].x
                                                  ) * dt[j];

            state7->pos[j].y = state1->pos[j].y + ( a7[0] * state1->vel[j].y
                                                  + a7[1] * state2->vel[j].y 
                                                  + a7[2] * state3->vel[j].y
                                                  + a7[3] * state4->vel[j].y
                                                  + a7[4] * state5->vel[j].y
                                                  + a7[5] * state6->vel[j].y
                                                  ) * dt[j];

            state7->pos[j].z = state1->pos[j].z + ( a7[0] * state1->vel[j].z
                                                  + a7[1] * state2->vel[j].z 
                                                  + a7[2] * state3->vel[j].z
                                                  + a7[3] * state4->vel[j].z
                                                  + a7[4] * state5->vel[j].z
                                                  + a7[5] * state6->vel[j].z
                                                  ) * dt[j];

            state7->vel[j].x = state1->vel[j].x + ( a7[0] * state1->acc[j].x
                                                  + a7[1] * state2->acc[j].x 
                                                  + a7[2] * state3->acc[j].x
                                                  + a7[3] * state4->acc[j].x
                                                  + a7[4] * state5->acc[j].x
                                                  + a7[5] * state6->acc[j].x
                                                  ) * dt[j];

            state7->vel[j].y = state1->vel[j].y + ( a7[0] * state1->acc[j].y
                                                  + a7[1] * state2->acc[j].y 
                                                  + a7[2] * state3->acc[j].y
                                                  + a7[3] * state4->acc[j].y
                                                  + a7[4] * state5->acc[j].y
                                                  + a7[5] * state6->acc[j].y
                                                  ) * dt[j];

            state7->vel[j].z = state1->vel[j].z + ( a7[0] * state1->acc[j].z
                                                  + a7[1] * state2->acc[j].z 
                                                  + a7[2] * state3->acc[j].z
                                                  + a7[3] * state4->acc[j].z
                                                  + a7[4] * state5->acc[j].z
                                                  + a7[5] * state6->acc[j].z
                                                  ) * dt[j];
        }
    }

    accelFunc(state7, field); // kv7 = f(x7, v7)

    double minErr = 0, maxErr = 0, sumErr = 0;
    long numFail = 0, numStopped = 0;

    #pragma omp parallel for reduction(+:numFail, numStopped, sumErr) reduction(max:maxErr) reduction(min:minErr)
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            // Absolute difference in acceleration between 4th- and 5th-order.
            double xErr = (b4[0]-b5[0]) * state1->acc[j].x
                        + (b4[1]-b5[1]) * state2->acc[j].x
                        + (b4[2]-b5[2]) * state3->acc[j].x
                        + (b4[3]-b5[3]) * state4->acc[j].x
                        + (b4[4]-b5[4]) * state5->acc[j].x
                        + (b4[5]-b5[5]) * state6->acc[j].x
                        + (b4[6]-b5[6]) * state7->acc[j].x;

            double yErr = (b4[0]-b5[0]) * state1->acc[j].y
                        + (b4[1]-b5[1]) * state2->acc[j].y
                        + (b4[2]-b5[2]) * state3->acc[j].y
                        + (b4[3]-b5[3]) * state4->acc[j].y
                        + (b4[4]-b5[4]) * state5->acc[j].y
                        + (b4[5]-b5[5]) * state6->acc[j].y
                        + (b4[6]-b5[6]) * state7->acc[j].y;

            double zErr = (b4[0]-b5[0]) * state1->acc[j].z
                        + (b4[1]-b5[1]) * state2->acc[j].z
                        + (b4[2]-b5[2]) * state3->acc[j].z
                        + (b4[3]-b5[3]) * state4->acc[j].z
                        + (b4[4]-b5[4]) * state5->acc[j].z
                        + (b4[5]-b5[5]) * state6->acc[j].z
                        + (b4[6]-b5[6]) * state7->acc[j].z;


            double err = dt[j] * sqrt((xErr*xErr + yErr*yErr + zErr*zErr)
                                     /(state1->vel[j].x*state1->vel[j].x
                                      +state1->vel[j].y*state1->vel[j].y
                                      +state1->vel[j].z*state1->vel[j].z));

            bool succ = err < rtol;

            if(lastIterationSuccess[j]) {
                double step_divisor = 1.25 * pow((err / rtol), error_scale_power);
                if( succ && step_divisor < 1/maxLengthen) step_divisor = 1/maxLengthen;
                if(!succ && step_divisor > maxFirstShorten ) step_divisor = maxFirstShorten;
                dt[j] /= step_divisor;
            }
            else if(!succ) {
                dt[j] /= maxOtherShorten;
            }
                        
            lastIterationSuccess[j] = succ;
            bool force;
            if(dt[j] < dt_min) {
                if(this->kill_failed_particles) {
                    state->running[j] = false;
                    numStopped++;
                } else {
                    dt[j] = dt_min;
                    force = true;
                }
            }
            if(succ || force) {
                
                state->pos[j].x = state7->pos[j].x;
                state->pos[j].y = state7->pos[j].y;
                state->pos[j].z = state7->pos[j].z;

                state->vel[j].x = state7->vel[j].x;
                state->vel[j].y = state7->vel[j].y;
                state->vel[j].z = state7->vel[j].z;

                state->acc[j].x = state7->acc[j].x;
                state->acc[j].y = state7->acc[j].y;
                state->acc[j].z = state7->acc[j].z;
            }

            minErr = fmin(minErr, err);
            maxErr = fmax(maxErr, err);
            sumErr += err;
            numFail += succ ? 0 : 1;
        }
    }
    state->N_running -= numStopped;
    
    // Proportion error in |vel| over the timestep, per proton.
    if(verbose) {
        printf("        Errors: min: %.1e, max: %.1e, mean: %.1e\n", minErr, maxErr, sumErr / N);

        if(numFail > 0)
            printf("            Failed for %li particles, will reduce their stepsize to try again.\n", numFail);
        
        if(numStopped > 0)
            printf("            %li particles failed to converge at miminum stepsize; stopping them.\n", numStopped);
    }
}
