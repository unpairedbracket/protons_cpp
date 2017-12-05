#include "main.h"

#include <cstdio>
#include <cassert>

#include <chrono>
#include <fstream>
#include <iostream>

#include <omp.h>

#include "util/math.h"
#include "util/physical_constants.h"
#include "detectors/detector.h"
#include "integrators/integrator.h"
#include "config/config_parser.h"

#ifdef USE_GL
#include "graphics/window.h"
#endif

int main(int argc, char *argv[]) {
    load_config("configs/config.yml", "defaults/config.yml");
    ParticleInfo* particleType = getParticleInfo();

    //Source
    ParticleSource* source = getSourceInfo();

    ParticleState* state = source->genParticleState(particleType);

    //Object
    FieldStructure* field = getFieldsInfo();
    initFieldArrays(field, state->N);
    field->initFields();

    //Detector
    ParticleDetector* detector = getDetectorInfo();

    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        state->pos[j].x -= state->pos[j].z * state->vel[j].x / state->vel[j].z; 
        state->pos[j].y -= state->pos[j].z * state->vel[j].y / state->vel[j].z;
        state->pos[j].z -= state->pos[j].z;
    }

    const long steps = 10000;

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    Integrator* integrator = getIntegratorInfo();
    integrator->initFromState(state);

    double sumTimes = 0;
    double sumSqTimes = 0;

    #ifdef USE_GL
        openWindow("shaders/shader.vert", "shaders/shader.frag");
        setupMatrix(source->distance, detector->distance * 1.01, source->divergence);
        setupBuffers(&state->pos[0].x, state->running, state->N);
        updateBuffers(&state->pos[0].x, state->running, state->N);
        draw(state->N);
        printf("Window opened\n");
    #endif

    long i;

    for(i = 0; i < steps && state->N_running > 0; i++) {
        begin = std::chrono::steady_clock::now();

        integrator->step(state, field);
        invalidateStates(state);

        #ifdef USE_GL
            updateBuffers(&state->pos[0].x, state->running, state->N);
            draw(state->N);
        #endif

        end = std::chrono::steady_clock::now();
        double nanoTaken = std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count();
        sumTimes += nanoTaken;
        sumSqTimes += nanoTaken*nanoTaken;
        printf("Iteration %li finished, %li particles still running, took %f us\n", i, state->N_running, nanoTaken/1000.0);
    }
    printf("%li iterations taken \n", i);
    
    detector->finalPush(state);

    #ifdef USE_GL
        updateBuffers(&state->pos[0].x, state->running, state->N);
        draw(state->N, true);
    #endif
    
    std::cout << "Average time taken: " << sumTimes / (1000.0 * i) << " us (" << sumTimes / (1000.0 * i * state->N) << " us per particle)" << std::endl;
    std::cout << "Standard Deviation: " << sqrt((sumSqTimes/i) - (sumTimes*sumTimes)/(i*i)) / 1000.0 << " us (" << sqrt((sumSqTimes/i) - (sumTimes*sumTimes)/(i*i)) / (1000.0 * state->N) << " us per particle)" << std::endl;

    detector->output(state);

    integrator->deinit();
    return 0;
}

void invalidateStates(ParticleState* particles) {
    #pragma omp parallel for
    for(long j = 0; j < particles->N; j++) {
        if(particles->running[j])
        {
            if(particles->pos[j].z > 10 * 0.001 * 100
            || particles->pos[j].x*particles->pos[j].x + particles->pos[j].y*particles->pos[j].y > 0.002*0.002
            || particles->vel[j].x*particles->vel[j].x + particles->vel[j].y*particles->vel[j].y + particles->vel[j].z*particles->vel[j].z > c*c
            ) {
                particles->running[j] = false;
                #pragma omp atomic
                particles->N_running--;
            }
        }
    }
}
