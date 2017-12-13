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
    field->initFieldArrays(state->N);
    field->initFields();

    //Detector
    ParticleDetector* detector = getDetectorInfo();

    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        double distance = state->pos[j].z - field->min_z;
        state->pos[j].x -= distance * state->vel[j].x / state->vel[j].z; 
        state->pos[j].y -= distance * state->vel[j].y / state->vel[j].z;
        state->pos[j].z -= distance;
    }

    const long steps = 1500;

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    Integrator* integrator = getIntegratorInfo();
    integrator->initFromState(state);

    double sumTimes = 0;
    double sumSqTimes = 0;

    field->orientBeam(state);
    
    #ifdef USE_GL
        openWindow("shaders/shader.vert", "shaders/shader.frag");
        setupMatrix(source, field, detector);
        setupBuffers(&state->pos[0].x, state->running, state->N);
        updateBuffers(&state->pos[0].x, state->running, state->N);
        draw(state->N, 0, nullptr);
        printf("Window opened\n");
    #endif

    long i;


    for(i = 0; i < steps && state->N_running > 0; i++) {
        begin = std::chrono::steady_clock::now();

        integrator->step(state, field);
        field->invalidatePositions(state);
        invalidateStates(state);

        #ifdef USE_GL
            updateBuffers(&state->pos[0].x, state->running, state->N);
            draw(state->N, 0, nullptr);
        #endif

        end = std::chrono::steady_clock::now();
        double nanoTaken = std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count();
        sumTimes += nanoTaken;
        sumSqTimes += nanoTaken*nanoTaken;
        printf("Iteration %li finished, %li particles still running, took %f us\n", i, state->N_running, nanoTaken/1000.0);
    }
    printf("%li iterations taken \n", i);

    field->deorientBeam(state);
    
    detector->finalPush(state);

    #ifdef USE_GL
    field->orientBeam(state);
        char name[11] = "output.png";
        updateBuffers(&state->pos[0].x, state->running, state->N);
        draw(state->N, 10, name, false);
    field->deorientBeam(state);
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
            if(particles->vel[j].x*particles->vel[j].x + particles->vel[j].y*particles->vel[j].y + particles->vel[j].z*particles->vel[j].z > c*c 
            ) {
                particles->running[j] = false;
                #pragma omp atomic
                particles->N_running--;
            }
        }
    }
}
