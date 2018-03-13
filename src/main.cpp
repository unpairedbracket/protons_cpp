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
    ParticleInfo particleType = getParticleInfo();

    //Source
    ParticleSource* source = getSourceInfo();

    ParticleState* state = source->createParticleState(particleType);

    Interpolator* interp = getInterpolatorInfo();
    bool shouldInterpolate = interp;

    //Object
    FieldStructure* field = getFieldsInfo();
    field->initFieldArrays(state->N);
    field->initFields();

    //Detector
    ParticleDetector* detector = getDetectorInfo();

    Integrator* integrator = getIntegratorInfo();
    integrator->initFromState(state);


    #ifdef USE_GL
        int drawType = 2;

        if(drawType == 0) {
            openWindow("shaders/shader.vert", "shaders/shader.frag");
            setupMatrix(source, field, detector);
            setupBuffers(&state->pos[0].x, state->running, state->N);
            updateBuffers(&state->pos[0].x, state->running, state->N);
            draw(state->N, 0, nullptr);
        } else if(drawType == 1) {
            openWindow("shaders/shader.vert", "shaders/shader.frag");
            setupMatrix(0.001, field);
            setupBuffers(&state->vel[0].x, state->running, state->N);
            updateBuffers(&state->vel[0].x, state->running, state->N);
            draw(state->N, 0, nullptr);
        } else if(drawType == 2) {
            printf("Setting up detector rendering\n");
            detector->setupGraphics();
        }
        printf("Window opened\n");
    #endif

    int number_runs = getNumberRuns();

    for(int iteration = 0; iteration<number_runs; iteration++) {
    source->setParticleState(state);
    setAllRunning(state);
    integrator->reset();
    detector->detectUndeviated(state);

    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        double distance = state->pos[j].z - field->min_z;
        state->pos[j].x -= distance * state->vel[j].x / state->vel[j].z;
        state->pos[j].y -= distance * state->vel[j].y / state->vel[j].z;
        state->pos[j].z -= distance;
    }

    if(shouldInterpolate)
        interp->setSamplePoints(state);


    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    double sumTimes = 0;
    double sumSqTimes = 0;

    field->orientBeam(state);

    long i = 0;
    const long steps = 15000;

    for(i = 0; i < steps && state->N_running > 0; i++) {
        begin = std::chrono::steady_clock::now();

        integrator->step(state, field);
        field->invalidatePositions(state);
        invalidateStates(state);

        #ifdef USE_GL
            if(drawType == 0) {
                updateBuffers(&state->pos[0].x, state->running, state->N);
                draw(state->N, 0, nullptr);
            } else if(drawType == 1) {
                updateBuffers(&state->pos[0].x, state->running, state->N);
                draw(state->N, 0, nullptr);
            }   
        #endif

        end = std::chrono::steady_clock::now();
        double nanoTaken = std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count();
        sumTimes += nanoTaken;
        sumSqTimes += nanoTaken*nanoTaken;
        //printf("Iteration %li finished, %li particles still running, took %f us\n", i, state->N_running, nanoTaken/1000.0);
    }
    printf("%d: %li iterations taken \n", iteration, i+1);

    std::cout << "Average time taken: " << sumTimes / (1000.0 * i) << " us (" << sumTimes / (1000.0 * i * state->N) << " us per particle)" << std::endl;
    std::cout << "Standard Deviation: " << sqrt((sumSqTimes/i) - (sumTimes*sumTimes)/(i*i)) / 1000.0 << " us (" << sqrt((sumSqTimes/i) - (sumTimes*sumTimes)/(i*i)) / (1000.0 * state->N) << " us per particle)" << std::endl;

    field->deorientBeam(state);
    
    double density = state->N /( 2*pi() * (1 - cos(source->divergence)) * (source->distance + detector->distance) * (source->distance + detector->distance));

    detector->setExpectedFluence((iteration+1) * density);

    if(!shouldInterpolate) {
        detector->finalPush(state);
        detector->detect(state);
        #ifdef USE_GL
            if(drawType == 0) {
                field->orientBeam(state);
                char name[11] = "output.png";
                updateBuffers(&state->pos[0].x, state->running, state->N);
                draw(state->N, 10, name, false);
                field->deorientBeam(state);
            } else if(drawType == 1) {
                field->orientBeam(state);
                char name[11] = "output.png";
                updateBuffers(&state->pos[0].x, state->running, state->N);
                draw(state->N, 10, name, false);
                field->deorientBeam(state);
            } else if(drawType == 2) {
                detector->draw();
            }
        #endif
    } else {
        begin = std::chrono::steady_clock::now();
        interp->setSampleValues(state);
        for(int i = 0; i < interp->iterations; i++) {
            printf("Initialising particles for iteration %d... ", i);
            interp->initState(particleType);
            if(i==69) {printf("NICE\n");} else {printf("DONE\n");}
            printf("Interpolating %ld particles... ", interp->interpParticles->N);
            interp->interpolate();
            printf("DONE\n");
            detector->finalPush(interp->interpParticles);
            detector->detect(interp->interpParticles);
            #ifdef USE_GL
                if(drawType == 0) {
                    field->orientBeam(state);
                    char name[11] = "output.png";
                    updateBuffers(&state->pos[0].x, state->running, state->N);
                    draw(state->N, 10, name, false);
                    field->deorientBeam(state);
                } else if(drawType == 1) {
                    field->orientBeam(state);
                    char name[11] = "output.png";
                    updateBuffers(&state->pos[0].x, state->running, state->N);
                    draw(state->N, 10, name, false);
                    field->deorientBeam(state);
                } else if(drawType == 2) {
                    detector->draw();
                }
            #endif
        }
        end = std::chrono::steady_clock::now();
        double nanoTaken = std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count();
        printf("Interpolation took %f microsec per particle\n", nanoTaken/(1000.0 * interp->iterations * interp->interpParticles->N));
    }
    if((iteration+1) % 1000 == 0) detector->output();
    }

    detector->performInversion();

    detector->output();

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
