#include "main.h"

int main(int argc, char *argv[]) {
    const double mass = 1.6726217E-27;
    const double charge = 1.6021766208E-19;

    load_config("configs/config.yml");
    ParticleInfo* particleType = getParticleInfo();

    //Source
    ParticleSource* source = getSourceInfo();

    ParticleState* state = source->genParticleState(particleType);
    print_status(state);

    //Object
    FieldStructure* field = getFieldsInfo();
    initFieldArrays(field, state->N);

    //Detector
    double detectorDistance = 1; //Distance from back object plane to image plane
    ParticleDetector detector;
    initDetector(&detector, particleType, detectorDistance);

    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        state->pos[j].x -= state->pos[j].z * state->vel[j].x / state->vel[j].z; 
        state->pos[j].y -= state->pos[j].z * state->vel[j].y / state->vel[j].z;
        state->pos[j].z -= state->pos[j].z;
    }

    print_status(state);

    printf("\n");
    const long steps = 10000;

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    Integrator* integrator = getIntegrator();
    integrator->initFromState(state);

    double sumTimes = 0;
    double sumSqTimes = 0;

    if(USE_GL) {
        openWindow("shaders/shader.vert", "shaders/shader.frag");
        //setupMatrix(source->distance, detector.distance * 1.01, source.divergence);
        setupMatrix(0.02, detector.distance * 1.01, 0.05);
        setupBuffers(&state->pos[0].x, state->running, state->N);
        updateBuffers(&state->pos[0].x, state->running, state->N);
        draw(state->N);
    }

    long i;

    for(i = 0; i < steps && state->N_running > 0; i++) {
        begin = std::chrono::steady_clock::now();

        integrator->step(state, field);
        invalidateStates(state);

        if(USE_GL) {
            updateBuffers(&state->pos[0].x, state->running, state->N);
            draw(state->N);
        }

        end = std::chrono::steady_clock::now();
        double nanoTaken = std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count();
        sumTimes += nanoTaken;
        sumSqTimes += nanoTaken*nanoTaken;
        printf("Iteration %li finished, %li particles still running, took %f us\n", i, state->N_running, nanoTaken/1000.0);
    }
    printf("%li iterations taken \n", i);
    
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        state->pos[j].x += (detector.distance - state->pos[j].z) * state->vel[j].x / state->vel[j].z;
        state->pos[j].y += (detector.distance - state->pos[j].z) * state->vel[j].y / state->vel[j].z;
        state->pos[j].z += (detector.distance - state->pos[j].z);
    }
    
    if(USE_GL) {
        updateBuffers(&state->pos[0].x, state->running, state->N);
        draw(state->N, true);
    }
    
    std::cout << "Average time taken: " << sumTimes / (1000.0 * i) << " us (" << sumTimes / (1000.0 * i * state->N) << " us per particle)" << std::endl;
    std::cout << "Standard Deviation: " << sqrt((sumSqTimes/i) - (sumTimes*sumTimes)/(i*i)) / 1000.0 << " us (" << sqrt((sumSqTimes/i) - (sumTimes*sumTimes)/(i*i)) / (1000.0 * state->N) << " us per particle)" << std::endl;

    print_status(state);
    print_status_raw(state);

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
            || particles->vel[j].x*particles->vel[j].x + particles->vel[j].y*particles->vel[j].y + particles->vel[j].z*particles->vel[j].z > 3E8*3E8
            ) {
                particles->running[j] = false;
                #pragma omp atomic
                particles->N_running--;
            }
        }
    }
}

void print_status(ParticleState* state) {
    for(long i = 0; i < 3; i++) {//state->N; i++) {
        printf("pos = [%f, %f, %f]; vel = [%f, %f, %f]\n", state->pos[i].x, state->pos[i].y, state->pos[i].z, state->vel[i].x, state->vel[i].y, state->vel[i].z);
    }
}

void print_status_raw(ParticleState* state) {
    std::ofstream posfile;
    posfile.open("pos.txt");
    for(long i = 0; i < state->N; i++) {
        posfile << state->pos[i].x << "," << state->pos[i].y << "," << state->pos[i].z << std::endl;
    }
    posfile.close();
    
    std::ofstream velfile;
    velfile.open("vel.txt");
    for(long i = 0; i < state->N; i++) {
        velfile << state->vel[i].x << "," << state->vel[i].y << "," << state->vel[i].z << std::endl;
    }
    velfile.close();
}

