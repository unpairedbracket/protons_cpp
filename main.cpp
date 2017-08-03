#include "main.h"

int main(int argc, char *argv[]) {
    const double mass = 1.6726217E-27;
    const double charge = 1.6021766208E-19;
    ParticleInfo protons;
    initParticle(&protons, mass, charge);

    long n_x;
    long n_y;

    if(argc != 3) {
        n_x = 100;
        n_y = 100;
    } else {
        n_x = atol(argv[1]);
        n_y = atol(argv[2]);
    }

    long N = n_x * n_y;

    //Source
    double sourceDistance = 0.02; //Distance from proton source to front object plane, metres
    double sourceDivergence = 0.04; //Radians from axis, i.e. pi/2 for hemisphere-ish. Must currently be smaller than pi/4 due to dodgy implementation.
    double sourceEnergy = 40E6 * charge; //Source energy, joules
    ParticleSource source;
    initSource(&source, &protons, sourceDistance, sourceDivergence, sourceEnergy, n_x, n_y);

    ParticleState state;
    initParticleState(&state, &source);
    initPos(&source, &state);
    initVel(&source, &state);

    print_status(&state);

    //Object
    double objectLength = 0.001; //Distance from front object plane to back object plane, metres
    double objectRadius = 0.0005; //Transverse size of object
    double fieldStrength = 25;
    CocoonField* field = initCocoonField(fieldStrength, objectRadius, objectLength, state.N);
    field->initFields();

    //Detector
    double detectorDistance = 1; //Distance from back object plane to image plane
    ParticleDetector detector;
    initDetector(&detector, &protons, detectorDistance);

    #pragma omp parallel for
    for(long j = 0; j < N; j++) {
        state.pos[j].x += source.distance * state.vel[j].x / state.vel[j].z; 
        state.pos[j].y += source.distance * state.vel[j].y / state.vel[j].z;
        state.pos[j].z += source.distance;
    }

    print_status(&state);

    printf("\n");
    const long steps = 1000;

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    integrator* RKDPIntegrator = makeIntegrator(N);

    double sumTimes, sumSqTimes;

    double dt_0 = 0.0001 / 87538948.0;

    #pragma omp parallel for
    for(long j = 0; j < N; j++) {
        RKDPIntegrator->dt[j] = dt_0;
    }

    long i;

    openWindow("shaders/shader.vert", "shaders/shader.frag");
    setupMatrix(source.distance, detector.distance, source.divergence);
    setupBuffers(&state.pos[0].x, state.running, state.N);
    updateBuffers(&state.pos[0].x, state.running, state.N);
    draw(state.N);

    for(i = 0; i < steps && state.N_running > 0; i++) {
        begin = std::chrono::steady_clock::now();

        step(
                RKDPIntegrator,
                &state, field,
                accel
            );
        invalidateStates(&state);
        updateBuffers(&state.pos[0].x, state.running, state.N);
        draw(state.N);
        end = std::chrono::steady_clock::now();
        double nanoTaken = std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count();
        sumTimes += nanoTaken;
        sumSqTimes += nanoTaken*nanoTaken;
        printf("Iteration %li finished, %li particles still running\n", i, state.N_running);
    }
    printf("%li iterations taken \n", i);
    /*
    #pragma omp parallel for
    for(long j = 0; j < N; j++) {
        state.posX[j] += state.velX[j] / state.velZ[j] * detector.distance;
        state.posY[j] += state.velY[j] / state.velZ[j] * detector.distance;
        state.posZ[j] += detector.distance;
    }
    */
    std::cout << "Average time taken: " << sumTimes / (1000.0 * i) << " us (" << sumTimes / (1000.0 * i * state.N) << " us per particle)" << std::endl;
    std::cout << "Standard Deviation: " << sqrt((sumSqTimes/i) - (sumTimes*sumTimes)/(i*i)) / 1000.0 << " us (" << sqrt((sumSqTimes/i) - (sumTimes*sumTimes)/(i*i)) / (1000.0 * state.N) << " us per particle)" << std::endl;

    print_status(&state);
    print_status_raw(&state);
    return 0;
}

void accel(ParticleState* state, FieldStructure* fields, Vector3* acc) {
    fields->getFields(state);
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            acc[j].x = state->particleInfo->qmratio * ( fields->E[j].x + state->vel[j].y * fields->B[j].z - state->vel[j].z * fields->B[j].y ); // aX = q/m ( vy Bz - vz By );
            acc[j].y = state->particleInfo->qmratio * ( fields->E[j].y + state->vel[j].z * fields->B[j].x - state->vel[j].x * fields->B[j].z ); // aY = q/m ( vz Bx - vx Bz );
            acc[j].z = state->particleInfo->qmratio * ( fields->E[j].z + state->vel[j].x * fields->B[j].y - state->vel[j].y * fields->B[j].x ); // aZ = q/m ( vx By - vy Bx );
        }
    }
}

void invalidateStates(ParticleState* particles) {
    #pragma omp parallel for
    for(long j = 0; j < particles->N; j++) {
        if(particles->running[j])
        {
            if(particles->pos[j].z > 10 * 0.001
            || particles->pos[j].x*particles->pos[j].x + particles->pos[j].y*particles->pos[j].y > 0.001*0.001
            || particles->vel[j].x*particles->vel[j].x + particles->vel[j].y*particles->vel[j].y + particles->vel[j].z*particles->vel[j].z > 3E8*3E8
            ) {
                particles->running[j] = false;
                #pragma omp atomic
                particles->N_running--;
            }
        }
    }
}

void initPos(ParticleSource* sourceInfo, ParticleState* particles) {
    for(long i = 0; i < particles->N; i++) {
        particles->pos[i].x = 0;
        particles->pos[i].y = 0;
        particles->pos[i].z = -sourceInfo->distance;
    }
}

void initVel(ParticleSource* sourceInfo, ParticleState* particles) {
    assert(sourceInfo->x_extent*sourceInfo->y_extent == particles->N);
    assert(sourceInfo->divergence < pi() / 4);

    double speed = sqrt(2 * sourceInfo->energy / particles->particleInfo->mass);
    printf("speed: %f\n", speed);
    for(long i = 0; i < particles->N; i++) {
        long j = i / sourceInfo->y_extent;
        long k = i % sourceInfo->y_extent;
        //printf("(%lu, %lu)", j, k);
        double angleX = (2*j / ((double)sourceInfo->x_extent-1) - 1) * sourceInfo->divergence;
        double angleY = (2*k / ((double)sourceInfo->y_extent-1) - 1) * sourceInfo->divergence;

        particles->vel[i].x = speed * sin(angleX);
        particles->vel[i].y = speed * sin(angleY);
        particles->vel[i].z = speed * sqrt(cos(angleX + angleY) * cos(angleX - angleY));
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

