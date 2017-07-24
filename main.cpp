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
        state.posX[j] += source.distance * state.velX[j] / state.velZ[j]; 
        state.posY[j] += source.distance * state.velY[j] / state.velZ[j];
        state.posZ[j] += source.distance;
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

    for(i = 0; i < steps && state.N_running > 0; i++) {
        begin = std::chrono::steady_clock::now();

        step(
                RKDPIntegrator,
                &state, field,
                accel
            );
        invalidateStates(&state);
        end = std::chrono::steady_clock::now();
        //printf("After iteration %ld\n", i+1);
        //print_status(state);
        //printf("\n");
        double nanoTaken = std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count();
        //std::cout << "Time difference = " << nanoTaken/1000.0 << " us" <<std::endl;
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

void accel(ParticleState* state, FieldStructure* fields, double* accx, double* accy, double* accz) {
    fields->getFields(state);
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            accx[j] = state->particleInfo->qmratio * ( fields->Ex[j] + state->velY[j] * fields->Bz[j] - state->velZ[j] * fields->By[j] ); // aX = q/m ( vy Bz - vz By );
            accy[j] = state->particleInfo->qmratio * ( fields->Ey[j] + state->velZ[j] * fields->Bx[j] - state->velX[j] * fields->Bz[j] ); // aY = q/m ( vz Bx - vx Bz );
            accz[j] = state->particleInfo->qmratio * ( fields->Ez[j] + state->velX[j] * fields->By[j] - state->velY[j] * fields->Bx[j] ); // aZ = q/m ( vx By - vy Bx );
        }
    }
}

void invalidateStates(ParticleState* particles) {
    #pragma omp parallel for
    for(long j = 0; j < particles->N; j++) {
        if(particles->running[j])
        {
            if(particles->posZ[j] > 10 * 0.001) {
                particles->running[j] = false;
                #pragma omp atomic
                particles->N_running--;
            }
        }
    }
}

void initPos(ParticleSource* sourceInfo, ParticleState* particles) {
    for(long i = 0; i < particles->N; i++) {
        particles->posX[i] = 0;
        particles->posY[i] = 0;
        particles->posZ[i] = -sourceInfo->distance;
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

        particles->velX[i] = speed * sin(angleX);
        particles->velY[i] = speed * sin(angleY);
        particles->velZ[i] = speed * sqrt(cos(angleX + angleY) * cos(angleX - angleY));
    }
}

void print_status(ParticleState* state) {
    for(long i = 0; i < 3; i++) {//state->N; i++) {
        printf("pos = [%f, %f, %f]; vel = [%f, %f, %f]\n", state->posX[i], state->posY[i], state->posZ[i], state->velX[i], state->velY[i], state->velZ[i]);
    }
}

void print_status_raw(ParticleState* state) {
    std::ofstream posfile;
    posfile.open("pos.txt");
    for(long i = 0; i < state->N; i++) {
        posfile << state->posX[i] << "," << state->posY[i] << "," << state->posZ[i] << std::endl;
    }
    posfile.close();
    
    std::ofstream velfile;
    velfile.open("vel.txt");
    for(long i = 0; i < state->N; i++) {
        velfile << state->velX[i] << "," << state->velY[i] << "," << state->velZ[i] << std::endl;
    }
    velfile.close();
}

