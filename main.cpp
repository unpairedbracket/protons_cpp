#include "main.h"

int main(int argc, char *argv[]) {
    const double mass = 1;
    const double charge = 1;
    ParticleInfo* protons = makeParticle(mass, charge);

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
    double sourceDistance = 1.0; //Distance from proton source to front object plane, metres
    double sourceDivergence = pi()/8.0; //Radians from axis, i.e. pi/2 for hemisphere-ish. Must currently be smaller than pi/4 due to dodgy implementation.
    double sourceEnergy = 1.0; //Source energy, joules
    ParticleSource* source = makeSource(protons, sourceDistance, sourceDivergence, sourceEnergy, n_x, n_y);

    ParticleState* state = makeParticleState(source);
    initPos(state);
    initVel(source, state);

    print_status(state);

    //Object
    double objectLength = 0.01; //Distance from front object plane to back object plane, metres
    FieldStructure* field = makeFieldStructure(N);
    double objectRadius = 0.5; //Transverse size of object

    //Detector
    double detectorDistance; //Distance from back object plane to image plane
    ParticleDetector* detector = makeDetector(protons, detectorDistance);

    #pragma omp parallel for
    for(long j = 0; j < N; j++) {
        state->posX[j] += state->velX[j] / state->velZ[j] * source->distance;
        state->posY[j] += state->velY[j] / state->velZ[j] * source->distance;
        state->posZ[j] += sourceDistance;
    }

    print_status(state);

    printf("\n");
    const long steps = 1000;

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    integrator* RKDPIntegrator = makeIntegrator(N);

    double sumTimes, sumSqTimes;

    #pragma omp parallel for
    for(long j = 0; j < N; j++) {
        RKDPIntegrator->dt[j] = 0.1;
    }

    for(long i = 0; i < steps; i++) {
        begin = std::chrono::steady_clock::now();

        step(
                RKDPIntegrator,
                state, field,
                accel
            );
        end = std::chrono::steady_clock::now();
        //printf("After iteration %ld\n", i+1);
        //print_status(state);
        //printf("\n");
        double nanoTaken = std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count();
        //std::cout << "Time difference = " << nanoTaken/1000.0 << " us" <<std::endl;
        sumTimes += nanoTaken;
        sumSqTimes += nanoTaken*nanoTaken;
    }
    std::cout << "Average time taken: " << sumTimes / (1000.0 * steps) << " us (" << sumTimes / (1000.0 * steps * state->N) << " us per particle)" << std::endl;
    std::cout << "Standard Deviation: " << sqrt((sumSqTimes/steps) - (sumTimes*sumTimes)/(steps*steps)) / 1000.0 << " us (" << sqrt((sumSqTimes/steps) - (sumTimes*sumTimes)/(steps*steps)) / (1000.0 * state->N) << " us per particle)" << std::endl;

//    print_status(positionX, positionY, positionZ, velocityX, velocityY, velocityZ, N);
    return 0;
}

void accel(ParticleInfo* particle, FieldStructure* fields, double* x, double* y, double* z, double* vx, double* vy, double* vz, double* accx, double* accy, double* accz, long N) {
    getFields(fields, x, y, z, N);
    #pragma omp parallel for
    for(long j = 0; j < N; j++) {
        accx[j] = particle->qmratio * ( fields->Ex[j] + vy[j] * fields->Bz[j] - vz[j] * fields->By[j] ); // aX = q/m ( vy Bz - vz By );
        accy[j] = particle->qmratio * ( fields->Ey[j] + vz[j] * fields->Bx[j] - vx[j] * fields->Bz[j] ); // aY = q/m ( vz Bx - vx Bz );
        accz[j] = particle->qmratio * ( fields->Ez[j] + vx[j] * fields->By[j] - vy[j] * fields->Bx[j] ); // aZ = q/m ( vx By - vy Bx );
    }
}

void initPos(ParticleState* particles) {
    for(long i = 0; i < particles->N; i++) {
        particles->posX[i] = 0;
        particles->posY[i] = 0;
        particles->posZ[i] = 0;
    }
}

void initVel(ParticleSource* sourceInfo, ParticleState* particles) {
    assert(sourceInfo->x_extent*sourceInfo->y_extent == particles->N);
    assert(sourceInfo->divergence < pi() / 4);

    double speed = sqrt(2 * sourceInfo->energy / particles->particleInfo->mass);
    printf("speed: %f", speed);
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

