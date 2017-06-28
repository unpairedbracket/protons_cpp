#include "main.h"

int main() {
    const double mass = 1;
    const double charge = 1;
    ParticleInfo* protons = makeParticle(mass, charge);

    long n = 1000;
    long m = 1000;
    long N = n * m;
    double* positionX = new double[N];
    double* positionY = new double[N];
    double* positionZ = new double[N];

    double* velocityX = new double[N];
    double* velocityY = new double[N];
    double* velocityZ = new double[N];

    double* accelX = new double[N];
    double* accelY = new double[N];
    double* accelZ = new double[N];

    double* Bx = new double[N];
    double* By = new double[N];
    double* Bz = new double[N];

    double* Ex = new double[N];
    double* Ey = new double[N];
    double* Ez = new double[N];

    //Source
    double sourceDistance = 1.0; //Distance from proton source to front object plane, metres
    double sourceDivergence = pi()/8.0; //Radians from axis, i.e. pi/2 for hemisphere
    double sourceEnergy = 1.0; //Source energy, joules
    ParticleSource* source = makeSource(protons, sourceDistance, sourceDivergence, sourceEnergy);

    //Object
    double objectLength = 0.01; //Distance from front object plane to back object plane, metres
    FieldStructure* field = makeFieldStructure(N);
    double objectRadius = 0.5; //Transverse size of object

    //Detector
    double detectorDistance; //Distance from back object plane to image plane
    ParticleDetector* detector = makeDetector(protons, detectorDistance);

    double* dt = new double[N];

    initPos(positionX, positionY, positionZ, N);
    initVel(source, velocityX, velocityY, velocityZ, n, m);

    #pragma omp parallel for
    for(long j = 0; j < N; j++) {
        positionX[j] += velocityX[j] / velocityZ[j] * sourceDistance;
        positionY[j] += velocityY[j] / velocityZ[j] * sourceDistance;
        positionZ[j] += sourceDistance;
    }

    printf("\n");
    const long steps = 1000;

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    integrator* RKDPIntegrator = makeIntegrator(N);

    for(long i = 0; i < steps; i++) {
        begin = std::chrono::steady_clock::now();

        step(
                RKDPIntegrator,
                protons, field,
                positionX, positionY, positionZ,
                velocityX, velocityY, velocityZ,
                dt, N,
                accel
            );
       
        end = std::chrono::steady_clock::now();
        //printf("After iteration %ld\n", i+1);
        //print_status(positionX, positionY, positionZ, velocityX, velocityY, velocityZ, N);
        //printf("\n");
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count()/1000.0 << " us" <<std::endl;
    }

//    print_status(positionX, positionY, positionZ, velocityX, velocityY, velocityZ, N);

    free(positionX);
    free(positionY);
    free(positionZ);
    free(velocityX);
    free(velocityY);
    free(velocityZ);
    free(accelX);
    free(accelY);
    free(accelZ);

    free(Bx);
    free(By);
    free(Bz);
    free(Ex);
    free(Ey);
    free(Ez);

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

void initPos(double* X, double* Y, double* Z, long N) {
    for(long i = 0; i < N; i++) {
        X[i] = 0;
        Y[i] = 0;
        Z[i] = 0;
    }
}

void initVel(ParticleSource* sourceInfo, double* vX, double* vY, double* vZ, long n, long m) {
    for(long i = 0; i < m*n; i++) {
        long j = i / n;
        long k = i % n;
        double angleX = (2*j / ((double)m-1) - 1) * sourceInfo->divergence;
        double angleY = (2*k / ((double)n-1) - 1) * sourceInfo->divergence;
        double speed = sqrt(2 * sourceInfo->energy / sourceInfo->particleInfo->mass);

        vX[i] = speed * sin(angleX);
        vY[i] = speed * sin(angleY);
        vZ[i] = speed * sqrt(cos(angleX + angleY) * cos(angleX - angleY));
    }
}

void print_status(double* X, double* Y, double* Z, double* vX, double* vY, double* vZ, long N) {
    for(long i = 0; i < N; i++) {
        printf("pos = [%f, %f, %f]; vel = [%f, %f, %f]\n", X[i], Y[i], Z[i], vX[i], vY[i], vZ[i]);
    }
}

