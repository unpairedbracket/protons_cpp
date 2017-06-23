#include "main.h"

#include <cmath>
#include <cstdio>

#include <chrono>
#include <iostream>

#include <omp.h>

int main() {
    const double mass = 1;
    const double charge = 1;

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

    double B = 0.01;
    double dt = 1;

    initPos(positionX, positionY, positionZ, N);
    initVel(velocityX, velocityY, velocityZ, n, m);

    //print_status(positionX, positionY, positionZ, velocityX, velocityY, velocityZ, N);
    printf("\n");
    const long steps = 1000;

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    for(long i = 0; i < steps; i++) {
        begin = std::chrono::steady_clock::now();
        #pragma omp parallel for
        for(long j = 0; j < N; j++) {
            accelX[j] = charge * B / mass * velocityY[j]; // aX = q/m ( vy Bz - vz By );
            accelY[j] = charge * B / mass *-velocityX[j]; // aY = q/m ( vz Bx - vx Bz );
            accelZ[j] = 0;                                // aZ = q/m ( vx By - vy Bx );
        //}

        //for(long j = 0; j < N; j++) {
            positionX[j] += velocityX[j] * dt;
            positionY[j] += velocityY[j] * dt;
            positionZ[j] += velocityZ[j] * dt;

            velocityX[j] += accelX[j] * dt;
            velocityY[j] += accelY[j] * dt;
            velocityZ[j] += accelZ[j] * dt;
        }
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

    return 0;
}

void initPos(double* X, double* Y, double* Z, long N) {
    for(long i = 0; i < N; i++) {
        X[i] = 0;
        Y[i] = 0;
        Z[i] = 0;
    }
}

void initVel(double* vX, double* vY, double* vZ, long n, long m) {
    for(long i = 0; i < m*n; i++) {
        long j = i / n;
        long k = i % n;
        double angleX = (2*j / ((double)m-1) - 1) * pi() / 8.0;
        double angleY = (2*k / ((double)n-1) - 1) * pi() / 8.0;

        vX[i] = sin(angleX);
        vY[i] = sin(angleY);
        vZ[i] = sqrt(cos(angleX + angleY) * cos(angleX - angleY));
    }
}

void print_status(double* X, double* Y, double* Z, double* vX, double* vY, double* vZ, long N) {
    for(long i = 0; i < N; i++) {
        printf("pos = [%f, %f, %f]; vel = [%f, %f, %f]\n", X[i], Y[i], Z[i], vX[i], vY[i], vZ[i]);
    }
}

