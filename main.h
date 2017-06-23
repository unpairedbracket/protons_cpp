#include <cmath>
#include <cstdio>

#include <chrono>
#include <iostream>

#include <omp.h>

constexpr double pi() { return std::atan(1)*4; }

void initPos(double* X, double* Y, double* Z, long N);

void initVel(double* vX, double* vY, double* vZ, long n, long m);

void print_status(double* X, double* Y, double* Z, double* vX, double* vY, double* vZ, long N);

int main();
