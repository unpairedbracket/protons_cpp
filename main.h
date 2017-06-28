#pragma once

#include <cmath>
#include <cstdio>

#include <chrono>
#include <iostream>

#include <omp.h>

#include "structs.h"
#include "integrator.h"
#include "fields.h"

constexpr double pi() { return std::atan(1)*4; }

void initPos(double* X, double* Y, double* Z, long N);

void initVel(ParticleSource* source, double* vX, double* vY, double* vZ, long n, long m);

void print_status(double* X, double* Y, double* Z, double* vX, double* vY, double* vZ, long N);

int main();
void accel(ParticleInfo* particle, FieldStructure* fields, double* x, double* y, double* z, double* vx, double* vy, double* vz, double* accx, double* accy, double* accz, long N);
