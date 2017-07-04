#pragma once

#include <cmath>
#include <cstdio>
#include <cassert>

#include <chrono>
#include <iostream>

#include <omp.h>

#include "structs.h"
#include "integrator.h"
#include "fields.h"

constexpr double pi() { return std::atan(1)*4; }

void initPos(ParticleState* particles);

void initVel(ParticleSource* sourceInfo, ParticleState* particles);

void print_status(ParticleState* state);

int main();

void accel(ParticleInfo* particle, FieldStructure* fields, double* x, double* y, double* z, double* vx, double* vy, double* vz, double* accx, double* accy, double* accz, long N);
