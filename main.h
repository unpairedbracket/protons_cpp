#pragma once

#include <cmath>
#include <cstdio>
#include <cassert>

#include <chrono>
#include <fstream>
#include <iostream>

#include <omp.h>

#include "structs.h"
#include "integrator.h"
#include "fields_cocoon.h"

#include "graphics/window.h"

constexpr double pi() { return std::atan(1)*4; }
const bool USE_GL = true;

void initPos(ParticleSource* sourceInfo, ParticleState* particles);

void initVel(ParticleSource* sourceInfo, ParticleState* particles);

void print_status(ParticleState* state);
void print_status_raw(ParticleState* state);

int main(int argc, char *argv[]);

void accel(ParticleState* state, FieldStructure* fields, Vector3* acc);

void invalidateStates(ParticleState* state);
