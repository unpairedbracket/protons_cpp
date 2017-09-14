#pragma once

#include <cstdio>
#include <cassert>

#include <chrono>
#include <fstream>
#include <iostream>

#include <omp.h>

#include "math.h"
#include "structs.h"
#include "integrator.h"
#include "fields_cocoon.h"
#include "config_parser.h"

#include "graphics/window.h"

const bool USE_GL = true;

void initPos(ParticleSource* sourceInfo, ParticleState* particles);

void initVel(ParticleSource* sourceInfo, ParticleState* particles);

void print_status(ParticleState* state);
void print_status_raw(ParticleState* state);

int main(int argc, char *argv[]);

void accel(ParticleState* state, FieldStructure* fields);

void invalidateStates(ParticleState* state);
