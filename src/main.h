#pragma once

#include "particles/structs.h"
#include "fields/fields.h"
#include "sources/source.h"

void initPos(ParticleSource* sourceInfo, ParticleState* particles);

void initVel(ParticleSource* sourceInfo, ParticleState* particles);

void print_status(ParticleState* state);

int main(int argc, char *argv[]);

void accel(ParticleState* state, FieldStructure* fields);

void invalidateStates(ParticleState* state);
