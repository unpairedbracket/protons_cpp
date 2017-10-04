#pragma once

#include "../particles/structs.h"
#include "../fields/fields.h"

void accel(ParticleState* state, FieldStructure* fields);
void accel_relativistic(ParticleState* state, FieldStructure* fields);
