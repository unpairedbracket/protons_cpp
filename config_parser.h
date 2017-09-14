#pragma once

#include <cstdio>
#include <iostream>

#include <yaml-cpp/yaml.h>

#include "math.h"
#include "source.h"
#include "structs.h"
#include "fields.h"
#include "fields_cocoon.h"
#include "integrator.h"

static YAML::Node config;

void load_config(const std::string& filename);

ParticleInfo* getParticleInfo();
ParticleSource* getSourceInfo();
FieldStructure* getFieldsInfo();
ParticleDetector* getParticleDetector();
Integrator* getIntegrator();

