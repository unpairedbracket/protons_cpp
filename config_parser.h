#pragma once

#include <cstdio>
#include <iostream>

#include <yaml-cpp/yaml.h>

#include "structs.h"
#include "fields.h"
#include "fields_cocoon.h"

static YAML::Node config;

void load_config(const std::string& filename);

ParticleInfo* getParticleInfo();
ParticleSource* getSourceInfo(ParticleInfo* particleType);
FieldStructure* getFieldsInfo();
ParticleDetector* getParticleDetector();

