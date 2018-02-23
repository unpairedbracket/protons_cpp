#pragma once

#include <cstdio>
#include <iostream>

#include <yaml-cpp/yaml.h>

#include "../util/math.h"
#include "../sources/source.h"
#include "../particles/structs.h"
#include "../fields/fields.h"
#include "../integrators/integrator.h"
#include "../detectors/detector.h"
#include "../interpolation/interpolator.h"

static YAML::Node config;
static YAML::Node defaults;

void load_config(const std::string& filename, const std::string& defaults_filename);

void resolve_defaults(YAML::Node node);
void replace_with_file(YAML::Node node, std::string prop);

int getNumberRuns();
ParticleInfo getParticleInfo();
ParticleSource* getSourceInfo();
FieldStructure* getFieldsInfo();
ParticleDetector* getDetectorInfo();
Integrator* getIntegratorInfo();
Interpolator* getInterpolatorInfo();

