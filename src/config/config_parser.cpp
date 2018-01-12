#include "config_parser.h"

#include "../fields/fields_cocoon.h"
#include "../fields/fields_flash.h"
#include "../interpolation/natural.h"
#include "../interpolation/linear.h"
#include "../util/physical_constants.h"

void load_config(const std::string& filename, const std::string& defaults_filename) {
    config = YAML::LoadFile(filename);
    defaults = YAML::LoadFile(defaults_filename);

    if(!config["relativistic"]) config["relativistic"] = defaults["relativistic"];
    if(!config["particleType"]) config["particleType"] = defaults["particleType"];
    if(!config["source"]) config["source"] = defaults["source"];
    if(!config["field"]) config["field"] = defaults["field"];
    if(!config["integrator"]) config["integrator"] = defaults["integrator"];
    if(!config["detector"]) config["detector"] = defaults["detector"];
    if(!config["middleware"]) config["middleware"] = defaults["middleware"];
    if(!config["interpolation"]) config["interpolation"] = defaults["interpolation"];

    resolve_defaults(config);

    std::cout << config << std::endl;
}

void resolve_defaults(YAML::Node node) {
    replace_with_file(node, "source");
    replace_with_file(node, "field");
    replace_with_file(node, "integrator");
    replace_with_file(node, "detector");
    replace_with_file(node, "interpolation");
}

void replace_with_file(YAML::Node node, std::string prop) {
    std::string newfile;
    YAML::Node newnode;
    if(node[prop].IsScalar()) {
        newfile = "defaults/" + prop + "/" + node[prop].as<std::string>() + ".yml";
    } else if(node[prop].IsMap()) {
        newfile = "defaults/" + prop + "/" + node[prop]["type"].as<std::string>() + ".yml";
    }
    newnode = YAML::LoadFile(newfile);

    if(node[prop].IsMap()) {
        newnode = YAML::LoadFile(newfile);
        for(auto setting : node[prop]) {
            newnode[setting.first.as<std::string>()] = setting.second;
        }
    }
    node[prop] = newnode;
}

ParticleInfo* getParticleInfo() {
    YAML::Node particleNode;
    double mass, charge;
    if(!config["particleType"]) {
        std::cout << "No particle type specified. Defaulting to protons." << std::endl;
        config["particleType"] = "proton";
    }

    particleNode = config["particleType"];

    if(!(particleNode.IsScalar() || particleNode.IsMap())) {
        std::cout << "Particle specified as an unsupported type. Defaulting to protons" << std::endl;
        config["particleType"] = "proton";
    }

    if(particleNode.IsScalar()) {
       //We have a named particle
        std::string name = particleNode.as<std::string>();
        if(name.compare("proton") == 0) {
            mass = m_p;
            charge = + e;
        } else if(name.compare("electron") == 0) {
            mass = m_e;
            charge = - e;
        } else {
            std::cout << "Particle " << name << " not implemented. Using protons as default." << std::endl;
            mass = m_p;
            charge = + e;
        }
    }

    if(particleNode.IsMap()) {
        if(!particleNode["mass"]) {
            std::cout << "Particle mass not specified. Using proton mass as default." << std::endl;
            particleNode["mass"] = m_p;
        }
        if(!particleNode["charge"]) {
            std::cout << "Particle charge not specified. Using elemental charge as default." << std::endl;
            particleNode["charge"] = e;
        }
        mass = particleNode["mass"].as<double>();
        charge = particleNode["charge"].as<double>();
    }

    ParticleInfo* particle = new ParticleInfo({mass, charge});
    initParticle(particle);
    return particle;
}

ParticleSource* getSourceInfo() {
    return getSourceInfo(config["source"]);
}

ParticleSource* getSourceInfo(YAML::Node sourceNode) {
    ParticleSource* source = nullptr;

    if(sourceNode.IsMap()) {
        std::string sourceType = sourceNode["type"].as<std::string>();

        if(sourceType.compare("square") == 0) {
            SquareSource* squareSource = new SquareSource();
            source = squareSource;

            squareSource->x_extent = sourceNode["N"][0].as<long>();
            squareSource->y_extent = sourceNode["N"][1].as<long>();
            squareSource->distance = sourceNode["distance"].as<double>();
            squareSource->divergence = sourceNode["divergence"].as<double>();
            squareSource->energy = sourceNode["energy"].as<double>();
        }
        if(sourceType.compare("helix") == 0) {
            if(!sourceNode["pitch"]) {
                std::cout << "No source pitch specified. Defaulting to golden ratio * 2 pi" << std::endl;
                sourceNode["pitch"] = (1 + sqrt(5)) * pi();
            }

            HelixSource* helixSource = new HelixSource();
            source = helixSource;

            helixSource->N = sourceNode["N"].as<long>();
            helixSource->distance = sourceNode["distance"].as<double>();
            helixSource->divergence = sourceNode["divergence"].as<double>();
            helixSource->energy = sourceNode["energy"].as<double>();
            helixSource->dphi = sourceNode["pitch"].as<double>();
        }
        if(sourceType.compare("scatter") == 0) {
            ScatterSource* scatterSource = new ScatterSource();
            source = scatterSource;

            scatterSource->N = sourceNode["N"].as<long>();
            scatterSource->distance = sourceNode["distance"].as<double>();
            scatterSource->divergence = sourceNode["divergence"].as<double>();
            scatterSource->energy = sourceNode["energy"].as<double>();
            scatterSource->phirand = std::uniform_real_distribution<double>(-pi(), pi());
            scatterSource->zrand   = std::uniform_real_distribution<double>(cos(scatterSource->divergence), 1);
        }
    }

    return source;
}

FieldStructure* getFieldsInfo() {
    FieldStructure* field = nullptr;
    YAML::Node fieldNode;

    fieldNode = config["field"];

    if(fieldNode.IsMap()) {
        if(!fieldNode["type"]) {
            std::cout << "No field type specified. Defaulting to cocoon." << std::endl;
            fieldNode["type"] = "cocoon";
        }

        std::string fieldType = fieldNode["type"].as<std::string>();

        if(!(fieldType.compare("cocoon") == 0
          || fieldType.compare("flash") == 0
          )) {
            std::cout << "Unsupported field type '" << fieldType << "'. Defaulting to cocoon." << std::endl;
            fieldType = std::string("cocoon");
        }

        if(fieldType.compare("cocoon") == 0) {
            CocoonField* cocoonField = new CocoonField();
            field = cocoonField;

            cocoonField->r_scale = fieldNode["radial_scale"].as<double>();
            cocoonField->z_scale = fieldNode["length_scale"].as<double>();
            cocoonField->B_strength = fieldNode["B_scale"].as<double>();
        }

        if(fieldType.compare("flash") == 0) {
            FlashField* flashField = new FlashField();
            field = flashField;

            if(!fieldNode["origin"].IsSequence() || fieldNode["origin"].size() != 3) {
                std::cout << "FLASH Field origin must be a 3d vector. Defaulting to [0, 0, 0]." << std::endl;
                fieldNode["origin"] = YAML::Load("[0, 0, 0]");
            }

            flashField->filename = fieldNode["filename"].as<std::string>();
            flashField->origin = {fieldNode["origin"][0].as<double>(), fieldNode["origin"][1].as<double>(), fieldNode["origin"][2].as<double>()};
        }

        if(!fieldNode["axis"].IsSequence() || fieldNode["axis"].size() != 3) {
            std::cout << "Field principal axis must be a 3d vector. Defaulting to [0, 0, 1]." << std::endl;
            fieldNode["axis"] = YAML::Load("[0, 0, 1]");
        }

        field->zaxis = {fieldNode["axis"][0].as<double>(), fieldNode["axis"][1].as<double>(), fieldNode["axis"][2].as<double>()};
        field->theta = fieldNode["theta"].as<double>() * pi() / 180.0;
        field->phi = fieldNode["phi"].as<double>();
    }

    return field;
}

ParticleDetector* getDetectorInfo() {
    ParticleDetector* detector = nullptr;
    YAML::Node detectorNode = config["detector"];

    std::string detectorType = detectorNode["type"].as<std::string>();

    if(detectorType.compare("none") == 0) {
        detector = new DetectorNoop();
    } else if(detectorType.compare("text") == 0) {
        detector = new DetectorTextFile();
    } else if(detectorType.compare("hdf5") == 0) {
        detector = new DetectorHDF5();
    } else if(detectorType.compare("fluence") == 0) {
        DetectorFluence* fluence = new DetectorFluence();
        fluence->detectorPixels[0] = detectorNode["pixels"][0].as<int>();
        fluence->detectorPixels[1] = detectorNode["pixels"][1].as<int>();
        fluence->detectorSize[0] = detectorNode["size"][0].as<double>();
        fluence->detectorSize[1] = detectorNode["size"][1].as<double>();
        fluence->detectorArray = new double[fluence->detectorPixels[0]*fluence->detectorPixels[1]];
        detector = fluence;
    }

    detector->distance = detectorNode["distance"].as<double>();

    return detector;
}

Integrator* getIntegratorInfo() {
    Integrator* integrator = nullptr;
    YAML::Node integratorNode;

    if(!config["integrator"]) {
        std::cout << "No integrator specified. Defaulting to RK4." << std::endl;
        config["integrator"] = "RK4";
    }

    integratorNode = config["integrator"];

    if(integratorNode.IsScalar()) {
        std::cout << "Integrator only specified by name. Will use defaults for that type" << std::endl;
        YAML::Node parentNode;
        parentNode["type"] = integratorNode;
        integratorNode = parentNode;
    }

    if(integratorNode.IsMap()) {
        if(!integratorNode["type"]) {
            std::cout << "No integrator type specified. Defaulting to RK4." << std::endl;
            integratorNode["type"] = "RK4";
        }

        std::string integratorType = integratorNode["type"].as<std::string>();

        if(!(integratorType.compare("euler") == 0
          || integratorType.compare("Euler") == 0
          || integratorType.compare("rk4") == 0
          || integratorType.compare("RK4") == 0
          || integratorType.compare("RKDP") == 0
          )) {
            std::cout << "Unsupported integrator type '" << integratorType << "'. Defaulting to RK4." << std::endl;
            integratorType = std::string("RK4");
        }

        if(integratorType.compare("euler") == 0 || integratorType.compare("Euler") == 0) {
            integrator = new EulerIntegrator();
        } else if (integratorType.compare("rk4") == 0 || integratorType.compare("RK4") == 0) {
            integrator = new RK4Integrator();
        } else if (integratorType.compare("RKDP") == 0) {
            RKDPIntegrator* rkdp = new RKDPIntegrator();
            integrator = rkdp;
            rkdp->dt_min = integratorNode["dt_min"].as<double>();
            rkdp->rtol = integratorNode["tol"].as<double>();
            rkdp->maxLengthen = integratorNode["max_lengthen"].as<double>();
            rkdp->maxFirstShorten = integratorNode["max_first_shorten"].as<double>();
            rkdp->maxOtherShorten = integratorNode["max_shorten"].as<double>();
            rkdp->kill_failed_particles = integratorNode["remove_failed_particles"].as<bool>();
            rkdp->verbose = integratorNode["verbose"].as<bool>();
        }

        if(!integratorNode["relativistic"]) {
            if(config["relativistic"]) {
                std::cout << "Integrator inheriting relativistic preference from global setting." << std::endl;
                integratorNode["relativistic"] = config["relativistic"];
            } else {
                std::cout << "Defaulting to non-relativistic algorithm." << std::endl;
                integratorNode["relativistic"] = false;
            }
        }

        if(!integratorNode["time_step"]) {
            std::cout << "No time step specified." << std::endl;
            integratorNode["time_step"] = 1E-10;
        }

        integrator->setInitTimestep(integratorNode["time_step"].as<double>());
        integrator->setRelativistic(integratorNode["relativistic"].as<bool>());
    }

    return integrator;

}

Interpolator* getInterpolatorInfo() {
    Interpolator* interpolator = nullptr;
    YAML::Node interpolatorNode = config["interpolation"];

    if(interpolatorNode.IsMap()) {
        std::string interpolatorType = interpolatorNode["type"].as<std::string>();

        if (interpolatorType.compare("linear") == 0) {
            interpolator = new LinearInterpolator();

            replace_with_file(interpolatorNode, "source");
            interpolator->interpSource = getSourceInfo(interpolatorNode["source"]);
            interpolator->iterations = interpolatorNode["iterations"].as<int>();
        } else if (interpolatorType.compare("natural_neighbor") == 0) {
            interpolator = new NaturalInterpolator();

            replace_with_file(interpolatorNode, "source");
            interpolator->interpSource = getSourceInfo(interpolatorNode["source"]);
            interpolator->iterations = interpolatorNode["iterations"].as<int>();
        }
    }

    return interpolator;
}
