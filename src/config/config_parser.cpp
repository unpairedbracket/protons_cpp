#include "config_parser.h"

#include "../fields/fields_cocoon.h"
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

    resolve_defaults(config);

    std::cout << config << std::endl;
}

void resolve_defaults(YAML::Node node) {
    replace_with_file(node, "source");
    replace_with_file(node, "field");
    replace_with_file(node, "integrator");
    replace_with_file(node, "detector");
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
    
    ParticleSource* source = nullptr;
    YAML::Node sourceNode;
    
    if(!config["source"]) {
        std::cout << "No source type specified. Defaulting to helix." << std::endl;
        config["source"] = "helix";
    }

    sourceNode = config["source"];

    if(sourceNode.IsScalar()) {
        std::cout << "Source type only specified by name. Will use defaults for that type" << std::endl;
        YAML::Node parentNode;
        parentNode["type"] = sourceNode;
        sourceNode = parentNode;
    }

    if(sourceNode.IsMap()) {
        if(!sourceNode["type"]) {
            std::cout << "No source type specified. Defaulting to cocoon." << std::endl;
            sourceNode["type"] = "cocoon";
        }

        std::string sourceType = sourceNode["type"].as<std::string>();

        if(!(sourceType.compare("square") == 0
          || sourceType.compare("helix") == 0
          )) {
            std::cout << "Unsupported source type '" << sourceType << "'. Defaulting to helix." << std::endl;
            sourceType = std::string("helix");
        }

        if(sourceType.compare("square") == 0) {
            if(!sourceNode["N"]) {
                std::cout << "No number of particles specified for source. Defaulting to 100x100" << std::endl;
                sourceNode["N"] = YAML::Load("[100, 100]");
            }
            if(!sourceNode["N"].IsSequence() || sourceNode["N"].size() != 2) {
                std::cout << "Wrong number of components of particle number specified: " << sourceNode["N"].size() << ". Expected 2: [x_number, y_number], defaulting to [100, 100]." << std::endl;
                sourceNode["N"] = YAML::Load("[100, 100]");
            }

            if(!sourceNode["distance"]) {
                std::cout << "No distance specified for source. Defaulting to 0.01" << std::endl;
                sourceNode["distance"] = 0.01;
            }

            if(!sourceNode["divergence"]) {
                std::cout << "No source divergence specified. Defaulting to 0.01 radian" << std::endl;
                sourceNode["divergence"] = 0.01;
            }

            if(!sourceNode["energy"]) {
                std::cout << "No source energy specified. Defaulting to 1 MeV" << std::endl;
                sourceNode["energy"] = 1E6 * 1.6E-19;
            }

            SquareSource* squareSource = new SquareSource();
            source = squareSource;
 
            squareSource->x_extent = sourceNode["N"][0].as<long>();
            squareSource->y_extent = sourceNode["N"][1].as<long>();
            squareSource->distance = sourceNode["distance"].as<double>();
            squareSource->divergence = sourceNode["divergence"].as<double>();
            squareSource->energy = sourceNode["energy"].as<double>();
        }
        if(sourceType.compare("helix") == 0) {
            if(!sourceNode["N"]) {
                std::cout << "No number of particles specified for source. Defaulting to 10000" << std::endl;
                sourceNode["N"] = 10000;
            }

            if(!sourceNode["distance"]) {
                std::cout << "No distance specified for source. Defaulting to 0.01" << std::endl;
                sourceNode["distance"] = 0.01;
            }

            if(!sourceNode["divergence"]) {
                std::cout << "No source divergence specified. Defaulting to 0.01 radian" << std::endl;
                sourceNode["divergence"] = 0.01;
            }

            if(!sourceNode["energy"]) {
                std::cout << "No source energy specified. Defaulting to 1 MeV" << std::endl;
                sourceNode["energy"] = 1E6 * 1.6E-19;
            }

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
    }

    return source;
}

FieldStructure* getFieldsInfo() {
    FieldStructure* field = nullptr;
    YAML::Node fieldNode;
    
    if(!config["field"]) {
        std::cout << "No field type specified. Defaulting to cocoon." << std::endl;
        config["field"] = "cocoon";
    }

    fieldNode = config["field"];

    if(fieldNode.IsScalar()) {
        std::cout << "Field type only specified by name. Will use defaults for that type" << std::endl;
        YAML::Node parentNode;
        parentNode["type"] = fieldNode;
        fieldNode = parentNode;
    }

    if(fieldNode.IsMap()) {
        if(!fieldNode["type"]) {
            std::cout << "No field type specified. Defaulting to cocoon." << std::endl;
            fieldNode["type"] = "cocoon";
        }

        std::string fieldType = fieldNode["type"].as<std::string>();

        if(!(fieldType.compare("cocoon") == 0
          || fieldType.compare("otherthing") == 0
          )) {
            std::cout << "Unsupported field type '" << fieldType << "'. Defaulting to cocoon." << std::endl;
            fieldType = std::string("cocoon");
        }

        if(fieldType.compare("cocoon") == 0) {
            if(!fieldNode["radial_scale"]) {
                std::cout << "No radial scale specified for cocoon. Defaulting to 0.001" << std::endl;
                fieldNode["radial_scale"] = 0.001;
            }

            if(!fieldNode["length_scale"]) {
                std::cout << "No length scale specified for cocoon. Defaulting to 0.0005" << std::endl;
                fieldNode["length_scale"] = 0.0005;
            }

            if(!fieldNode["B_scale"]) {
                std::cout << "No magentic field strength scale specified for cocoon. Defaulting to 1 Tesla" << std::endl;
                fieldNode["B_scale"] = 1;
            }

            CocoonField* cocoonField = new CocoonField();
            field = cocoonField;
 
            cocoonField->r_scale = fieldNode["radial_scale"].as<double>();
            cocoonField->z_scale = fieldNode["length_scale"].as<double>();
            cocoonField->B_strength = fieldNode["B_scale"].as<double>();
        }
    }

    return field;
}

ParticleDetector* getDetectorInfo() {
    ParticleDetector* detector = nullptr;
    YAML::Node detectorNode = config["detector"];

    std::string detectorType = detectorNode["type"].as<std::string>();

    if(detectorType.compare("none") == 0) {
        detector = new DetectorNoop();
        detector->distance = 1;
    } else if(detectorType.compare("text") == 0) {
        detector = new DetectorTextFile();
        detector->distance = detectorNode["distance"].as<double>();
    } else if(detectorType.compare("hdf5") == 0) {
        detector = new DetectorHDF5();
        detector->distance = detectorNode["distance"].as<double>();
    }

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
            integrator = new RKDPIntegrator();
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

