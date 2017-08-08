#include "config_parser.h"

void load_config(const std::string& filename) {
    config = YAML::LoadFile(filename);
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
            mass = 1.6726217E-27;
            charge = 1.6021766208E-19;
        } else if(name.compare("electron") == 0) {
            mass = 9.10938356E-31;
            charge = -1.6021766208E-19;
        } else {
            std::cout << "Particle " << name << " not implemented. Using protons as default." << std::endl;
            mass = 1.6726217E-27;
            charge = 1.6021766208E-19;
        }
    }
    
    if(particleNode.IsMap()) {
        if(!particleNode["mass"]) {
            std::cout << "Particle mass not specified. Using proton mass as default." << std::endl;
            particleNode["mass"] = 1.6726217E-27;
        }
        if(!particleNode["charge"]) {
            std::cout << "Particle charge not specified. Using elemental charge as default." << std::endl;
            particleNode["charge"] = 1.6021766208E-19;
        }
        mass = particleNode["mass"].as<double>();
        charge = particleNode["charge"].as<double>();
    }

    ParticleInfo* particle = new ParticleInfo({mass, charge});
    initParticle(particle);
    return particle;
}

ParticleSource* getSourceInfo(ParticleInfo* particleType);
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

ParticleDetector* getParticleDetector();

