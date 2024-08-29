#ifndef __DRIVER_CXX__
#define __DRIVER_CXX__

#include "Driver.h"
#include "supera/algorithm/BBoxInteraction.h"
#include "supera/algorithm/LArTPCMLReco3D.h"

namespace supera {

    std::string Driver::DumpConfig(const YAML::Node& cfg)
    {
        std::stringstream ss;
        ss << "Dump of the configurations (IsDefined, IsNull) ... (IsScalar,IsList,IsMap)"
        << std::endl;

        for (auto const &item: cfg) {
            ss << item.first.as<std::string>() << " : ("
            << item.second.IsDefined() << ","
            << item.second.IsNull() << ") ... ("
            << item.second.IsScalar() << ","
            << item.second.IsSequence() << ","
            << item.second.IsMap() << ")" 
            << std::endl;
        }

        return ss.str();
    }

    void Driver::Configure(const YAML::Node& cfg) 
    {
        LOG_INFO() << "starting"<<std::endl;
        if(cfg["DriverConfig"]) {
            auto driver_cfg = cfg["DriverConfig"];
            if(!driver_cfg.IsMap()) {
                throw meatloaf("DriverConfig must be a dictionary type");
            }
            if(driver_cfg["LogLevel"]) {
                auto log_level = cfg["LogLevel"].as<std::string>();
                this->SetLogConfig(supera::msg::parseStringThresh(log_level));
                LOG_INFO() << "Set the logger level to " << log_level << std::endl;
            }
            _assert_inout_voxel_count = cfg["AssertInOutVoxelCount"].as<bool>();
        }

        LOG_WARNING()<<"BBox config..."<<std::endl;
        if(cfg["BBoxAlgorithm"] && cfg["BBoxConfig"] && cfg["BBoxConfig"].IsMap())
        {
            std::string name = cfg["BBoxAlgorithm"].as<std::string>();
            if(name == "BBoxInteraction") {
                _algo_bbox = new BBoxInteraction();
                _algo_bbox->Configure(cfg["BBoxConfig"]);
            }
            else{
                std::string msg = name + " is not known to Supera...";
                throw meatloaf(msg);
            }
        }
        else
        {
            LOG_FATAL() << "BBoxAlgorithm or BBoxConfig not found or invalid in the configuration."
            << std::endl
            << this->DumpConfig(cfg)
            << std::endl;
            throw meatloaf("Failed to configure");
        }

        LOG_WARNING()<<"Label config..."<<std::endl;
        if(cfg["LabelAlgorithm"] && cfg["LabelConfig"] && cfg["LabelConfig"].IsMap())
        {
            std::string name = cfg["LabelAlgorithm"].as<std::string>();
            if(name == "LArTPCMLReco3D") {
                _algo_label = new LArTPCMLReco3D();
                _algo_label->Configure(cfg["LabelConfig"]);
            }
            else{
                std::string msg = name + " is not known to Supera...";
                throw meatloaf(msg);
            }
        }
        else
        {
            LOG_FATAL() << "LabelAlgorithm or LabelConfig not found or invalid in the configuration."
            << std::endl
            << this->DumpConfig(cfg)
            << std::endl;
            throw meatloaf("Failed to configure");
        }
    }

    /*
    void Driver::ConfigureBBoxAlgorithm(const std::string& name,
        const std::map<std::string,std::string>& params)
    {
        if(name == "BBoxInteraction") {
            _algo_bbox = new BBoxInteraction();
            PSet cfg;
            cfg.data = params;
            _algo_bbox->Configure(cfg);

            if (cfg.exists("LogLevel") > 0)
                _algo_bbox->SetLogConfig(supera::msg::parseStringThresh(cfg.get<std::string>("LogLevel")));
        }
        else{
            std::string msg = name + " is not known to Supera...";
            throw meatloaf(msg);
        }

    }

    void Driver::ConfigureLabelAlgorithm(const std::string& name,
        const std::map<std::string,std::string>& params)
    {
        if(name == "LArTPCMLReco3D") {
            _algo_label = new LArTPCMLReco3D();
            PSet cfg;
            cfg.data = params;
            _algo_label->Configure(cfg);

            if (cfg.exists("LogLevel") > 0)
                _algo_label->SetLogConfig(supera::msg::parseStringThresh(cfg.get<std::string>("LogLevel")));

        }
        else{
            std::string msg = name + " is not known to Supera...";
            throw meatloaf(msg);
        }
    }
    */

    void Driver::GenerateImageMeta(const EventInput& data)
    {
        if(!_algo_bbox) 
            throw meatloaf("BBoxAlgorithm is not configured yet!");

        _meta.clear();
        _label = EventOutput();
        _meta  = _algo_bbox->Generate(data);
    }


    void Driver::GenerateLabel(const EventInput& data)
    {
        if(!_algo_label) 
            throw meatloaf("LabelAlgorithm is not configured yet!");

        if(!_meta.valid())
            throw meatloaf("BBoxAlgorithm must be run first");
        _label = _algo_label->Generate(data, _meta);

        if(_assert_inout_voxel_count) {
        // Make sure all energy depositions are accounted for
            std::vector<supera::EDep> edeps;
            edeps.reserve(_label._energies.size());
            for(auto const& p : data) {
                for(auto const& edep : p.pcloud) {
                    edeps.push_back(edep);
                }
            }
            for(auto const& unass_edeps : data.unassociated_edeps) {
                for(auto const& edep : unass_edeps) {
                    edeps.push_back(edep);
                }
            }

            auto vs = _meta.edep2voxelset(edeps);

            if(vs.size() != _label._energies.size()) {
                LOG_FATAL() << "Unexpected count of output voxels (EventOuptut::_energies) "
                << _label._energies.size()
                << " for the input with "
                << vs.size() 
                << " voxels" << std::endl;
                throw meatloaf("Detected potentially unaccounted input voxels");
            }

            for(size_t i=0; i<vs.size(); ++i) {
                auto const& id_a = vs.as_vector()[i].id();
                auto const& id_b = _label._energies.as_vector()[i].id();

                if(id_a == id_b)
                    continue;

                if(id_a < id_b)
                    LOG_FATAL() << "Detected the input voxel ID " << id_a << " is missing in the output voxel sets." << std::endl;
                else
                    LOG_FATAL() << "Detected the output voxel ID " << id_b << " did not exist in the input voxel set." << std::endl;
                throw meatloaf("Detected inconsistent voxel sets between the input and output of Supera");
            }

        }
    }

    void Driver::Generate(const EventInput& data)
    {
        this->Reset();
        this->GenerateImageMeta(data);
        this->GenerateLabel(data);
    }

}
#endif
