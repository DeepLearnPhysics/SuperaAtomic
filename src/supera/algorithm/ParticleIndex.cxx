#ifndef __PARTICLEINDEX_CXX__
#define __PARTICLEINDEX_CXX__

#include "ParticleIndex.h"
#include "supera/base/meatloaf.h"
#include "supera/data/Event.h"

namespace supera{ 

  void ParticleIndex::_configure(const YAML::Node& cfg)
  {
        for (auto const &item: cfg) {
            std::cout << item.first.as<std::string>() << std::endl
            << item.second.IsDefined() << " "
            << item.second.IsNull() << " ... "
            << item.second.IsScalar() << " "
            << item.second.IsSequence() << " "
            << item.second.IsMap()
            << std::endl;
        } 
  }

  void ParticleIndex::InferParentage(const EventInput& larmcp_v)
  {
    /*
      Assumptions
        - self and parent id are valid for all particles
        - for primary particles, its parent and own id should be identical
    */
    _pdgcode_v.resize(larmcp_v.size());
    _parent_id_v.resize(larmcp_v.size());
    _parent_pdg_v.resize(larmcp_v.size());
    _ancestor_id_v.resize(larmcp_v.size());
    _ancestor_pdg_v.resize(larmcp_v.size());
    _parent_history_v.resize(larmcp_v.size());

    for(size_t i=0; i<larmcp_v.size(); ++i) {
      _pdgcode_v[i] = _parent_pdg_v[i] = _ancestor_pdg_v[i] = supera::kINVALID_PDG;
      _parent_id_v[i] = _ancestor_id_v[i] = supera::kINVALID_INSTANCEID;
    }

    for(auto& history : _parent_history_v) history.clear();

    // fill in the ParticleIndex's working structures.

    // first: create the mapping between GEANT4 trackid <-> index in the particle array
    // (FYI: the larmcp objects correspond each to one GEANT4 particle
    //       + any associated energy deposits)
    for(size_t index=0; index<larmcp_v.size(); ++index) {
      auto const& larmcp = larmcp_v[index];  // pull off the GEANT4 track information component
      if(larmcp.id == supera::kINVALID_INSTANCEID) {
        LOG_FATAL() << "Instance ID cannot be invalid\n";
        throw supera::meatloaf();
      }
      _pdgcode_v[index] = abs(larmcp.part.pdg);

      // parents
      _parent_id_v[index] = larmcp.parent_id;
      _parent_pdg_v[index] = larmcp_v.at(larmcp.parent_id).part.pdg;

      // ancestors
      InstanceID_t current_id = larmcp.id;
      InstanceID_t parent_id  = larmcp.parent_id;
      while(parent_id != current_id) {
        _parent_history_v[index].push_back(parent_id);
        current_id = parent_id;
        parent_id  = larmcp_v.at(parent_id).parent_id;
      }

      auto const& ancestor_id = larmcp_v.at(parent_id).ancestor_id;
      if(ancestor_id != kINVALID_INSTANCEID && ancestor_id != parent_id) {
        LOG_FATAL() << "For the input particle ID=" << larmcp.id
        << " and ancestor ID=" << ancestor_id 
        << " the parentage trace yielded the ancestor ID=" << parent_id
        << std::endl;
        throw supera::meatloaf();
      }
      _ancestor_id_v[index] = parent_id;
      _ancestor_pdg_v[index] = larmcp_v.at(parent_id).part.pdg;

    }
  }

  const std::vector< InstanceID_t >& 
  ParticleIndex::ParentIdArray(const InstanceID_t id) const
  {
    if(id >= _parent_history_v.size()) {
      LOG_ERROR() << "Instance ID " << id << " is not valid. "
      << "Registered range 0 => " << _pdgcode_v.size() 
      << " ... Returning an empty list.\n";
      return _empty_trackid_v;
    }
    return _parent_history_v[id];
  }


}
#endif
