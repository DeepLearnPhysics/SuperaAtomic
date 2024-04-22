#include "Particle.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>
#include "supera/base/meatloaf.h"

namespace supera {


  std::string Particle::dump() const
  {
    std::stringstream ss;
    std::stringstream buf;
    ss  << "      \033[95m" << "Particle " << " (PdgCode,TrackID) = (" << StringifyPdgCode(pdg) << "," << StringifyTrackID(trackid) << ")\033[00m " << std::endl
        << "      ... with Parent (" << StringifyPdgCode(parent_pdg) << "," << StringifyTrackID(parent_trackid) << ")" << std::endl
        << "      ... + Ancestor  (" << StringifyPdgCode(ancestor_pdg) << "," << StringifyTrackID(ancestor_trackid) << ")" << std::endl;

    buf << "      ";

    ss<< buf.str() << "Vertex   (x, y, z, t) = (" << StringifyPoint3D(vtx.pos) << " T=" << vtx.time << std::endl
      << buf.str() << "Momentum (px, py, pz) = (" << StringifyDouble(px) << ","
      << buf.str() << "Momentum (end_px, end_py, end_pz) = (" << end_px << "," << end_py << "," << end_pz << ")" << std::endl
      << StringifyDouble(py) << "," 
      << StringifyDouble(pz) << ")" << std::endl
      << buf.str() << "Initial  Energy  = " << StringifyDouble(energy_init) << std::endl
      << buf.str() << "Deposit  Energy  = " << StringifyDouble(energy_deposit) << std::endl
      << buf.str() << "Creation Process = " << process << std::endl
      << buf.str() << "Instance ID      = " << StringifyInstanceID(id) << std::endl
      << buf.str() << "Group ID         = " << StringifyInstanceID(group_id) << std::endl
      << buf.str() << "Interaction ID   = " << StringifyInstanceID(interaction_id) << std::endl
      << buf.str() << "Type     = " << type << std::endl
      << buf.str() << "Shape    = " << shape << std::endl;
    ss << buf.str() << "Children = ";
    for (const auto & child :  children_id)
      ss << " " << StringifyInstanceID(child);
    ss << std::endl;

    return ss.str();

  }


  std::string Particle::dump2cpp(const std::string &instanceName) const
  {
    std::stringstream ss;
    
    ss << "supera::Particle " << instanceName << ";\n";

    ss << instanceName << ".id = " << id << ";\n";
    ss << instanceName << ".shape = static_cast<supera::SemanticType_t>(" << shape << ");\n";
    ss << instanceName << ".trackid = " << StringifyTrackID(trackid) << ";\n";
    ss << instanceName << ".pdg = " << StringifyPdgCode(pdg) << ";\n";
    ss << instanceName << ".px = " << StringifyDouble(px) << ";\n";
    ss << instanceName << ".py = " << StringifyDouble(py) << ";\n";
    ss << instanceName << ".pz = " << StringifyDouble(pz) << ";\n";
    ss << instanceName << ".end_px = " << end_px << ";\n";
    ss << instanceName << ".end_py = " << end_py << ";\n";
    ss << instanceName << ".end_pz = " << end_pz << ";\n";
    ss << instanceName << ".vtx = {" << StringifyDouble(vtx.pos.x) << ", "
                                     << StringifyDouble(vtx.pos.y) << ", "
                                     << StringifyDouble(vtx.pos.z) << ", "
                                     << StringifyDouble(vtx.time ) << "};\n";
    ss << instanceName << ".end_pt = {" << StringifyDouble(end_pt.pos.x) << ", "
                                        << StringifyDouble(end_pt.pos.y) << ", "
                                        << StringifyDouble(end_pt.pos.z) << ", "
                                        << StringifyDouble(end_pt.time ) << "};\n";
    ss << instanceName << ".first_step = {" << StringifyDouble(first_step.pos.x) << ", "
                                            << StringifyDouble(first_step.pos.y) << ", "
                                            << StringifyDouble(first_step.pos.z) << ", "
                                            << StringifyDouble(first_step.time ) << "};\n";
    ss << instanceName << ".last_step = {" << StringifyDouble(last_step.pos.x) << ", "
                                           << StringifyDouble(last_step.pos.y) << ", "
                                           << StringifyDouble(last_step.pos.z) << ", "
                                           << StringifyDouble(last_step.time ) << "};\n";
    ss << instanceName << ".dist_travel = " << StringifyDouble(dist_travel) << ";\n";
    ss << instanceName << ".energy_init = " << StringifyDouble(energy_init) << ";\n";
    ss << instanceName << ".energy_deposit = " << StringifyDouble(energy_deposit) << ";\n";
    ss << instanceName << ".process = \"" << process << "\";\n";

    ss << instanceName << ".parent_trackid = " << StringifyTrackID(parent_trackid) << ";\n";
    ss << instanceName << ".parent_pdg = " << StringifyPdgCode(parent_pdg) << ";\n";
    ss << instanceName << ".parent_vtx = {" << StringifyDouble(parent_vtx.pos.x) << ", "
                                            << StringifyDouble(parent_vtx.pos.y) << ", "
                                            << StringifyDouble(parent_vtx.pos.z) << ", "
                                            << StringifyDouble(parent_vtx.time ) << "};\n";

    // this particle is at the top so it's its own ancestor
    ss << instanceName << ".ancestor_trackid = " << StringifyTrackID(ancestor_trackid) << ";\n";
    ss << instanceName << ".ancestor_pdg = " << StringifyPdgCode(ancestor_pdg) << ";\n";
    ss << instanceName << ".ancestor_vtx = {" << StringifyDouble(ancestor_vtx.pos.x) << ", "
                                              << StringifyDouble(ancestor_vtx.pos.y) << ", "
                                              << StringifyDouble(ancestor_vtx.pos.z) << ", "
                                              << StringifyDouble(ancestor_vtx.time ) << "};\n";
    ss << instanceName << ".ancestor_process = \"" << ancestor_process << "\";\n";

    ss << instanceName << ".parent_process = \"" << parent_process << "\";\n";
    ss << instanceName << ".parent_id = " << parent_id << ";\n";

    ss << instanceName << ".children_id = { ";
    for (const supera::InstanceID_t & chid : children_id)
      ss << StringifyInstanceID(chid) << (chid != children_id.back() ? ", " : "");
    ss << " };\n";

    ss << instanceName << ".group_id = " << StringifyInstanceID(group_id) << ";\n";
    ss << instanceName << ".interaction_id = " << StringifyInstanceID(interaction_id) << ";\n";

    return ss.str();
  }


  bool Particle::operator==(const Particle &rhs) const
  {
    // Particles are the same if their data is the same...
    return (id == rhs.id) &&
           (shape == rhs.shape) &&
           (trackid == rhs.trackid) &&
           (pdg == rhs.pdg) &&
           (px == rhs.px) &&
           (py == rhs.py) &&
           (pz == rhs.pz) &&
           (end_px == rhs.end_px) &&
           (end_py == rhs.end_py) &&
           (end_pz == rhs.end_pz) &&
           (vtx == rhs.vtx) &&
           (end_pt == rhs.end_pt) &&
           (first_step == rhs.first_step) &&
           (last_step == rhs.last_step) &&
           (dist_travel == rhs.dist_travel) &&
           (energy_init == rhs.energy_init) &&
           (energy_deposit == rhs.energy_deposit) &&
           (process == rhs.process) &&
           (parent_trackid == rhs.parent_trackid) &&
           (parent_pdg == rhs.parent_pdg) &&
           (parent_vtx == rhs.parent_vtx) &&
           (ancestor_trackid == rhs.ancestor_trackid) &&
           (ancestor_pdg == rhs.ancestor_pdg) &&
           (ancestor_vtx == rhs.ancestor_vtx) &&
           (ancestor_process == rhs.ancestor_process) &&
           (parent_process == rhs.parent_process) &&
           (parent_id == rhs.parent_id) &&
           (children_id == rhs.children_id) &&
           (group_id == rhs.group_id) &&
           (interaction_id == rhs.interaction_id);
  }

  // --------------------------------------------------------

  std::string ParticleInput::dump2cpp(const std::string &instanceName) const
  {
    std::stringstream ss;

    ss << "supera::ParticleInput " << instanceName << ";\n";

    ss << instanceName << ".pcloud.reserve(" << pcloud.size() << ");\n";
    for (std::size_t idx = 0; idx < pcloud.size(); idx++)
    {
      std::string edepInstance = instanceName + "_edep" + std::to_string(idx);
      std::string edepCode = pcloud[idx].dump2cpp(edepInstance);
      ss << edepCode;
      ss << instanceName << ".pcloud.emplace_back(std::move(" << edepInstance << "));\n";
    }

    ss << instanceName << ".valid = " << valid << ";\n";

    std::string partInstance = instanceName + "_particle";
    ss << part.dump2cpp(partInstance);
    ss << instanceName << ".part = std::move(" << partInstance << ");\n";

    return ss.str();
  }


  // --------------------------------------------------------

  ParticleLabel::ParticleLabel()
  : valid(false)
  , merge_id(supera::kINVALID_TRACKID)
  {}

  void ParticleLabel::UpdateFirstPoint(const EDep& pt)
  { 
    if(pt.x == supera::kINVALID_DOUBLE) return; 
    if(pt.t < first_pt.t || first_pt.t == supera::kINVALID_DOUBLE ) 
      { 
        /*
        std::cout << "Updating track " << part.trackid << " first pt " << first_pt.t << " => " << pt.t << std::endl;
        std::cout << "    " << first_pt.x << "," << first_pt.y << "," << first_pt.z 
        << " => "
        << pt.x << "," << pt.y << "," << pt.z 
        << std::endl;
        */
        first_pt = pt; 
      } 
  }

  void ParticleLabel::UpdateLastPoint(const EDep& pt)
  { 
    if(pt.x == supera::kINVALID_DOUBLE) return; 
    if(pt.t > last_pt.t  || last_pt.t  == supera::kINVALID_DOUBLE ) 
      { 
        /*
        std::cout << "Updating track " << part.trackid << " last  pt " << last_pt.t << " => " << pt.t << std::endl;
        std::cout << "    " << last_pt.x << "," << last_pt.y << "," << last_pt.z 
        << " => "
        << pt.x << "," << pt.y << "," << pt.z 
        << std::endl;
        */
        last_pt  = pt; 
      }
  }

  void ParticleLabel::SizeCheck() const
  {

    if(dedx.size() && energy.size() != dedx.size()) {
      std::cerr << "Size mismatch: " << energy.size() << " v.s. " << dedx.size() << std::endl;
      throw std::exception();
    }

  }


  size_t ParticleLabel::Size() const
  { 
    SizeCheck();
    return energy.size();
  }

  void ParticleLabel::Merge(ParticleLabel& child) {
    //std::cout<<"Merging child track " << child.part.trackid << " into parent track " << part.trackid << std::endl;
    if(!(this->valid)) {
      std::cerr<<"Cannot merge into an invalid parent!" << std::endl;
      std::cerr<<"Parent info..."<<std::endl<<this->dump()<<std::endl;
      std::cerr<<"Child info..."<<std::endl<<child.dump()<<std::endl;
      throw meatloaf();
    }
    
    for(auto const& vox : child.energy.as_vector())
      this->energy.emplace(vox.id(),vox.value(),true);
    for(auto const& vox : child.dedx.as_vector())
      this->dedx.emplace(vox.id(),vox.value(),true);
    

    if(child.part.shape == kShapeTrack || child.part.shape == kShapeShower || first_pt.x == kINVALID_DOUBLE)
      this->UpdateFirstPoint(child.first_pt);
    if(child.part.shape == kShapeTrack || first_pt.x == kINVALID_DOUBLE) {
      //std::cout<<child.dump2cpp()<<std::endl;
      this->UpdateLastPoint(child.last_pt);
    }
    this->merged_v.push_back(child.part.trackid);
    for(auto const& trackid : child.merged_v)
      this->merged_v.push_back(trackid);

    child.energy.clear_data();
    child.dedx.clear_data();
    child.valid=false;
    child.merge_id=this->part.trackid;
  }

  // semantic classification (supera::SemanticType_t)
  
  /*
  supera::SemanticType_t ParticleLabel::shape() const
  {

    if(type == kInvalidProcess) return supera::kShapeUnknown;
    if(type == kDelta) return supera::kShapeDelta;
    if(type == kNeutron) //return supera::kShapeUnknown;

    if(part.pdg == 11 || part.pdg == 22 || part.pdg == -11) {
      if(type == kCompton || type == kPhoton || type == kPrimary || type == kConversion || type==kOtherShower)
        return supera::kShapeShower;
      if(type == kDecay) {
        if(part.parent_pdg == 13 || part.parent_pdg == -13)
          return supera::kShapeMichel;
        else
          return supera::kShapeShower;
      }
      return supera::kShapeLEScatter;
    }else
    return supera::kShapeTrack;
  }
*/

  std::string ParticleLabel::dump2cpp(const std::string &instanceName) const
  {
    std::stringstream ss;

    ss << "supera::ParticleLabel " << instanceName << ";\n";

    std::string partInstance = instanceName + "_part";
    ss << part.dump2cpp(partInstance);
    ss << instanceName << ".part = std::move(" << partInstance << ");\n";
    ss << instanceName << ".valid = " << valid << ";\n";
    ss << instanceName << ".merged_v = { ";
    for (std::size_t tkid : merged_v)
      ss << tkid << (tkid == merged_v.back() ? "" : ", ");
    ss << " };\n";

    std::string energyInstance = instanceName + "_energyVoxSet";
    ss << energy.dump2cpp(energyInstance);
    ss << instanceName << ".energy = std::move(" << energyInstance << ");\n";

    std::string dedxInstance = instanceName + "_dedxVoxSet";
    ss << dedx.dump2cpp(dedxInstance);
    ss << instanceName << ".dedx = std::move(" << dedxInstance << ");\n";

    std::string firstptInstance = instanceName + "_firstEdep";
    ss << first_pt.dump2cpp(firstptInstance);
    ss << instanceName << ".first_pt = std::move(" << firstptInstance << ");\n";

    std::string lastptInstance = instanceName + "_lastEdep";
    ss << last_pt.dump2cpp(lastptInstance);
    ss << instanceName << ".last_pt = std::move(" << lastptInstance << ");\n";

    return ss.str();
  }

  bool ParticleLabel::operator==(const ParticleLabel &rhs) const
  {
    // ParticleLabels equivalent if all their data match...
    return (part == rhs.part) &&
           (valid == rhs.valid) &&
           (merged_v == rhs.merged_v) &&
           (energy == rhs.energy) &&
           (dedx == rhs.dedx) &&
           (first_pt == rhs.first_pt) &&
           (last_pt == rhs.last_pt);

  }

  std::string ParticleLabel::dump() const
  {
    std::stringstream st;

    st << "ParticleLabel object:\n";
    st << "  valid = " << valid << "\n";
    st << "  num voxels: " << energy.size() << "\n";
    st << "  first step: " << StringifyEDep(first_pt) << "\n";
    st << "  last  step: " << StringifyEDep(last_pt) << "\n";
    st << "  merged into track ID: " << StringifyTrackID(merge_id) << "\n";
    st << "  merged track IDs: ";
    for (const TrackID_t tk : merged_v)
      st << tk << " ";
    st << "\n";
    st << "  Particle:\n";
    st << part.dump() << "\n";

    return st.str();
  }

  // --------------------------------------------------------
  /*
  const supera::VoxelSet &EventOutput::VoxelDeDxs() const
  {
    // recompute only if particle list has changed under us
    if (IsDirty(DIRTY_FLAG::kDeDx))
    {
      // we want an energy-weighted mean here.
      // first add all the dEdXs up, weighted by their energies...
      supera::VoxelSet dEdXs;
      for (const supera::ParticleLabel & part : Particles())
      {
        for (const supera::Voxel & dedx : part.dedx.as_vector())
          dEdXs.emplace(dedx.id(), dedx.value() * VoxelEnergies().find(dedx.id()).value(), true);
      }

      // now renormalize them...
      for (const supera::Voxel & dedx : dEdXs.as_vector())
      {
        auto energy = VoxelEnergies().find(dedx.id()).value();
        if (energy > 0)
          dEdXs.emplace(dedx.id(), dedx.value() / energy, false);
      }

      _dEdXs = std::move(dEdXs);
    } // if (IsDirty(...))

    return _dEdXs;
  }
  
  // --------------------------------------------------------

  const supera::VoxelSet &EventOutput::VoxelEnergies() const
  {
    // recompute only if particle list has changed under us
    if (IsDirty(DIRTY_FLAG::kEnergy))
    {
      supera::VoxelSet energies;
      for (const supera::ParticleLabel & part : Particles())
        energies.emplace(part.energy, true);
      _energies = std::move(energies);
    }

    return _energies;
  }

  // --------------------------------------------------------
  
  const supera::VoxelSet &
  EventOutput::VoxelLabels(const std::vector<supera::SemanticType_t> &semanticPriority) const
  {
    // recompute only if particle list has changed under us
    if (IsDirty(DIRTY_FLAG::kLabel))
    {
      supera::VoxelSet semantics;
      for (const supera::ParticleLabel & part : Particles())
      {
        auto const &vs = part.energy;
        SemanticType_t semantic = part.part.shape;
        for (auto const &vox : vs.as_vector())
        {
          auto const &prev = semantics.find(vox.id());
          if (prev.id() == supera::kINVALID_VOXELID)
            semantics.emplace(vox.id(), semantic, false);
          else
          {
            // todo: what if the new voxel has 10x the energy??
            SemanticType_t prioritized_semantic = EventOutput::_SemanticPriority(static_cast<SemanticType_t>(prev.value()), semantic, semanticPriority);
            if (prioritized_semantic != static_cast<SemanticType_t>(prev.value()))
              semantics.emplace(vox.id(), semantic, false);
          }
        } // for (vox)
      } // for (part)

      _semanticLabels = std::move(semantics);
    } // if (IsDirty(...))

    return _semanticLabels;
  }
*/
  // --------------------------------------------------------



  const ParticleLabel& EventOutput::Particle(InstanceID_t id) const
  {
    // check if it's already ordered. Then we're done
    if (_particles.size() > id && _particles[id].part.id == id)
      return _particles[id];

    for(auto const& p : _particles) {
      if(p.part.id != id)
        continue;
      return p;
    }
    std::cerr<<"No particle found with an id "<<id<<std::endl;
    throw meatloaf();
  }

  supera::SemanticType_t EventOutput::_SemanticPriority(supera::SemanticType_t a, supera::SemanticType_t b,
                                                        const std::vector<supera::SemanticType_t> & semanticPriority)
  {
    if (a == b)
      return a;
    for (auto const &semantic : semanticPriority)
    {
      if (a == semantic)
        return a;
      if (b == semantic)
        return b;
    }
    return a;
  }

  // --------------------------------------------------------

  bool EventOutput::operator==(const EventOutput &rhs) const
  {
    // the event outputs are the same if their ParticleLabels are the same.

    // definitely different if different lengths
    if (_particles.size() != rhs._particles.size())
      return false;

    // the ParticleLabels might not be in the same order,
    // but if the *contents* are the same, they're still equivalent.
    // we'll sort them by GEANT trackID using a map from track id -> vector index
    // (since maps are inherently sorted).
    // then we can walk through the track indices, check those first,
    // and only if they are equal do the actual comparison between ParticleList objects.
    std::map<supera::TrackID_t, std::size_t> lhs_tracks, rhs_tracks;
    for (std::size_t idx = 0; idx < _particles.size(); idx++)
    {
      lhs_tracks[_particles[idx].part.trackid] = idx;
      rhs_tracks[rhs._particles[idx].part.trackid] = idx;
    }

    auto it_lhs = lhs_tracks.begin();
    auto it_rhs = rhs_tracks.begin();
    for (; it_lhs != lhs_tracks.end() && it_rhs != rhs_tracks.end(); it_lhs++, it_rhs++ )
    {
      // since the map is sorted, any mismatch in track IDs already means they're not equal
      if (it_lhs->first != it_rhs->first)
        return false;

      // if we have to, check the objects themselves
      if (_particles[it_lhs->second] != rhs._particles[it_rhs->second])
        return false;
    }

    return true;

  }

  // --------------------------------------------------------

  std::string EventOutput::dump2cpp(const std::string &instanceName) const
  {
    std::stringstream ss;

    ss << "supera::EventOutput " << instanceName << ";\n";

    ss << instanceName << ".Particles().reserve(" << _particles.size() << ");\n";
    for (std::size_t idx = 0; idx < _particles.size(); idx++)
    {
      const supera::ParticleLabel & part = _particles[idx];
      std::string partInstance = instanceName + "_part" + std::to_string(idx) + "_label";
      ss << part.dump2cpp(partInstance);
      ss << instanceName << ".Particles().push_back(std::move(" << partInstance << "));\n";
    }

    return ss.str();
  }

  /// Helper function that generates a vector of vector of voxel id and value
  void EventOutput::FillClustersEnergy(std::vector<std::vector<supera::VoxelID_t> >& ids,
    std::vector<std::vector<float> >& values,
    bool fill_unassociated) const
  {
    ids.resize(_particles.size());
    values.resize(_particles.size());

    for(size_t i=0; i<_particles.size(); ++i) {
      auto const& target = _particles[i].energy;
      auto& id_v = ids[i];
      auto& value_v = values[i];
      id_v.clear();
      value_v.clear();
      id_v.reserve(target.size());
      value_v.reserve(target.size());
      for(auto const& vox : target.as_vector()) {
        id_v.push_back(vox.id());
        value_v.push_back(vox.value());
      }
    }

    if(fill_unassociated) {
      auto const& unass = _unassociated_voxels.as_vector();
      std::vector<supera::VoxelID_t> id_v;
      std::vector<float> value_v;
      id_v.resize(unass.size());
      value_v.resize(unass.size());
      for(size_t i=0; i<unass.size(); ++i){
        id_v[i] = unass[i].id();
        value_v[i] = unass[i].value();
      }
      ids.push_back(id_v);
      values.push_back(value_v);
    }
  }

  void EventOutput::FillClustersdEdX(std::vector<std::vector<supera::VoxelID_t> >& ids,
    std::vector<std::vector<float> >& values,
    bool fill_unassociated) const
  {
    ids.resize(_particles.size());
    values.resize(_particles.size());

    for(size_t i=0; i<_particles.size(); ++i) {
      auto const& target = _particles[i].dedx;
      auto& id_v = ids[i];
      auto& value_v = values[i];
      id_v.clear();
      value_v.clear();
      id_v.reserve(target.size());
      value_v.reserve(target.size());
      for(auto const& vox : target.as_vector()) {
        id_v.push_back(vox.id());
        value_v.push_back(vox.value());
      }
    }
    if(fill_unassociated) {
      auto const& unass = _unassociated_voxels.as_vector();
      std::vector<supera::VoxelID_t> id_v;
      std::vector<float> value_v;
      id_v.resize(unass.size());
      value_v.resize(unass.size(),0.);
      for(size_t i=0; i<unass.size(); ++i)
        id_v[i] = unass[i].id();
      ids.push_back(id_v);
      values.push_back(value_v);
    }
  }

  void EventOutput::FillTensorSemantic(std::vector<VoxelID_t>& ids,
    std::vector<float>& values) const
  {
    ids.clear();
    values.clear();
    ids.reserve(_semanticLabels.size());
    values.reserve(_semanticLabels.size());
    for(auto const& vox : _semanticLabels.as_vector())
    {
      ids.push_back(vox.id());
      values.push_back(vox.value());
    }
  }

  void EventOutput::FillTensorEnergy(std::vector<VoxelID_t>& ids,
    std::vector<float>& values) const
  {
    ids.clear();
    values.clear();
    ids.reserve(_energies.size());
    values.reserve(_energies.size());
    for(auto const& vox : _energies.as_vector())
    {
      ids.push_back(vox.id());
      values.push_back(vox.value());
    }
  }

} // namespace supera
