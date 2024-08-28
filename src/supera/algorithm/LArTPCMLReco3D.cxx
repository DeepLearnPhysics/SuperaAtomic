#ifndef __LARTPCMLRECO3D_CXX__
#define __LARTPCMLRECO3D_CXX__

#include "LArTPCMLReco3D.h"
#include <algorithm>
#include <cassert>
#include <set>
#include <string>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <numeric>

namespace supera {

    LArTPCMLReco3D::LArTPCMLReco3D(std::string name)
    : LabelAlgorithm(name)
    {}

    // --------------------------------------------------------------------
    void LArTPCMLReco3D::_configure(const YAML::Node& cfg)
    {
        LOG_INFO() << "starting" << std::endl;

        _semantic_priority.clear();
        if(cfg["SemanticPriority"])
            _semantic_priority = cfg["SemanticPriority"].as<std::vector<size_t> >();
        this->SetSemanticPriority(_semantic_priority);

        _touch_threshold = 1;
        if(cfg["TouchDistance"]) 
            _touch_threshold = cfg["TouchDistance"].as<size_t>();

        _edep_threshold = 0.01;
        if(cfg["EnergyDepositThreshold"]) 
            _edep_threshold = cfg["EnergyDepositThreshold"].as<double>();

        _delta_size = 3;
        if(cfg["DeltaSize"])
            _delta_size = cfg["DeltaSize"].as<size_t>();

        _compton_size = 10;
        if(cfg["ComptonSize"])
            _compton_size = cfg["ComptonSize"].as<size_t>();

        _lescatter_size = 2;
        if(cfg["LEScatterSize"])
            _lescatter_size = cfg["LEScatterSize"].as<size_t>();

        _store_lescatter = true;
        if(cfg["StoreLEScatter"])
            _store_lescatter = cfg["StoreLEScatter"].as<bool>();

        _rewrite_interactionid = true;
        if(cfg["RewriteInteractionID"])
            _rewrite_interactionid = cfg["RewriteInteractionID"].as<bool>();

        _merge_unassociated_edeps = true;
        if(cfg["MergeUnassociatedEDeps"])
            _merge_unassociated_edeps = cfg["MergeUnassociatedEDeps"].as<bool>();

        std::vector<double> min_coords(3,std::numeric_limits<double>::lowest());
        std::vector<double> max_coords(3,std::numeric_limits<double>::max());
        if(cfg["WorldBoundMin"])
            min_coords = cfg["WorldBoundMin"].as<std::vector<double> >();
        if(cfg["WorldBoundMax"])
            max_coords = cfg["WorldBoundMax"].as<std::vector<double> >();

        _world_bounds.update(min_coords.at(0),min_coords.at(1),min_coords.at(2),
            max_coords.at(0),max_coords.at(1),max_coords.at(2));

        LOG_INFO() << "done" << std::endl;

    }
    // --------------------------------------------------------------------

    void LArTPCMLReco3D::SetSemanticPriority(std::vector<size_t>& order)
    {
        LOG_INFO() << "starting" << std::endl;

        std::vector<size_t> result;
        std::vector<bool> assigned((size_t)(supera::kShapeUnknown),false);
        for(auto const& type : order) {
            if(type >= supera::kShapeUnknown) {
                LOG_FATAL() << "SemanticPriority received an unsupported semantic type " << type << "\n";
                throw meatloaf(std::to_string(__LINE__));
            }
            bool ignore = false;
            for(auto const& used : result) {
                if(used != type) continue;
                ignore = true;
            }
            if(ignore) {
                LOG_FATAL() << "Duplicate SemanticPriority received for type " << type << "\n";
                throw meatloaf(std::to_string(__LINE__));
            }
            result.push_back(type);
            assigned[type]=true;
        }

        // Now add other types to make sure
        for(size_t i=0; i<assigned.size(); ++i)
        {
            if(assigned[i]) continue;
            result.push_back(i);
        }
        if(result.size() != (size_t)(kShapeUnknown)) {
            LOG_FATAL() << "Logic error!\n";
            throw meatloaf(std::to_string(__LINE__));
        }
        order = result;
        LOG_INFO() << "done" << std::endl;
    }

    EventOutput LArTPCMLReco3D::Generate(const EventInput& data, const ImageMeta3D& meta)
    {
        LOG_INFO() << "starting" << std::endl;

        EventOutput result;

        // fill in the working structures that link the list of particles and its genealogy
        _mcpl.InferParentage(data);

        // Assign the initial labels for each particle.
        // They will be grouped together in various ways in the subsequent steps.
        supera::VoxelSet unassociated_voxels;
        std::vector<supera::ParticleLabel> labels = this->InitializeLabels(data, meta, unassociated_voxels);

        //// Now group the labels together in certain cases
        //// (e.g.: electromagnetic showers, neutron clusters, ...)
        //// There are lots of edge cases so the logic is spread out over many methods.
        ////this->MergeShowerIonizations(labels); // merge supera::kIonization = too small delta rays into parents
        //// ** TODO identify and merge too-small shower fragments to other touching showers **
        ////this->MergeShowerTouchingElectron(meta, labels); // merge larcv::kShapeLEScatter to touching shower
        //// Apply energy threshold (may drop some pixels)


        this->ApplyEnergyThreshold(labels);
        this->SetSemanticType(labels);
        this->MergeShowerConversion(labels); // merge supera::kConversion a photon merged to a parent photon
        this->MergeTouchingLEScatter(meta,labels);
        this->MergeShowerFamilyTouching(meta, labels); // merge supera::kShapeShower to touching parent shower/delta/michel
        this->MergeShowerTouching(meta, labels); // merge supera::kShapeShower to touching shower in the same family tree

        //this->MergeTouchingLEScatter(meta,labels);

        // Re-classify small photons into ShapeLEScatter
        for(auto& label : labels) {
            if(!label.valid) continue;
            if(label.part.type != supera::kPhoton) continue;
            if(label.energy.size() < _compton_size)
                label.part.shape = supera::kShapeLEScatter;
        }
        this->MergeTouchingLEScatter(meta,labels);
    	// Re-classify large ShapeLEScatter into shower
        for(auto& label : labels) {
         if(!label.valid) continue;
         if(label.part.shape != supera::kShapeLEScatter) continue;
         if(label.energy.size() < _compton_size)
           continue;
       label.part.shape = supera::kShapeShower;
   }

        // ** TODO consider this separate from MergeShowerIonizations?? **
   this->MergeDeltas(labels); // merge supera::kDelta to a parent if too small

        // After all merging is 

        // Now that we have grouped the true particles together,
        // at this point we're ready to build a new set of labels
        // which contain only the top particle of each merged group.

   auto output2input = this->RegisterOutputParticles(labels);

   this->SetGroupID(labels,output2input);

   this->SetAncestorAttributes(labels,output2input);

        // maybe the user has already set these upstream
        // (like in DUNE ND-LAr case)
   if (_rewrite_interactionid)
      this->SetInteractionID(labels);

        // We're finally to fill in the output container.
        // There are two things we need:
        //  (1) labels for each voxel (what semantic type is each one?)
        //  (2) labels for particle groups.
        // The output format is an object containing a collection of supera::ParticleLabel,
        // each of which has voxels attached to it, so that covers both things.
        // EventOutput computes VoxelSets with the sum across all particles
        // for voxel energies and semantic labels
  this->BuildOutputLabels(labels,result,output2input,unassociated_voxels);

  LOG_INFO() << "done" << std::endl;

  return result;
}

    // --------------------------------------------------------------------

void LArTPCMLReco3D::BuildOutputLabels(std::vector<supera::ParticleLabel>& labels,
    supera::EventOutput& result, 
    const std::vector<InstanceID_t>& output2input,
    const supera::VoxelSet& unass) const
{
    LOG_INFO() << "starting" << std::endl;
        // Build the outupt
    std::vector<supera::ParticleLabel> output_particles;
    output_particles.reserve(output2input.size());
    for(auto const& index : output2input) {
        output_particles.emplace_back(std::move(labels[index]));
        labels[index].valid=false;
    }

        // Fill the energy and semantic label tensor for unassociated 3D points
    result._energies.reserve(unass.size());
    result._semanticLabels.reserve(unass.size());
    for(auto const& vox : unass.as_vector()) {
        result._energies.emplace(vox.id(),vox.value(),true);
        result._semanticLabels.emplace(vox.id(),supera::kShapeLEScatter,false);
    }

        // Semantic label
        // dedx energy semanticlabels
    for (auto rit = _semantic_priority.crbegin(); rit != _semantic_priority.crend(); ++rit)
    {
        auto stype = supera::SemanticType_t((*rit));
        for(auto& label : output_particles) {
            if(label.part.shape != stype)
                continue;
                // Contribute to the output
            assert(label.energy.size() == label.dedx.size());
                //result._dEdXs.reserve(label.energy.size()+result._dEdXs.size());
            result._energies.reserve(label.energy.size()+result._energies.size());
            result._semanticLabels.reserve(label.energy.size()+result._semanticLabels.size());
            auto const& input_dedx   = label.dedx.as_vector();
            auto const& input_energy = label.energy.as_vector();
            for(size_t i=0; i<input_dedx.size(); ++i) {
                    //result._dEdXs.emplace(input_dedx[i].id(),input_dedx[i].value(),true);
                result._energies.emplace(input_energy[i].id(),input_energy[i].value(),true);
                result._semanticLabels.emplace(input_energy[i].id(),(float)(stype),false);
                if(std::isnan(input_energy[i].value())) {
                    LOG_ERROR() << "NAN found (HE)" << std::endl;
                    LOG_ERROR() << label.dump() << std::endl;
                    throw meatloaf(std::to_string(__LINE__));
                }
            }
        }
            // If this is LEScatter type, and if _store_lescatter == false, make sure to add here
        if(stype == supera::kShapeLEScatter && !_store_lescatter) 
        {
            for(auto& label : labels){
                
                if(!label.valid) continue;

                if(label.part.shape != supera::kShapeUnknown) {
                    LOG_FATAL() << "Unexpected (logic error): valid particle remaining that is not kShapeUnknown shape...\n"
                    << label.dump() << "\n";
                    throw meatloaf(std::to_string(__LINE__));
                }
                    // Contribute to the output
                assert(label.energy.size() == label.dedx.size());
                    //result._dEdXs.reserve(label.energy.size()+result._dEdXs.size());
                result._energies.reserve(label.energy.size()+result._energies.size());
                result._semanticLabels.reserve(label.energy.size()+result._semanticLabels.size());
                auto const& input_dedx   = label.dedx.as_vector();
                auto const& input_energy = label.energy.as_vector();
                for(size_t i=0; i<input_dedx.size(); ++i) {
                        //result._dEdXs.emplace(input_dedx[i].id(),input_dedx[i].value(),true);
                    result._energies.emplace(input_energy[i].id(),input_energy[i].value(),true);
                    result._semanticLabels.emplace(input_energy[i].id(),(float)(stype),false);
                    if(std::isnan(input_energy[i].value())) {
                        LOG_ERROR() << "NAN found (LE)" << std::endl;
                        LOG_ERROR() << label.dump() << std::endl;
                        throw meatloaf(std::to_string(__LINE__));
                    }
                }
            }
        }
    }

    result = std::move(output_particles);
    result._unassociated_voxels = unass;
    LOG_INFO() << "done" << std::endl;

}

    // --------------------------------------------------------------------

    // ------------------------------------------------------
void LArTPCMLReco3D::MergeParticleLabel(std::vector<supera::ParticleLabel>& labels,
    supera::InstanceID_t dest_id,
    supera::InstanceID_t target_id) const 
{
    auto& dest   = labels.at(dest_id);
    auto& target = labels.at(target_id);
    dest.Merge(target);
    for(auto const& id : target.merged_v)
        labels.at(id).merge_id = dest.id;
}
    // ------------------------------------------------------


    // ------------------------------------------------------
void LArTPCMLReco3D::SetGroupID(std::vector<supera::ParticleLabel>& labels,
    const std::vector<InstanceID_t>& output2input) const
{
    LOG_INFO() << "starting" << std::endl;
    for(auto& label : labels){
        if(!label.valid) continue;

        auto& part = label.part;

            // Primary particles get its own group ID
        if( label.id == label.parent_id) {
            part.group_id = label.id;
        }else{
            auto parent_index = output2input.at(label.parent_id);

            switch(part.shape)
            {
                case kShapeTrack:
                case kShapeMichel:
                part.group_id = part.id;
                break;

                case kShapeDelta:
                if(parent_index == kINVALID_INDEX || !labels[parent_index].valid) {
                    LOG_FATAL() << "Delta ray with an invalid parent is not allowed!\n";
                    throw meatloaf(std::to_string(__LINE__));
                }
                        //part.group_id = labels[parent_index].part.id;
                part.group_id = part.id;
                break;

                case kShapeShower:

                part.group_id = part.id;
                for(auto const& parent_index : _mcpl.ParentIdArray(output2input.at(label.id)))
                {
                    auto parent_shape = labels[parent_index].part.shape;
                    auto parent_id    = labels[parent_index].id;

                    if(parent_index == kINVALID_INDEX)
                        continue;
                    if(!labels[parent_index].valid)
                        continue;
                    if(parent_shape == kShapeLEScatter)
                        continue;
                    if(parent_shape != kShapeShower && parent_shape != kShapeMichel && parent_shape != kShapeDelta)
                        break;
                    if(parent_id == kINVALID_INSTANCEID)
                        continue;
                    part.group_id = labels[parent_index].id;
                }

                break;

                case kShapeLEScatter:
                break;

                default:
                LOG_FATAL() << " Unexpected shape type " << part.shape << "\n";
                throw meatloaf(std::to_string(__LINE__));
                break;
            }
        }
    }
    LOG_INFO() << "done" << std::endl;
}

    // ------------------------------------------------------

void LArTPCMLReco3D::SetAncestorAttributes(std::vector<supera::ParticleLabel>& labels,
    const std::vector<InstanceID_t> &output2input) const
{
    LOG_INFO() << "starting" << std::endl;
    for(auto& label : labels){

        if(!label.valid) continue;
        auto label_id = output2input.at(label.id);
        auto parent_id = output2input.at(label.parent_id);
        auto ancestor_id = output2input.at(label.ancestor_id);
        auto const& parent_input_id_v = _mcpl.ParentIdArray(label_id);

            // Now parent_trackid must be filled unless the input data was insufficient
        if(parent_id == kINVALID_INSTANCEID){
            LOG_FATAL() << "Parent instance ID missing for the instance ID " 
            << label.id << "\n"
            << "Check the input data and make sure all particles have a parent track ID\n";
            throw meatloaf(std::to_string(__LINE__));
        }

            // If ancestor_trackid is invalid, set it to the parent.
        if(ancestor_id == kINVALID_INSTANCEID) {
            LOG_INFO() << "Ancestor instance ID not set for the instance ID "
            << label.id << "\n"
            << "Setting it to the parent instance ID " << parent_id << "\n";
            ancestor_id = parent_id;
        }

            // Consistency check
        if(!parent_input_id_v.empty() && parent_input_id_v.front() != parent_id) {
            LOG_FATAL() << "Logic error: the parent instance ID " << parent_id
            << " != the first in the ancestory instance IDs " << parent_input_id_v.front() << "\n";
            throw meatloaf(std::to_string(__LINE__));
        }

        if(!parent_input_id_v.empty() && parent_input_id_v.back() != ancestor_id) {
            LOG_FATAL() << "Logic error: the ancestor instance ID " << ancestor_id
            << " != the most distant parent ID " << parent_input_id_v.back() << "\n";
            throw meatloaf(std::to_string(__LINE__));
        }


            // Attempt to fill parent info. 
        if(parent_id != kINVALID_INSTANCEID) {
            auto const& parent = labels[parent_id];
            label.part.parent_trackid = parent.part.trackid;
            label.part.parent_id  = parent.part.id;
            label.part.parent_pdg = parent.part.pdg;
            label.part.parent_vtx = parent.part.vtx;
            label.part.parent_process = parent.part.process; 
        }

            // Attempt to fill ancestor info
        if(ancestor_id != kINVALID_INSTANCEID) {
            auto const& ancestor = labels[ancestor_id];
            label.part.ancestor_trackid = ancestor.part.trackid;
            label.part.ancestor_id  = ancestor.part.id;
            label.part.ancestor_pdg = ancestor.part.pdg;
            label.part.ancestor_vtx = ancestor.part.vtx;
            label.part.ancestor_process = ancestor.part.process; 
        }
    }
    LOG_INFO() << "done" << std::endl;
}


    // ------------------------------------------------------


void LArTPCMLReco3D::SetInteractionID(std::vector<supera::ParticleLabel>& labels) const
{
    LOG_INFO() << "starting" << std::endl;
    std::vector<supera::Vertex> int2vtx;
    for(auto& label : labels) {
        if(!label.valid) continue;

        InstanceID_t iid = kINVALID_INSTANCEID;
        for(InstanceID_t cand_iid=0; cand_iid<int2vtx.size(); ++cand_iid)
        {
            auto const& vtx = int2vtx[cand_iid];
            if(vtx == label.part.ancestor_vtx){
                iid = cand_iid;
                break;
            }
        }
        if(iid == kINVALID_INSTANCEID){
            iid = int2vtx.size();
            int2vtx.push_back(label.part.ancestor_vtx);
        }
        label.part.interaction_id = iid;
    }
    LOG_INFO() << "done" << std::endl;
}
    // ------------------------------------------------------


    // ------------------------------------------------------
void LArTPCMLReco3D::ApplyEnergyThreshold(std::vector<supera::ParticleLabel>& labels) const
{
    LOG_INFO() << "starting" << std::endl;
        // Loop again and eliminate voxels that has energy below threshold
    for (auto &label : labels)
    {
        supera::VoxelSet energies, dEdXs;
        energies.reserve (label.energy.size() );
        dEdXs.reserve    (label.dedx.size()   );

        const auto energy_vec =  label.energy.as_vector();
        const auto dedx_vec   =  label.dedx.as_vector();
        for (std::size_t idx = 0; idx < energy_vec.size(); idx++)
        {
            const auto & vox = energy_vec[idx];
            if (vox.value() < _edep_threshold)
                continue;
            if (dedx_vec[idx].id() != vox.id()) {
                LOG_FATAL() << "Unmatched voxel ID between dE/dX and energy voxels \n";
                throw meatloaf(std::to_string(__LINE__));
            }
            energies.emplace (vox.id(), vox.value(),           true);
            dEdXs.emplace    (vox.id(), dedx_vec[idx].value(), true);
        }
        label.energy = std::move(energies);
        label.dedx = std::move(dEdXs);
    }
    LOG_INFO() << "done" << std::endl;
} // LArTPCMLReco3D::ApplyEnergyThreshold()

    // ------------------------------------------------------
void LArTPCMLReco3D::SetSemanticType(std::vector<supera::ParticleLabel>& labels) const
{
    LOG_INFO() << "starting" << std::endl;
    for(auto& label : labels) {
        if(!label.valid) continue;

        switch(label.part.type) {
            case kInvalidProcess:
            label.part.shape = supera::kShapeUnknown;
            throw meatloaf("'kInvalidProcess' particle process encountered:\n" + label.dump());

            case kTrack:
            label.part.shape = supera::kShapeTrack;
            break;

            case kPrimary:
            if(std::abs(label.part.pdg) != 11 && label.part.pdg != 22) 
            {
                label.part.shape = supera::kShapeTrack;
                break;
            }
            label.part.shape = supera::kShapeShower;
            break;

            case kDelta:
            if(label.energy.size() < _delta_size)
                label.part.shape = kShapeLEScatter;
            else
                label.part.shape = kShapeDelta;
            break;

            case kDecay:
            if(std::abs(label.part.pdg) == 11 && std::abs(label.part.parent_pdg) == 13)
                label.part.shape = kShapeMichel;
            else if(std::abs(label.part.pdg)==11 || label.part.pdg == 22) {
                if(label.energy.size() > _compton_size)
                    label.part.shape = kShapeShower;
                else {
                    label.part.shape = kShapeLEScatter;
                    LOG_WARNING() << "Assigned to kShapeLEScatter " << std::endl << label.dump() << std::endl;
                }
            }else{
                label.part.shape = kShapeTrack;
            }
            break;

            case kIonization:
            case kPhotoElectron:
            case kNeutron:
            label.part.shape = kShapeLEScatter;
                    //LOG_WARNING() << "\"Neutron\" type assigned to kShapeLEScatter " << std::endl << label.dump() << std::endl;
            break;

            case kPhoton:
            label.part.shape = kShapeShower;
            break;

            case kConversion:
            case kCompton:
            case kOtherShower:
            if(std::abs(label.part.pdg)==11 || label.part.pdg == 22) {
                if(label.energy.size() > _compton_size)
                    label.part.shape = kShapeShower;
                else {
                    label.part.shape = kShapeLEScatter;
                            //LOG_WARNING() << "Assigned to kShapeLEScatter " << std::endl << label.dump() << std::endl;
                }
            }else{
                label.part.shape = kShapeTrack;
            }
            break;

            case kNucleus:
            if(label.energy.size() > _compton_size)
                label.part.shape = kShapeTrack;
            else
                label.part.shape = kShapeLEScatter;
            break;

        }
    }
    LOG_INFO() << "done" << std::endl;
}


    // ------------------------------------------------------

std::vector<InstanceID_t>
LArTPCMLReco3D::RegisterOutputParticles(std::vector<supera::ParticleLabel> &inputLabels) const 
{
    LOG_INFO() << "starting" << std::endl;
        /*
            This function registers the output particles to be stored.
            It entails 2 actions
                1. Provide output ID = index number in the output vector of particles
                    1.1. First treat particles that are not kShapeLEScatter
                    1.2. Second, provide output ID to kShapeLEScatter particles
                2. Set merged particle track ID to also point to the output ID of its superset particle
                3. Record track ID => output index mapping
            In addition, for those particles to be stored:
                4. Set the first and the last step 
        */
    std::vector<InstanceID_t> output2input, input2output, lescatter_id_v;
    input2output.resize(inputLabels.size(),kINVALID_INSTANCEID);
    output2input.reserve(inputLabels.size());
    lescatter_id_v.reserve(inputLabels.size());

        // first create the track id list
    output2input.clear();
    output2input.reserve(input2output.size());

        // assign particle group ID numbers and make sure they have all info set
    LOG_VERBOSE() << "Considering incoming particles:\n";
    for (size_t label_index=0; label_index<inputLabels.size(); ++label_index)
    {
        auto& inputLabel = inputLabels[label_index];
        auto input_id = inputLabel.id;
        
        if (!inputLabel.valid)
        {
            LOG_VERBOSE() << "   --> invalid particle (i.e. already merged), skipping \n";
            continue;
        }
        
        if (inputLabel.part.shape == supera::kShapeLEScatter) {
            LOG_VERBOSE() << "   --> LEScatter, skipping in the first loop\n";
            lescatter_id_v.push_back(label_index);
            continue;
        }
        if (inputLabel.part.shape == supera::kShapeUnknown) {
            LOG_FATAL()   << "   --> ShapeUnknown found and unexpected!\n"
            << inputLabel.dump() << "\n";
            throw meatloaf(std::to_string(__LINE__));
        }

            // Register this particle
        if(input2output[input_id] == kINVALID_INSTANCEID) {
            input2output[input_id]=output2input.size();
            output2input.push_back(input_id);
        }

            // Register all parents if not yet
        for(auto const& parent_id : _mcpl.ParentIdArray(input_id)){
            if(input2output[parent_id] != kINVALID_INSTANCEID)
                continue;
            input2output[parent_id]=output2input.size();
            output2input.push_back(parent_id);
            inputLabels[parent_id].valid = true;
        }
    }

    if(_store_lescatter) {
        for(auto const& input_id : lescatter_id_v) {
            if(input2output[input_id] == kINVALID_INSTANCEID) {
                input2output[input_id] = output2input.size();
                output2input.push_back(input_id);
            }
            for(auto const& parent_id : _mcpl.ParentIdArray(input_id)) {
                if(input2output[parent_id] != kINVALID_INSTANCEID)
                    continue;
                input2output[parent_id]=output2input.size();
                output2input.push_back(parent_id);
                inputLabels[parent_id].valid = true;
            }
        }
    }

        // Set particle's id to invalid if not to be stored in the output.
    for(InstanceID_t input_id=0; input_id<input2output.size(); ++input_id) {

        auto& label = inputLabels[input_id];

        if(input2output[input_id] != kINVALID_INSTANCEID) {
            if(!(label.valid)) {
                LOG_FATAL() << "Unexpected: invalid particle instance ID " << input_id
                << " registered to be stored" << std::endl;
                throw meatloaf();
            }
            continue;
        }
        if(label.valid) {
            LOG_FATAL() << "Unexpected: valid particle instance ID " << input_id
            << " not registered to be stored" << std::endl;
            throw meatloaf();
        }

        label.id = label.parent_id = label.ancestor_id = kINVALID_INSTANCEID;
    }

        // Now update the labels for those particles to be stored
    for(InstanceID_t output_id=0; output_id<output2input.size(); ++output_id) {

        auto const& input_id = output2input[output_id];
        auto& label = inputLabels[input_id];
        auto& part  = label.part;

        label.id = output_id;
        label.parent_id = input2output[label.parent_id];
        label.ancestor_id = input2output[label.ancestor_id];

        part.id = label.id;
        part.parent_id = label.parent_id;
        part.ancestor_id = label.ancestor_id;
        part.energy_deposit = label.energy.size() ? label.energy.sum() : 0.;

        auto const &first_pt = label.first_pt;
        auto const &last_pt = label.last_pt;
        if (first_pt.t != kINVALID_DOUBLE) {
            part.first_step = supera::Vertex(first_pt.x, first_pt.y, first_pt.z, first_pt.t);
        }else if (label.energy.size()){
            LOG_WARNING() << "Non-empty particle yet the first step is not filled!" << std::endl << label.dump() << std::endl;
        }

        if (last_pt.t != kINVALID_DOUBLE){
            part.last_step = supera::Vertex(last_pt.x, last_pt.y, last_pt.z, last_pt.t);
        }else if (label.energy.size() && part.shape == supera::kShapeTrack){
            LOG_WARNING() << "Non-empty particle yet the last step is not filled!" << std::endl << label.dump() << std::endl;
        }

    }
    LOG_INFO() << "done" << std::endl;
    return output2input;

} // LArTPCMLReco3D::RegisterOutputParticles()

    // ------------------------------------------------------

std::vector<size_t> sort_by_voxel_count(const std::vector<supera::ParticleLabel> &parts) {

    std::vector<size_t> idx_v(parts.size());
    std::iota(idx_v.begin(), idx_v.end(), 0);

    std::stable_sort(idx_v.begin(), idx_v.end(), 
        [&parts](size_t i1, size_t i2) {return parts[i1].energy.size() < parts[i2].energy.size();});

    return idx_v;
}


std::vector<supera::ParticleLabel>
LArTPCMLReco3D::InitializeLabels(const EventInput &evtInput, 
    const supera::ImageMeta3D &meta, 
    VoxelSet &unassociated_voxels) const
{
    LOG_INFO() << "starting" << std::endl;
        // this default-constructs the whole lot of them, which fills their values with defaults/invalid values
    std::vector<supera::ParticleLabel> labels(evtInput.size());

    LOG_DEBUG() << "Initializing labels with incoming particles...\n";
    for (std::size_t idx = 0; idx < evtInput.size(); idx++)
    {
        auto& label = labels[idx];
        label.id = evtInput[idx].id;
        label.parent_id = evtInput[idx].parent_id;
        label.ancestor_id = evtInput[idx].ancestor_id;
        label.part  = evtInput[idx].part;
        label.part.parent_pdg = _mcpl.ParentPdgCode()[idx];

        if(label.part.parent_pdg != supera::kINVALID_PDG)
            label.valid = true;

        for (const supera::EDep & edep : evtInput[idx].pcloud)
        {
            auto vox_id = meta.id(edep);
            if(vox_id == supera::kINVALID_VOXELID || !_world_bounds.contains(edep)) {
                LOG_VERBOSE() << "Skipping EDep from instance ID " << label.id
                << " E=" << edep.e
                << " pos=" << edep.x << "," << edep.y << "," << edep.z << ")\n";
                continue;
            }

            label.energy.emplace (vox_id, edep.e,    true);
            label.dedx.emplace   (vox_id, edep.dedx, true);
            label.UpdateFirstPoint(edep);
            label.UpdateLastPoint(edep);
        }

        LOG_VERBOSE() << label.dump() << "\n";

    }  // for (idx)

        //
        // Handle unassociated_edeps in EventInput
        //

    std::vector<supera::VoxelSet> unassociated_dedx_v;
    std::vector<supera::VoxelSet> unassociated_vox_v;

    unassociated_dedx_v.resize(evtInput.unassociated_edeps.size());
    unassociated_vox_v.resize(evtInput.unassociated_edeps.size());

    size_t invalid_unass_ctr=0;
    size_t valid_unass_ctr=0;
    for(size_t cidx=0; cidx<evtInput.unassociated_edeps.size(); ++cidx) {

        auto const& input_cluster = evtInput.unassociated_edeps[cidx];
        auto& dedx_cluster = unassociated_dedx_v[cidx];
        auto& vox_cluster  = unassociated_vox_v[cidx];

        dedx_cluster.reserve(input_cluster.size());
        vox_cluster.reserve(input_cluster.size());

        for(auto const& edep : input_cluster) {

            auto vox_id = meta.id(edep);
            if(vox_id == supera::kINVALID_VOXELID) {
                invalid_unass_ctr++;
                continue;
            }
            vox_cluster.emplace(vox_id, edep.e, true);
            dedx_cluster.emplace(vox_id, edep.dedx, true);
            valid_unass_ctr++;
        }
    }
    if(invalid_unass_ctr){
        size_t total_unass_count=0;
        for(auto const& vs : evtInput.unassociated_edeps)
            total_unass_count += vs.size();
        LOG_WARNING() << invalid_unass_ctr << "/" << total_unass_count
        << " unassociated packets are ignored (outside BBox)" << std::endl;
    }


    if(!_merge_unassociated_edeps) {
        unassociated_voxels.reserve(valid_unass_ctr);
        for(auto const& cluster : unassociated_vox_v) {
            for(auto const& vox : cluster.as_vector()) {
                unassociated_voxels.add(vox);
            }
        }

    }else{
        VoxelSet unass_after;
        unass_after.reserve(valid_unass_ctr);

        for(size_t cidx=0; cidx<unassociated_vox_v.size(); cidx++) {

            auto const& vox_vs  = unassociated_vox_v[cidx];
            auto const& dedx_vs = unassociated_dedx_v[cidx];
            bool merged=false;
            for(auto const& lidx : sort_by_voxel_count(labels)) {
                auto& label = labels[lidx];
                if(this->IsTouching(meta, vox_vs,label.energy)) {
                        //merge
                    auto const& vox_v = vox_vs.as_vector();
                    auto const& dedx_v = dedx_vs.as_vector();
                    for(size_t ivox=0; ivox<vox_v.size(); ++ivox) {
                        label.energy.add(vox_v[ivox]);
                        label.dedx.add(dedx_v[ivox]);
                    }
                    merged=true;
                    break;
                }
            }
            if(!merged) {
                for(auto const& vox : vox_vs.as_vector()) {
                    unass_after.add(vox);
                }
            }
        }
        LOG_WARNING() << "Unassociated voxels " << valid_unass_ctr << " => " << unass_after.size() << std::endl;
        unassociated_voxels=unass_after;
    }
    LOG_INFO() << "done" << std::endl;
    return labels;
}  // LArTPCMLReco3D::InitializeLabels()


    // ------------------------------------------------------

void LArTPCMLReco3D::MergeShowerConversion(std::vector<supera::ParticleLabel>& labels) const
{
    LOG_INFO() << "starting" << std::endl;
    int merge_ctr = 0;
        //int invalid_ctr = 0;
    do
    {
        merge_ctr = 0;
        for (auto &label : labels)
        {
            if (!label.valid) continue;
                //if(grp.part.type != supera::kIonization && grp.part.type != supera::kConversion) continue;
            if (label.part.type != supera::kConversion) continue;
            if (std::abs(label.part.pdg) != 11) {
                LOG_FATAL() << "Unexpected: type kConversion for a particle that is not electron!\n";
                throw meatloaf(std::to_string(__LINE__));
            }

            auto const& parent_id_v = _mcpl.ParentIdArray(label.id);
            InstanceID_t found_id = kINVALID_INSTANCEID;
            for(auto const& parent_id : parent_id_v) 
            {
                LOG_DEBUG() << "Inspecting: instance id " << label.id << " => parent instance id " << parent_id << "\n";
                if (parent_id == supera::kINVALID_INSTANCEID || !labels[parent_id].valid)
                {
                    LOG_VERBOSE() << "Missing/Invalid parent particle with a instance id " << parent_id << "\n"
                    << "Could not find a parent for instance id " << label.id 
                    << " Track ID " << label.part.trackid
                    << " PDG " << label.part.pdg
                    << " " << label.part.process << " E = " << label.part.energy_init
                    << " (" << label.part.energy_deposit << ") MeV\n";
                    continue;
                }
                found_id = parent_id;
                break;
            }
            if (found_id != kINVALID_INSTANCEID) {
                this->MergeParticleLabel(labels,found_id,label.id);
                merge_ctr++;
            }
        }

        LOG_INFO() << "Merge counter: " << merge_ctr << "\n";
    } while (merge_ctr > 0);
    LOG_INFO() << "done" << std::endl;
}  // LArTPCMLReco3D::MergeShowerConversion()

    // ------------------------------------------------------

void LArTPCMLReco3D::MergeDeltas(std::vector<supera::ParticleLabel>& labels) const
{
    LOG_INFO() << "starting" << std::endl;
    for (auto &label : labels)
    {
            //if(label.part.type != supera::kDelta) continue;
        if (label.part.shape != supera::kShapeDelta) continue;
        auto parent_id = label.parent_id;
        if(parent_id == kINVALID_INSTANCEID) continue;
        auto &parent = labels[parent_id];
        if (!parent.valid) continue;

            // allows the test on unique voxels to be put in the if() below
            // and only used if needed due to short-circuiting.
        const auto UniqueVoxelCount = [](const supera::ParticleLabel & grp, const supera::ParticleLabel & parent)
        {
            size_t unique_voxel_count = 0;
            for (auto const &vox : grp.energy.as_vector())
            {
                if (parent.energy.find(vox.id()).id() == supera::kINVALID_VOXELID)
                    ++unique_voxel_count;
            }
            return unique_voxel_count;
        };

            // if voxel count is smaller than delta ray requirement, simply merge
        if (label.energy.size() < _delta_size || UniqueVoxelCount(label, parent) < _delta_size)
        {
                // if parent is found, merge
            LOG_INFO() << "Merging delta id " << label.part.id << " PDG " << label.part.pdg
            << " " << label.part.process << " vox count " << label.energy.size() 
            << " (unique " << UniqueVoxelCount(label, parent) << ")\n"
            << " ... parent found with instance id " << parent.id
            << " PDG " << parent.part.pdg << " " << parent.part.process << "\n";
            LOG_INFO() << "Time difference: " << label.part.first_step.time - parent.part.first_step.time << "\n";
            this->MergeParticleLabel(labels, parent.id, label.id);
        }
        else
        {
            LOG_INFO() << "NOT merging delta id " << label.id << " PDG " << label.part.pdg
            << " " << label.part.process << " vox count " << label.energy.size() << "\n"
            <<" ... parent found with instance id " << parent.id
            << " PDG " << parent.part.pdg << " " << parent.part.process << "\n";

        }

    }
    LOG_INFO() << "done" << std::endl;
} // LArTPCMLReco3D::MergeShowerDeltas()
    // ------------------------------------------------------

void LArTPCMLReco3D::MergeShowerFamilyTouching(const supera::ImageMeta3D& meta,
 std::vector<supera::ParticleLabel>& labels) const
{
    LOG_INFO() << "starting" << std::endl;
        // Merge touching shower fragments
        // Direct parentage between kShapeShower => kShapeShower/kShapeDelta/kShapeMichel
    int merge_ctr = 0;
    int invalid_ctr = 0;
    do {
        merge_ctr = 0;
        for (auto& label : labels) {
            if (!label.valid) continue;
            if (label.part.shape != supera::kShapeShower) continue;
            if (label.parent_id == supera::kINVALID_INSTANCEID) continue;  // primaries can't have parents
                // search for a possible parent
            auto parent_id = kINVALID_INSTANCEID;
            LOG_VERBOSE() << "   Found particle group with shape 'shower', PDG=" << label.part.pdg
            << "\n    instance id=" << label.id
            << ", and alleged parent id=" << label.parent_id << "\n";
                // a direct parent ?
            if (parent_id != kINVALID_INSTANCEID && labels[parent_id].valid)
                parent_id = label.parent_id;
            else
            {
                for (size_t shower_index = 0; shower_index < labels.size(); ++shower_index)
                {
                    auto const &candidate_grp = labels[shower_index];
                    if (candidate_grp.id == label.parent_id || !candidate_grp.valid)
                        continue;
                    for (auto const &id : candidate_grp.merged_v)
                    {
                        if (id != label.parent_id)
                            continue;
                        parent_id = static_cast<InstanceID_t>(candidate_grp.id);
                        break;
                    }
                    if (parent_id != kINVALID_INSTANCEID)
                        break;
                }
            }
            if (parent_id == kINVALID_INSTANCEID || parent_id == label.id) continue;
            auto& parent = labels[parent_id];
                //auto parent_type = labels[parent_trackid].part.type;
                //if(parent_type == supera::kTrack || parent_type == supera::kNeutron) continue;
            if (parent.part.shape != supera::kShapeShower && 
                parent.part.shape != supera::kShapeDelta && 
                parent.part.shape != supera::kShapeMichel)
                continue;
            if (!parent.valid) continue;
            if (this->IsTouching(meta, label.energy, parent.energy)) {
                    // if parent is found, merge
                this->MergeParticleLabel(labels, parent.id, label.id);
                LOG_VERBOSE() << "   Merged to group w/ id=" << parent.id << "\n";
                merge_ctr++;
            }
        }
        LOG_DEBUG() << "Merge counter: " << merge_ctr << " invalid counter: " << invalid_ctr << "\n";
    } while (merge_ctr>0);
    LOG_INFO() << "done" << std::endl;
} // LArTPCMLReco3D::MergeShowerFamilyTouching()


    // ------------------------------------------------------

void LArTPCMLReco3D::MergeShowerIonizations(std::vector<supera::ParticleLabel>& labels) const
{
    LOG_INFO() << "starting" << std::endl;
        // Loop over particles of a type kIonization (=touching to its parent physically by definition)
        // If a parent is found, merge to the parent
    int merge_ctr = 0;
    int invalid_ctr = 0;
    do
    {
        merge_ctr = 0;
        for (auto &label : labels)
        {
            if (!label.valid) continue;
            if (label.part.type != supera::kIonization) continue;
                // merge to a valid "parent"

            bool parent_found = false;
            auto parent_id = kINVALID_INSTANCEID;
            for(auto const& id : _mcpl.ParentIdArray(label.id))
            {
                parent_id = id;
                if(parent_id == kINVALID_INSTANCEID) continue;
                if(!labels[parent_id].valid) continue;
                parent_found = true;
                break;
            }
                // if parent is found, merge
            if (parent_found)
            {
                this->MergeParticleLabel(labels,parent_id,label.id); 
                merge_ctr++;
            }
        } // for (grp)
        LOG_DEBUG() << "Ionization merge counter: " << merge_ctr << " invalid counter: " << invalid_ctr << "\n";
    } while (merge_ctr > 0);
    LOG_INFO() << "done" << std::endl;
} // LArTPCMLReco3D::MergeShowerIonizations()


    // ------------------------------------------------------
void LArTPCMLReco3D::MergeShowerTouching(const supera::ImageMeta3D& meta,
   std::vector<supera::ParticleLabel>& labels) const
{
    LOG_INFO() << "starting" << std::endl;
        // Go over all pair-wise combination of two shower instances
        // For each shower, find all consecutive parents of shower/michel/delta type (break if track found)
        // If there is a common parent in two list AND if two showers are physically touching, merge
    int merge_ctr = 0;
    do
    {
        merge_ctr = 0;
        for (size_t i = 0; i < labels.size(); ++i)
        {
            auto &label_a = labels[i];
            if (!label_a.valid) continue;
            if (label_a.part.shape != supera::kShapeShower) continue;

            for (size_t j = i+1; j < labels.size(); ++j)
            {
                if (i == j) continue;
                auto &label_b = labels[j];
                if (!label_b.valid) continue;
                if (label_b.part.shape != supera::kShapeShower) continue;

                    // check if these showers share the parentage
                    // list a's parents
                std::set<supera::InstanceID_t> parent_list_a;
                std::set<supera::InstanceID_t> parent_list_b;

                auto parents_a = this->ParentShowerIDs(label_a.id, labels);
                for (auto const &parent_id : parents_a) parent_list_a.insert(parent_id);
                    parent_list_a.insert(label_a.id);

                auto parents_b = this->ParentShowerIDs(label_b.id, labels);
                for (auto const &parent_id : parents_b) parent_list_b.insert(parent_id);
                    parent_list_b.insert(label_b.id);

                bool same_family = false;
                for (auto const &parent_id : parent_list_a)
                {
                    if (parent_list_b.find(parent_id) != parent_list_b.end())
                        same_family = true;
                    if (same_family) break;
                }
                for (auto const &parent_id : parent_list_b)
                {
                    if (parent_list_a.find(parent_id) != parent_list_a.end())
                        same_family = true;
                    if (same_family) break;
                }


                if (same_family && this->IsTouching(meta, label_a.energy, label_b.energy))
                {

                    if (label_a.energy.size() > label_b.energy.size())
                        this->MergeParticleLabel(labels, label_a.id, label_b.id);
                    else
                        this->MergeParticleLabel(labels, label_b.id, label_a.id);
                    merge_ctr++;
                }
            }
        }
        LOG_INFO() << "Merge counter: " << merge_ctr << "\n";
    } while (merge_ctr > 0);
    LOG_INFO() << "done" << std::endl;
} // LArTPCMLReco3D::MergeShowerTouching()

    // ------------------------------------------------------

void LArTPCMLReco3D::MergeShowerTouchingElectron(const supera::ImageMeta3D& meta,
    std::vector<supera::ParticleLabel>& labels) const
{
    LOG_INFO() << "starting" << std::endl;
    size_t merge_ctr = 1;
    while (merge_ctr)
    {
        merge_ctr = 0;
        for (auto &label : labels)
        {
                //if (!label.valid || label.energy.size() < 1 || label.shape() != supera::kShapeLEScatter) continue;
            if( !label.valid || label.energy.size()<1 || 
                label.energy.size()>_compton_size ||
                std::abs(label.part.pdg) != 11)
                continue;
            if( label.part.type != kPhotoElectron && 
                label.part.type != kIonization && 
                label.part.type != kCompton &&
                label.part.type != kConversion)
                continue;

            auto const &parents = _mcpl.ParentIdArray(label.id);

            LOG_VERBOSE() << "Inspecting LEScatter Instance ID " << label.id
            << " PDG " << label.part.pdg
            << " " << label.part.process << "\n";
            LOG_VERBOSE() << "  ... parents:\n";
            for(auto const& parent_id : parents)
                LOG_VERBOSE() << "     "<< parent_id << "\n";

            for (auto const &parent_id : parents)
            {
                if(parent_id == kINVALID_INDEX) continue;
                auto &parent = labels[parent_id];
                if (!parent.valid || parent.energy.size() < 1) continue;
                if (this->IsTouching(meta, label.energy, parent.energy))
                {
                    LOG_VERBOSE() << "Merging LEScatter instance id = " << label.id
                    << " into touching parent shower group (id=" << parent.part.group_id << ")"
                    << " with instance id = " << parent.id << "\n";
                    this->MergeParticleLabel(labels,parent_id,label.id);
                    merge_ctr++;
                    break;
                }
            } // for (parent_trackid)
        } // for (grp)
    } // while (merge_ctr)
    LOG_INFO() << "done" << std::endl;
} // LArTPCMLReco3D::MergeShowerTouchingElectron()

    // ------------------------------------------------------

void LArTPCMLReco3D::MergeTouchingLEScatter(const supera::ImageMeta3D& meta,
    std::vector<supera::ParticleLabel>& labels) const
{
    LOG_INFO() << "starting" << std::endl;
    std::vector<supera::SemanticType_t> priorities = {supera::kShapeShower, supera::kShapeMichel, supera::kShapeDelta,
        supera::kShapeTrack, supera::kShapeGhost};

        size_t merge_ctr = 1;
        while (merge_ctr)
        {
            merge_ctr = 0;
            for (auto &label : labels)
            {
                //if (!label.valid || label.energy.size() < 1 || label.shape() != supera::kShapeLEScatter) continue;
                if( !label.valid || label.energy.size()<1 || 
                    //label.energy.size()>_lescatter_size ||
                    label.part.shape != supera::kShapeLEScatter)
                    continue;

                if( label.part.type == supera::kNeutron || label.part.type == supera::kNucleus)
                    continue;

                bool merged = false;
                // Try shower first
                for(auto const& shape : priorities) {
                    for(auto &dest : labels) {
                        if(!dest.valid || dest.part.shape != shape)
                            continue;

                        if(this->IsTouching(meta, label.energy, dest.energy))
                        {
                            this->MergeParticleLabel(labels,dest.id,label.id);
                            merge_ctr++;
                            merged=true;
                            break;
                        }
                    }
                    if(merged) break;
                }
            } // for (grp)
        } // while (merge_ctr)
        LOG_INFO() << "done" << std::endl;
    } // LArTPCMLReco3D::MergeTouchingLEScatter()

    // ------------------------------------------------------

    bool LArTPCMLReco3D::IsTouching(const ImageMeta3D& meta, const VoxelSet& vs1, const VoxelSet& vs2, bool verbose) const
    {
        LOG_DEBUG() << "starting" << std::endl;
        bool touching = false;
        size_t ix1, iy1, iz1;
        size_t ix2, iy2, iz2;
        size_t diffx, diffy, diffz;

        // Test1: is there an overlapping voxel?
        if(vs1.size()<vs2.size()) {
            for(auto const& vox : vs1.as_vector()) 
            {
                auto const& overlap = vs2.find(vox.id());
                if(overlap.id() != kINVALID_VOXELID)
                    return true;
            }
        }else{
            for(auto const& vox : vs2.as_vector()) 
            {
                auto const& overlap = vs1.find(vox.id());
                if(overlap.id() != kINVALID_VOXELID)
                    return true;
            }
        }


        // Test2: is there a voxel close enough? (distance calculation)
        for (auto const &vox1 : vs1.as_vector())
        {
            meta.id_to_xyz_index(vox1.id(), ix1, iy1, iz1);
            for (auto const &vox2 : vs2.as_vector())
            {
                meta.id_to_xyz_index(vox2.id(), ix2, iy2, iz2);
                if (ix1 > ix2) diffx = ix1 - ix2; else diffx = ix2 - ix1;
                if (iy1 > iy2) diffy = iy1 - iy2; else diffy = iy2 - iy1;
                if (iz1 > iz2) diffz = iz1 - iz2; else diffz = iz2 - iz1;
                touching = diffx <= _touch_threshold && diffy <= _touch_threshold && diffz <= _touch_threshold;
                if (touching && verbose)
                {
                    LOG_INFO()<<"Touching ("<<ix1<<","<<iy1<<","<<iz1<<") ("<<ix2<<","<<iy2<<","<<iz2<<")\n";
                    LOG_INFO()<<"    Dist ("<<diffx<<","<<diffy<<","<<diffz<<")\n";
                    break;
                }
            }
            if (touching) break;
        }
        LOG_DEBUG() << "done" << std::endl;
        return touching;
    } // LArTPCMLReco3D::IsTouching()


    // ------------------------------------------------------

    std::vector<supera::InstanceID_t>
    LArTPCMLReco3D::ParentShowerIDs(InstanceID_t id,
        const std::vector<supera::ParticleLabel>& labels,
        bool include_lescatter) const
    {
        LOG_DEBUG() << "starting" << std::endl;
        std::vector<supera::InstanceID_t> result;
        if( id >= labels.size() ) {
            LOG_FATAL() << "Cannot infer the input id " << id << " which is larger than the max possible " << labels.size() << std::endl;
            throw meatloaf();
        }
        if( id == kINVALID_INSTANCEID )
            return result;
        auto const& parents = _mcpl.ParentIdArray(id);
        result.reserve(parents.size());

        for(auto const& parent_id : parents) {

            if(parent_id == kINVALID_INDEX) continue;

            auto const& grp = labels[parent_id];

            if(grp.part.shape == supera::kShapeTrack ||
             grp.part.shape == supera::kShapeUnknown)
                break;

            if(!grp.valid) continue;

            if(grp.part.shape == supera::kShapeMichel ||
             grp.part.shape == supera::kShapeShower ||
             grp.part.shape == supera::kShapeDelta  ||
             (grp.part.shape == supera::kShapeLEScatter && include_lescatter))
                result.push_back(parent_id);
        }
        LOG_DEBUG() << "done" << std::endl;
        return result;
    } // LArTPCMLReco3D::ParentShowerTrackIDs()

    // ------------------------------------------------------




}
#endif
