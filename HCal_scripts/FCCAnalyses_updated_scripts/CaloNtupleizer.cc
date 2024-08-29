#include "FCCAnalyses/CaloNtupleizer.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "edm4hep/MCParticleData.h"

#include <math.h>

#include "DD4hep/Detector.h"

namespace FCCAnalyses{

namespace CaloNtupleizer{

dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
dd4hep::DDSegmentation::BitFieldCoder* m_decoderECal;
dd4hep::DDSegmentation::BitFieldCoder* m_decoderHCal;
dd4hep::DDSegmentation::BitFieldCoder* m_decoderBoth = new dd4hep::DDSegmentation::BitFieldCoder("system:4");

void loadGeometry(std::string xmlGeometryPath, std::string readoutName){
  dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
  dd4hepgeo->fromCompact(xmlGeometryPath);
  dd4hepgeo->volumeManager();
  dd4hepgeo->apply("DD4hepVolumeManager", 0, 0);
  m_decoder = dd4hepgeo->readout(readoutName).idSpec().decoder();
}

void loadGeometryBoth(std::string xmlGeometryPath, std::string readoutName, std::string readoutECal, std::string readoutHCal){
  dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
  dd4hepgeo->fromCompact(xmlGeometryPath);
  dd4hepgeo->volumeManager();
  dd4hepgeo->apply("DD4hepVolumeManager", 0, 0);
  m_decoder = dd4hepgeo->readout(readoutName).idSpec().decoder();
  m_decoderECal = dd4hepgeo->readout(readoutECal).idSpec().decoder();
  m_decoderHCal = dd4hepgeo->readout(readoutHCal).idSpec().decoder();
}


sel_layers::sel_layers(int arg_min, int arg_max) : _min(arg_min), _max(arg_max) {};
ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>  sel_layers::operator() (const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& in) {

  ROOT::VecOps::RVec<edm4hep::CalorimeterHitData> res;
  for (auto & p: in){
    dd4hep::DDSegmentation::CellID cellId = p.cellID;
    int layer = m_decoder->get(cellId, "layer");
    if (layer>_min && layer<_max)res.emplace_back(p);
  }
  return res;
}
// SIM calo hit
ROOT::VecOps::RVec<float> getSimCaloHit_r (const ROOT::VecOps::RVec<edm4hep::SimCalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(std::sqrt(p.position.x * p.position.x + p.position.y * p.position.y));
  }
  return result;
}

ROOT::VecOps::RVec<float> getSimCaloHit_x (const ROOT::VecOps::RVec<edm4hep::SimCalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.position.x);
  }
  return result;
}

ROOT::VecOps::RVec<float> getSimCaloHit_y (const ROOT::VecOps::RVec<edm4hep::SimCalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.position.y);
  }
  return result;
}

ROOT::VecOps::RVec<float> getSimCaloHit_z (const ROOT::VecOps::RVec<edm4hep::SimCalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.position.z);
  }
  return result;
}

ROOT::VecOps::RVec<float> getSimCaloHit_phi (const ROOT::VecOps::RVec<edm4hep::SimCalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    TVector3 t3;
    t3.SetXYZ(p.position.x, p.position.y, p.position.z);
    result.push_back(t3.Phi());
  }
  return result;
}

ROOT::VecOps::RVec<float> getSimCaloHit_theta (const ROOT::VecOps::RVec<edm4hep::SimCalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    TVector3 t3;
    t3.SetXYZ(p.position.x, p.position.y, p.position.z);
    result.push_back(t3.Theta());
  }
  return result;
}

ROOT::VecOps::RVec<float> getSimCaloHit_eta (const ROOT::VecOps::RVec<edm4hep::SimCalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    TVector3 t3;
    t3.SetXYZ(p.position.x, p.position.y, p.position.z);
    result.push_back(t3.Eta());
  }
  return result;
}

ROOT::VecOps::RVec<float> getSimCellID (const ROOT::VecOps::RVec<edm4hep::SimCalorimeterHitData>& in){
  ROOT::VecOps::RVec<int> result;
  for (auto & p: in){
    result.push_back(p.cellID);
  }
  return result;
}

ROOT::VecOps::RVec<float> getSimCaloHit_energy (const ROOT::VecOps::RVec<edm4hep::SimCalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.energy);
  }
  return result;
}

ROOT::VecOps::RVec<int> getSimCaloHit_depth (const ROOT::VecOps::RVec<edm4hep::SimCalorimeterHitData>& in,const int decodingVal){
  ROOT::VecOps::RVec<int> result;
  for (auto & p: in){
    result.push_back(p.cellID >> decodingVal & (8-1) );
  }
  return result;
}

ROOT::VecOps::RVec<TVector3> getSimCaloHit_positionVector3 (const ROOT::VecOps::RVec<edm4hep::SimCalorimeterHitData>& in){
  ROOT::VecOps::RVec<TVector3> result;
  for (auto & p: in){
    TVector3 t3;
    t3.SetXYZ(p.position.x, p.position.y, p.position.z);
    result.push_back(t3);
  }
  return result;
}

// calo hit
ROOT::VecOps::RVec<float> getCaloHit_x (const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.position.x);
  }
  return result;
}

ROOT::VecOps::RVec<float> getCaloHit_y (const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.position.y);
  }
  return result;
}

ROOT::VecOps::RVec<float> getCaloHit_z (const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.position.z);
  }
  return result;
}

ROOT::VecOps::RVec<float> getCaloHit_phi (const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    TVector3 t3;
    t3.SetXYZ(p.position.x, p.position.y, p.position.z);
    result.push_back(t3.Phi());
  }
  return result;
}

ROOT::VecOps::RVec<int>
getCaloHit_phiIdx(const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData> &in) {
  ROOT::VecOps::RVec<int> result;
  for (auto & p: in){
    dd4hep::DDSegmentation::CellID cellId = p.cellID;
    result.push_back(m_decoder->get(cellId, "phi"));
  }
  return result;
}

ROOT::VecOps::RVec<int> getCaloHit_moduleIdx(
    const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData> &in) {
  ROOT::VecOps::RVec<int> result;
  for (auto &p : in) {
    dd4hep::DDSegmentation::CellID cellId = p.cellID;
    result.push_back(m_decoder->get(cellId, "module"));
  }
  return result;
}

ROOT::VecOps::RVec<int>
getCaloHit_thetaIdx(const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData> &in) {
  ROOT::VecOps::RVec<int> result;
  for (auto &p : in) {
    dd4hep::DDSegmentation::CellID cellId = p.cellID;
    result.push_back(m_decoder->get(cellId, "theta"));
  }
  return result;
}

ROOT::VecOps::RVec<float> getCaloHit_theta (const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    TVector3 t3;
    t3.SetXYZ(p.position.x, p.position.y, p.position.z);
    result.push_back(t3.Theta());
  }
  return result;
}

ROOT::VecOps::RVec<float> getCaloHit_eta (const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    TVector3 t3;
    t3.SetXYZ(p.position.x, p.position.y, p.position.z);
    result.push_back(t3.Eta());
  }
  return result;
}

ROOT::VecOps::RVec<int>
getCaloHit_etaIdx(const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData> &in) {
  ROOT::VecOps::RVec<int> result;
  for (auto & p: in){
    dd4hep::DDSegmentation::CellID cellId = p.cellID;
    result.push_back(m_decoder->get(cellId, "eta"));
  }
  return result;
}

ROOT::VecOps::RVec<float> getCaloHit_energy (const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.energy);
  }
  return result;
}

ROOT::VecOps::RVec<int> getCaloHit_layer (const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& in){
  ROOT::VecOps::RVec<int> result;
  for (auto & p: in){
    dd4hep::DDSegmentation::CellID cellId = p.cellID;
    result.push_back(m_decoder->get(cellId, "layer"));
  }
  return result;
}

ROOT::VecOps::RVec<TVector3> getCaloHit_positionVector3 (const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& in){
  ROOT::VecOps::RVec<TVector3> result;
  for (auto & p: in){
    TVector3 t3;
    t3.SetXYZ(p.position.x, p.position.y, p.position.z);
    result.push_back(t3);
  }
  return result;
}

// calo cluster
ROOT::VecOps::RVec<float> getCaloCluster_x (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.position.x);
  }
  return result;
}

ROOT::VecOps::RVec<float> getCaloCluster_y (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.position.y);
  }
  return result;
}

ROOT::VecOps::RVec<float> getCaloCluster_z (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.position.z);
  }
  return result;
}

ROOT::VecOps::RVec<float> getCaloCluster_phi (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    TVector3 t3;
    t3.SetXYZ(p.position.x, p.position.y, p.position.z);
    result.push_back(t3.Phi());
  }
  return result;
}

ROOT::VecOps::RVec<float> getCaloCluster_theta (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    TVector3 t3;
    t3.SetXYZ(p.position.x, p.position.y, p.position.z);
    result.push_back(t3.Theta());
  }
  return result;
}

ROOT::VecOps::RVec<float> getCaloCluster_eta (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    TVector3 t3;
    t3.SetXYZ(p.position.x, p.position.y, p.position.z);
    result.push_back(t3.Eta());
  }
  return result;
}

ROOT::VecOps::RVec<float> getCaloCluster_energy (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.energy);
  }
  return result;
}

ROOT::VecOps::RVec<TVector3> getCaloCluster_positionVector3 (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in){
  ROOT::VecOps::RVec<TVector3> result;
  for (auto & p: in){
    TVector3 t3;
    t3.SetXYZ(p.position.x, p.position.y, p.position.z);
    result.push_back(t3);
  }
  return result;
}

ROOT::VecOps::RVec<int> getCaloCluster_firstCell (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in){
  ROOT::VecOps::RVec<int> result;
  for (auto & p: in){
    result.push_back(p.hits_begin);
  }
  return result;
}

ROOT::VecOps::RVec<int> getCaloCluster_lastCell (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in){
  ROOT::VecOps::RVec<int> result;
  for (auto & p: in){
    result.push_back(p.hits_end);
  }
  return result;
}

ROOT::VecOps::RVec<std::vector<float>>
getCaloCluster_energyInLayers (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in,
                               const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& cells,
                               const int nLayers) {
  static const int layer_idx = m_decoder->index("layer");
  static const int cryo_idx = m_decoder->index("cryo");
  ROOT::VecOps::RVec<std::vector<float>> result;
  result.reserve(in.size());

  for (const auto & c: in) {
    std::vector<float> energies(nLayers, 0);
    for (auto i = c.hits_begin; i < c.hits_end; i++) {
      int layer = m_decoder->get(cells[i].cellID, layer_idx);
      int cryoID = m_decoder->get(cells[i].cellID, cryo_idx);
      if(cryoID == 0) {
        energies[layer] += cells[i].energy;
      }
    }
    result.push_back(energies);
  }
  return result;
}

ROOT::VecOps::RVec<std::vector<float>>
getCaloCluster_energyInLayersHCal (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in,
                               const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& cells,
                               const int nLayers) {
  static const int layer_idx = m_decoder->index("layer");
  //static const int cryo_idx = m_decoder->index("cryo");
  ROOT::VecOps::RVec<std::vector<float>> result;
  result.reserve(in.size());

  for (const auto & c: in) {
    std::vector<float> energies(nLayers, 0);
    for (auto i = c.hits_begin; i < c.hits_end; i++) {
      int layer = m_decoder->get(cells[i].cellID, layer_idx);
      //int cryoID = m_decoder->get(cells[i].cellID, cryo_idx);
      //if(cryoID == 0) {
        energies[layer] += cells[i].energy;
        if(layer>13){ 
          std::cout << "standalone layer " << layer << " \n ";     
          std::cout << "standalone ene " << cells[i].energy << " \n "; 
        }
    }
    result.push_back(energies);
  }
  return result;
}

ROOT::VecOps::RVec<std::vector<float>>
getCaloCluster_energyInLayersBoth (const ROOT::VecOps::RVec<edm4hep::ClusterData>& in,
                               const ROOT::VecOps::RVec<edm4hep::CalorimeterHitData>& cells,
                               const int nLayers, const int nLayersECal) {

  ROOT::VecOps::RVec<std::vector<float>> result;
  result.reserve(in.size());

  // loop over clusters
  for (const auto & c: in) {
    std::vector<float> energies(nLayers, 0);

    // loop over cluster hits = cells 
    for (auto i = c.hits_begin; i < c.hits_end; i++) {

      // identify calo system
      auto systemId = m_decoderBoth->get(cells[i].cellID, "system");
      //std::cout << systemId << " \n ";

      if (systemId == 4) { // ECAL BARREL system id
        static const int layer_idx = m_decoderECal->index("layer");
        if (layer_idx != 4){ 
        std::cout << "ECal layer_idx " << layer_idx << " \n "; 
        } 
	static const int cryo_idx = m_decoderECal->index("cryo");
        int layerECal = m_decoderECal->get(cells[i].cellID, layer_idx);
        //std::cout << "ECal layer" << layerECal << " \n "; 
        int cryoID = m_decoderECal->get(cells[i].cellID, cryo_idx);
        if(cryoID == 0) {
          energies[layerECal] += cells[i].energy;        
        }
      }

      if (systemId == 8) { // HCAL BARREL system id
        static const int layer_idx = m_decoderHCal->index("layer"); 
        if (layer_idx != 1){  
        std::cout << "HCal layer_idx " << layer_idx << " \n "; 
      }
       	int layerHCal = m_decoderHCal->get(cells[i].cellID, layer_idx);
        //std::cout << "HCal layer " << layerHCal << " \n ";
        energies[nLayersECal+layerHCal] += cells[i].energy;
        
       	//std::cout << "HCal cells ene " << cells[i].energy << " \n "; 
      }
    }

    result.push_back(energies);
  }
  //std::cout << "Mici si supis!";
  return result;
}


ROOT::VecOps::RVec<float> getSimParticleSecondaries_x (const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in){
    result.push_back(p.vertex.x);
  }
  return result;
}


ROOT::VecOps::RVec<float> getSimParticleSecondaries_y (const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in){
  ROOT::VecOps::RVec<float> result;
for (auto & p: in) {
  result.push_back(p.vertex.y);
}
return result;
}


ROOT::VecOps::RVec<float> getSimParticleSecondaries_z (const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in){
  ROOT::VecOps::RVec<float> result;
for (auto & p: in) {
  result.push_back(p.vertex.z);
}
return result;
}


  ROOT::VecOps::RVec<float> getSimParticleSecondaries_PDG (const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in){
    ROOT::VecOps::RVec<float> result;
    for (auto & p: in) {
      result.push_back(p.PDG);
    }
    return result;
  }

ROOT::VecOps::RVec<float> getSimParticleSecondaries_phi (const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Phi());
  }
  return result;
}


ROOT::VecOps::RVec<float> getSimParticleSecondaries_theta (const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Theta());
  }
  return result;
}

ROOT::VecOps::RVec<float> getSimParticleSecondaries_eta (const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Eta());
  }
  return result;
}

ROOT::VecOps::RVec<float> getSimParticleSecondaries_energy (const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.E());
  }
  return result;
}

getFloatAt::getFloatAt(size_t pos) {m_pos = pos;}
ROOT::RVecF getFloatAt::operator()(const ROOT::VecOps::RVec<std::vector<float>>& in) {
  ROOT::RVecF result;
  for (auto & v : in) {
    result.push_back(v[m_pos]);
  }
  return result;
}

}//end NS CaloNtupleizer

}//end NS FCCAnalyses
