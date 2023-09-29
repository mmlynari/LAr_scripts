/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// Run with
// root
// .L displayV3.C+
// display()
//
/******************************************************************************/


/******************************************************************************/
// dependencies
/******************************************************************************/
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLine.h>
#include <TArc.h>
#include <TCrown.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TBox.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TGLUtil.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TMaterial.h>
#include <TView.h>
#include <TGeoTube.h>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>
#include <TTUBE.h>
#include <TEvePointSet.h>
#include <TEveBoxSet.h>
#include <TEveStraightLineSet.h>
#include <TEveLine.h>
#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include <TEveQuadSet.h>
#include <TEveGeoShape.h>
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEveGeoNode.h>
#include <TEveProjectionAxes.h>
#include <TEveProjectionManager.h>
#include <TEveBrowser.h>
#include <TEveRGBAPaletteOverlay.h>
#include <TGLViewer.h>
#include <TGTextEntry.h>
#include <TGTab.h>
#include <TGButton.h>
#include <TGLAnnotation.h>
#include <TPRegexp.h>
#include <TSystem.h>
#include <iostream>
#include <cmath>
using namespace std;


/******************************************************************************/
// SETTINGS FOR THE DISPLAY
// - geometry file (and corresponding merging of theta cells/modules vs layer
// - data file
// - flags to draw or not elements of the event (particles, hits, ..)
// - minimum energy of particles, hits, cells, clusters
/******************************************************************************/

// G4 geometry file
std::string geomFile = "ECalBarrel.root";
  
// merging along theta and module directions and drawing options
//std::string evtFile = "output_fullCalo_SimAndDigi_withTopoCluster_MagneticField_False_pMin_10000_MeV_ThetaMinMax_40_140_pdgId_11_pythiaFalse_NoiseFalse.root";
//std::string evtFile = "test.root";
//const std::vector<int> mergedCells_Theta = {4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
//const std::vector<int> mergedModules = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
//bool drawMergedCells = false;

//std::string evtFile = "output_fullCalo_SimAndDigi_withTopoCluster_MagneticField_False_pMin_10000_MeV_ThetaMinMax_40_140_pdgId_11_pythiaFalse_NoiseFalse_testmerge.root";
//const std::vector<int> mergedCells_Theta = {4, 8, 2, 1, 8, 4, 2, 1, 4, 2, 1, 8};
//const std::vector<int> mergedModules = {2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1};
//bool drawMergedCells = true;

// 50 electrons 100 GeV
// std::string evtFile = "output_evts_50_100_GeV_ThetaMinMax_40_140_PhiMinMax_0_6.28318.root";
// 1 electron 50 GeV
// std::string evtFile = "output_evts_1_pdg_11_50_GeV_ThetaMinMax_90_90_PhiMinMax_1.570795_1.570795_MagneticField_False_NoiseFalse.root";
// 1 pi0  50 GeV
std::string evtFile = "output_evts_1_pdg_111_50_GeV_ThetaMinMax_90_90_PhiMinMax_1.570795_1.570795_MagneticField_False_NoiseFalse.root";
// 1 photon 50 GeV
// std::string evtFile = "output_evts_1_pdg_22_50_GeV_ThetaMinMax_90_90_PhiMinMax_1.570795_1.570795_MagneticField_False_NoiseFalse.root";
// 1 pi0  10 GeV
// std::string evtFile = "output_evts_1_pdg_111_10_GeV_ThetaMinMax_90_90_PhiMinMax_1.570795_1.570795_MagneticField_False_NoiseFalse.root";
// 1 photon 10 GeV
// std::string evtFile = "output_evts_1_pdg_22_10_GeV_ThetaMinMax_90_90_PhiMinMax_1.570795_1.570795_MagneticField_False_NoiseFalse.root";

const std::vector<int> mergedCells_Theta = {4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
const std::vector<int> mergedModules = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

bool drawParticles = true;
bool drawHits = true;
bool drawCells = true;
bool drawMergedCells = false;
bool drawTopoClusters = true;

// min particle energy (GeV)
float ParticleEnergyThreshold = 1.0;
// min hit and cell energy (GeV)
float HitEnergyThreshold = 0.0;
float CellEnergyThreshold = 0.0;
// min cluster energy (in GeV)
float TopoClusterEnergyThreshold = 2.;


/******************************************************************************/
// GEOMETRICAL CONSTANTS
// They should match the geometry used in the geometry file and the readout
// used to produce the event file
/******************************************************************************/

// units
// G4 uses mm but Root uses cm..
const double cm = 1.0;
const double mm = 0.1;
const std::string units = "cm";
  
// z extension of barrel
const double zMin = -3100.*mm;
const double zMax =  3100.*mm;

// radial extension of barrel
const double rMin = 217.28*cm;
const double rMax = 257.33*cm;

// nominal radial thickness of layers 
const std::vector<double> drNom(
  {
    15.*mm, 35.*mm, 35.*mm, 35.*mm, 35.*mm, 35.*mm, 
    35.*mm, 35.*mm, 35.*mm, 35.*mm, 35.*mm, 35.*mm
  }
);

// number of layers
const int nLayers = drNom.size();

// number of electrodes
const int nModules = 1545;

// inclination angle of electrodes
const double alpha = 50*TMath::Pi()/180.;

// grid in theta
// - size of cell
const double thetaGrid = 0.00981748/4;
// - theta of edge of first cell
const double thetaMin = 0.59027850 - thetaGrid/2.;
// - n(bins)
const int nThetaBins = 800;



/******************************************************************************/
// GEOMETRY HELPER FUNCTIONS AND DERIVED GEOMETRY PARAMETERS
/******************************************************************************/

// calculate length along electrode at r=r_out for electrode
// starting at r=r_in and inclined in phi by alpha
// for r_in>r_out, when the solution exists, it's always > 0
double _L(double alpha, double r_in, double r_out) {
  if (r_out > r_in)
    return sqrt(r_out*r_out - r_in*r_in*sin(alpha)*sin(alpha)) - r_in*cos(alpha);
  else
    return r_in*cos(alpha) - sqrt(r_out*r_out - r_in*r_in*sin(alpha)*sin(alpha));
}

// calculate phi offset between point at r=r_out with respect 
// to point at r=r_in if they are on the same electrode 
// inclined in phi by alpha
double _phi(double alpha, double r_in, double r_out) {
  double L = _L(alpha, r_in, r_out);
  if (r_out > r_in)
    return TMath::ASin(L/r_out * sin(alpha));
  else
    return -TMath::ASin(L/r_out * sin(alpha));
}

// calculate radius of point on electrode starting at 
// r=r_in, inclined in phi by alpha, and at distance L 
// from beginning of electrode
double _r(double alpha, double r_in, double L) {
  return sqrt( (r_in+L*cos(alpha)) * (r_in+L*cos(alpha))
	       + (L*sin(alpha)) * (L*sin(alpha)) );
}

// length of electrodes
static const double Ltot = _L(alpha, rMin, rMax);

// phi separation of modules
const double gridPhi = TMath::TwoPi()/nModules;

// r for electrode length = L/2
static const double rAvg = _r(alpha, rMin, Ltot/2.0);

// delta phi of point at L/2
static const double dPhiAvg = _phi(alpha, rMin, rAvg);

// phi edge of module 0
static const double phiMin = -alpha - gridPhi/2.0;

// other quantities, calculated by calcGeom
// radial position of each layer
std::vector<double> r;
std::vector<double> dr;
// length of electrode along each layer (dRnom/cos(alpha))
std::vector<double> dL;

// calculate derived parameters of geometry depending on the main ones
// (such as radial/electrode length of each layer)
// and print them to screen
void calcGeom() {

  // print initial information
  cout << "r(min) = " << rMin << " " << units << endl;
  cout << "r(max) = " << rMax << " " << units << endl;
  cout << "total thickness = " << rMax-rMin << " " << units << endl;

  // calculate total length
  cout << "electrode length = " << Ltot << " " << units << endl << endl;

  // calculate length along electrode of each layer  
  cout << "n(layers) = " << nLayers << endl << endl;
  std::vector<double> L;
  double sumL(0.0);
  dL.clear();
  for (int i=0; i<nLayers; i++) {
    dL.push_back(drNom[i]/cos(alpha));
    sumL+=dL[i];
  }
  cout << "electrode length/layer: " << endl;
  L.push_back(0.0);
  for (int i=0; i<nLayers; i++) {
    dL[i] = dL[i]*Ltot/sumL;
    cout <<  i  <<  " " << dL[i] << " " << units << endl;
    L.push_back(L[i] + dL[i]);
  }
  cout << endl;
  
  // calculate r and dr of each layer
  cout << "radial position of each layer: " << endl;
  r.resize(nLayers+1);
  dr.resize(nLayers);
  r[0] = rMin;
  cout <<  "0 " << r[0] << " " << units << endl;
  for (int iLayer=1; iLayer<=nLayers; iLayer++) {
    r[iLayer] = _r(alpha, r[0], L[iLayer]);
    dr[iLayer-1] = r[iLayer]-r[iLayer-1];
    cout <<  iLayer  <<  " " << r[iLayer] << " " << units << endl;
  }
  cout << endl;

  // print calculated thickness of each layer
  cout << "radial thickness of each layer: " << endl;
  for (int iLayer=0; iLayer<nLayers; iLayer++) {
    cout <<  iLayer  <<  " " << dr[iLayer] << " " << units << endl;
  }
  cout << endl;

  // calculate radial position of points at middle length of
  // electrode in each layer
  cout << "r of center " << endl;
  std::vector<double> rMed;
  for (int i=0; i<nLayers; i++) {
    rMed.push_back(_r(alpha, r[0], (L[i]+L[i+1])/2.0));
    cout <<  i  <<  " " << rMed[i] << " " << units << endl;
  }
  cout << "r(avg) = " << rAvg << " " << units << endl;

  cout << "phi grid : " << gridPhi << endl;
  cout << "phi offset of electrode mid point wrt initial point : " << dPhiAvg << endl;
  cout << "initial phi of edge of module 0 : " << phiMin << endl;
}


/******************************************************************************/
// HELPER FUNCTIONS related to the readout
/******************************************************************************/

// extract layer number from cellID
ULong_t Layer(ULong_t cellID) {
  const ULong_t mask = (1<<8) -1;
  return (cellID >> 11) & mask;
}

// extract module number from cellID
ULong_t Module(ULong_t cellID) {
  const ULong_t mask = (1<<11) -1;
  return (cellID >> 19) & mask;
}

// extract theta bin from cellID
ULong_t ThetaBin(ULong_t cellID) {
  const ULong_t mask = (1<<10) -1;
  return (cellID >> 30) & mask;
}

// return the sign of a float
int sgn(float val) {
  return (val > 0.) - (val < 0.);
}


/******************************************************************************/
// HELPER FUNCTIONS related to the graphic system
/******************************************************************************/

void makeGui();

// derived TGLAnnotation class that overrides MouseEnter method
// to avoid editing of annotation
class TGLConstAnnotation : public TGLAnnotation
{
 public:
  TGLConstAnnotation(TGLViewerBase *parent, const char *text, Float_t posx, Float_t posy) :
  TGLAnnotation(parent, text, posx, posy)
    {
      ;
    }
  Bool_t MouseEnter(TGLOvlSelectRecord& /*rec*/)
    {
      fActive = kFALSE;
      return kTRUE;
    }
};


/******************************************************************************/
// ELEMENTS OF THE EVENT DISPLAY
/******************************************************************************/

Int_t eventId = 0; // Current event id
Int_t nEvents = 0; // Number of events in file
TEveTrackList* particles = nullptr;
TEvePointSet* hits = nullptr;
TEvePointSet* cells = nullptr;
TEvePointSet* cells_merged = nullptr;
TEvePointSet* clusters = nullptr;
std::vector<TEveQuadSet*> qs_rhoz;
std::vector<TEveQuadSet*> qs_rhophi;
TEveElementList* topoclusters_rhoz = nullptr;
TEveElementList* topoclusters_rhophi = nullptr;
std::vector<TEveBoxSet*> bs;
TEveElementList* topoclusters_3D = nullptr;

TEveGeoShape* barrel  = nullptr;
TEveViewer* rhoPhiView = nullptr;
TEveViewer* rhoZView = nullptr;
TGLViewer* rhoPhiGLView = nullptr;
TGLViewer* rhoZGLView = nullptr;
TEveScene* rhoPhiScene = nullptr;
TEveScene* rhoPhiEventScene = nullptr;
TEveScene* rhoPhiEventSceneManual = nullptr;
TEveScene* rhoZScene = nullptr;
TEveScene* rhoZEventScene = nullptr;
TEveScene* rhoZEventSceneManual = nullptr;
TEveProjectionManager* rhoPhiProjManager = nullptr;
TEveProjectionManager* rhoZProjManager = nullptr;

TGLConstAnnotation* eventLabel = nullptr;
TGTextEntry *textEntry   = nullptr;

/******************************************************************************/
// CLASS TO READ AN EVENT FROM THE EVENT FILE
// AND FILL THE EVENT DISPLAY OBJECTS
/******************************************************************************/
class EventReader
{
 private:
  TTreeReader* fReader = nullptr;
  // 1a. - primary particles
  TTreeReaderArray<Int_t> *genParticles_PDG = nullptr;
  TTreeReaderArray<Int_t> *genParticles_generatorStatus = nullptr;
  TTreeReaderArray<Int_t> *genParticles_simulatorStatus = nullptr;
  TTreeReaderArray<Float_t> *genParticles_charge = nullptr;
  TTreeReaderArray<Float_t> *genParticles_time = nullptr;
  TTreeReaderArray<Double_t> *genParticles_mass = nullptr;
  TTreeReaderArray<Double_t> *genParticles_vertex_x = nullptr;
  TTreeReaderArray<Double_t> *genParticles_vertex_y = nullptr;
  TTreeReaderArray<Double_t> *genParticles_vertex_z = nullptr;
  //TTreeReaderArray<Double_t> *genParticles_endpoint_x = nullptr;
  //TTreeReaderArray<Double_t> *genParticles_endpoint_y = nullptr;
  //TTreeReaderArray<Double_t> *genParticles_endpoint_z = nullptr;
  TTreeReaderArray<Float_t> *genParticles_momentum_x = nullptr;
  TTreeReaderArray<Float_t> *genParticles_momentum_y = nullptr;
  TTreeReaderArray<Float_t> *genParticles_momentum_z = nullptr;

  // 1a. - secondary particles
  TTreeReaderArray<Int_t> *SimParticleSecondaries_PDG = nullptr;
  //TTreeReaderArray<Int_t> *SimParticleSecondaries_generatorStatus = nullptr;
  //TTreeReaderArray<Int_t> *SimParticleSecondaries_simulatorStatus = nullptr;
  //TTreeReaderArray<Float_t> *SimParticleSecondaries_charge = nullptr;
  TTreeReaderArray<Float_t> *SimParticleSecondaries_time = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_mass = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_vertex_x = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_vertex_y = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_vertex_z = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_endpoint_x = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_endpoint_y = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_endpoint_z = nullptr;
  TTreeReaderArray<Float_t> *SimParticleSecondaries_momentum_x = nullptr;
  TTreeReaderArray<Float_t> *SimParticleSecondaries_momentum_y = nullptr;
  TTreeReaderArray<Float_t> *SimParticleSecondaries_momentum_z = nullptr;

  // 1b. hits
  TTreeReaderArray<ULong_t> *ECalBarrelPositionedHits_cellID = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedHits_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedHits_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedHits_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedHits_position_z = nullptr;

  // 1c. - cells
  TTreeReaderArray<ULong_t> *ECalBarrelPositionedCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells_position_z = nullptr;

  // 1d. - cells with coarser merging
  TTreeReaderArray<ULong_t> *ECalBarrelPositionedCells2_cellID = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells2_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells2_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells2_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells2_position_z = nullptr;

  // 1e. - the corrected calo topo clusters
  TTreeReaderArray<Float_t> *CorrectedCaloTopoClusters_energy = nullptr;
  TTreeReaderArray<Float_t> *CorrectedCaloTopoClusters_position_x = nullptr;
  TTreeReaderArray<Float_t> *CorrectedCaloTopoClusters_position_y = nullptr;
  TTreeReaderArray<Float_t> *CorrectedCaloTopoClusters_position_z = nullptr;
  TTreeReaderArray<UInt_t> *CorrectedCaloTopoClusters_hits_begin = nullptr;
  TTreeReaderArray<UInt_t> *CorrectedCaloTopoClusters_hits_end = nullptr;  
  
  // 1f. - cells in the topo clusters
  TTreeReaderArray<ULong_t> *PositionedCaloTopoClusterCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *PositionedCaloTopoClusterCells_energy = nullptr;
  TTreeReaderArray<Float_t> *PositionedCaloTopoClusterCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *PositionedCaloTopoClusterCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *PositionedCaloTopoClusterCells_position_z = nullptr;

 public:
  EventReader(TFile* f)
    {
      fReader = new TTreeReader("events", f);
      nEvents = fReader->GetEntries();
      
      // 1a. - primary particles
      genParticles_PDG = new TTreeReaderArray<Int_t>(*fReader, "genParticles.PDG");
      genParticles_generatorStatus = new TTreeReaderArray<Int_t>(*fReader, "genParticles.generatorStatus");
      genParticles_simulatorStatus = new TTreeReaderArray<Int_t>(*fReader, "genParticles.simulatorStatus");
      genParticles_charge = new TTreeReaderArray<Float_t>(*fReader, "genParticles.charge");
      genParticles_time = new TTreeReaderArray<Float_t>(*fReader, "genParticles.time");
      genParticles_mass = new TTreeReaderArray<Double_t>(*fReader, "genParticles.mass");
      genParticles_vertex_x = new TTreeReaderArray<Double_t>(*fReader, "genParticles.vertex.x");
      genParticles_vertex_y = new TTreeReaderArray<Double_t>(*fReader, "genParticles.vertex.y");
      genParticles_vertex_z = new TTreeReaderArray<Double_t>(*fReader, "genParticles.vertex.z");
      //genParticles_endpoint_x = new TTreeReaderArray<Double_t>(*fReader, "genParticles.endpoint.x");
      //genParticles_endpoint_y = new TTreeReaderArray<Double_t>(*fReader, "genParticles.endpoint.y");
      //genParticles_endpoint_z = new TTreeReaderArray<Double_t>(*fReader, "genParticles.endpoint.z");
      genParticles_momentum_x = new TTreeReaderArray<Float_t>(*fReader, "genParticles.momentum.x");
      genParticles_momentum_y = new TTreeReaderArray<Float_t>(*fReader, "genParticles.momentum.y");
      genParticles_momentum_z = new TTreeReaderArray<Float_t>(*fReader, "genParticles.momentum.z");

      // 1a. - secondary particles
      if (fReader->GetTree()->FindBranch("SimParticleSecondaries.PDG")) {
	SimParticleSecondaries_PDG = new TTreeReaderArray<Int_t>(*fReader, "SimParticleSecondaries.PDG");
	//SimParticleSecondaries_generatorStatus = new TTreeReaderArray<Int_t>(*fReader, "SimParticleSecondaries.generatorStatus");
	//SimParticleSecondaries_simulatorStatus = new TTreeReaderArray<Int_t>(*fReader, "SimParticleSecondaries.simulatorStatus");
	//SimParticleSecondaries_charge = new TTreeReaderArray<Float_t>(*fReader, "SimParticleSecondaries.charge");
	SimParticleSecondaries_time = new TTreeReaderArray<Float_t>(*fReader, "SimParticleSecondaries.time");
	SimParticleSecondaries_mass = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.mass");
	SimParticleSecondaries_vertex_x = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.vertex.x");
	SimParticleSecondaries_vertex_y = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.vertex.y");
	SimParticleSecondaries_vertex_z = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.vertex.z");
	SimParticleSecondaries_endpoint_x = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.endpoint.x");
	SimParticleSecondaries_endpoint_y = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.endpoint.y");
	SimParticleSecondaries_endpoint_z = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.endpoint.z");
	SimParticleSecondaries_momentum_x = new TTreeReaderArray<Float_t>(*fReader, "SimParticleSecondaries.momentum.x");
	SimParticleSecondaries_momentum_y = new TTreeReaderArray<Float_t>(*fReader, "SimParticleSecondaries.momentum.y");
	SimParticleSecondaries_momentum_z = new TTreeReaderArray<Float_t>(*fReader, "SimParticleSecondaries.momentum.z");
      }
	
      // 1b. - hits
      ECalBarrelPositionedHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, "ECalBarrelPositionedHits.cellID");
      ECalBarrelPositionedHits_energy     = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedHits.energy");
      ECalBarrelPositionedHits_position_x = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedHits.position.x");
      ECalBarrelPositionedHits_position_y = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedHits.position.y");
      ECalBarrelPositionedHits_position_z = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedHits.position.z");
      
      // 1c. - cells
      ECalBarrelPositionedCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, "ECalBarrelPositionedCells.cellID");
      ECalBarrelPositionedCells_energy     = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells.energy");
      ECalBarrelPositionedCells_position_x = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells.position.x");
      ECalBarrelPositionedCells_position_y = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells.position.y");
      ECalBarrelPositionedCells_position_z = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells.position.z");
      
      // 1d. - cells with coarser merging
      ECalBarrelPositionedCells2_cellID     = new TTreeReaderArray<ULong_t>(*fReader, "ECalBarrelPositionedCells2.cellID");
      ECalBarrelPositionedCells2_energy     = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells2.energy");
      ECalBarrelPositionedCells2_position_x = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells2.position.x");
      ECalBarrelPositionedCells2_position_y = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells2.position.y");
      ECalBarrelPositionedCells2_position_z = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells2.position.z");
      
      // 1e. - the corrected calo topo clusters
      CorrectedCaloTopoClusters_energy     = new TTreeReaderArray<Float_t>(*fReader, "CorrectedCaloTopoClusters.energy");
      CorrectedCaloTopoClusters_position_x = new TTreeReaderArray<Float_t>(*fReader, "CorrectedCaloTopoClusters.position.x");
      CorrectedCaloTopoClusters_position_y = new TTreeReaderArray<Float_t>(*fReader, "CorrectedCaloTopoClusters.position.y");
      CorrectedCaloTopoClusters_position_z = new TTreeReaderArray<Float_t>(*fReader, "CorrectedCaloTopoClusters.position.z");
      CorrectedCaloTopoClusters_hits_begin = new TTreeReaderArray<UInt_t>(*fReader, "CorrectedCaloTopoClusters.hits_begin");
      CorrectedCaloTopoClusters_hits_end   = new TTreeReaderArray<UInt_t>(*fReader, "CorrectedCaloTopoClusters.hits_end");

      // 1f. - cells in the topo clusters
      PositionedCaloTopoClusterCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, "PositionedCaloTopoClusterCells.cellID");
      PositionedCaloTopoClusterCells_energy     = new TTreeReaderArray<Float_t>(*fReader, "PositionedCaloTopoClusterCells.energy");
      PositionedCaloTopoClusterCells_position_x = new TTreeReaderArray<Float_t>(*fReader, "PositionedCaloTopoClusterCells.position.x");
      PositionedCaloTopoClusterCells_position_y = new TTreeReaderArray<Float_t>(*fReader, "PositionedCaloTopoClusterCells.position.y");
      PositionedCaloTopoClusterCells_position_z = new TTreeReaderArray<Float_t>(*fReader, "PositionedCaloTopoClusterCells.position.z");

      qs_rhoz.clear();
      qs_rhophi.clear();      
    }

  ~EventReader() {
    delete genParticles_PDG;
    delete genParticles_generatorStatus;
    delete genParticles_simulatorStatus;
    delete genParticles_charge;
    delete genParticles_time;
    delete genParticles_mass;
    delete genParticles_vertex_x;
    delete genParticles_vertex_y;
    delete genParticles_vertex_z;
    //delete genParticles_endpoint_x;
    //delete genParticles_endpoint_y;
    //delete genParticles_endpoint_z;
    delete genParticles_momentum_x;
    delete genParticles_momentum_y;
    delete genParticles_momentum_z;

    if (SimParticleSecondaries_PDG) {
      delete SimParticleSecondaries_PDG;
      // delete SimParticleSecondaries_generatorStatus;
      // delete SimParticleSecondaries_simulatorStatus;
      // delete SimParticleSecondaries_charge;
      delete SimParticleSecondaries_time;
      delete SimParticleSecondaries_mass;
      delete SimParticleSecondaries_vertex_x;
      delete SimParticleSecondaries_vertex_y;
      delete SimParticleSecondaries_vertex_z;
      delete SimParticleSecondaries_endpoint_x;
      delete SimParticleSecondaries_endpoint_y;
      delete SimParticleSecondaries_endpoint_z;
      delete SimParticleSecondaries_momentum_x;
      delete SimParticleSecondaries_momentum_y;
      delete SimParticleSecondaries_momentum_z;
    }
    delete ECalBarrelPositionedHits_cellID;
    delete ECalBarrelPositionedHits_energy;
    delete ECalBarrelPositionedHits_position_x;
    delete ECalBarrelPositionedHits_position_y;
    delete ECalBarrelPositionedHits_position_z;
    delete ECalBarrelPositionedCells_cellID;
    delete ECalBarrelPositionedCells_energy;
    delete ECalBarrelPositionedCells_position_x;
    delete ECalBarrelPositionedCells_position_y;
    delete ECalBarrelPositionedCells_position_z;
    delete ECalBarrelPositionedCells2_cellID;
    delete ECalBarrelPositionedCells2_energy;
    delete ECalBarrelPositionedCells2_position_x;
    delete ECalBarrelPositionedCells2_position_y;
    delete ECalBarrelPositionedCells2_position_z;
    delete CorrectedCaloTopoClusters_energy;
    delete CorrectedCaloTopoClusters_position_x;
    delete CorrectedCaloTopoClusters_position_y;
    delete CorrectedCaloTopoClusters_position_z;
    delete CorrectedCaloTopoClusters_hits_begin;
    delete CorrectedCaloTopoClusters_hits_end;
    delete PositionedCaloTopoClusterCells_cellID;
    delete PositionedCaloTopoClusterCells_energy;
    delete PositionedCaloTopoClusterCells_position_x;
    delete PositionedCaloTopoClusterCells_position_y;
    delete PositionedCaloTopoClusterCells_position_z;
    delete fReader;
  }

  void loadEvent(int event = -1) {

    if (event != -1) eventId = event;

    printf("Loading event %d.\n", eventId);
    textEntry->SetTextColor(0xff0000);
    textEntry->SetText(Form("Loading event %d...",eventId));

    fReader->SetEntry(eventId);

    TString partType;
    double pmax=0.0;
    int ipmax=-1;
    for (unsigned int i = 0; i < genParticles_generatorStatus->GetSize(); i ++) {
      if ( (*genParticles_generatorStatus)[i] != 1 ) continue;
      float _px = (*genParticles_momentum_x)[i];
      float _py = (*genParticles_momentum_y)[i];
      float _pz = (*genParticles_momentum_z)[i];
      float _p = sqrt( _px*_px + _py*_py + _pz*_pz );
      //cout << p << endl;
      if (_p==0.) continue;
      if (_p>pmax) {
	pmax = _p;
	ipmax = i;
      }
    }
    int pdgID = (*genParticles_PDG)[ipmax];
    if (pdgID == 111) partType = "pi0";
    else if (pdgID == 22) partType = "y";
    else if (pdgID == 11) partType = "e-";
    else if (pdgID == -11) partType = "e+";
    else partType = "unknown";
    
    //
    // particles
    //
    if (drawParticles) {
      cout << "Creating particles" << endl;
      if (particles == nullptr) {
	particles = new TEveTrackList("particles");
	TEveTrackPropagator* trkProp = particles->GetPropagator();
	trkProp->SetMagField( 0.01 );
	particles->SetMainColor(kYellow);
	particles->SetLineWidth(2);
	gEve->AddElement(particles);
      }
      else
	particles->DestroyElements();

      // handle differently the e/gamma vs pi0 particle guns
      // (for pi0, need to look for the two photons in secondary particle list
      // unfortunately cannot just use the latter to show particles since
      // info about endpoint is broken (identical to vertex)
      float m = (*genParticles_mass)[ipmax];
      float px = (*genParticles_momentum_x)[ipmax];
      float py = (*genParticles_momentum_y)[ipmax];
      float pz = (*genParticles_momentum_z)[ipmax];
      float p = sqrt( px*px + py*py + pz*pz );
	
      double t = (*genParticles_time)[ipmax];
      double x1 = (*genParticles_vertex_x)[ipmax] * mm;
      double y1 = (*genParticles_vertex_y)[ipmax] * mm;
      double z1 = (*genParticles_vertex_z)[ipmax] * mm;
      double r1 = sqrt(x1*x1+y1*y1);
      double sintheta = sqrt(px*px + py*py)/p;
      double x2 = x1 + px/p * rMax / sintheta;
      double y2 = y1 + py/p * rMax / sintheta;
      double z2 = z1 + pz/p * rMax / sintheta;
      double r2 = sqrt(x2*x2+y2*y2);

      TEveMCTrack mct;
      mct.SetPdgCode( pdgID );
      mct.SetMomentum( px, py, pz, sqrt(p*p + m*m) );
      mct.SetProductionVertex( x1, y1, z1, t );
      TEveTrack* track = new TEveTrack(&mct, particles->GetPropagator());
      track->SetAttLineAttMarker(particles);
      track->SetElementTitle(Form("p = %.3f GeV\ntheta = %f\nphi = %f\nx = %f cm\ny = %f cm\nz= %f cm",
				  p, acos(pz/p), atan2(py,px),
				  x1/cm, y1/cm, z1/cm));
      particles->AddElement(track);

      // if the particle is a pi0, also draw the two photons, and set the endpoint
      // of the pi0 track
      if (pdgID == 111) {
	bool decayVtxSet = false;
	for (unsigned int i = 0; i < SimParticleSecondaries_PDG->GetSize(); i ++) {
	  pdgID = (*SimParticleSecondaries_PDG)[i];
	  // keep only photons
	  if (pdgID!=22) continue;
	  //if ( (*SimParticleSecondaries_generatorStatus)[i] != 1 ) continue;
	  px = (*SimParticleSecondaries_momentum_x)[i];
	  py = (*SimParticleSecondaries_momentum_y)[i];
	  pz = (*SimParticleSecondaries_momentum_z)[i];
	  m = (*SimParticleSecondaries_mass)[i];
	  p = sqrt( px*px + py*py + pz*pz );
	  float e = sqrt(p*p + m*m);
	  // cout << "p = "<< p << endl;
	  if (p<ParticleEnergyThreshold) continue;

	  // cout << "PDG = "<< pdgID << endl;
	  t = (*SimParticleSecondaries_time)[i];
	  x1 = (*SimParticleSecondaries_vertex_x)[i] * mm;
	  y1 = (*SimParticleSecondaries_vertex_y)[i] * mm;
	  z1 = (*SimParticleSecondaries_vertex_z)[i] * mm;
	  r1 = sqrt(x1*x1+y1*y1);
	  // the two photons from a pi0 in the origin must come from small R
	  if (r1>1.) continue;
	  sintheta = sqrt(px*px + py*py)/p;
	  x2 = x1 + px/p * rMax / sintheta;
	  y2 = y1 + py/p * rMax / sintheta;
	  z2 = z1 + pz/p * rMax / sintheta;
	  //double x2 = (*SimParticleSecondaries_endpoint_x)[i] * mm;
	  //double y2 = (*SimParticleSecondaries_endpoint_y)[i] * mm;
	  //double z2 = (*SimParticleSecondaries_endpoint_z)[i] * mm;
	  r2 = sqrt(x2*x2+y2*y2);
	  // cout << "x1 y1 z1 x2 y2 z2 = "
	  //      << x1 << " " << y1 << " " << z1 << " "
	  //      << x2 << " " << y2 << " " << z2 << endl;
	  // set pi0 decay point
	  if (!decayVtxSet) {
	    TEveVectorF v; v[0] = x1; v[1]=y1; v[2]=z1;
	    TEvePathMark mark(TEvePathMark::kDecay, v);
	    track->AddPathMark(mark);
	    decayVtxSet = true;
	  }
	  TEveMCTrack mct;
	  mct.SetPdgCode( pdgID );
	  mct.SetMomentum( px, py, pz, sqrt(p*p + m*m) );
	  mct.SetProductionVertex( x1, y1, z1, t );
	  TEveTrack* track = new TEveTrack(&mct, particles->GetPropagator());
	  track->SetAttLineAttMarker(particles);
	  track->SetElementTitle(Form("p = .3%f GeV\ntheta = %f\nphi = %f\nx = %f cm\ny = %f cm\nz= %f cm",
				      p, acos(pz/p), atan2(py,px),
				      x1/cm, y1/cm, z1/cm));
	  
	  particles->AddElement(track);
	}
      }
      particles->MakeTracks();
    }
    
    //
    // hits
    //
    if (drawHits) {
      cout << "Creating hits" << endl;
      if (hits == nullptr) {
	hits = new TEvePointSet();
	hits->SetName(Form("hits (E>%.1f GeV)",CellEnergyThreshold));
	hits->SetMarkerStyle(4);
	hits->SetMarkerSize(1);
	hits->SetMarkerColor(kRed);
	gEve->AddElement(hits);
      }
      else
	hits->Reset();
      for (unsigned int i = 0; i < ECalBarrelPositionedHits_position_x->GetSize(); i ++) {
	float E = (*ECalBarrelPositionedHits_energy)[i];
	if (E<HitEnergyThreshold) continue;
	// ULong_t cellID = (*ECalBarrelPositionedHits_cellID)[i];
	// ULong_t layer = Layer(cellID);
	hits->SetNextPoint( (*ECalBarrelPositionedHits_position_x)[i] * mm ,
			    (*ECalBarrelPositionedHits_position_y)[i] * mm ,
			    (*ECalBarrelPositionedHits_position_z)[i] * mm );
      }
    }
    
    //
    // cells
    //
    if (drawCells) {
      cout << "Creating cells" << endl; 
      if (cells == nullptr) {
	cells = new TEvePointSet();
	cells->SetName(Form("cells (E>%.1f GeV)",CellEnergyThreshold));
	cells->SetMarkerStyle(4);
	cells->SetMarkerSize(3);
	cells->SetMarkerColor(kBlue);
	gEve->AddElement(cells);
      }
      else
	cells->Reset();
      for (unsigned int i = 0; i < ECalBarrelPositionedCells_position_x->GetSize(); i ++) {
	float E = (*ECalBarrelPositionedCells_energy)[i];
	if (E<CellEnergyThreshold) continue;
	// ULong_t cellID = (*ECalBarrelPositionedCells_cellID)[i];
	// ULong_t layer = Layer(cellID);
	cells->SetNextPoint( (*ECalBarrelPositionedCells_position_x)[i] * mm ,
			     (*ECalBarrelPositionedCells_position_y)[i] * mm ,
			     (*ECalBarrelPositionedCells_position_z)[i] * mm );
      }
    }
    
    //
    // cells merged 
    //
    if (drawMergedCells) {
      cout << "Creating merged cells" << endl;
      if (cells_merged == nullptr) {
	cells_merged = new TEvePointSet();
	cells_merged->SetName("cells_merged");
	cells_merged->SetMarkerStyle(4);
	cells_merged->SetMarkerSize(3);
	cells_merged->SetMarkerColor(kBlue);
	gEve->AddElement(cells_merged);
      }
      else
	cells_merged->Reset();
      for (unsigned int i = 0; i < ECalBarrelPositionedCells2_position_x->GetSize(); i ++) {
	float E = (*ECalBarrelPositionedCells2_energy)[i];
	// if (E<minCellE) continue;
	// ULong_t cellID = (*ECalBarrelPositionedCells_cellID)[i];
	// ULong_t layer = Layer(cellID);
	cells_merged->SetNextPoint( (*ECalBarrelPositionedCells2_position_x)[i] * mm ,
				    (*ECalBarrelPositionedCells2_position_y)[i] * mm ,
				    (*ECalBarrelPositionedCells2_position_z)[i] * mm );
      }
    }


    //
    // clusters
    //
    if (drawTopoClusters) {
      cout << "Creating clusters" << endl;

      // centers of the clusters
      if (clusters == nullptr) {
	clusters = new TEvePointSet();
	clusters->SetName(Form("clusters (E>%.1f GeV)",TopoClusterEnergyThreshold));
	clusters->SetMarkerStyle(4);
	clusters->SetMarkerSize(6);
	clusters->SetMarkerColor(kGreen);
	gEve->AddElement(clusters);
      }
      else
	clusters->Reset();
      
      for (unsigned int i = 0; i < CorrectedCaloTopoClusters_position_x->GetSize(); i ++) {
	float E = (*CorrectedCaloTopoClusters_energy)[i];
	if (E < TopoClusterEnergyThreshold) continue;
	// cluster positions are in cm and hits/cells in mm ...
	clusters->SetNextPoint( (*CorrectedCaloTopoClusters_position_x)[i] * cm ,
				(*CorrectedCaloTopoClusters_position_y)[i] * cm ,
				(*CorrectedCaloTopoClusters_position_z)[i] * cm ); 
      }

      
      // the DestroyElements method does not work for TExeQuadSet as for the TEvePointSet..
      // so I have to destroy and recreate the quad and box sets
      for (auto qs : qs_rhoz) {
	if (qs)
	  qs->Destroy();
      }
      qs_rhoz.clear();
      
      for (auto qs : qs_rhophi) {
	if (qs)
	  qs->Destroy();
      }
      qs_rhophi.clear();

      for (auto _bs : bs) {
	if (_bs)
	  _bs->Destroy();
      }
      bs.clear();

      TEveRGBAPalette *pal = new TEveRGBAPalette(0, 1000);
      
      // clusters in 3D
      if (topoclusters_3D==nullptr) {
	topoclusters_3D = new TEveElementList("Clusters in 3D (no E cut)");
	gEve->AddElement(topoclusters_3D);
      }
      else
	topoclusters_3D->DestroyElements();
      
      // clusters in 2D
      if (topoclusters_rhoz==nullptr) {
	topoclusters_rhoz = new TEveElementList(Form("Clusters in rho-z (E>%.1f GeV)",
						     TopoClusterEnergyThreshold));
	// add to scene or manager?
	// rhoZProjManager->AddElement(topoclusters_rhoz);
	rhoZEventSceneManual->AddElement(topoclusters_rhoz);
	gEve->AddToListTree(topoclusters_rhoz,false);
      }
      else
	topoclusters_rhoz->DestroyElements();

      if (topoclusters_rhophi==nullptr) {
	topoclusters_rhophi = new TEveElementList(Form("Clusters in rho-phi (E>%.1f GeV)",
						       TopoClusterEnergyThreshold));
	// add to scene or manager? -- scene that is not auto-projected!
	// rhoPhiProjManager->AddElement(topoclusters_rhophi);
	rhoPhiEventSceneManual->AddElement(topoclusters_rhophi);
	gEve->AddToListTree(topoclusters_rhophi,false);
      }
      else
	topoclusters_rhophi->DestroyElements();

      for (unsigned int i = 0; i < CorrectedCaloTopoClusters_energy->GetSize(); i ++) {
	float energy = (*CorrectedCaloTopoClusters_energy)[i];
	float xcl = (*CorrectedCaloTopoClusters_position_x)[i];
	float ycl = (*CorrectedCaloTopoClusters_position_y)[i];
	float zcl = (*CorrectedCaloTopoClusters_position_z)[i];
	float rcl = sqrt(xcl*xcl + ycl*ycl);
	float phicl = atan2(ycl, xcl);
	float thetacl = atan2(rcl, zcl);

	if (energy < TopoClusterEnergyThreshold) {
	  qs_rhoz.push_back(nullptr);
	  qs_rhophi.push_back(nullptr);
	}
	else {
	  TEveQuadSet* aqs = new TEveQuadSet(TEveQuadSet::kQT_FreeQuad, false, 32,
					     Form("cluster %d", (int) i));
	  aqs->SetMainTransparency(80);
	  // by calling SetOwnIds(kTRUE) the digit-set becomes
	  // the owner of the assigned objects and deletes
	  // them on destruction.
	  aqs->SetOwnIds(kTRUE);
	  aqs->SetPalette(pal);
	  aqs->SetTitle(Form("E = %f GeV\nR = %f cm\ntheta = %f\nphi = %f",
			     energy,
			     rcl,
			     thetacl,
			     phicl));
	  qs_rhoz.push_back(aqs);
	  topoclusters_rhoz->AddElement(aqs);
	  
	  TEveQuadSet* aqs2 = new TEveQuadSet(TEveQuadSet::kQT_FreeQuad, false, 32,
					      Form("cluster %d", (int) i));
	  aqs2->SetMainTransparency(80);
	  aqs2->SetOwnIds(kTRUE);
	  aqs2->SetPalette(pal);
	  aqs2->SetTitle(Form("E = %f GeV\nR = %f cm\ntheta = %f\nphi = %f",			    
			      energy,
			      rcl,
			      thetacl,
			      phicl));
	  qs_rhophi.push_back(aqs2);
	  topoclusters_rhophi->AddElement(aqs2);
	}
	
	TEveBoxSet* _bs = new TEveBoxSet(Form("cluster %d", (int) i),"");
	_bs->Reset(TEveBoxSet::kBT_FreeBox, false, 32);
	_bs->SetMainTransparency(80);
	_bs->SetOwnIds(kTRUE);
	_bs->SetPalette(pal);
	_bs->SetTitle(Form("E = %f GeV\nR = %f cm\ntheta = %f\nphi = %f",
			   energy,
			   rcl,
			   thetacl,
			   phicl));
	bs.push_back(_bs);
	topoclusters_3D->AddElement(_bs);
      }
      
      // loop over cells and attach them to clusters
      for (unsigned int i = 0; i < PositionedCaloTopoClusterCells_energy->GetSize(); i ++) {
	int icl = -1;
	for (unsigned int j = 0; j < CorrectedCaloTopoClusters_energy->GetSize(); j ++) {
	  if (i >= (*CorrectedCaloTopoClusters_hits_begin)[j] &&
	      i < (*CorrectedCaloTopoClusters_hits_end)[j]) {
	    icl = j;
	    break;
	  }
	}
	if (icl==-1) continue; // should never happen..
	ULong_t cellID = (*PositionedCaloTopoClusterCells_cellID)[i];
	int layer = (int) Layer(cellID);
	float x_center = (*PositionedCaloTopoClusterCells_position_x)[i] * mm;
	float y_center = (*PositionedCaloTopoClusterCells_position_y)[i] * mm;
	float z_center = (*PositionedCaloTopoClusterCells_position_z)[i] * mm;
	float r_center = sqrt(x_center*x_center + y_center*y_center);
	float r_in = r[layer];
	float r_out = r[layer+1];
	float theta_center = atan2(r_center, z_center);
	float phi_center = atan2(y_center, x_center);
	float energy = (*PositionedCaloTopoClusterCells_energy)[i];
	
	// cluster cells in rho-z projection
	float verts[12];
	verts[0] = r_in/tan(theta_center - thetaGrid*mergedCells_Theta[layer]/2.);
	verts[1] = r_in*sgn(y_center);
	verts[3] = r_in/tan(theta_center + thetaGrid*mergedCells_Theta[layer]/2.);
	verts[4] = r_in*sgn(y_center);
	verts[6] = r_out/tan(theta_center + thetaGrid*mergedCells_Theta[layer]/2.);
	verts[7] = r_out*sgn(y_center);
	verts[9] = r_out/tan(theta_center - thetaGrid*mergedCells_Theta[layer]/2.);
	verts[10] = r_out*sgn(y_center);
	verts[2] = verts[5] = verts[8] = verts[11] = 0.;
	if (qs_rhoz[icl]!=nullptr) {
	  qs_rhoz[icl]->AddQuad(verts);
	  qs_rhoz[icl]->QuadValue( (int) (1000 * energy) );
	  qs_rhoz[icl]->QuadId( new TNamed(Form("Cell %lu", cellID), "Dong!") );
	}
	// cluster cells in rho-phi projection
	int module = (int) Module(cellID);
	double Lin = _L(alpha, rMin, r[layer]);
	double Lout = _L(alpha, rMin, r[layer+1]);
	double deltaL = rMin*sin(gridPhi/2.0)*sin(alpha);
	for (int j=0; j<mergedModules[layer]; j++) {
	  int iModule = module + j;
	  double phi0 = phiMin + iModule*gridPhi;
	  verts[0] = rMin*cos(phi0) + (Lin+deltaL)*cos(phi0+alpha);
	  verts[1] = rMin*sin(phi0) + (Lin+deltaL)*sin(phi0+alpha);
	  verts[3] = rMin*cos(phi0) + (Lout+deltaL)*cos(phi0+alpha);
	  verts[4] = rMin*sin(phi0) + (Lout+deltaL)*sin(phi0+alpha);
	  verts[6] = rMin*cos(phi0+gridPhi) + (Lout-deltaL)*cos(phi0+gridPhi+alpha);
	  verts[7] = rMin*sin(phi0+gridPhi) + (Lout-deltaL)*sin(phi0+gridPhi+alpha);
	  verts[9] = rMin*cos(phi0+gridPhi) + (Lin-deltaL)*cos(phi0+gridPhi+alpha);
	  verts[10] = rMin*sin(phi0+gridPhi) + (Lin-deltaL)*sin(phi0+gridPhi+alpha);
	  verts[2] = verts[5] = verts[8] = verts[11] = 0.;
	  if (qs_rhophi[icl]!=nullptr) {
	    qs_rhophi[icl]->AddQuad(verts);
	    qs_rhophi[icl]->QuadValue( (int) (1000 * energy) );
	    qs_rhophi[icl]->QuadId(new TNamed(Form("Cell %lu", cellID), "Dong!"));
	  }
	}

	// cluster cells in 3D
	float verts3D[24];
	for (int j=0; j<mergedModules[layer]; j++) {
	  int iModule = module + j;
	  double phi0 = phiMin + iModule*gridPhi;
	  
	  verts3D[0] = rMin*cos(phi0) + (Lin+deltaL)*cos(phi0+alpha);
	  verts3D[1] = rMin*sin(phi0) + (Lin+deltaL)*sin(phi0+alpha);
	  verts3D[2] = r_in/tan(theta_center - thetaGrid*mergedCells_Theta[layer]/2.);
	  
	  verts3D[3] = rMin*cos(phi0+gridPhi) + (Lin-deltaL)*cos(phi0+gridPhi+alpha);
	  verts3D[4] = rMin*sin(phi0+gridPhi) + (Lin-deltaL)*sin(phi0+gridPhi+alpha);
	  verts3D[5] = r_in/tan(theta_center - thetaGrid*mergedCells_Theta[layer]/2.);
	  
	  verts3D[6] = rMin*cos(phi0+gridPhi) + (Lin-deltaL)*cos(phi0+gridPhi+alpha);
	  verts3D[7] = rMin*sin(phi0+gridPhi) + (Lin-deltaL)*sin(phi0+gridPhi+alpha);
	  verts3D[8] = r_in/tan(theta_center + thetaGrid*mergedCells_Theta[layer]/2.);
	  
	  verts3D[9] = rMin*cos(phi0) + (Lin+deltaL)*cos(phi0+alpha);
	  verts3D[10] = rMin*sin(phi0) + (Lin+deltaL)*sin(phi0+alpha);
	  verts3D[11] = r_in/tan(theta_center + thetaGrid*mergedCells_Theta[layer]/2.);
	  
	  verts3D[12] = rMin*cos(phi0) + (Lout+deltaL)*cos(phi0+alpha);
	  verts3D[13] = rMin*sin(phi0) + (Lout+deltaL)*sin(phi0+alpha);
	  verts3D[14] = r_out/tan(theta_center - thetaGrid*mergedCells_Theta[layer]/2.);
	  
	  verts3D[15] = rMin*cos(phi0+gridPhi) + (Lout-deltaL)*cos(phi0+gridPhi+alpha);
	  verts3D[16] = rMin*sin(phi0+gridPhi) + (Lout-deltaL)*sin(phi0+gridPhi+alpha);
	  verts3D[17] = r_out/tan(theta_center - thetaGrid*mergedCells_Theta[layer]/2.);
	  
	  verts3D[18] = rMin*cos(phi0+gridPhi) + (Lout-deltaL)*cos(phi0+gridPhi+alpha);
	  verts3D[19] = rMin*sin(phi0+gridPhi) + (Lout-deltaL)*sin(phi0+gridPhi+alpha);
	  verts3D[20] = r_out/tan(theta_center + thetaGrid*mergedCells_Theta[layer]/2.);
	  
	  verts3D[21] = rMin*cos(phi0) + (Lout+deltaL)*cos(phi0+alpha);
	  verts3D[22] = rMin*sin(phi0) + (Lout+deltaL)*sin(phi0+alpha);
	  verts3D[23] = r_out/tan(theta_center + thetaGrid*mergedCells_Theta[layer]/2.);
	  
	  bs[icl]->AddBox(verts3D);
	  bs[icl]->DigitValue( (int) (1000 * energy) );
	  //bs->BoxId(new TNamed(Form("Cell %lu", cellID), "Dong!"));
	}
      }
    }

    if (eventLabel == nullptr) {
      eventLabel = new TGLConstAnnotation(gEve->GetDefaultGLViewer(),
					  Form("%s, %.1f GeV\nEvent %d",
					       partType.Data(), pmax, eventId), 0.1, 0.9);
      eventLabel->SetTextSize(0.05);// % of window diagonal
      eventLabel->SetAllowClose(false);
    }
    else  {
      eventLabel->SetText(Form("%s, %.1f GeV\nEvent %d", partType.Data(), pmax, eventId));
    }
    
    TEveElement* top = (TEveElement*) gEve->GetCurrentEvent();
    // if nothing is drawn, top is null
    if (top) {
      rhoPhiEventScene->DestroyElements();
      rhoPhiProjManager->ImportElements(top, rhoPhiEventScene);

    // slow??
    //TEveRGBAPaletteOverlay *po = new TEveRGBAPaletteOverlay(pal, 0.55, 0.1, 0.4, 0.05);
    //rhoPhiGLView->AddOverlayElement(po);

      rhoZEventScene->DestroyElements();
      rhoZProjManager->ImportElements(top, rhoZEventScene);
      
      //TEveRGBAPaletteOverlay *po2 = new TEveRGBAPaletteOverlay(pal, 0.55, 0.1, 0.4, 0.05);
      //rhoZGLView->AddOverlayElement(po2);
    }
    
    cout << "Done" << endl;

    textEntry->SetTextColor((Pixel_t)0x000000);
    textEntry->SetText(Form("Event %d loaded", eventId));

  }
};



/******************************************************************************/
// THE MAIN FUNCTION
/******************************************************************************/

EventReader* eventReader = nullptr;

void display(int evt = 0) {

  // calculate the geometry parameters
  calcGeom();

  // create the eve manageer
  TEveManager::Create();

  // see palettes here: https://root.cern.ch/doc/master/classTColor.html
  // gStyle->SetPalette(kAvocado);
  gStyle->SetPalette(kSienna);
  
  // first tab
  gEve->GetDefaultGLViewer()->SetGuideState(TGLUtil::kAxesOrigin, false, false, 0);
  gEve->GetDefaultGLViewer()->DrawGuides();
  gEve->GetDefaultViewer()->SetElementName("3D view");
  gEve->GetDefaultGLViewer()->CurrentCamera().RotateRad(-.7, 0.5);

  // Create the geometry and the readout
  TEveElementList* geom = new TEveElementList("Geometry");
  TEveElementList* PCBs = new TEveElementList("PCBs");
  TEveElementList* actives = new TEveElementList("Active elements");
  TEveElementList* passives = new TEveElementList("Passive elements");
  TEveElementList* readout = new TEveElementList("Readout");
  bool useG4geom = true;
  if (useG4geom) {
    //auto fGeom = TFile::Open(geomFile.c_str(), "CACHEREAD");
    auto fGeom = TFile::Open(geomFile.c_str(), "READ");
    if (!fGeom) return;
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) fGeom->Get("ECalBarrel_vol");
    barrel = TEveGeoShape::ImportShapeExtract(gse, 0);
    barrel->SetMainTransparency(70);
    fGeom->Close();
    delete fGeom;
    barrel->SetPickableRecursively(kTRUE);
    geom->AddElement(barrel);
    TPRegexp re;
    // set transparency of the subvolumes of the bath
    re = TPRegexp("LAr_bath*");
    TEveElement* bath = barrel->FindChild(re);
    TEveElement::List_t matches;
    re = TPRegexp("ECAL_Cryo*");
    barrel->FindChildren(matches, re);
    for (auto a : matches) a->SetMainTransparency(70);
    re = TPRegexp("services*");
    barrel->FindChildren(matches, re);
    for (auto a : matches) a->SetMainTransparency(70);
    // make lists of elements inside bath to turn on/off simultaneously
    if (bath) {
      TEveElementList* newbath = new TEveElementList("LAr_bath");
      barrel->AddElement(newbath);
      TEveElement::List_t matches;
      re = TPRegexp("PCB*");
      bath->FindChildren(matches, re);
      for (auto a : matches) PCBs->AddElement(a);
      newbath->AddElement(PCBs);
      TEveElement::List_t matches2;
      re = TPRegexp("active*");
      bath->FindChildren(matches2, re);
      for (auto a : matches2) actives->AddElement(a);
      newbath->AddElement(actives); 
      TEveElement::List_t matches3;     
      re = TPRegexp("passive*");
      bath->FindChildren(matches3, re);
      for (auto a : matches3) passives->AddElement(a);
      newbath->AddElement(passives);
      barrel->RemoveElement(bath);
      // hide elements inside bath by default because they are slow in 3D
      newbath->SetRnrSelfChildren(true, false);
    }
    gEve->AddGlobalElement(geom);
    gEve->AddToListTree(geom, true);
  }
  else {
    // the barrel envelope
    barrel = new TEveGeoShape("Barrel");
    barrel->SetShape(new TGeoTube(rMin, rMax, zMax));
    barrel->SetMainColor(kCyan);
    barrel->SetMainTransparency(90);
    barrel->SetNSegments(128);
    // geom->AddElement(barrel);
    gEve->AddGlobalElement(barrel);
    // the barrel layers
    TEveElementList* layers = new TEveElementList("layers");
    barrel->AddElement(layers);
    TEveGeoShape* b;
    for (int iLayer=0; iLayer<nLayers; iLayer++) {
      b = new TEveGeoShape(Form("Barrel layer %d",iLayer));
      b->SetShape(new TGeoTube(r[iLayer], r[iLayer+1], zMax));
      b->SetMainColor(kCyan);
      b->SetMainTransparency(90);
      b->SetNSegments(128);
      layers->AddElement(b);
    }    
    // the electrodes
    TEveElementList* modules = new TEveElementList("modules");
    barrel->AddElement(modules);
    for (int iModule=0; iModule<nModules; iModule++) {
      double phi0 = iModule*gridPhi-gridPhi/12.; // small extra shift is due to finite width of element (?)
      double phi = phi0 + dPhiAvg;
      b = new TEveGeoShape(Form("Module %d", iModule));
      b->SetShape(new TGeoBBox(Ltot/2, 0.01, (zMax-zMin)/2.0));
      b->SetMainColor(kGray);
      b->SetMainTransparency(95);
      TGeoRotation* rot = new TGeoRotation();
      rot->SetAngles((alpha+iModule*gridPhi)*180./TMath::Pi(), 0., 0.);
      TGeoCombiTrans* c1 = new TGeoCombiTrans(rAvg*cos(phi), rAvg*sin(phi), 0.0, rot);
      b->SetTransMatrix(*c1);
      modules->AddElement(b);
    }
    gEve->AddToListTree(barrel, true);
  }

  gEve->AddToListTree(readout,true);

  // create second tab (R-phi view)
  rhoPhiView = gEve->SpawnNewViewer("Projection Rho-Phi");
  // two scenes, for geometry and event
  rhoPhiScene = gEve->SpawnNewScene("Rho-Phi geometry",
				    "Scene holding projected geometry data for the RhoPhi view.");
  rhoPhiView->AddScene(rhoPhiScene);
  rhoPhiEventScene = gEve->SpawnNewScene("RhoPhi Event Data",
					 "Scene holding projected event-data for the RhoPhi view.");
  rhoPhiView->AddScene(rhoPhiEventScene);
  rhoPhiEventSceneManual = gEve->SpawnNewScene("RhoPhi Event Data 2",
					       "Scene holding hand-crafted event-data for the RhoPhi view.");
  rhoPhiView->AddScene(rhoPhiEventSceneManual);
  rhoPhiGLView = rhoPhiView->GetGLViewer();
  // set camera orientation
  rhoPhiGLView->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  // create 3D->2D projection manager for rho-phi
  rhoPhiProjManager = new TEveProjectionManager();
  rhoPhiProjManager->SetProjection(TEveProjection::kPT_RPhi);
  auto axes = new TEveProjectionAxes(rhoPhiProjManager);
  axes->SetElementName("Rho-Phi projection axes");
  rhoPhiScene->AddElement(axes);
  if (useG4geom)
    rhoPhiProjManager->ImportElements(geom, rhoPhiScene);
  else
    rhoPhiProjManager->ImportElements(barrel, rhoPhiScene);

  // the merged module grid
  TEveStraightLineSet* gridmod = new TEveStraightLineSet("phi readout merged");
  gridmod->SetLineColor(kViolet+2);
  gridmod->SetLineWidth(8);
  for (int iLayer=0; iLayer<nLayers; iLayer++) {
    double Lin = _L(alpha, rMin, r[iLayer]);
    double Lout = _L(alpha, rMin, r[iLayer+1]);
    for (int iModule=0; iModule<nModules; iModule++) {
      if (iModule % mergedModules[iLayer] != 0) continue;
      double phi0 = phiMin + iModule*gridPhi;
      double x1 = rMin*cos(phi0) + Lin*cos(phi0+alpha);
      double y1 = rMin*sin(phi0) + Lin*sin(phi0+alpha);
      double x2 = rMin*cos(phi0) + Lout*cos(phi0+alpha);
      double y2 = rMin*sin(phi0) + Lout*sin(phi0+alpha);
      gridmod->AddLine(x1, y1, 0., x2, y2, 0.);
    }
  }
  // add to scene?
  rhoPhiScene->AddElement(gridmod);
  readout->AddElement(gridmod);
  
  // third tab (R-z view)
  rhoZView = gEve->SpawnNewViewer("Projection Rho-Z");
  rhoZScene = gEve->SpawnNewScene("Rho-Z geometry",
				  "Scene holding projected geometry data for the RhoZ view.");
  rhoZView->AddScene(rhoZScene);
  rhoZEventScene = gEve->SpawnNewScene("RhoZ Event Data",
				       "Scene holding projected event-data for the RhoZ view.");
  rhoZView->AddScene(rhoZEventScene);
  rhoZEventSceneManual = gEve->SpawnNewScene("RhoZ Event Data 2",
					     "Scene holding hand-crafted event-data for the RhoZ view.");
  rhoZView->AddScene(rhoZEventSceneManual);
  rhoZGLView = rhoZView->GetGLViewer();
  rhoZGLView->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);

  rhoZProjManager = new TEveProjectionManager();
  rhoZProjManager->SetProjection(TEveProjection::kPT_RhoZ);
  auto axes2 = new TEveProjectionAxes(rhoZProjManager);
  axes2->SetElementName("Rho-Z projection axes");
  rhoZScene->AddElement(axes2);
  if (useG4geom)
    rhoZProjManager->ImportElements(geom, rhoZScene);
  else
    rhoZProjManager->ImportElements(barrel, rhoZScene);

  // the theta readout grid
  TEveStraightLineSet* grid = new TEveStraightLineSet("theta readout");
  grid->SetLineColor(kViolet);
  for (int iTheta=0; iTheta<=nThetaBins; iTheta++) {
    double theta = thetaMin + iTheta*thetaGrid;
    double r1 = rMin;
    double r2 = rMax;
    double z1 = r1*cos(theta)/sin(theta);
    double z2 = r2*cos(theta)/sin(theta);
    if (z1<zMax && z1>zMin) {
      if (z2>zMax) {
	z2 = zMax;
	r2 = z2*sin(theta)/cos(theta);
      }
      if (z2<zMin) {
	z2 = zMin;
	r2 = z2*sin(theta)/cos(theta);
      }
      grid->AddLine(z1,  r1, 0., z2,  r2, 0.);
      grid->AddLine(z1, -r1, 0., z2, -r2, 0.);
    }
  }
  rhoZScene->AddElement(grid);
  readout->AddElement(grid);
  
  // the merged grid
  TEveStraightLineSet* grid2 = new TEveStraightLineSet("theta readout merged");
  grid2->SetLineColor(kViolet+2);
  grid2->SetLineWidth(8);
  for (int iLayer=0; iLayer<nLayers; iLayer++) {
    double r1 = r[iLayer];
    double r2 = r[iLayer+1];
    for (int iTheta=0; iTheta<=nThetaBins; iTheta++) {
      if (iTheta % mergedCells_Theta[iLayer] != 0) continue;
      double theta = thetaMin + iTheta*thetaGrid;
      double z1 = r1*cos(theta)/sin(theta);
      double z2 = r2*cos(theta)/sin(theta);
      double r2tmp = r2;
      if (z1<zMax && z1>zMin) {
	if (z2>zMax) {
	  z2 = zMax;
	  r2tmp = z2*sin(theta)/cos(theta);
	}
	if (z2<zMin) {
	  z2 = zMin;
	  r2tmp = z2*sin(theta)/cos(theta);
	}
	grid2->AddLine(z1, r1, 0.,
		       z2, r2tmp, 0.);
	grid2->AddLine(z1, -r1, 0.,
		       z2, -r2tmp, 0.);
      }
    }
  }
  rhoZScene->AddElement(grid2);
  readout->AddElement(grid2);


  gEve->Redraw3D(true);

  
  makeGui();

  //
  // data
  //

  // setup the reader
  TFile* f = TFile::Open(evtFile.c_str(), "READ");
  eventReader = new EventReader(f);

  // load and display the requested event
  eventReader->loadEvent(evt);

  // Set the 3D view as the active tab and rotate the camera
  gEve->GetBrowser()->GetTabRight()->SetTab(0);

  // Draw
  gEve->Redraw3D(true);
}


/******************************************************************************/
// GUI
/******************************************************************************/

//______________________________________________________________________________
//
// EvNavHandler class is needed to connect GUI signals.

class EvNavHandler
{
public:
  void Fwd()
  {
    if (eventId < nEvents - 1) {
      ++eventId;
      eventReader->loadEvent();
    } else {
      textEntry->SetTextColor(0xff0000);
      textEntry->SetText("Already at last event");
      printf("Already at last event.\n");
    }
  }
  void Bck()
  {
    if (eventId > 0) {
      --eventId;
      eventReader->loadEvent();
    } else {
      textEntry->SetTextColor(0xff0000);
      textEntry->SetText("Already at first event");
      printf("Already at first event.\n");
    }
  }
};

void makeGui()
{
  // Create minimal GUI for event navigation.

  TEveBrowser* browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft);

  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("GUI");
  frmMain->SetCleanup(kDeepCleanup);
  
  TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
  {
    
    TString icondir( Form("%s/icons/", gSystem->Getenv("ROOTSYS")) );
    TGPictureButton* b = 0;
    EvNavHandler    *fh = new EvNavHandler;
    
    b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 2, 10, 10));
    b->Connect("Clicked()", "EvNavHandler", fh, "Bck()");
    
    b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 10, 10, 10));
    b->Connect("Clicked()", "EvNavHandler", fh, "Fwd()");
    
    textEntry = new TGTextEntry(hf);
    textEntry->SetEnabled(kFALSE);
    hf->AddFrame(textEntry, new TGLayoutHints(kLHintsLeft | kLHintsCenterY  |
					      kLHintsExpandX, 2, 10, 10, 10));
    
  }
  frmMain->AddFrame(hf);
  
  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();
  
  browser->StopEmbedding();
  browser->SetTabTitle("Event Control", 0);
}
