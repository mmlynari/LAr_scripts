#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include <iostream>
using namespace std;

TTree *T = nullptr;
//const std::string filename = "neighbours_map_barrel_thetamodulemerged.root";
const std::string filename = "neighbours_map_HCalBarrel.root";
const std::string treename = "neighbours";

TTree *Tnoise = nullptr;
const std::string filenameNoise = "cellNoise_map_electronicsNoiseLevel_ecalB_thetamodulemerged_hcalB_thetaphi.root";
const std::string treenameNoise = "noisyCells";

ULong64_t cID;
std::vector<unsigned long> *neighbours = 0;
ULong64_t cIDNoise;
double noiseLevel;
double noiseOffset;

bool useHCalReadoutWithRows = false;

// HELPER FUNCTIONS
ULong_t ReadNbitsAtPositionMFromCellID(int n, int m, ULong_t cellID)
{
  const ULong_t mask = (1<<n) - 1;
  return (cellID >> m) & mask;
}

// ECalBarrel: system:4,cryo:1,type:3,subtype:3,layer:8,module:11,theta:10
// HCalBarrel with row: system:4,layer:5,row:9,theta:9,phi:10
// HCalBarrel without row: system:4,layer:5,theta:9,phi:10

// extract system ID from cellID
ULong_t SystemID(ULong_t cellID)
{
  return ReadNbitsAtPositionMFromCellID(4, 0, cellID);
}
// extract layer number from cellID
ULong_t Layer(ULong_t cellID)
{
  if (SystemID(cellID)==4)
    return ReadNbitsAtPositionMFromCellID(8, 11, cellID);
  if (SystemID(cellID)==8)
    return ReadNbitsAtPositionMFromCellID(5, 4, cellID);
  return 999999999;
}

// extract module number from cellID
ULong_t ECalBarrelModule(ULong_t cellID)
{
  if (SystemID(cellID)==4)
    return ReadNbitsAtPositionMFromCellID(11, 19, cellID);
  return 999999999;
}

// extract theta bin from cellID
ULong_t ThetaBin(ULong_t cellID)
{
  if (SystemID(cellID)==4)
    return ReadNbitsAtPositionMFromCellID(10, 30, cellID);
  if (SystemID(cellID)==8)
  {
    if (useHCalReadoutWithRows)
      return ReadNbitsAtPositionMFromCellID(9, 18, cellID);
    else
      return ReadNbitsAtPositionMFromCellID(9, 9, cellID);
  }
  return 999999999;
}

// extract row number from cellID
ULong_t HCalBarrelRow(ULong_t cellID)
{
  if (SystemID(cellID)==8)
  {
    if (useHCalReadoutWithRows)
      return ReadNbitsAtPositionMFromCellID(9, 9, cellID);
  }
  return 999999999;
}

// extract phi number from cellID
ULong_t HCalBarrelPhiBin(ULong_t cellID)
{
  if (SystemID(cellID)==8) {
    if (useHCalReadoutWithRows)
      return ReadNbitsAtPositionMFromCellID(10, 27, cellID);
    else
      return ReadNbitsAtPositionMFromCellID(10, 18, cellID);
  }
  return 999999999;   
}

void LoadNeighboursMap()
{
  if (T == nullptr && filename!="")
  {
    TFile *f = TFile::Open(filename.c_str(), "READ");
    T = (TTree *)f->Get(treename.c_str());
    T->SetBranchAddress("cellId", &cID);
    T->SetBranchAddress("neighbours", &neighbours);
  }
}

void LoadNoiseMap()
{
  if (Tnoise == nullptr && filenameNoise!="")
  {
    TFile *fNoise = TFile::Open(filenameNoise.c_str(), "READ");
    Tnoise = (TTree *)fNoise->Get(treenameNoise.c_str());
    Tnoise->SetBranchAddress("cellId", &cIDNoise);
    Tnoise->SetBranchAddress("noiseLevel", &noiseLevel);
    Tnoise->SetBranchAddress("noiseOffset", &noiseOffset);
  }
}

void printCell(ULong_t cellID)
{
  cout << "cellID: " << cellID << endl;
  cout << "System: " << SystemID(cellID) << endl;
  if (SystemID(cellID)==4) {  
    cout << "Layer: " << Layer(cellID) << endl;
    cout << "Theta bin: " << ThetaBin(cellID) << endl;
    cout << "Module: " << ECalBarrelModule(cellID) << endl;
  }
  if (SystemID(cellID)==8) {  
    cout << "Layer: " << Layer(cellID) << endl;
    cout << "Row: " << HCalBarrelRow(cellID) << endl;
    cout << "Theta bin: " << ThetaBin(cellID) << endl;
    cout << "Phi bin: " << HCalBarrelPhiBin(cellID) << endl;
  }
  cout << endl;
}

void printCell(ULong64_t iEntry, bool showNeighbours, bool showNoise)
{
  T->GetEntry(iEntry);
  cout << "=================================================" << endl;
  cout << endl;
  printCell(cID);
  if (showNeighbours)
  {
    cout << "Neighbours: " << endl
         << endl;
    for (unsigned int i = 0; i < neighbours->size(); i++)
    {
      printCell(neighbours->at(i));
    }
  }
  if (showNoise)
  {
    for (ULong64_t jEntry = 0; jEntry < Tnoise->GetEntries(); jEntry++)
    {
      Tnoise->GetEntry(jEntry);
      if (cIDNoise == cID)
      {
        cout << "Noise: level = " << noiseLevel << " , offset = " << noiseOffset
            << endl << endl;
            break;
      }
    }
  }
  cout << "=================================================" << endl;
}

void printNeighboursOfCell(ULong_t cellID)
{
  LoadNeighboursMap();
  for (ULong64_t iEntry = 0; iEntry < T->GetEntries(); iEntry++)
  {
    T->GetEntry(iEntry);
    if (cID == cellID)
    {
      printCell(iEntry, true, false);
      return;
    }
  }
  cout << "CellID not found" << endl;
}

void printNeighbours(int n = 10)
{
  LoadNeighboursMap();
  for (int i = 0; i < n; i++)
  {
    int entry = (int)gRandom->Uniform(T->GetEntries());
    printCell(entry, true, false);
  }
}
