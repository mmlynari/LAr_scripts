TTree* T = nullptr;
const std::string filename = "neighbours_map_barrel_thetamodulemerged.root";
const std::string treename = "neighbours";
ULong64_t cID;
std::vector<unsigned long> *neighbours=0;

// HELPER FUNCTIONS

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

void printCell(ULong_t cellID) {
  cout << "cellID: " << cellID << endl;
  cout << "Layer: " << Layer(cellID) << endl;
  cout << "Theta bin: " << ThetaBin(cellID) << endl;
  cout << "Module: " << Module(cellID) << endl;
  cout << endl;
}

void LoadNeighboursMap() {
  if (T==nullptr) {
    TFile* f = TFile::Open(filename.c_str(),"READ");
    T = (TTree*) f->Get(treename.c_str());
    T->SetBranchAddress("cellId", &cID);
    T->SetBranchAddress("neighbours", &neighbours);
  }
}

void printCellAndNeighbours(ULong64_t iEntry) {
  T->GetEntry(iEntry);
  cout << "=================================================" << endl;
  cout << endl;
  printCell(cID);
  cout << "Neighbours: " << endl << endl;
  for (unsigned int i=0; i<neighbours->size(); i++) {
    printCell(neighbours->at(i));
  }
  cout << "=================================================" << endl;
}

void printNeighboursOfCell(ULong_t cellID) {
  LoadNeighboursMap();
  for (ULong64_t iEntry=0; iEntry<T->GetEntries(); iEntry++) {
    T->GetEntry(iEntry);
    if (cID == cellID) {
      printCellAndNeighbours(iEntry);
      return;
    }
  }
  cout << "CellID not found" << endl;
}


void printNeighbours(int n=10) {
  LoadNeighboursMap();
  for (int i=0; i<n; i++) {
    int entry = (int) gRandom->Uniform(T->GetEntries());
    printCellAndNeighbours(entry);
  }
}
