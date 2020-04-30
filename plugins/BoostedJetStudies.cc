/**
   
*/

// system include files
#include <memory>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "L1Trigger/L1TCaloLayer1/src/UCTParameters.hh"

#include "L1Trigger/L1TCaloLayer1/src/UCTLayer1.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTCrate.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTCard.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTTower.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"

#include "L1Trigger/L1TCaloSummary/src/UCTObject.hh"
#include "L1Trigger/L1TCaloSummary/src/UCTSummaryCard.hh"
#include "L1Trigger/L1TCaloSummary/src/UCTGeometryExtended.hh"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "L1Trigger/Run3Ntuplizer/plugins/helpers.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// GCT and RCT data formats
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
//#include "L1Trigger/L1TCaloLayer1/src/L1UCTCollections.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Muon.h"


#include "L1Trigger/L1TCaloLayer1/src/L1TCaloLayer1FetchLUTs.hh"

#pragma extra_include "TLorentzVector.h";
#pragma link C++ class std::vector<TLorentzVector>;

using namespace l1tcalo;
using namespace l1extra;
using namespace std;

bool compareByPt (l1extra::L1JetParticle i, l1extra::L1JetParticle j) { return(i.pt()>j.pt()); };

float towerEtaMap[28]= {
    0.0435,
    0.1305, 0.2175, 0.3045, 0.3915, 0.4785,
    0.5655, 0.6525, 0.7395, 0.8265, 0.9135,
    1.0005, 1.0875, 1.1745, 1.2615, 1.3485,
    1.4355, 1.5225, 1.6095, 1.6965, 1.7835,
    1.8705, 1.9575, 2.0445, 2.217, 2.391,
    2.565, //2.739,
    2.913,
    };

float towerPhiMap[72]=
    {-0.131, -0.044, 0.044, 0.131, 0.218, 0.305, 0.393, 0.480, 0.567, 0.654, 0.742, 0.829, 0.916, 1.004, 1.091, 1.178, 1.265, 1.353, 1.440, 1.527, 1.614, 1.702, 1.789, 1.876, 1.963, 2.051, 2.138, 2.225, 2.313, 2.400, 2.487, 2.574, 2.662, 2.749, 2.836, 2.923, 3.011, 3.098,
      -3.098, -3.011, -2.923, -2.836, -2.749, -2.662, -2.574, -2.487, -2.400, -2.313, -2.225, -2.138, -2.051, -1.963, -1.876, -1.789, -1.702, -1.614, -1.527, -1.440, -1.353, -1.265, -1.178, -1.091, -1.004, -0.916, -0.829, -0.742, -0.654, -0.567, -0.480, -0.393, -0.305, -0.218};

float getRecoEta(int ieta, short zside){
  float eta = -999;
  if(ieta<0 || ieta>(28*2)){
    std::cout<<"Error!!! towereta out of bounds in triggerGeometryTools.h "<<std::endl;
    std::cout<<"ieta "<<ieta<<std::endl;
    exit(0);
  }
  if(zside == 1)
    eta = towerEtaMap[ieta];
  else if(zside == -1)
    eta = towerEtaMap[ieta];
  else{
    std::cout<<"Error!!! zside out of bounds in triggerGeometryTools.h "<<std::endl;
    std::cout<<"zside "<<zside<<std::endl;
    exit(0);
  }
  return eta;
};

float getRecoEtaNew(int caloEta){
  float eta = -999.;
  static bool first = true;
  static double twrEtaValues[42];
  if(first) {
    twrEtaValues[0] = 0;
    for(unsigned int i = 0; i < 20; i++) {
      twrEtaValues[i + 1] = 0.0436 + i * 0.0872;
    }
    twrEtaValues[21] = 1.785;
    twrEtaValues[22] = 1.880;
    twrEtaValues[23] = 1.9865;
    twrEtaValues[24] = 2.1075;
    twrEtaValues[25] = 2.247;
    twrEtaValues[26] = 2.411;
    twrEtaValues[27] = 2.575;
    twrEtaValues[28] = 2.825;
    twrEtaValues[29] = 999.;
    twrEtaValues[30] = (3.15+2.98)/2.;
    twrEtaValues[31] = (3.33+3.15)/2.;
    twrEtaValues[32] = (3.50+3.33)/2.;
    twrEtaValues[33] = (3.68+3.50)/2.;
    twrEtaValues[34] = (3.68+3.85)/2.;
    twrEtaValues[35] = (3.85+4.03)/2.;
    twrEtaValues[36] = (4.03+4.20)/2.;
    twrEtaValues[37] = (4.20+4.38)/2.;
    twrEtaValues[38] = (4.74+4.38*3)/4.;
    twrEtaValues[39] = (4.38+4.74*3)/4.;
    twrEtaValues[40] = (5.21+4.74*3)/4.;
    twrEtaValues[41] = (4.74+5.21*3)/4.;
    first = false;
  }
  uint32_t absCaloEta = abs(caloEta);
  if(absCaloEta <= 41) {
    if(caloEta < 0)
      eta =  -twrEtaValues[absCaloEta];
    else
      eta = +twrEtaValues[absCaloEta];
  }
  return eta;
};

float getRecoPhi(int iphi){
  return towerPhiMap[iphi-1];
};

float getRecoPhiNew(int caloPhi){
  float phi = -999.;
  if(caloPhi > 72) phi = +999.;
  uint32_t absCaloPhi = std::abs(caloPhi) - 1;
  if(absCaloPhi < 36)
    phi = (((double) absCaloPhi + 0.5) * 0.0872);
  else
    phi = (-(71.5 - (double) absCaloPhi) * 0.0872);
  return phi;
};

int TPGEtaRange(int ieta){
  int iEta = 0;
  if(ieta < 0)
    iEta = ieta + 28;
  else if(ieta > 0)
    iEta = ieta + 27;
  return iEta;
};

//
// class declaration
//

class BoostedJetStudies : public edm::EDAnalyzer {
public:
  explicit BoostedJetStudies(const edm::ParameterSet&);
  ~BoostedJetStudies();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void zeroOutAllVariables();

private:
  void analyze(const edm::Event& evt, const edm::EventSetup& es);      
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  void print();

  // ----------member data ---------------------------
  edm::InputTag genSrc_;
  edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalTPSource;
  std::string ecalTPSourceLabel;
  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalTPSource;
  std::string hcalTPSourceLabel;

  edm::EDGetTokenT<vector<pat::Jet> > jetSrc_;
  edm::EDGetTokenT<vector<pat::Jet> > jetSrcAK8_;

  std::vector< std::vector< std::vector < uint32_t > > > ecalLUT;
  std::vector< std::vector< std::vector < uint32_t > > > hcalLUT;
  std::vector< std::vector< uint32_t > > hfLUT;

  uint32_t nPumBins;

  std::vector< std::vector< std::vector < uint32_t > > > pumLUT;

  std::vector< UCTTower* > twrList;

  bool useLSB;
  bool useCalib;
  bool useECALLUT;
  bool useHCALLUT;
  bool useHFLUT;

  double caloScaleFactor;

  uint32_t jetSeed;
  uint32_t tauSeed;
  float tauIsolationFactor;
  uint32_t eGammaSeed;
  double eGammaIsolationFactor;

  bool verbose;

  UCTParameters uctParameters;
  UCTLayer1 *layer1;
  UCTSummaryCard *summaryCard;

  TH1F* nEvents;
  TH1F* recoJet_pt;
  TH1F* recoJet_eta;
  TH1F* recoJet_phi;

  TH1F* recoJetAK8_pt;
  TH1F* recoJetAK8_eta;
  TH1F* recoJetAK8_phi;

  TTree* l1Tree;
  int run, lumi, event;

  double genPt, genEta, genPhi;
  double recoPt, recoEta, recoPhi;
  double l1Pt, l1Eta, l1Phi;

  double genPt_1, genEta_1, genPhi_1;
  double recoPt_1, recoEta_1, recoPhi_1;
  double l1Pt_1, l1Eta_1, l1Phi_1;
  
  double genPt_2, genEta_2, genPhi_2;
  double recoPt_2, recoEta_2, recoPhi_2;
  double l1Pt_2, l1Eta_2, l1Phi_2;

  double genDeltaEta, genDeltaPhi, genDeltaR, genMass;
  double recoDeltaEta, recoDeltaPhi, recoDeltaR, recoMass;
  double l1DeltaEta, l1DeltaPhi, l1DeltaR, l1Mass;

  int l1NthJet_1, l1NthJet_2;
  int l1NTau_1, l1NTau_2;
  int recoNthJet_1, recoNthJet_2;

  double vbfBDT;
  double recoPt_;
  std::vector<int> nSubJets, nBHadrons, HFlav, nTausInfo;
  std::vector<std::vector<int>> subJetHFlav;
  std::vector<float> tau1, tau2, tau3;

  std::vector<TLorentzVector> *allRegions  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *allEcalTPGs  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *allHcalTPGs  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *caloClusters  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *l1Jets  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *ak8Jets  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *subJets  = new std::vector<TLorentzVector>;

  int nGenJets, nRecoJets, nL1Jets;
  int l1Matched_1, l1Matched_2;
  void createBranches(TTree *tree);
  TTree* efficiencyTree;
  edm::Service<TFileService> tfs_;  

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BoostedJetStudies::BoostedJetStudies(const edm::ParameterSet& iConfig) :
  ecalTPSource(consumes<EcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("ecalToken"))),
  ecalTPSourceLabel(iConfig.getParameter<edm::InputTag>("ecalToken").label()),
  hcalTPSource(consumes<HcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("hcalToken"))),
  hcalTPSourceLabel(iConfig.getParameter<edm::InputTag>("hcalToken").label()),
  ecalLUT(28, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(256))),
  hcalLUT(28, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(256))),
  hfLUT(12, std::vector < uint32_t >(256)),
  nPumBins(iConfig.getParameter<unsigned int>("nPumBins")),
  pumLUT(nPumBins, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(13))),
  useLSB(iConfig.getParameter<bool>("useLSB")),
  useCalib(iConfig.getParameter<bool>("useCalib")),
  useECALLUT(iConfig.getParameter<bool>("useECALLUT")),
  useHCALLUT(iConfig.getParameter<bool>("useHCALLUT")),
  useHFLUT(iConfig.getParameter<bool>("useHFLUT")),
  caloScaleFactor(iConfig.getParameter<double>("caloScaleFactor")),
  jetSeed(iConfig.getParameter<unsigned int>("jetSeed")),
  tauSeed(iConfig.getParameter<unsigned int>("tauSeed")),
  tauIsolationFactor(iConfig.getParameter<double>("tauIsolationFactor")),
  eGammaSeed(iConfig.getParameter<unsigned int>("eGammaSeed")),
  eGammaIsolationFactor(iConfig.getParameter<double>("eGammaIsolationFactor")),
  verbose(iConfig.getParameter<bool>("verbose")),
  uctParameters(iConfig.getParameter<double>("activityFraction"), 
		iConfig.getParameter<double>("ecalActivityFraction"), 
		iConfig.getParameter<double>("miscActivityFraction")),
  jetSrc_(    consumes<vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("recoJets"))),
  jetSrcAK8_( consumes<vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("recoJetsAK8"))),
  genSrc_((        iConfig.getParameter<edm::InputTag>( "genParticles")))
{
  std::vector<double> pumLUTData;
  char pumLUTString[10];
  for(uint32_t pumBin = 0; pumBin < nPumBins; pumBin++) {
    for(uint32_t side = 0; side < 2; side++) {
      if(side == 0) sprintf(pumLUTString, "pumLUT%2.2dp", pumBin);
      else sprintf(pumLUTString, "pumLUT%2.2dn", pumBin);
      pumLUTData = iConfig.getParameter<std::vector < double > >(pumLUTString);
      for(uint32_t iEta = 0; iEta < std::max((uint32_t) pumLUTData.size(), MaxUCTRegionsEta); iEta++) {
	pumLUT[pumBin][side][iEta] = (uint32_t) round(pumLUTData[iEta] / caloScaleFactor);
      }
      if(pumLUTData.size() != (MaxUCTRegionsEta))
	std::cerr << "PUM LUT Data size integrity check failed; Expected size = " << MaxUCTRegionsEta
		  << "; Provided size = " << pumLUTData.size()
		  << "; Will use what is provided :(" << std::endl;
    }
  }
  /*
  produces< L1CaloRegionCollection >();
  produces< L1EmParticleCollection >( "Isolated" ) ;
  produces< L1EmParticleCollection >( "NonIsolated" ) ;
  produces< L1JetParticleCollection >( "Central" ) ;
  produces< L1JetParticleCollection >( "Forward" ) ;
  produces< L1JetParticleCollection >( "Tau" ) ;
  produces< L1JetParticleCollection >( "IsoTau" ) ;
  produces< L1EtMissParticleCollection >( "MET" ) ;
  produces< L1EtMissParticleCollection >( "MHT" ) ;
  */
  layer1 = new UCTLayer1(&uctParameters);
  summaryCard = new UCTSummaryCard(layer1, &pumLUT, jetSeed, tauSeed, tauIsolationFactor, eGammaSeed, eGammaIsolationFactor);
  vector<UCTCrate*> crates = layer1->getCrates();
  for(uint32_t crt = 0; crt < crates.size(); crt++) {
    vector<UCTCard*> cards = crates[crt]->getCards();
    for(uint32_t crd = 0; crd < cards.size(); crd++) {
      vector<UCTRegion*> regions = cards[crd]->getRegions();
      for(uint32_t rgn = 0; rgn < regions.size(); rgn++) {
	vector<UCTTower*> towers = regions[rgn]->getTowers();
	for(uint32_t twr = 0; twr < towers.size(); twr++) {
	  twrList.push_back(towers[twr]);
	}
      }
    }
  }
  // Initialize the Trees

  recoPt_      = iConfig.getParameter<double>("recoPtCut");
  nEvents      = tfs_->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );
  recoJet_pt   = tfs_->make<TH1F>( "recoJet_pt" , "p_{t}", 300,  0., 300. );
  recoJet_eta  = tfs_->make<TH1F>( "recoJet_eta"  , "eta", 100,  -3, 3. );
  recoJet_phi  = tfs_->make<TH1F>( "recoJet_phi"  , "phi", 100,  -4, 4. );
  
  recoJetAK8_pt   = tfs_->make<TH1F>( "recoJetAK8_pt" , "p_{t}", 300,  0., 300. );
  recoJetAK8_eta  = tfs_->make<TH1F>( "recoJetAK8_eta"  , "eta", 100,  -3, 3. );
  recoJetAK8_phi  = tfs_->make<TH1F>( "recoJetAK8_phi"  , "phi", 100,  -4, 4. );

  efficiencyTree = tfs_->make<TTree>("efficiencyTree", "Gen Matched Jet Tree ");
  createBranches(efficiencyTree);
  
}

BoostedJetStudies::~BoostedJetStudies() {
  if(layer1 != 0) delete layer1;
  if(summaryCard != 0) delete summaryCard;
}

//
// member functions
//

// ------------ method called to produce the data  ------------

void BoostedJetStudies::analyze( const edm::Event& evt, const edm::EventSetup& es )
{
  using namespace edm;

  nEvents->Fill(1);
  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();
  Handle<L1CaloRegionCollection> regions;
   
  std::vector<pat::Jet> goodJets;
  std::vector<pat::Jet> goodJetsAK8;

  allRegions->clear();
  allEcalTPGs->clear();
  allHcalTPGs->clear();
  caloClusters->clear();
  l1Jets->clear();
  ak8Jets->clear();
  subJets->clear();
  nTausInfo.clear();

  // Start Running Layer 1
  edm::Handle<EcalTrigPrimDigiCollection> ecalTPs;
  evt.getByToken(ecalTPSource, ecalTPs);
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPs;
  evt.getByToken(hcalTPSource, hcalTPs);

  std::unique_ptr<L1CaloRegionCollection> rgnCollection (new L1CaloRegionCollection);
  std::unique_ptr<L1EmParticleCollection> iEGCands(new L1EmParticleCollection);
  std::unique_ptr<L1EmParticleCollection> nEGCands(new L1EmParticleCollection);
  std::unique_ptr<L1JetParticleCollection> iTauCands(new L1JetParticleCollection);
  std::unique_ptr<L1JetParticleCollection> nTauCands(new L1JetParticleCollection);
  std::unique_ptr<L1JetParticleCollection> cJetCands(new L1JetParticleCollection);
  std::unique_ptr<L1JetParticleCollection> fJetCands(new L1JetParticleCollection);
  std::unique_ptr<L1JetParticleCollection> bJetCands(new L1JetParticleCollection);
  std::unique_ptr<L1EtMissParticleCollection> metCands(new L1EtMissParticleCollection);
  std::unique_ptr<L1EtMissParticleCollection> mhtCands(new L1EtMissParticleCollection);

  uint32_t expectedTotalET = 0;

  if(!layer1->clearEvent()) {
    std::cerr << "UCT: Failed to clear event" << std::endl;
    exit(1);
  }

  for ( const auto& ecalTp : *ecalTPs ) {
    int caloEta = ecalTp.id().ieta();
    int caloPhi = ecalTp.id().iphi();
    int et = ecalTp.compressedEt();
    bool fgVeto = ecalTp.fineGrain();
    if(et != 0) {
      UCTTowerIndex t = UCTTowerIndex(caloEta, caloPhi);
      if(!layer1->setECALData(t,fgVeto,et)) {
	std::cerr << "UCT: Failed loading an ECAL tower" << std::endl;
	return;
      }
      expectedTotalET += et;
      int ieta = TPGEtaRange(caloEta);
      int zside = ecalTp.id().zside();
      //float eta = getRecoEta(ieta, zside);
      float eta = getRecoEtaNew(caloEta);
      std::cout<<"ECAL   "<<"eta: "<<eta<<std::endl;
      //float phi = getRecoPhi(caloPhi);
      float phi = getRecoPhiNew(caloPhi);
      TLorentzVector temp;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      allEcalTPGs->push_back(temp);
    }
  }

  for ( const auto& hcalTp : *hcalTPs ) {
    int caloEta = hcalTp.id().ieta();
    uint32_t absCaloEta = abs(caloEta);
    // Tower 29 is not used by Layer-1
    if(absCaloEta == 29) {
      continue;
    }
    // Prevent usage of HF TPs with Layer-1 emulator if HCAL TPs are old style
    else if(hcalTp.id().version() == 0 && absCaloEta > 29) {
      continue;
    }
    else if(absCaloEta <= 41) {
      int caloPhi = hcalTp.id().iphi();
      if(caloPhi <= 72) {
	int et = hcalTp.SOI_compressedEt();
	bool fg = hcalTp.SOI_fineGrain();
	if(et != 0) {
	  UCTTowerIndex t = UCTTowerIndex(caloEta, caloPhi);
	  uint32_t featureBits = 0;
	  if(fg) featureBits = 0x1F; // Set all five feature bits for the moment - they are not defined in HW / FW yet!
	  if(!layer1->setHCALData(t, featureBits, et)) {
	    std::cerr << "caloEta = " << caloEta << "; caloPhi =" << caloPhi << std::endl;
	    std::cerr << "UCT: Failed loading an HCAL tower" << std::endl;
	    return; 
	  }
	  expectedTotalET += et;
          int ieta = TPGEtaRange(caloEta);
          int zside = hcalTp.id().zside();
          //float eta = getRecoEta(ieta, zside);
          float eta = getRecoEtaNew(caloEta);
          //float phi = getRecoPhi(caloPhi);
          float phi = getRecoPhiNew(caloPhi);
          TLorentzVector temp;
          temp.SetPtEtaPhiE(et,eta,phi,et);
          allHcalTPGs->push_back(temp);
	}
      }
      else {
	std::cerr << "Illegal Tower: caloEta = " << caloEta << "; caloPhi =" << caloPhi << std::endl;	
      }
    }
    else {
      std::cerr << "Illegal Tower: caloEta = " << caloEta << std::endl;
    }
  }

  if(!layer1->process()) {
    std::cerr << "UCT: Failed to process layer 1" << std::endl;
    exit(1);
  }

  // Crude check if total ET is approximately OK!
  // We can't expect exact match as there is region level saturation to 10-bits
  // 1% is good enough
  int diff = abs((int)layer1->et() - (int)expectedTotalET);
  if(verbose && diff > 0.01 * expectedTotalET ) {
    print();
    std::cout << "Expected " 
	      << std::showbase << std::internal << std::setfill('0') << std::setw(10) << std::hex
	      << expectedTotalET << std::dec << std::endl;
  }
 
  UCTGeometry g;

  vector<UCTCrate*> crates = layer1->getCrates();
  for(uint32_t crt = 0; crt < crates.size(); crt++) {
    vector<UCTCard*> cards = crates[crt]->getCards();
    for(uint32_t crd = 0; crd < cards.size(); crd++) {
      vector<UCTRegion*> regions = cards[crd]->getRegions();
      for(uint32_t rgn = 0; rgn < regions.size(); rgn++) {
	uint32_t rawData = regions[rgn]->rawData();
	uint32_t regionData = rawData & 0x0000FFFF;
	// uint32_t regionLoc = rawData >> LocationShift;
	// uint32_t regionET = rawData & 0x3FF;
	uint32_t crate = regions[rgn]->getCrate();
	uint32_t card = regions[rgn]->getCard();
	uint32_t region = regions[rgn]->getRegion();
	bool negativeEta = regions[rgn]->isNegativeEta();
	uint32_t rPhi = g.getUCTRegionPhiIndex(crate, card);
	// We want to reuse L1CaloRegion and L1CaloRegionDetID
	// We do not want to change those classes too much
	// We want comparison to legacy for Barrel and Endcap to work transparently
	// Noting that rEta is packed in 5 bits of L1CaloRegionDetID, we have a scheme!
	// We store the Barrel and Endcap regions in the same location as done for RCT
	// HF has changed in the upgrade, 6x2 HF regions instead of 4x2 in case of RCT
	// Note that for upgrade region numbers range 0-6 for Barrel/Endcap and 7-12 for HF
	// So, the scheme used for rEta for upgrade is:
	// rEta= 0- 3 for -HF regions 7-10
	// rEta= 4-10 for -B/E regions 0-6
	// rEta=11-17 for +B/E regions 0-6
	// rEta=18-23 for +HF regions 7-12
	// rEta=30 for -HF region 11
	// rEta=31 for -HF region 12
	uint32_t rEta = 10 - region;
	if(negativeEta && region == 11) rEta = 30;
	if(negativeEta && region == 12) rEta = 31;
	if(!negativeEta) rEta = 11 + region; // Positive eta portion is offset by 11
	rgnCollection->push_back(L1CaloRegion((uint16_t) regionData, (unsigned) rEta, (unsigned) rPhi, (int16_t) 0));   
      }
    }
  }  

  for(vector<L1CaloRegion>::const_iterator testRegion = rgnCollection->begin(); testRegion != rgnCollection->end(); ++testRegion){
    uint16_t test_raw = testRegion->raw();
    uint32_t test_et = testRegion->et();
    uint32_t test_rEta = testRegion->id().ieta();
    uint32_t test_rPhi = testRegion->id().iphi();
    UCTRegionIndex test_rIndex = g.getUCTRegionIndexFromL1CaloRegion(test_rEta, test_rPhi);
    UCTTowerIndex test_tIndex = g.getUCTTowerIndexFromL1CaloRegion(test_rIndex, test_raw);
    int test_cEta = test_tIndex.first;
    int test_cPhi = test_tIndex.second;
    bool test_negativeEta = g.getNegativeSide(test_cEta);

    if(testRegion->et()>0 && fabs(test_cEta)<28 && test_cPhi < 72){
      float pt = test_et;
      float eta = 0;
      if(fabs(test_cEta)<28){
        //eta = towerEtaMap[(int)fabs(test_cEta)];
        //eta = eta*fabs(test_cEta)/test_cEta;
        eta = getRecoEtaNew(test_cEta);
      }
      //float phi = towerPhiMap[test_cPhi];
      float phi = getRecoPhiNew(test_cPhi);
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(pt,eta,phi,pt);
      allRegions->push_back(temp);
    }
  }

  //evt.put(std::move(rgnCollection), "");

  if(!summaryCard->process()) {
    std::cerr << "UCT: Failed to process summary card" << std::endl;
    exit(1);      
  }

  double pt = 0;
  double eta = -999.;
  double phi = -999.;
  double mass = 0;
  double caloScaleFactor = 0.5;
  /* Do not bother with all of these things...  
  std::list<UCTObject*> emObjs = summaryCard->getEMObjs();
  for(std::list<UCTObject*>::const_iterator i = emObjs.begin(); i != emObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    nEGCands->push_back(L1EmParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1EmParticle::kNonIsolated));
  }
  std::list<UCTObject*> isoEMObjs = summaryCard->getIsoEMObjs();
  for(std::list<UCTObject*>::const_iterator i = isoEMObjs.begin(); i != isoEMObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    iEGCands->push_back(L1EmParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1EmParticle::kIsolated));
  }
  std::list<UCTObject*> tauObjs = summaryCard->getTauObjs();
  for(std::list<UCTObject*>::const_iterator i = tauObjs.begin(); i != tauObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    nTauCands->push_back(L1JetParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1JetParticle::kTau));
  }
  std::list<UCTObject*> isoTauObjs = summaryCard->getIsoTauObjs();
  for(std::list<UCTObject*>::const_iterator i = isoTauObjs.begin(); i != isoTauObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    iTauCands->push_back(L1JetParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1JetParticle::kTau));
  }
  std::list<UCTObject*> centralJetObjs = summaryCard->getCentralJetObjs();
  for(std::list<UCTObject*>::const_iterator i = centralJetObjs.begin(); i != centralJetObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    cJetCands->push_back(L1JetParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1JetParticle::kCentral));
    if(pt > 150.) {
      std::cout << "Jet: pt = " << pt << " eta = " << eta << " phi = " << phi << std::endl;
    }
  }
  std::list<UCTObject*> forwardJetObjs = summaryCard->getForwardJetObjs();
  for(std::list<UCTObject*>::const_iterator i = forwardJetObjs.begin(); i != forwardJetObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    fJetCands->push_back(L1JetParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1JetParticle::kForward));
  }
*/
  std::list<UCTObject*> boostedJetObjs = summaryCard->getBoostedJetObjs();
  for(std::list<UCTObject*>::const_iterator i = boostedJetObjs.begin(); i != boostedJetObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    bJetCands->push_back(L1JetParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1JetParticle::kCentral));// using kCentral for now, need a new type
    nTausInfo.push_back(object->nTaus());
  }

  /*
  const UCTObject* et = summaryCard->getET();
  pt = ((double) et->et()) * caloScaleFactor;
  double totET = pt;
  const UCTObject* met = summaryCard->getMET();
  pt = ((double) met->et()) * caloScaleFactor;
  eta = g.getUCTTowerEta(met->iEta());
  phi = g.getUCTTowerPhi(met->iPhi());
  metCands->push_back(L1EtMissParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1EtMissParticle::kMET, totET));
  */
  
  /*
  const UCTObject* ht = summaryCard->getHT();
  pt = ((double) ht->et()) * caloScaleFactor;
  double totHT = pt;
  const UCTObject* mht = summaryCard->getMHT();
  pt = ((double) mht->et()) * caloScaleFactor;
  eta = g.getUCTTowerEta(mht->iEta());
  phi = g.getUCTTowerPhi(mht->iPhi());
  mhtCands->push_back(L1EtMissParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1EtMissParticle::kMHT, totHT));
  */

  /*
  evt.put(std::move(iEGCands), "Isolated");
  evt.put(std::move(nEGCands), "NonIsolated");
  evt.put(std::move(iTauCands), "IsoTau");
  evt.put(std::move(nTauCands), "Tau");
  evt.put(std::move(cJetCands), "Central");
  evt.put(std::move(fJetCands), "Forward");
  evt.put(std::move(bJetCands), "Boosted");
  evt.put(std::move(metCands), "MET");
  evt.put(std::move(mhtCands), "MHT");
  */

  // Finish Running Layer 1

  // Start Runing Analysis
  Handle<vector<pat::Jet> > jets;
  if(evt.getByToken(jetSrc_, jets)){//Begin Getting Reco Taus
    for (const pat::Jet &jet : *jets) {
      recoJet_pt->Fill( jet.pt() );
      recoJet_eta->Fill( jet.eta() );
      recoJet_phi->Fill( jet.phi() );
      //get rid of the low pt stuff for analysis to save disk space
      if(jet.pt() > recoPt_ ) {
	goodJets.push_back(jet);
      }
    }
  }
  else
    cout<<"Error getting reco jets"<<std::endl;

  Handle<vector<pat::Jet> > jetsAK8;

  if(evt.getByToken(jetSrcAK8_, jetsAK8)){//Begin Getting Reco Jets
    for (const pat::Jet &jetAK8 : *jetsAK8) {
      //recoJetAK8_pt->Fill( jetAK8.pt() );
      //recoJetAK8_eta->Fill( jetAK8.eta() );
      //recoJetAK8_phi->Fill( jetAK8.phi() );
      //get rid of the cruft for analysis to save disk space
      if(jetAK8.pt() > recoPt_ ) {
	goodJetsAK8.push_back(jetAK8);
        TLorentzVector temp ;
        temp.SetPtEtaPhiE(jetAK8.pt(),jetAK8.eta(),jetAK8.phi(),jetAK8.et());
        ak8Jets->push_back(temp);
      }
    }
  }
  else
    cout<<"Error getting AK8 jets"<<std::endl;
  std::cout<<"AK8 jets size: "<<jetsAK8->size()<<std::endl;

  zeroOutAllVariables();

  if(goodJetsAK8.size()>0){

    for(auto jet:goodJetsAK8){
      tau1.push_back(jet.userFloat("NjettinessAK8Puppi:tau1"));
      tau2.push_back(jet.userFloat("NjettinessAK8Puppi:tau2"));
      tau3.push_back(jet.userFloat("NjettinessAK8Puppi:tau3"));
      nSubJets.push_back(jet.subjets("SoftDropPuppi").size());
      HFlav.clear();
      for(unsigned int isub=0; isub<((jet.subjets("SoftDropPuppi")).size()); isub++){
        HFlav.push_back(jet.subjets("SoftDropPuppi")[isub]->hadronFlavour());
        TLorentzVector temp;
        temp.SetPtEtaPhiE(jet.subjets("SoftDropPuppi")[isub]->pt(),jet.subjets("SoftDropPuppi")[isub]->eta(),jet.subjets("SoftDropPuppi")[isub]->phi(),jet.subjets("SoftDropPuppi")[isub]->et());
        subJets->push_back(temp);
      }
      subJetHFlav.push_back(HFlav);
      nBHadrons.push_back(jet.jetFlavourInfo().getbHadrons().size());
      std::cout<<"N subjets: "<< jet.subjets("SoftDropPuppi").size()<<std::endl;
      std::cout<<"N BHadrons: "<< jet.jetFlavourInfo().getbHadrons().size()<<std::endl;
      //take more variables from here: https://github.com/gouskos/HiggsToBBNtupleProducerTool/blob/opendata_80X/NtupleAK8/src/FatJetInfoFiller.cc#L215-L217
      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
    }

    //Match to boosted jets and see if we can match subjettiness functions...
    vector<l1extra::L1JetParticle> l1JetsSorted;
    for( vector<l1extra::L1JetParticle>::const_iterator l1Jet = bJetCands->begin(); l1Jet != bJetCands->end(); l1Jet++ ){
      l1JetsSorted.push_back(*l1Jet);
    }
    if(l1JetsSorted.size() > 1){  std::sort(l1JetsSorted.begin(),l1JetsSorted.end(),compareByPt);}
    pat::Jet recoJet_1;
    pat::Jet recoJet_2;

    recoPt_1  = goodJetsAK8.at(0).pt();
    recoEta_1 = goodJetsAK8.at(0).eta();
    recoPhi_1 = goodJetsAK8.at(0).phi();
    recoJet_1 = goodJetsAK8.at(0);

    if(goodJetsAK8.size()>1){
      recoJet_2 = goodJetsAK8.at(1);
      recoPt_2  = goodJetsAK8.at(1).pt();
      recoEta_2 = goodJetsAK8.at(1).eta();
      recoPhi_2 = goodJetsAK8.at(1).phi();
      recoDeltaEta = recoEta_1 - recoEta_2;
      recoDeltaPhi = recoPhi_1 - recoPhi_2;
      recoDeltaR = reco::deltaR(recoJet_1.p4(), recoJet_1.p4() );
      recoMass = (recoJet_1.p4() + recoJet_2.p4()).mass();
    }

    int i = 0;
    int foundL1Jet_1 = 0;
    int foundL1Jet_2 = 0;
    l1extra::L1JetParticle l1Jet_1;
    l1extra::L1JetParticle l1Jet_2;
    //for(std::list<UCTObject*>::const_iterator i = boostedJetObjs.begin(); i != boostedJetObjs.end(); i++) {
      //const UCTObject test = *i;
    if(l1JetsSorted.size() > 0){
      for(auto jet : l1JetsSorted){
        TLorentzVector temp;
        temp.SetPtEtaPhiE(jet.pt(),jet.eta(),jet.phi(),jet.et());
        l1Jets->push_back(temp);
        if(reco::deltaR(jet, recoJet_1)<0.4 && foundL1Jet_1 == 0 ){
          l1Jet_1 = jet;
          l1Pt_1  = jet.pt();
          l1Eta_1 = jet.eta();
          l1Phi_1 = jet.phi();
          l1NthJet_1 = i;
          l1NTau_1 = nTausInfo[i];
          foundL1Jet_1 = 1;
        }
        if(l1JetsSorted.size() > 1){
          if(recoPt_2 > 0 && reco::deltaR(jet, recoJet_2)<0.4 && foundL1Jet_2 == 0 ){
            l1Jet_2 = jet;
            l1Pt_2  = jet.pt();
            l1Eta_2 = jet.eta();
            l1Phi_2 = jet.phi();
            l1NthJet_2 = i;
            l1NTau_2 = nTausInfo[i];
            foundL1Jet_2 = 1;
          }
        }
        i++;
      }
    }
  }  

  efficiencyTree->Fill();
}

void BoostedJetStudies::zeroOutAllVariables(){
  genPt=-99; genEta=-99; genPhi=-99;
  recoPt=-99; recoEta=-99; recoPhi=-99;
  l1Pt=-99; l1Eta=-99; l1Phi=-99;
  genPt_1=-99; genEta_1=-99; genPhi_1=-99;
  recoPt_1=-99; recoEta_1=-99; recoPhi_1=-99;
  l1Pt_1=-99; l1Eta_1=-99; l1Phi_1=-99;
  genPt_2=-99; genEta_2=-99; genPhi_2=-99;
  recoPt_2=-99; recoEta_2=-99; recoPhi_2=-99;
  l1Pt_2=-99; l1Eta_2=-99; l1Phi_2=-99;
  genDeltaEta=-99; genDeltaPhi=-99; genDeltaR=-99; genMass=-99;
  recoDeltaEta=-99; recoDeltaPhi=-99; recoDeltaR=-99; recoMass=-99;
  l1DeltaEta=-99; l1DeltaPhi=-99; l1DeltaR=-99; l1Mass=-99;
  l1NthJet_1=-99; l1NthJet_2=-99;
  l1NTau_1=-99; l1NTau_2=-99;
  recoNthJet_1=-99; recoNthJet_2=-99;
  vbfBDT=-99; recoPt_=-99;
  nGenJets=-99; nRecoJets=-99; nL1Jets=-99;
  l1Matched_1=-99; l1Matched_2=-99;
 
  nSubJets.clear(); nBHadrons.clear(); subJetHFlav.clear();
  tau1.clear(); tau2.clear(); tau3.clear();
}

void BoostedJetStudies::print() {
  vector<UCTCrate*> crates = layer1->getCrates();
  for(uint32_t crt = 0; crt < crates.size(); crt++) {
    vector<UCTCard*> cards = crates[crt]->getCards();
    for(uint32_t crd = 0; crd < cards.size(); crd++) {
      vector<UCTRegion*> regions = cards[crd]->getRegions();
      for(uint32_t rgn = 0; rgn < regions.size(); rgn++) {
	if(regions[rgn]->et() > 10) {
	  int hitEta = regions[rgn]->hitCaloEta();
	  int hitPhi = regions[rgn]->hitCaloPhi();
	  vector<UCTTower*> towers = regions[rgn]->getTowers();
	  for(uint32_t twr = 0; twr < towers.size(); twr++) {
	    if(towers[twr]->caloPhi() == hitPhi && towers[twr]->caloEta() == hitEta) {
	      std::cout << "*";
	    }
	    if(towers[twr]->et() > 10) std::cout << *towers[twr];
	  }
	  std::cout << *regions[rgn];
	}
      }
      std::cout << *cards[crd];
    }
    std::cout << *crates[crt];
  }
  std::cout << *layer1;
}


void BoostedJetStudies::createBranches(TTree *tree){
    tree->Branch("run",     &run,     "run/I");
    tree->Branch("lumi",    &lumi,    "lumi/I");
    tree->Branch("event",   &event,   "event/I");

    tree->Branch("recoPt_1",      &recoPt_1,     "recoPt_1/D");
    tree->Branch("recoEta_1",     &recoEta_1,    "recoEta_1/D");
    tree->Branch("recoPhi_1",     &recoPhi_1,    "recoPhi_1/D");
    tree->Branch("recoNthJet_1",  &recoNthJet_1, "recoNthJet_1/I");

    tree->Branch("recoPt_2",      &recoPt_2,      "recoPt_2/D");
    tree->Branch("recoEta_2",     &recoEta_2,     "recoEta_2/D");
    tree->Branch("recoPhi_2",     &recoPhi_2,     "recoPhi_2/D");
    tree->Branch("recoNthJet_2",  &recoNthJet_2,  "recoNthJet_2/I");

    tree->Branch("recoDeltaEta",  &recoDeltaEta, "recoDeltaEta/D");
    tree->Branch("recoDeltaPhi",  &recoDeltaPhi, "recoDeltaPhi/D");
    tree->Branch("recoDeltaR",    &recoDeltaR,   "recoDeltaR/D");
    tree->Branch("recoMass",      &recoMass,     "recoMass/D");
      
    tree->Branch("l1Pt_1",        &l1Pt_1,       "l1Pt_1/D"); 
    tree->Branch("l1Eta_1",       &l1Eta_1,      "l1Eta_1/D");
    tree->Branch("l1Phi_1",       &l1Phi_1,      "l1Phi_1/D");
    tree->Branch("l1NthJet_1",    &l1NthJet_1,   "l1NthJet_1/I");
    tree->Branch("l1NTau_1",      &l1NTau_1,     "l1NTau_1/I");

    tree->Branch("l1Pt_2",        &l1Pt_2,       "l1Pt_2/D"); 
    tree->Branch("l1Eta_2",       &l1Eta_2,      "l1Eta_2/D");
    tree->Branch("l1Phi_2",       &l1Phi_2,      "l1Phi_2/D");
    tree->Branch("l1NthJet_2",    &l1NthJet_2,   "l1NthJet_2/I");
    tree->Branch("l1NTau_2",      &l1NTau_2,     "l1NTau_2/I");

    tree->Branch("l1DeltaEta",    &l1DeltaEta,   "l1DeltaEta/D");
    tree->Branch("l1DeltaPhi",    &l1DeltaPhi,   "l1DeltaPhi/D");
    tree->Branch("l1DeltaR",      &l1DeltaR,     "l1DeltaR/D");
    tree->Branch("l1Mass",        &l1Mass,       "l1Mass/D");

    tree->Branch("l1Matched_1",   &l1Matched_1, "l1Matched_1/I");
    tree->Branch("l1Matched_2",   &l1Matched_2, "l1Matched_2/I");
    tree->Branch("nRecoJets",     &nRecoJets,    "nRecoJets/I");
    tree->Branch("nL1Jets",       &nL1Jets,      "nL1Jets/I");
    tree->Branch("vbfBDT",        &vbfBDT,       "vbfBDT/D");

    tree->Branch("tau1",          &tau1);
    tree->Branch("tau2",          &tau2);
    tree->Branch("tau3",          &tau3);
    tree->Branch("nSubJets",      &nSubJets);
    tree->Branch("subJetHFlav",   &subJetHFlav);
    tree->Branch("nBHadrons",     &nBHadrons);

    tree->Branch("allRegions", "vector<TLorentzVector>", &allRegions, 32000, 0);
    tree->Branch("hcalTPGs", "vector<TLorentzVector>", &allHcalTPGs, 32000, 0);
    tree->Branch("ecalTPGs", "vector<TLorentzVector>", &allEcalTPGs, 32000, 0);
    tree->Branch("caloClusters", "vector<TLorentzVector>", &caloClusters, 32000, 0);
    tree->Branch("l1Jets", "vector<TLorentzVector>", &l1Jets, 32000, 0);
    tree->Branch("ak8Jets", "vector<TLorentzVector>", &ak8Jets, 32000, 0);
    tree->Branch("subJets", "vector<TLorentzVector>", &subJets, 32000, 0);
  }


// ------------ method called once each job just before starting event loop  ------------
void 
BoostedJetStudies::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BoostedJetStudies::endJob() {
}

// ------------ method called when starting to processes a run  ------------

void
BoostedJetStudies::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  if(!L1TCaloLayer1FetchLUTs(iSetup, ecalLUT, hcalLUT, hfLUT, useLSB, useCalib, useECALLUT, useHCALLUT, useHFLUT)) {
    std::cerr << "L1TCaloLayer1::beginRun: failed to fetch LUTS - using unity" << std::endl;
  }
  for(uint32_t twr = 0; twr < twrList.size(); twr++) {
    twrList[twr]->setECALLUT(&ecalLUT);
    twrList[twr]->setHCALLUT(&hcalLUT);
    twrList[twr]->setHFLUT(&hfLUT);
  }
}
 
// ------------ method called when ending the processing of a run  ------------
/*
  void
  BoostedJetStudies::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  BoostedJetStudies::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  BoostedJetStudies::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BoostedJetStudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BoostedJetStudies);
