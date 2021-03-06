#ifndef L1TRegionNtupleProducer_H
#define L1TRegionNtupleProducer_H

// system include files
#include <memory>
#include <unistd.h>
#include <iostream>
#include <fstream>

#include <iostream>
#include <fstream>
#include <vector>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

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


#include <memory>
#include <math.h>
#include <vector>
#include <list>

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctJetCand.h"
#include "Math/LorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "L1Trigger/L1TCaloLayer1/src/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"

//
// class declaration
//
using std::vector;

class L1TRegionNtupleProducer : public edm::EDAnalyzer {

 public:
  
  // Constructor
  L1TRegionNtupleProducer(const edm::ParameterSet& ps);
  
  // Destructor
  virtual ~L1TRegionNtupleProducer();

  edm::Service<TFileService> tfs_;

  std::ofstream file0, file1, file10;
  
  std::vector<float> vRegionEt;
  std::vector<float> vRegionEta;
  std::vector<float> vRegionPhi;
  std::vector<float> vRegionTau;
  std::vector<float> vRegionEG;

  TH1F* nEvents;

  TH1F* regionEta;
  TH1F* regionPhi;
  TH1F* regionPt;
  TH1F* regionEtaFine;
  TH1F* regionPhiFine;
  TH1F* regionTotal;

  TH1F* regionHitEta;
  TH1F* regionHitPhi;
  TTree* regionTree;
  TFileDirectory folder;

  int run, lumi, event;
  float nvtx;
  void antikt();
  void initializeHCALTPGMap(const edm::Handle<HcalTrigPrimDigiCollection> hcal, const  edm::ESHandle<L1CaloHcalScale> hcalScale, double hTowerETMap[73][57], bool testMode = false);
  void initializeECALTPGMap(edm::Handle<EcalTrigPrimDigiCollection> ecal, double eTowerETMap[73][57], bool testMode = false);

 protected:
  // Analyze
  void analyze(const edm::Event& evt, const edm::EventSetup& es);
  
  // BeginJob
  void beginJob(const edm::EventSetup &es);
  
  // EndJob
  void endJob(void);

  
 private:
  // ----------member data ---------------------------

  int nev_; // Number of events processed
  bool verbose_;
  std::ofstream logFile_;
  edm::InputTag rctSource_; 

  edm::EDGetTokenT<vector<pat::PackedCandidate> > pfCandsToken_;  
  edm::EDGetTokenT<L1CaloRegionCollection> L1RegionCollection;
  edm::EDGetTokenT<L1CaloEmCollection> L1EMCollection_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;
  edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalSrc_; 
  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalSrc_;
  //edm::EDGetTokenT<double> recoPt_;
  //edm::EDGetTokenT<std::string> folderName_;
  edm::EDGetTokenT<reco::VertexCollection> vtxLabel_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> discriminatorMu_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> discriminatorIso_;
  edm::EDGetTokenT<vector<pat::Tau> > tauSrc_;
  //edm::EDGetTokenT<L1GctJetCandCollection> gctIsoTauJetsSource_;
  //edm::EDGetTokenT<L1GctJetCandCollection> gctTauJetsSource_;
  edm::EDGetTokenT<vector <l1extra::L1JetParticle> > l1ExtraIsoTauSource_;
  edm::EDGetTokenT<vector <l1extra::L1JetParticle> > l1ExtraTauSource_;
  edm::EDGetTokenT<BXVector <l1t::Tau> > l1Stage2TauSource_;
  edm::EDGetTokenT<BXVector <l1t::Tau> > l1Stage1TauSource_;
  edm::EDGetTokenT<BXVector <l1t::Tau> > l1Stage1IsoTauSource_;
  edm::EDGetTokenT<vector <L1GctJetCand> > l1GctTauSource_;
  edm::EDGetTokenT<vector <L1CaloRegion> > regionSource_;

  std::string folderName_;
  double recoPt_;

		 
 int TPGEtaRange(int ieta){
   int iEta = 0;
   // So here, -28 becomes 0.  -1 be comes 27.  +1 becomes 28. +28 becomes 55.
   // And we have mapped [-28, -1], [1, 28] onto [0, 55]   
   if(ieta < 0)
     iEta = ieta + 28;
   else if(ieta > 0)
     iEta = ieta + 27;
   return iEta;
 }

  int convertGenEta(double inputEta) {
    const double tpgEtaValues[27] = {
      0.087,      
      0.174, // HB and inner HE bins are 0.348 wide
      0.261,
      0.348,
      0.522,
      0.609,
      0.696,
      0.783,
      0.870,
      0.957,
      1.044,
      1.131,
      1.218,
      1.305,
      1.392,
      1.479,
      1.566,
      1.653,
      1.74,
      1.848,
      1.956, // Last two HE bins are 0.432 and 0.828 wide
      2.064,
      2.172,
      2.379,
      2.586,
      2.793,
      3
      //IGNORING HF
      //3.250, // HF bins are 0.5 wide
      //3.750,
      //4.250,
      //4.750
    };


    for (int n=1; n<29; n++){
      //std::cout<<"inputEta "<<inputEta<< " n "<< n <<" tpgEtaValues[n-1] "<< tpgEtaValues[n-1] << " abs(inputEta)<tpgEtaValues[n-1]"<<std::endl;
      if (std::fabs(inputEta)<tpgEtaValues[n-1]) {
	//std::cout<<"found to be true"<<std::endl;
	//int tpgEta = n;
	//Positive eta is >28
	//negative eta is 0 to 27
	if(inputEta>0){
	  //std::cout<<"returning input eta >0 so + 28"<<std::endl;
	  return n + 28;}
	else{
	  //std::cout<<"returning input eta <0 so n"<<std::endl;
	  return n;}
	break;
      }
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputeta: "<<inputEta<<std::endl;
    return -9;
  }

  //-pi < phi <= +pi,
  int convertGenPhi(double inputPhi){
    double posPhi[36];
    for(int n = 0; n < 36; n++)
      posPhi[n] = (0.087) * n + 0.0435;
    double negPhi[36];
    for(int n = 0; n < 36; n++)
      negPhi[n] = -3.14159 + 0.087 * n - 0.0435;

    //1 to 36 is 0 to pi
    if( 3.1416 > inputPhi && inputPhi >= 0){

      for(int n = 1; n < 36; n++){
	//std::cout<<"inputPhi "<<inputPhi<< " posPhi[n-1] "<< posPhi[n-1] << " n "<<n<<std::endl;
	if(inputPhi <= posPhi[n-1]){
	  int tpgPhi = n;
	  return tpgPhi;
	}
      }
    }

    //37 to 72 is -pi to 0
    else if(-3.1416 < inputPhi && inputPhi < 0){
      for(int n = 1; n < 36; n++)
	if(inputPhi < negPhi[n-1]){
	  int tpgPhi = n + 36;
	  return tpgPhi;
	}
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputphi: "<<inputPhi<<std::endl;
    return -9;
  }


  float convertRCTEta(uint32_t inputEta) {
    const double regionEtaValues[22] = {
      -4.75,
      -4.25,
      -3.75,
      -3.25,
      -2.5,
      -1.93,
      -1.566,
      -1.218,
      -0.87,
      -0.522,
      -0.174,
      0.174,
      0.522,
      0.87,
      1.218,
      1.566,
      1.93,
      2.5,
      3.25,
      3.75,
      4.25,
      4.75
    };
    return regionEtaValues[inputEta];
  };



  float convertRCTPhi(uint32_t inputPhi) {
    const double regionPhiValues[20] = {
      0.000,
      0.349,
      0.698,
      1.047,
      1.396,
      1.744,
      2.093,
      2.442,
      2.791,
      -3.14159,
      -2.791,
      -2.442,
      -2.093,
      -1.744,
      -1.396,
      -1.047,
      -0.698,
      -0.349
    };
    return regionPhiValues[inputPhi];
  };

float towerEtaMap[28]=   { 
  //-2.913, //-2.739, 
      //-2.565, -2.391, 2.217, 	//switch to smaller trigger towers here    
      //-2.0445, -1.9575, -1.8705, -1.7835, -1.6965, 
      //-1.6095, -1.5225, -1.4355, -1.3485, -1.2615, 
      //-1.1745, -1.0875, -1.0005, -0.9135, -0.8265, 
      //-0.7395, -0.6525, -0.5655, -0.4785, -0.3915, 
      //-0.3045, -0.2175, -0.1305, -0.0435, 
  0.0435, 
  0.1305, 0.2175, 0.3045, 0.3915, 0.4785, 
  0.5655, 0.6525, 0.7395, 0.8265, 0.9135, 
  1.0005, 1.0875, 1.1745, 1.2615, 1.3485, 
  1.4355, 1.5225, 1.6095, 1.6965, 1.7835, 
  1.8705, 1.9575, 2.0445, 2.217, 2.391, 
  2.565, //2.739,
  2.913,
  };

  //this is for mapping from RCT only
float towerPhiMap[72]=                        
  {-0.131, -0.044, 0.044, 0.131, 0.218, 0.305, 0.393, 0.480, 0.567, 0.654, 0.742, 0.829, 0.916, 1.004, 1.091, 1.178, 1.265, 1.353, 1.440, 1.527, 1.614, 1.702, 1.789, 1.876, 1.963, 2.051, 2.138, 2.225, 2.313, 2.400, 2.487, 2.574, 2.662, 2.749, 2.836, 2.923, 3.011, 3.098,
    -3.098, -3.011, -2.923, -2.836, -2.749, -2.662, -2.574, -2.487, -2.400, -2.313, -2.225, -2.138, -2.051, -1.963, -1.876, -1.789, -1.702, -1.614, -1.527, -1.440, -1.353, -1.265, -1.178, -1.091, -1.004, -0.916, -0.829, -0.742, -0.654, -0.567, -0.480, -0.393, -0.305, -0.218};


};



#endif
