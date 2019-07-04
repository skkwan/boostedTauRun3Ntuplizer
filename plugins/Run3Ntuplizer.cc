/*
 * \file Run3Ntuplizer.cc
 *
 * \author I. Ojalvo
 * Written for miniAOD
 */

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "L1Trigger/Run3Ntuplizer/interface/Run3Ntuplizer.h"
#include "L1Trigger/Run3Ntuplizer/interface/L1TRegionNtupleProducer.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Common/interface/RefToPtr.h"


#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <fstream>

using namespace edm;
using std::cout;
using std::endl;
using std::vector;

bool compareByPtJets (l1extra::L1JetParticle i,l1extra::L1JetParticle j) { return(i.pt()>j.pt()); };

Run3Ntuplizer::Run3Ntuplizer( const ParameterSet & cfg ) :
  ecalSrc_(consumes<EcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  jetSrc_(consumes<vector<pat::Jet> >(cfg.getParameter<edm::InputTag>("recoJets"))),
  jetSrcAK8_(consumes<vector<pat::Jet> >(cfg.getParameter<edm::InputTag>("recoJetsAK8"))),
  regionSource_(consumes<vector <L1CaloRegion> >(cfg.getParameter<edm::InputTag>("UCTRegion"))),
  centralJets_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1UCTCentralJets"))),
  forwardJets_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1UCTForwardJets")))
  {


    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    recoPt_              = cfg.getParameter<double>("recoPtCut");
    folder               = tfs_->mkdir(folderName_);
    //folder->cd();
    regionTree = folder.make<TTree>("EfficiencyTree", "Efficiency Tree");
    regionTree->Branch("run",        &run,     "run/I");
    regionTree->Branch("lumi",       &lumi,    "lumi/I");
    regionTree->Branch("event",      &event,   "event/I");
    regionTree->Branch("nvtx",       &nvtx,    "nvtx/D");
    regionTree->Branch("vRegionEt",  &vRegionEt  );
    regionTree->Branch("vRegionEta", &vRegionEta );
    regionTree->Branch("vRegionPhi", &vRegionPhi );
    regionTree->Branch("vRegionEG",  &vRegionEG  );
    regionTree->Branch("vRegionTau", &vRegionTau );

    nEvents       = folder.make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );

    regionHitEta  = folder.make<TH1F>( "regionHit_eta"  , "eta", 16, 1, 16. );
    regionHitPhi  = folder.make<TH1F>( "regionHit_phi"  , "phi", 16, 1, 16. );
    regionTotal   = folder.make<TH1F>( "regionHit_total"  , "fullmap", 16, 1, 16. );

    regionEta     = folder.make<TH1F>( "region_eta"  , "eta", 22, 1, 22. );
    regionPhi     = folder.make<TH1F>( "region_phi"  , "phi", 72, 1, 72. );
    regionPt      = folder.make<TH1F>( "region_pt"  , "pt", 100, 0, 100. );

    regionEtaFine   = folder.make<TH1F>( "region_eta_Fine"  , "eta", 88, 1, 88. );
    regionPhiFine   = folder.make<TH1F>( "region_phi_Fine"  , "phi", 72, 1, 72. );

    recoJet_pt   = folder.make<TH1F>( "recoJet_pt" , "p_{t}", 300,  0., 300. );
    recoJet_eta  = folder.make<TH1F>( "recoJet_eta"  , "eta", 100,  -3, 3. );
    recoJet_phi  = folder.make<TH1F>( "recoJet_phi"  , "phi", 100,  -4, 4. );

    recoJetAK8_pt   = folder.make<TH1F>( "recoJetAK8_pt" , "p_{t}", 300,  0., 300. );
    recoJetAK8_eta  = folder.make<TH1F>( "recoJetAK8_eta"  , "eta", 100,  -3, 3. );
    recoJetAK8_phi  = folder.make<TH1F>( "recoJetAK8_phi"  , "phi", 100,  -4, 4. );

    efficiencyTreeAK8 = folder.make<TTree>("EfficiencyTreeAK8", "Efficiency Tree AK8");
    efficiencyTreeAK8->Branch("run",    &run,     "run/I");
    efficiencyTreeAK8->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTreeAK8->Branch("event",  &event,   "event/I");
    efficiencyTreeAK8->Branch("nvtx",   &nvtx,     "nvtx/D");
    
    efficiencyTreeAK8->Branch("recoPt",    &recoPtAK8,   "recoPt/D");
    
    efficiencyTreeAK8->Branch("jetPt",     &jetPtAK8, "jetPt/D");

    efficiencyTreeAK8->Branch("recoEta",    &recoEtaAK8,   "recoEta/D");
    efficiencyTreeAK8->Branch("jetEta",     &jetEtaAK8, "jetEta/D");
    
    efficiencyTreeAK8->Branch("recoPhi",    &recoPhiAK8,   "recoPhi/D");
    efficiencyTreeAK8->Branch("jetPhi",     &jetPhiAK8, "jetPhi/D");

    efficiencyTreeAK8->Branch("l1Matched",  &l1MatchedAK8, "l1Matched/I");

    efficiencyTree = folder.make<TTree>("EfficiencyTree", "Efficiency Tree ");
    efficiencyTree->Branch("run",     &run,     "run/I");
    efficiencyTree->Branch("lumi",    &lumi,     "lumi/I");
    efficiencyTree->Branch("event",   &event,    "event/I");
    efficiencyTree->Branch("nvtx",    &nvtx,     "nvtx/D");
    
    efficiencyTree->Branch("recoPt",  &recoPt,   "recoJetPt/D");
    
    efficiencyTree->Branch("jetPt",   &jetPt,    "l1JetPt/D"); 

    efficiencyTree->Branch("recoEta", &recoEta,  "recoJetEta/D");
    efficiencyTree->Branch("jetEta",  &jetEta,   "l1JetEta/D");
    
    efficiencyTree->Branch("recoPhi", &recoPhi,   "recoJetPhi/D");
    efficiencyTree->Branch("jetPhi",  &jetPhi,    "l1JetPhi/D");

    efficiencyTree->Branch("l1Matched",  &l1Matched, "l1Matched/I");


  }

void Run3Ntuplizer::beginJob( const EventSetup & es) {
   std::cout<<"begin job..."<<std::endl;
}

void Run3Ntuplizer::analyze( const Event& evt, const EventSetup& es )
 {
   std::cout<<"Analyzing..."<<std::endl;
   nEvents->Fill(1);
   
   run = evt.id().run();
   lumi = evt.id().luminosityBlock();
   event = evt.id().event();
   Handle<L1CaloRegionCollection> regions;
   
   std::vector<pat::Jet> goodJets;
   std::vector<pat::Jet> goodJetsAK8;
  
   edm::Handle < vector<l1extra::L1JetParticle> > l1CentralJets;
   edm::Handle < vector<l1extra::L1JetParticle> > l1ForwardJets;

   edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs;
   edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;
   
   
  if(!evt.getByToken(ecalSrc_, ecalTPGs))
    std::cout<<"ERROR GETTING THE ECAL TPGS"<<std::endl;
  if(!evt.getByToken(hcalSrc_, hcalTPGs))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;

  if(!evt.getByToken(centralJets_, l1CentralJets))
    std::cout<<"ERROR GETTING THE CENTRAL JETS"<<std::endl;

  if(!evt.getByToken(forwardJets_, l1ForwardJets))
    std::cout<<"ERROR GETTING THE FORWARD JETS"<<std::endl;

  //sort the L1 jets
  vector<l1extra::L1JetParticle> l1JetsSorted;
  vector<l1extra::L1JetParticle> l1JetsSortedEtaRestricted2p4;
  for( vector<l1extra::L1JetParticle>::const_iterator l1Jet = l1CentralJets->begin(); l1Jet != l1CentralJets->end(); l1Jet++ ){
    l1JetsSorted.push_back(*l1Jet);
    if(abs(l1Jet->eta()) < 2.4) l1JetsSortedEtaRestricted2p4.push_back(*l1Jet);
  }

  for( vector<l1extra::L1JetParticle>::const_iterator l1Jet = l1ForwardJets->begin(); l1Jet != l1ForwardJets->end(); l1Jet++ ){
    l1JetsSorted.push_back(*l1Jet);
    if(abs(l1Jet->eta()) < 2.4) l1JetsSortedEtaRestricted2p4.push_back(*l1Jet);
  }

  std::sort(l1JetsSorted.begin(),l1JetsSorted.end(),compareByPtJets);
  std::sort(l1JetsSortedEtaRestricted2p4.begin(),l1JetsSortedEtaRestricted2p4.end(),compareByPtJets);

  ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

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
    std::cout<<"Error getting reco jets"<<std::endl;

  Handle<vector<pat::Jet> > jetsAK8;
  
  if(evt.getByToken(jetSrcAK8_, jetsAK8)){//Begin Getting Reco Jets
    for (const pat::Jet &jetAK8 : *jetsAK8) {
      recoJetAK8_pt->Fill( jetAK8.pt() );
      recoJetAK8_eta->Fill( jetAK8.eta() );
      recoJetAK8_phi->Fill( jetAK8.phi() );
      //get rid of the cruft for analysis to save disk space
      if(jetAK8.pt() > recoPt_ ) {
	goodJetsAK8.push_back(jetAK8);

      }
    }
  }
  else
    std::cout<<"Error getting AK8 jets"<<std::endl;


  //fill the jet variables here!
  //include delta eta between the two jets,
  // pt, eta, phi of each jet
  //delta phi between two jets
  //invariant mass of two jets
  //total number of jets in the event
  for(auto jet : goodJets){
    recoPt =  jet.pt();
    recoEta =  jet.eta();
    recoPhi =  jet.phi();
    jetPt = -99;
    jetEta = -99;
    jetPhi = -99;
    for(auto l1jet : l1JetsSorted){
      if(reco::deltaR(jet.eta(), jet.phi(), l1jet.eta(), l1jet.phi()< 0.2)){
	//now fill the tree!
	jetPt = l1jet.pt();
	jetEta = l1jet.eta();
	jetPhi = l1jet.phi();	
      }
    }

    efficiencyTree->Fill();
  }
  

  std::cout<<"making regions"<<std::endl;
  vRegionEt.clear();
  vRegionEta.clear();
  vRegionPhi.clear();
  vRegionTau.clear();
  vRegionEG.clear();

  UCTGeometry g;
  //************* Get Regions and make region plots
  if(!evt.getByToken(regionSource_,regions)){
    std::cout<<"ERROR GETTING THE REGIONS!!!"<<std::endl;}
  else{
    for(vector<L1CaloRegion>::const_iterator testRegion = regions->begin(); testRegion != regions->end(); ++testRegion){
      //UCTRegionProcess uctRegion(*region);
       uint16_t test_raw = testRegion->raw();
       uint32_t test_et = testRegion->et();
       //testRegionTotET += test_et;
       uint32_t test_rEta = testRegion->id().ieta();
       uint32_t test_rPhi = testRegion->id().iphi();
       UCTRegionIndex test_rIndex = g.getUCTRegionIndexFromL1CaloRegion(test_rEta, test_rPhi);
       UCTTowerIndex test_tIndex = g.getUCTTowerIndexFromL1CaloRegion(test_rIndex, test_raw);
       int test_cEta = test_tIndex.first;
       int test_cPhi = test_tIndex.second;
       bool test_negativeEta = g.getNegativeSide(test_cEta);
       //uint32_t test_crate = g.getCrate(test_cEta, test_cPhi);
       //uint32_t test_card = g.getCard(test_cEta, test_cPhi);
       //uint32_t test_region = g.getRegion(test_cEta, test_cPhi);
       //uint32_t test_iEta = g.getiEta(test_cEta);
       //uint32_t test_iPhi = g.getiPhi(test_cPhi);

       if(testRegion->et()>0 && fabs(test_cEta)<28 && test_cPhi < 72){
	 float pt = test_et;
	 float eta = 0;
	 if(fabs(test_cEta)<28){
	   eta = towerEtaMap[(int)fabs(test_cEta)];
	   eta = eta*fabs(test_cEta)/test_cEta;
	 }

	 float phi = towerPhiMap[test_cPhi];
	 float isEgammaLike = 0;
	 float isTauLike    = 0;
	 if(!((l1tcalo::RegionEGVeto & testRegion->raw()) == l1tcalo::RegionEGVeto))
	   isEgammaLike = 1;

	 if(!((l1tcalo::RegionTauVeto & testRegion->raw()) == l1tcalo::RegionTauVeto))
	   isTauLike = 1;
	 //std::cout<<"region eta,phi: "<<eta<<" , "<<phi<<std::endl;
	 vRegionEt.push_back(pt);
	 vRegionEta.push_back(eta);
	 vRegionPhi.push_back(phi);
	 vRegionEG.push_back(isEgammaLike);
	 vRegionEG.push_back(isTauLike);
	 
	 regionEta->Fill(eta);
	 regionPhi->Fill(phi);
	 regionPt->Fill(pt);
	
      }
    }
    regionTree->Fill();
  }


}


  
/*
 * Get the ECAL TPGS create a TPG map for the event
 *
 */
  
void Run3Ntuplizer::initializeECALTPGMap(Handle<EcalTrigPrimDigiCollection> ecal, double eTowerETMap[73][57], bool testMode){
  
  //std::cout << "ECAL TPGS" << std::endl;
  for (size_t i = 0; i < ecal->size(); ++i) {
    int cal_ieta = (*ecal)[i].id().ieta();
    int cal_iphi = (*ecal)[i].id().iphi();
    int iphi = cal_iphi-1;
    int ieta = TPGEtaRange(cal_ieta);
    // TPG iPhi starts at 1 and goes to 72.  Let's index starting at zero.
    // TPG ieta ideal goes from 0-55.
    double LSB = 0.5;
    double et= (*ecal)[i].compressedEt()*LSB;
    //if(et>0)std::cout<<"et "<< et<<std::endl;
    if(testMode && iphi == 34 && ieta == 11){
      et = 40;
    }

    //if(et>0)
    //cout<<"Before filling eTower"
    //<<"ECAL ieta:"<<ieta<<" cal_ieta:"<< cal_ieta<<" iphi:"<<iphi<<" et:"<<et<<endl;

    if (iphi >= 0 && iphi <= 72 &&
	ieta >= 0 && ieta <= 55) {
      eTowerETMap[iphi][ieta] = et; 
    }

  }

}

 void Run3Ntuplizer::initializeHCALTPGMap(const Handle<HcalTrigPrimDigiCollection> hcal, 
					 const ESHandle<L1CaloHcalScale>  hcalScale, 
					 double hTowerETMap[73][57], bool testMode){
  for (size_t i = 0; i < hcal->size(); ++i) {
    HcalTriggerPrimitiveDigi tpg = (*hcal)[i];
    int cal_ieta = tpg.id().ieta();
    int cal_iphi = tpg.id().iphi();
    int iphi = cal_iphi-1;
    int ieta = TPGEtaRange(cal_ieta);
    short absieta = std::abs(tpg.id().ieta());
    short zside = tpg.id().zside();
    double energy = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
    //if(energy>0)std::cout<<"energy "<< energy<<std::endl;

    if(testMode && iphi == 34 && ieta == 12){
      energy = 40;
    }

    if (iphi >= 0 && iphi <= 71 &&
	ieta >= 0 && ieta <= 55) {

      //(*hcal)[i].SOI_compressedEt(), absieta, zside)*LSB; //*LSB
      //if(energy>0)
      //std::cout<<"hcal iphi "<<iphi<<" ieta "<<ieta<<" energy "<<energy<<std::endl;
      hTowerETMap[iphi][ieta] = energy;
      //TPGSum_ +=energy;
      //TPGH_ += energy;
      //double alpha_h = TPGSFp_[cal_ieta]; //v3
      //hCorrTowerETMap[cal_iphi][cal_ieta] = alpha_h*energy;
      //cTPGH_ += alpha_h*energy;
      //if (energy > 0) {
      //std::cout << "hcal eta/phi=" << ieta << "/" << iphi
      //<< " = (" << getEtaTPG(ieta) << "/" << getPhiTPG(iphi) << ") "
      //<< " et=" << (*hcal)[i].SOI_compressedEt()
      //<< " energy=" << energy
      //<< " rctEta="<< twrEta2RegionEta(cal_ieta) << " rctPhi=" << twrPhi2RegionPhi(cal_iphi)
      //<< " fg=" << (*hcal)[i].SOI_fineGrain() << std::endl;
      //}
      //if (energy>maxTPGHPt){
      //maxTPGHPt=energy;
      //maxTPGHPt_phi = cal_iphi; //this one starts at 0-72
      //maxTPGHPt_eta = cal_ieta; //this one is 0-54
      //} 
    }
    //else
      //std::cout<<"HCAL failed checks iphi "<<iphi<<" ieta "<<ieta<<std::endl;
  }//end HCAL TPG
}

void Run3Ntuplizer::endJob() {
}

Run3Ntuplizer::~Run3Ntuplizer(){

}

DEFINE_FWK_MODULE(Run3Ntuplizer);
