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
  forwardJets_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1UCTForwardJets"))),
  genJets_(consumes<vector <reco::GenJet> >(cfg.getParameter<edm::InputTag>("genJets")))
  {


    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    recoPt_              = cfg.getParameter<double>("recoPtCut");
    isData_              = cfg.getParameter<bool>("isData");
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

    genTree = folder.make<TTree>("genJetTree", "gen Matched Jet Tree ");
    createBranches(genTree);
    createBranchesGen(genTree);

    efficiencyTree = folder.make<TTree>("efficiencyTree", "Reco Matched Jet Tree ");
    createBranches(efficiencyTree);

  }


void Run3Ntuplizer::createBranches(TTree *tree){
    tree->Branch("run",     &run,     "run/I");
    tree->Branch("lumi",    &lumi,     "lumi/I");
    tree->Branch("event",   &event,    "event/I");

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

    tree->Branch("l1Pt_2",        &l1Pt_2,       "l1Pt_2/D"); 
    tree->Branch("l1Eta_2",       &l1Eta_2,      "l1Eta_2/D");
    tree->Branch("l1Phi_2",       &l1Phi_2,      "l1Phi_2/D");
    tree->Branch("l1NthJet_2",    &l1NthJet_2,   "l1NthJet_2/I");

    tree->Branch("l1DeltaEta",    &l1DeltaEta,   "l1DeltaEta/D");
    tree->Branch("l1DeltaPhi",    &l1DeltaPhi,   "l1DeltaPhi/D");
    tree->Branch("l1DeltaR",      &l1DeltaR,     "l1DeltaR/D");
    tree->Branch("l1Mass",        &l1Mass,       "l1Mass/D");

    tree->Branch("l1Matched_1",   &l1Matched_1, "l1Matched_1/I");
    tree->Branch("l1Matched_2",   &l1Matched_2, "l1Matched_2/I");
    tree->Branch("nRecoJets",     &nRecoJets,    "nRecoJets/I");
    tree->Branch("nL1Jets",       &nL1Jets,      "nL1Jets/I");

  }

void Run3Ntuplizer::createBranchesGen(TTree *tree){
    tree->Branch("genPt_1",  &genPt_1,   "genPt_1/D");
    tree->Branch("genEta_1", &genEta_1,  "genEta_1/D");
    tree->Branch("genPhi_1", &genPhi_1,  "genPhi_1/D");

    tree->Branch("genPt_2",  &genPt_2,   "genPt_2/D");
    tree->Branch("genEta_2", &genEta_2,  "genEta_2/D");
    tree->Branch("genPhi_2", &genPhi_2,  "genPhi_2/D");

    tree->Branch("genDeltaEta",   &genDeltaEta,   "genDeltaEta/D");
    tree->Branch("genDeltaPhi",   &genDeltaPhi,   "genDeltaPhi/D");
    tree->Branch("genDeltaR",     &genDeltaR,     "genDeltaR/D");
    tree->Branch("genMass",       &genMass,       "genMass/D");
    
    tree->Branch("genMatched_1",  &genMatched_1, "genMatched_1/I");
    tree->Branch("genMatched_2",  &genMatched_2, "genMatched_2/I");

    tree->Branch("nGenJets",      &nGenJets,     "nGenJets/I");
	
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

   edm::Handle < vector<reco::GenJet> > genJets;
   
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

  if(!isData_)
    if(!evt.getByToken(genJets_, genJets))
      std::cout<<"ERROR GETTING THE GEN JETS"<<std::endl;

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

  zeroOutAllVariables();

  //fill the jet variables here!
  //include delta eta between the two jets,
  // pt, eta, phi of each jet
  //delta phi between two jets
  //invariant mass of two jets
  //total number of jets in the event
  reco::GenJet genJet_1;
  reco::GenJet genJet_2;
  if(!isData_)
    if(genJets->size()>0){
      genPt_1  = genJets->at(0).pt();
      genEta_1 = genJets->at(0).eta();
      genPhi_1 = genJets->at(0).phi();
      genJet_1 = genJets->at(0);
      
      if(genJets->size()>1){
	genPt_2  = genJets->at(1).pt();
	genEta_2 = genJets->at(1).eta();
	genPhi_2 = genJets->at(1).phi();
	genJet_2 = genJets->at(1);
	
	genDeltaEta = genEta_1 - genEta_2;
	genDeltaPhi = genPhi_1 - genPhi_2;
	genDeltaR = reco::deltaR(genJet_1,genJet_2 );
	genMass = (genJet_1.p4() + genJet_2.p4()).mass();
	
      }
      
      int i = 0;
      int foundRecoJet_1 = 0;
      int foundRecoJet_2 = 0;
      for(auto jet : goodJets){
	if(reco::deltaR(jet, genJet_1)<0.1 && foundRecoJet_1 == 0 ){
	  recoPt_1  = jet.pt();
	  recoEta_1 = jet.eta();
	  recoPhi_1 = jet.phi();
	  recoNthJet_1 = i;
	  foundRecoJet_1 = 1;
	}
	if(genPt_2 > 0 && reco::deltaR(jet, genJet_2)<0.1 && foundRecoJet_2 == 0 ){
	  recoPt_2  = jet.pt();
	  recoEta_2 = jet.eta();
	  recoPhi_2 = jet.phi();
	  recoNthJet_2 = i;
	  foundRecoJet_2 = 1; 
	}
	i++;
      }
      
      if(foundRecoJet_1>0 && foundRecoJet_2>0){
	recoDeltaEta = recoEta_1 - recoEta_2;
	recoDeltaPhi = recoPhi_1 - recoPhi_2;
	recoDeltaR = reco::deltaR(goodJets.at(recoNthJet_1), goodJets.at(recoNthJet_2) );
	recoMass = (goodJets.at(recoNthJet_1).p4() + goodJets.at(recoNthJet_2).p4()).mass();
      }
      
      
      i = 0;
      int foundL1Jet_1 = 0;
      int foundL1Jet_2 = 0;
      for(auto jet : l1JetsSorted){
      if(reco::deltaR(jet, genJet_1)<0.1 && foundL1Jet_1 == 0 ){
	l1Pt_1  = jet.pt();
	l1Eta_1 = jet.eta();
	l1Phi_1 = jet.phi();
	l1NthJet_1 = i;
	foundL1Jet_1 = 1;
      }
      if(genPt_2 > 0 && reco::deltaR(jet, genJet_2)<0.1 && foundL1Jet_2 == 0 ){
	l1Pt_2  = jet.pt();
	l1Eta_2 = jet.eta();
	l1Phi_2 = jet.phi();
	l1NthJet_2 = i;
	foundL1Jet_2 = 1;
      }
      i++;
      }
      
      if(foundL1Jet_1>0 && foundL1Jet_2>0){
	l1DeltaEta = l1Eta_1 - l1Eta_2;
	l1DeltaPhi = l1Phi_1 - l1Phi_2;
	l1DeltaR = reco::deltaR(l1JetsSorted.at(l1NthJet_1), l1JetsSorted.at(l1NthJet_2) );
	l1Mass = (l1JetsSorted.at(l1NthJet_1).p4() + l1JetsSorted.at(l1NthJet_2).p4()).mass();
      }
      
      nGenJets = genJets->size();
      nRecoJets = goodJets.size();
      nL1Jets = l1JetsSorted.size();
      genTree->Fill();
    }
  
 
  zeroOutAllVariables();

  //fill the jet variables here!
  //include delta eta between the two jets,
  // pt, eta, phi of each jet
  //delta phi between two jets
  //invariant mass of two jets
  //total number of jets in the event
  pat::Jet recoJet_1;
  pat::Jet recoJet_2;
  if(goodJets.size()>0){

    recoPt_1  = goodJets.at(0).pt();
    recoEta_1 = goodJets.at(0).eta();
    recoPhi_1 = goodJets.at(0).phi();
    recoJet_1 = goodJets.at(0);
    
    if(goodJets.size()>1){
      recoJet_2 = goodJets.at(1);
      recoPt_2  = goodJets.at(1).pt();
      recoEta_2 = goodJets.at(1).eta();
      recoPhi_2 = goodJets.at(1).phi();
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
    for(auto jet : l1JetsSorted){
      if(reco::deltaR(jet, recoJet_1)<0.1 && foundL1Jet_1 == 0 ){
	l1Jet_1 = jet;
	l1Pt_1  = jet.pt();
	l1Eta_1 = jet.eta();
	l1Phi_1 = jet.phi();
	l1NthJet_1 = i;
	foundL1Jet_1 = 1;
      }
      if(genPt_2 > 0 && reco::deltaR(jet, recoJet_2)<0.1 && foundL1Jet_2 == 0 ){
	l1Jet_2 = jet;
	l1Pt_2  = jet.pt();
	l1Eta_2 = jet.eta();
	l1Phi_2 = jet.phi();
	l1NthJet_2 = i;
	foundL1Jet_2 = 1;
      }
      i++;
    }

    if(foundL1Jet_1>0 && foundL1Jet_2>0){
      l1DeltaEta = l1Eta_1 - l1Eta_2;
      l1DeltaPhi = l1Phi_1 - l1Phi_2;
      l1DeltaR = reco::deltaR(l1Jet_1, l1Jet_2);
      l1Mass = (l1Jet_1.p4() + l1Jet_2.p4()).mass();
    }
    
    nRecoJets = goodJets.size();
    nL1Jets = l1JetsSorted.size();
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

void Run3Ntuplizer::zeroOutAllVariables(){
  genDeltaEta=-99;   genDeltaPhi=-99;   genDeltaR=-99; genMass=-99;
  recoPt_1=-99;  recoEta_1=-99;    recoPhi_1=-99;    recoNthJet_1=-99;    
  recoPt_2=-99;     recoEta_2=-99;    recoPhi_2=-99;    recoNthJet_2=-99;    
  recoDeltaEta=-99;    recoDeltaPhi=-99;      recoDeltaR=-99;   recoMass=-99;   
  l1Pt_1=-99;       l1Eta_1=-99;      l1Phi_1=-99;      l1NthJet_1=-99; 
  l1Pt_2=-99;       l1Eta_2=-99;      l1Phi_2=-99;      l1NthJet_2=-99; 
  l1DeltaEta=-99;   l1DeltaPhi=-99;   l1DeltaR=-99;     l1Mass=-99;     
  l1Matched_1=-99;   l1Matched_2=-99;               
  genMatched_1=-99;      genMatched_2=-99;    
  nGenJets=-99;   
  nRecoJets=-99;  
  nL1Jets=-99;    

};

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
