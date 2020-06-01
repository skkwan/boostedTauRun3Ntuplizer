#define rateplot_cxx
#include "rateplot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

bool compareByPt (TLorentzVector i, TLorentzVector j) { return(i.Pt()>j.Pt()); };

int main(int argc, char *argv[])
{
  
  if(argc > 1)
    { 
      rateplot t(argv[1], argv[2], argv[3], argv[4]);
      t.Loop(argv[3], argv[4]);
    }
  return 0;
}

using namespace std;

void rateplot::Loop(const char* pattern, const char* l1pt)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   std::stringstream sstrm1(pattern);
   std::string PATTERN;
   sstrm1 >> PATTERN;
   bool recoetacut;
   std::stringstream sstrm2(l1pt);
   double l1ptcut;
   sstrm2 >> l1ptcut;
   double SF = 1.2;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      vector<TLorentzVector> l1JetsSorted;
      for( vector<TLorentzVector>::const_iterator l1Jet = l1Jets->begin(); l1Jet != l1Jets->end(); l1Jet++ ){
         l1JetsSorted.push_back(*l1Jet);
      }
      if(l1JetsSorted.size() > 1){  std::sort(l1JetsSorted.begin(),l1JetsSorted.end(),compareByPt);}
      int m[l1JetsSorted.size()];
      for (size_t i = 0; i < l1JetsSorted.size(); i++){
         for (size_t j = 0; j < l1Jets->size(); j++){
            if(l1JetsSorted.at(i).Pt() == l1Jets->at(j).Pt()) m[i] = j;
         }
      }
      if(PATTERN == "none" && l1JetsSorted.size() > 0) {
         l1jetpt->Fill(l1JetsSorted.at(0).Pt()*SF);
         l1jeteta->Fill(l1JetsSorted.at(0).Eta());
         l1jetphi->Fill(l1JetsSorted.at(0).Phi());
      }
      if(PATTERN == "condensed" && l1JetsSorted.size() > 0) {
         for(size_t i = 0; i < l1JetsSorted.size(); i++){
            if(mEtaBits->at(m[i]) == "000000001010" || mPhiBits->at(m[i]) == "000000001010"){
               l1jetpt->Fill(l1JetsSorted.at(i).Pt()*SF);
               l1jeteta->Fill(l1JetsSorted.at(i).Eta());
               l1jetphi->Fill(l1JetsSorted.at(i).Phi());
               break;
            }
         }
      }
      if(PATTERN == "condensed12" && l1JetsSorted.size() > 0) {
         for(size_t i = 0; i < l1JetsSorted.size(); i++){
            if(mEtaBits12->at(m[i]) == "000000001010" || mPhiBits12->at(m[i]) == "000000001010"){
               l1jetpt->Fill(l1JetsSorted.at(i).Pt()*SF);
               l1jeteta->Fill(l1JetsSorted.at(i).Eta());
               l1jetphi->Fill(l1JetsSorted.at(i).Phi());
               break;
            }
         }
      }
      l1JetsSorted.clear();
   }
}

void rateplot::BookHistos(const char* file2){
   fileName = new TFile(file2, "RECREATE");
   fileName->cd();
   char name[100];
   Float_t binb[] = { -5., -4.5, -4., -3.5, -3., -2.5, -2.25, -2., -1.75, -1.5, -1.25, -1., -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 3., 3.5, 4., 4.5, 5.};
   Int_t  binnum = 30;
   
   sprintf(name, "l1jetpt");
   l1jetpt = new TH1F (name,"l1jetpt", 40, 0, 500);
   l1jetpt->SetTitle("Leading L1 jet");
   l1jetpt->GetXaxis()->SetTitle("p_{T} [GeV]");
   
   sprintf(name, "l1jeteta");
   l1jeteta = new TH1F (name,"l1jeteta", binnum, binb);
   l1jeteta->SetTitle("Leading L1 jet");
   l1jeteta->GetXaxis()->SetTitle("#eta");

   sprintf(name, "l1jetphi");
   l1jetphi = new TH1F (name,"l1jetphi", 20, -M_PI, M_PI);
   l1jetphi->SetTitle("Leading L1 jet");
   l1jetphi->GetXaxis()->SetTitle("#phi");
}

