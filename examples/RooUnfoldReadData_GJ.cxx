//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldExample.cxx 279 2011-02-11 18:23:44Z T.J.Adye $
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================


#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "TFile.h"
#include "TVectorD.h"

#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TPostScript.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TProfile.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;

//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================


//==============================================================================
// Example Unfolding
//==============================================================================
Double_t RelativePhi(Double_t mphi,Double_t vphi) {

  if (vphi < -1 * TMath::Pi())
    vphi += (2 * TMath::Pi());
  else if (vphi > TMath::Pi())
    vphi -= (2 * TMath::Pi());
  if (mphi < -1 * TMath::Pi())
    mphi += (2 * TMath::Pi());
  else if (mphi > TMath::Pi())
    mphi -= (2 * TMath::Pi());
  double dphi = mphi - vphi;
  if (dphi < -1 * TMath::Pi())
    dphi += (2 * TMath::Pi());
  else if (dphi > TMath::Pi())
    dphi -= (2 * TMath::Pi());
  return dphi; // dphi in [-Pi, Pi]                                                                                                                                                         
}


void RooUnfoldReadData_GJ(std::string cFiles2 = "",std::string tag = "", std::string date = "", int fileindex=0)
{


  //errors from the diagonal elements of the covariance matrix
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
 
 
  
   double Emin_gamma=100;
   const int bin_det_xj=6;

   const int bin_det_rg=6;


   Double_t xjmin_det=0.4;


   Double_t xjmax_det=1.5;


   Double_t rgmax_det=0.2;

  
    Double_t rgmin_det=0.;


   //***************************************************

   //det binning in Xj
   Double_t ybins[bin_det_xj];

   ybins[0]=xjmin_det;
   ybins[1]=0.6;
   ybins[2]=0.8;
   ybins[3]=1.;
   ybins[4]=1.2;
   ybins[5]=1.5;



   //det binning in Rg                                                                                                                                                 
   Double_t xbins[bin_det_rg];
   xbins[0]=-0.05;
   xbins[1]=0.;
   xbins[2]=0.05;
   xbins[3]=0.1;
   xbins[4]=0.15;
   xbins[5]=0.2;
  


   //the raw correlation
   TH2D *h2raw(0);
   h2raw=new TH2D("raw","raw",bin_det_rg-1,xbins,bin_det_xj-1,ybins);
   h2raw->Sumw2();


   //******************Fill the data input*******************//
  
  Int_t hiBin=0;  
  Int_t nref=0;
  Float_t jtpt[500];
  Float_t jteta[500];
  Float_t jtphi[500];
  Float_t refpt[500];
  Float_t refeta[500];
  Float_t refphi[500];
  Float_t jtsym[500];
  Float_t refsym[500];
  Float_t jtrg[500];
  Float_t refrg[500];
  Float_t jtdynkt[500];
  Float_t refdynkt[500];
  Float_t jtangu[500];
  Float_t weight=0;
  Float_t pthat=0;
  std::vector<float> *phoEt=0;
  std::vector<float> *phoEta=0;
  std::vector<float> *phoPhi=0; 

  std::vector<float> *elePt=0;
  std::vector<float> *eleEta=0;
  std::vector<float> *elePhi=0;
 
   std::vector<float> *phoE=0;
   std::vector<float> *pho_genMatchedIndex=0;
    std::vector<float> *phoSigmaEtaEta_2012=0;
    std::vector<float> *mcEt=0;
    std::vector<float> *mcPhi=0;
    std::vector<float> *mcPID=0;
    std::vector<float> *pfpIso4subUE=0;
    std::vector<float> *pfcIso4subUE=0;
    std::vector<float> *pfnIso4subUE=0;
    std::vector<float> *pho_ecalClusterIsoR4=0;
    std::vector<float> *pho_hcalRechitIsoR4=0;                                                                                          
    std::vector<float> *pho_trackIsoR4PtCut20=0;      


    std::vector<float> *phoE1x5_2012=0;
     std::vector<float> *phoE2x5_2012=0;
      std::vector<float> *phoE3x5_2012=0;
        std::vector<float> *phoE5x5_2012=0;
          std::vector<float> *phoE3x3_2012=0;
	  std::vector<float> *phoHoverE=0;

  Int_t nPho=0;



   ifstream infile;
   infile.open(cFiles2);
   char filename[300];
   int count=0;
    while(infile>>filename){
    TFile *input=TFile::Open(filename);
    if(!input || input->IsZombie()) { cout <<"The file could not be opened!"; delete input; continue;}
    TDirectoryFile *jetbranch;
    jetbranch=(TDirectoryFile*)input->Get("akCs2PFJetAnalyzer_substructure");
    TDirectoryFile *photonbranch=(TDirectoryFile*)input->Get("ggHiNtuplizerGED");
    TDirectoryFile *hiEvtAnalyzer=(TDirectoryFile*)input->Get("hiEvtAnalyzer");
    count=count+1; 
    cout<<"reading data file "<<count<<endl;
    if(count==100) break;
    TTree *jets=(TTree*)jetbranch->Get("t");
    TTree *evprop=(TTree*)hiEvtAnalyzer->Get("HiTree");
    TTree *gammas=(TTree*)photonbranch->Get("EventTree");
   
    jets->AddFriend(gammas);
    jets->AddFriend(evprop);
    Int_t nEv=jets->GetEntries();
 
  jets->SetBranchAddress("jtpt", &jtpt);
  jets->SetBranchAddress("jteta", &jteta);
  jets->SetBranchAddress("jtphi", &jtphi);
  jets->SetBranchAddress("nref",&nref);  
  jets->SetBranchAddress("hiBin",&hiBin);
  jets->SetBranchAddress("jtsym", &jtsym);
  jets->SetBranchAddress("jtrg", &jtrg);
  jets->SetBranchAddress("jtdynkt", &jtdynkt);
  jets->SetBranchAddress("jtangu", &jtangu); 

   jets->SetBranchAddress("nPho", &nPho);
   jets->SetBranchAddress("phoEt",&phoEt);
   jets->SetBranchAddress("phoEta",&phoEta);
   jets->SetBranchAddress("phoPhi",&phoPhi);
   
   jets->SetBranchAddress("elePt",&elePt);
   jets->SetBranchAddress("eleEta",&eleEta);
   jets->SetBranchAddress("elePhi",&elePhi);
   jets->SetBranchAddress("phoSigmaEtaEta_2012",&phoSigmaEtaEta_2012);

       //isolation variables
       jets->SetBranchAddress("pfcIso4subUE",&pfcIso4subUE);
       jets->SetBranchAddress("pfnIso4subUE",&pfnIso4subUE);
       jets->SetBranchAddress("pfpIso4subUE",&pfpIso4subUE);
       //iso variables for 2015 data 
       jets->SetBranchAddress("pho_ecalClusterIsoR4",&pho_ecalClusterIsoR4);
       jets->SetBranchAddress("pho_hcalRechitIsoR4",&pho_hcalRechitIsoR4);
       jets->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20);
       //cluster shape
       jets->SetBranchAddress("phoHoverE",&phoHoverE);
       
   for(int iEntry=0; iEntry< nEv; iEntry++){
   jets->GetEntry(iEntry);
   if(hiBin>60) continue; 
   
   int ing=-1;
   double ptmax=0;
   int genin=-1;
   for(Int_t m=0;m<phoEt->size();m++){
     if(TMath::Abs(phoEta->at(m))>1.44) continue;
      if(phoEt->at(m)>ptmax){
	ptmax=phoEt->at(m);
        ing=m;
       
      }
   }

   if(ing==-1) continue;
   if(phoEt->at(ing)<Emin_gamma) continue;
  

   double isolation=pfcIso4subUE->at(ing)+pfnIso4subUE->at(ing)+pfpIso4subUE->at(ing);
   if(isolation>1) continue; 
    int flagEle=0;
    for(Int_t n=0;n<elePt->size();n++){
      if(elePt->at(n)<10) continue;
      double dee=TMath::Abs(eleEta->at(n)-phoEta->at(ing));
      double dephi=RelativePhi(elePhi->at(n),phoPhi->at(ing));
      if(dee<0.02 && dephi<0.15) flagEle=1;}
    
   if(flagEle==1) continue;
   if(phoHoverE->at(ing)>=0.1) continue;
  
   
   if(phoSigmaEtaEta_2012->at(ing)>0.01) continue;
  
     int leadjet=-1;
     double ptlead=-200;
     for(Int_t k=0;k<nref;k++){
     

       double deltaphi=TMath::Abs(RelativePhi(phoPhi->at(ing),jtphi[k]));

       if(deltaphi>TMath::Pi()-0.6){

	 if(jtpt[k]>ptlead) {ptlead=jtpt[k];
	   leadjet=k;}
       }}
    if(leadjet==-1) continue;
    double xj=jtpt[leadjet]/phoEt->at(ing);
    double rg=jtrg[leadjet];
    if(xj>xjmax_det || xj<xjmin_det) continue; 
    
    if(rg>rgmax_det || rg<rgmin_det) continue;
   
    if(rg==0) rg=-0.025;
    h2raw->Fill(rg,xj);}
   delete jets;
   delete evprop;
   delete gammas;
   delete jetbranch;
   delete photonbranch;
   delete hiEvtAnalyzer;
    input->Close();
    delete input;
    }

   
 
    TFile *fout=new TFile (Form("/data_CMS/cms/cbaldene/UnfoldingTest/RooUnfold/test/total_%s_%i_%s.root",tag.c_str(),fileindex,date.c_str()),"RECREATE");
 fout->cd();
 


  h2raw->SetName("raw");
  h2raw->Write();
 
	  


	     
}

#ifndef __CINT__
int main () { RooUnfoldReadData_GJ(); return 0; }  // Main program when run stand-alone
#endif
