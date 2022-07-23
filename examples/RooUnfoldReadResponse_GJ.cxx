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
//#include "RooUnfoldTestHarness2D.h"
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


void RooUnfoldReadResponse_GJ(std::string cFiles2="",std::string tag = "", std::string date = "", int fileindex=0)
{

#ifdef __CINT__
  gSystem->Load("libRooUnfold.so");
#endif
  //errors from the diagonal elements of the covariance matrix
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
 
  
   double Emin_gamma=100;
   const int bin_det_xj=6;
   const int bin_true_xj=7;
   const int bin_det_rg=6;
   const int bin_true_rg=6;

   Double_t xjmin_det=0.4;
   Double_t xjmin_true=0;

   Double_t xjmax_det=1.5;
   Double_t xjmax_true=2;

   Double_t rgmax_det=0.2;
   Double_t rgmax_true=0.3;
  
    Double_t rgmin_det=0.;
   Double_t rgmin_true=0.;

   //***************************************************

   //det binning in xj
     Double_t ybins[bin_det_xj];
      
         ybins[0]=xjmin_det;
         ybins[1]=0.6;
	 ybins[2]=0.8;
 	 ybins[3]=1.;
         ybins[4]=1.2;
         ybins[5]=1.5;
       
     
       //true binning in xj
        Double_t ybins_true[bin_true_xj];
         ybins_true[0]=xjmin_true;
         ybins_true[1]=0.4;
	 ybins_true[2]=0.6;
 	
         ybins_true[3]=0.8;
         ybins_true[4]=1.;
         ybins_true[5]=1.5;
         ybins_true[6]=2.;

        
        

        //det binning in Rg
         Double_t xbins[bin_det_rg];
	 xbins[0]=-0.05;
	 xbins[1]=0.;
	 xbins[2]=0.05;
         xbins[3]=0.1;
	 xbins[4]=0.15;
 	 xbins[5]=0.2;
	
        
       
        //true binning in Rg
         Double_t xbins_true[bin_true_rg];
	 xbins_true[0]=-0.05;
	 xbins_true[1]=0.;
	 xbins_true[2]=0.05;
         xbins_true[3]=0.1;
	 xbins_true[4]=0.2;
 	 xbins_true[5]=0.3;
           
	 
  
   //MC detector correlation
   TH2D *h2smeared(0);
   h2smeared=new TH2D("smeared","smeared",bin_det_rg-1,xbins,bin_det_xj-1,ybins);
   //true correlations with det-level cuts
   TH2D *h2true(0);
    h2true=new TH2D("true","true",bin_true_rg-1,xbins_true,bin_true_xj-1,ybins_true);
     
    //fully-efficient true correlation (no det-level cuts)
    TH2D *h2fulleff(0);
    h2fulleff=new TH2D("truef","truef",bin_true_rg-1,xbins_true,bin_true_xj-1,ybins_true);
     TH2D *effnum=(TH2D*)h2fulleff->Clone("effnum");
     TH2D *effdenom=(TH2D*)h2fulleff->Clone("effdenom");

 
   effnum->Sumw2();
   effdenom->Sumw2();
   h2smeared->Sumw2();
   h2true->Sumw2();
  
   h2fulleff->Sumw2();
   //setup the binnings of the 4D response
    RooUnfoldResponse response;
    response.Setup(h2smeared,h2true);

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
  Float_t weight_pthat=0;
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



 
   //******************Fill the MCresponse*******************//
   int countmc=0;
   ifstream infile2;
   infile2.open(cFiles2);
   char filename2[300];
   
    while(infile2>>filename2){
    TFile *input2=TFile::Open(filename2);
    countmc=countmc+1;
    cout<<"reading mc file "<<countmc<<endl; 
    if(countmc==10) break;
  
    if(!input2 || input2->IsZombie()) { cout <<"The file could not be opened!"; delete input2; continue;}
    TDirectoryFile *jetbranchmc;
    jetbranchmc=(TDirectoryFile*)input2->Get("akCs2PFJetAnalyzer_substructure");
    TDirectoryFile *photonbranchmc=(TDirectoryFile*)input2->Get("ggHiNtuplizerGED");
    TDirectoryFile *hiEvtAnalyzermc=(TDirectoryFile*)input2->Get("hiEvtAnalyzer");
   
    TTree *jets_mc=(TTree*)jetbranchmc->Get("t");
    TTree *evprop_mc=(TTree*)hiEvtAnalyzermc->Get("HiTree");
    TTree *gammas_mc=(TTree*)photonbranchmc->Get("EventTree");
   
    jets_mc->AddFriend(gammas_mc);
    jets_mc->AddFriend(evprop_mc);
    Int_t nEvmc=jets_mc->GetEntries();
 
  jets_mc->SetBranchAddress("jtpt", &jtpt);
  jets_mc->SetBranchAddress("refpt", &refpt);
  jets_mc->SetBranchAddress("jteta", &jteta);
  jets_mc->SetBranchAddress("jtphi", &jtphi);
  jets_mc->SetBranchAddress("nref",&nref);  
  jets_mc->SetBranchAddress("hiBin",&hiBin);
  jets_mc->SetBranchAddress("pthat",&pthat); 
   jets_mc->SetBranchAddress("weight_pthat",&weight_pthat); 
  jets_mc->SetBranchAddress("jtsym", &jtsym);
  jets_mc->SetBranchAddress("refsym", &refsym);  
  jets_mc->SetBranchAddress("jtrg", &jtrg);
  jets_mc->SetBranchAddress("refrg", &refrg);  
  jets_mc->SetBranchAddress("jtdynkt", &jtdynkt);
  jets_mc->SetBranchAddress("refdynkt", &refdynkt);
  jets_mc->SetBranchAddress("nPho", &nPho);
  jets_mc->SetBranchAddress("phoEt",&phoEt);
  jets_mc->SetBranchAddress("phoEta",&phoEta);
  jets_mc->SetBranchAddress("phoPhi",&phoPhi);
  jets_mc->SetBranchAddress("pho_genMatchedIndex",&pho_genMatchedIndex);
  jets_mc->SetBranchAddress("mcEt",&mcEt);
  jets_mc->SetBranchAddress("mcPhi",&mcPhi);
  jets_mc->SetBranchAddress("mcPID",&mcPID); 

   for(int iEntry=0; iEntry< nEvmc; iEntry++){
   jets_mc->GetEntry(iEntry);
   if(hiBin>60) continue; 
   int ing=-1;
   int genin=-1;
   double ptmax=0;
   double scale=weight_pthat; 

   

  for(Int_t m=0;m<phoEt->size();m++){
     if(TMath::Abs(phoEta->at(m))>1.44) continue;
       if(phoEt->at(m)>ptmax){
	ptmax=phoEt->at(m);
        ing=m;
        genin=pho_genMatchedIndex->at(m);
}}

   if(ing==-1 || genin==-1) continue;
    if(mcPID->at(genin)!=22) continue;
   if(phoEt->at(ing)<Emin_gamma) continue;
    
    int leadjet=-1;
    double ptlead=-200;
     for(Int_t k=0;k<nref;k++){


       double deltaphi=TMath::Abs(RelativePhi(phoPhi->at(ing),jtphi[k]));

       if(deltaphi>TMath::Pi()-0.6){

	 if(jtpt[k]>ptlead) {ptlead=jtpt[k];
	   leadjet=k;}
       }}
    if(leadjet==-1) continue;
    double xj_det=jtpt[leadjet]/phoEt->at(ing);
    double rg_det=jtrg[leadjet];
    double xj_true=refpt[leadjet]/mcEt->at(genin);
    double rg_true=refrg[leadjet];

    // cout<<refpt[leadjet]<<" "<<jtpt[leadjet]<<" "<<phoEt->at(ing)<<" "<<mcEt->at(genin)<<endl;
    //cout<<xj_det<<" "<<xj_true<<endl;

    if(xj_true>xjmax_true || xj_true<xjmin_true) continue;
      if(rg_true>rgmax_true || rg_true<rgmin_true) continue;
      if(rg_true==0) rg_true=-0.025; 
      h2fulleff->Fill(rg_true,xj_true,scale);
    
      if(xj_det>xjmax_det || xj_det<xjmin_det) continue;
       if(rg_det>rgmax_det || rg_det<rgmin_det) continue; 
       if(rg_det == 0) rg_det=-0.025;
       //cout<<"passes "<<xj_det<<" "<<xj_true<<endl;
       h2smeared->Fill(rg_det,xj_det,scale);
       h2true->Fill(rg_true,xj_true,scale);
       response.Fill(rg_det,xj_det,rg_true,xj_true,scale);}

		      
   delete jets_mc;
   delete evprop_mc;
   delete gammas_mc;
   delete photonbranchmc;
   delete jetbranchmc;
   delete hiEvtAnalyzermc;
    input2->Close();
    delete input2;
    }

   



 


 
 TFile *fout=new TFile (Form("/data_CMS/cms/cbaldene/UnfoldingTest/RooUnfold/test/ResponseIsRead_total_%s_%i_%s.root",tag.c_str(),fileindex,date.c_str()),"RECREATE");
 fout->cd();

  h2smeared->SetName("smeared");
  h2smeared->Write();
  h2fulleff->Write();
  h2true->Write();  
  response.Write();
	  


	     
}

#ifndef __CINT__
int main () { RooUnfoldReadResponse_GJ(); return 0; }  // Main program when run stand-alone
#endif
