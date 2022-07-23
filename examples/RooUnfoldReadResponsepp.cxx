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

                                              
TH2D* CorrelationHist(const TMatrixD& cov, const char* name, const char* title, Int_t na, Int_t nb, Int_t kbin)
{

  TH2D* h=new TH2D (name, title,na,0,na,na,0,na);
  for(int l=0;l<na;l++){
    for(int n=0;n<na;n++){
      int index1=l+na*kbin;
      int index2=n+na*kbin;
      double vv=cov(index1,index1)*cov(index2,index2);
      if(vv>0.0) h->SetBinContent(l+1,n+1,cov(index1,index2)/sqrt(vv));
    }}
  return h;

}



void RooUnfoldReadResponsepp(std::string cFiles2="",std::string tag = "", std::string date = "", int fileindex=0)
{

#ifdef __CINT__
  gSystem->Load("libRooUnfold.so");
#endif
  //errors from the diagonal elements of the covariance matrix
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
 
  
   



    double Emin_gamma=100;
    double pHECut=0.0233;
    double pIsoCut=-0.479;
    double pSigmaCut=0.01;
    double purity=0.938;

   

   const int bin_det_xj=3;
   const int bin_true_xj=4;
   const int bin_det_rg=5;
   const int bin_true_rg=6;

   Double_t xjmin_det=0.;
   Double_t xjmin_true=0;

   Double_t xjmax_det=1.6;
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
       
        
         ybins[2]=1.6;

     
       //true binning in xj
         Double_t ybins_true[bin_true_xj];
         ybins_true[0]=xjmin_true;
         ybins_true[1]=0.6;
        
         ybins_true[2]=1.6;
         ybins_true[3]=2;
        
        

        //det binning in Rg
         Double_t xbins[bin_det_rg];
	 xbins[0]=-0.05;
	 xbins[1]=0.;
         xbins[2]=0.05;
 	 xbins[3]=0.1;
         xbins[4]=0.2;
        	
        
       
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
    Float_t refangu[500];
  
    Float_t pthat=0;
    Bool_t eleRej=0;
    Float_t SumCalIso=0;
    Float_t phoHoverE=0;
    Float_t phoEtCorrected=0;
    Float_t phoSigmaIEtaIEta_2012=0;
    Float_t phoPhi=0;
    Float_t phoEta=0;
    Float_t phoEt=0;
    Float_t weight=0;
    Float_t weight_pthat=0;
    Float_t weight_cent=0;
    Int_t mcPID=0;
    Float_t mcEt=0;
    
 
   //******************Fill the MCresponse*******************//
   int countmc=0;
   ifstream infile2;
   infile2.open(cFiles2);
   char filename2[300];
   

   ifstream infile;
   infile.open(cFiles2);
   char filename[300];
   int count=0;

   while(infile>>filename){
     TFile *input=TFile::Open(filename);
     if(!input || input->IsZombie()) { cout <<"The file could not be opened!"; delete input; continue;}

     TTree *jets=(TTree*)input->Get("jet_tree");
     Int_t nEv=jets->GetEntries();

     jets->SetBranchAddress("jtpt", &jtpt);
     jets->SetBranchAddress("jteta", &jteta);
     jets->SetBranchAddress("jtphi", &jtphi);
    
     
     jets->SetBranchAddress("jtsym", &jtsym);
     jets->SetBranchAddress("jtrg", &jtrg);
     jets->SetBranchAddress("jtdynkt", &jtdynkt);
     jets->SetBranchAddress("jtangu", &jtangu);

     jets->SetBranchAddress("refsym", &refsym);
     jets->SetBranchAddress("refrg", &refrg);
     jets->SetBranchAddress("refdynkt", &refdynkt);
     jets->SetBranchAddress("refangu", &refangu);
     jets->SetBranchAddress("refpt", &refpt); 
     jets->SetBranchAddress("weight",&weight);
     jets->SetBranchAddress("weight_pthat",&weight_pthat);
     

     jets->SetBranchAddress("phoEt",&phoEt);
     jets->SetBranchAddress("phoEta",&phoEta);
     jets->SetBranchAddress("phoPhi",&phoPhi);
     jets->SetBranchAddress("phoEtCorrected",&phoEtCorrected);
     jets->SetBranchAddress("mcEt",&mcEt);
     jets->SetBranchAddress("phoSigmaIEtaIEta_2012",&phoSigmaIEtaIEta_2012);
     jets->SetBranchAddress("mcPID",&mcPID);      
     jets->SetBranchAddress("nref",&nref);
     jets->SetBranchAddress("eleRej",&eleRej);
     jets->SetBranchAddress("SumCalIso",&SumCalIso);
     jets->SetBranchAddress("phoHoverE",&phoHoverE);


   for(int iEntry=0; iEntry< nEv;iEntry++){
   jets->GetEntry(iEntry);
  
   double scale=weight_pthat*weight; 
   //cout<<scale<<" "<<weight<<" "<<weight_cent<<endl;
   if(phoHoverE>pHECut) continue;
   if(phoSigmaIEtaIEta_2012>pSigmaCut)continue;
   if(eleRej==1) continue;
   if(mcPID!=22) continue;
   if(phoEtCorrected<Emin_gamma) continue;
   if(SumCalIso>pIsoCut) continue; 
    int leadjet=-1;
    double ptlead=-200;
    // cout<<nref<<"number of jets"<<endl;
      for(Int_t k=0;k<nref;k++){
       double deltaphi=TMath::Abs(RelativePhi(phoPhi,jtphi[k]));
       if(deltaphi>TMath::Pi()-1.5){
	 if(jtpt[k]>ptlead) {ptlead=jtpt[k];
	   leadjet=k;}
       }}
    if(leadjet==-1) continue;
    if(jtpt[leadjet]<40) continue;     
    double xj_det=jtpt[leadjet]/phoEtCorrected;
    double rg_det=jtrg[leadjet];
    double xj_true=refpt[leadjet]/mcEt;
    double rg_true=refrg[leadjet];

    
    

    if(xj_true>xjmax_true || xj_true<xjmin_true) continue;
      if(rg_true>rgmax_true || rg_true<rgmin_true) continue;
      if(rg_true==0) rg_true=-0.025; 
      h2fulleff->Fill(rg_true,xj_true,scale);
    
      if(xj_det>xjmax_det || xj_det<xjmin_det) continue;
       if(rg_det>rgmax_det || rg_det<rgmin_det) continue; 
       if(rg_det == 0) rg_det=-0.025;
    
       h2smeared->Fill(rg_det,xj_det,scale);
       h2true->Fill(rg_true,xj_true,scale);
       response.Fill(rg_det,xj_det,rg_true,xj_true,scale);}

		      
    delete jets;
 
 
 
    input->Close();
    delete input;
   }

   

   //********************************************************

   //i,j,k,l are indexes for rg^det,xj^det,rg^true,xj^true                                                                                                  
   //number of bins are (n1,n2,n3,n4) for (rg^det,xj^det,rg^true,xj^true)                                                                                                                        
   int n1=4;
   int n2=2;
   int n3=5;
   int n4=3;
   double sum=0;
   for(Int_t k=0;k<n3;k++){
     for(Int_t l=0;l<n4;l++){
       sum=0;
       for(Int_t i=0;i<n1;i++){
	 for(Int_t j=0;j<n2;j++){
	   int index_det=i+n1*j;
	   int index_true=k+n3*l;
     	   sum=sum+response(index_det,index_true);
	 }}

       cout<<"true indexes"<<k<<" "<<l<<" "<<"Pkl is prob.of true k,l causes "<<sum<<endl;

     }}

   //the roounfold object is a 2 D matrix response(detindex,trueindex)


    
   ////////////******************************************
   


 
  TFile *fout=new TFile (Form("/data_CMS/cms/lcunquei/UnfoldingSetup/OutputJuly2022/ResponseIsRead_pp_%s_%i_%s.root",tag.c_str(),fileindex,date.c_str()),"RECREATE");
 fout->cd();

  h2smeared->SetName("smeared");
  h2smeared->Write();
  h2fulleff->Write();
  h2true->Write();  
  response.Write();
	  



	     
}

#ifndef __CINT__
int main () { RooUnfoldReadResponsepp(); return 0; }  // Main program when run stand-alone
#endif
