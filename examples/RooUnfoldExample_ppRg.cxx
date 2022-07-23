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

void RooUnfoldExample_ppRg(Int_t flag=0)
{
  //flag=0 is default
  //flag=1 is lower truncation
  //flag=2 is higher truncation
  //flag=3 is change in raw binning
  //flag=4 is change in tracking efficiency
  //flag=9 unagged bin somwhere else
  //flag=10 is change in prior
  //flag=20 is dummy unfolding: raw input is det-level MC, just consistency check of binning etc  

  

  //errors from the diagonal elements of the covariance matrix
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
   //Get the tree for data
  TString fnamedata;
  fnamedata="/grid_mnt/data__data.polcms/cms/lcunquei/AliceRgMCandData/PP2021Rg/ppdatazcut04.root";
  TFile *inputdata;
  inputdata=TFile::Open(fnamedata);
  TTree *data=(TTree*)inputdata->Get("JetSubstructure_Jet_AKTChargedR040_tracks_pT0150_E_scheme_TCRawTree_Data_NoSub_Incl"); 
  cout<<data<<endl;

  //these are the MC trees, one for the default case, one with degraded tracking efficiency for the systematics
   TString cFiles2;
   cFiles2="filesmcpp_alice.txt";
   if(flag==4)cFiles2="filesmcpp_loweff_alice.txt";

   const int bin_det_pt=6;
   const int bin_true_pt=9;
   const int bin_det_rg=8;
   const int bin_true_rg=9;

   Double_t ptmin_det=20;
   Double_t ptmin_true=0;

   Double_t ptmax_det=80;
   Double_t ptmax_true=160;

   Double_t rgmax_det=0.35;
   Double_t rgmax_true=0.6;
  
   Double_t zcut=0.1;
    if(flag==1) ptmin_det=17;
    if(flag==2) ptmin_det=23;
    if(flag==11) ptmin_det=15;


   //***************************************************

   //this is the binning in jet pT
     Double_t xbins[bin_det_pt];
      
         xbins[0]=ptmin_det;
         xbins[1]=30;
	 xbins[2]=40;
 	 xbins[3]=50;
         xbins[4]=60;
         xbins[5]=80;
        
       
      	 //xbinsb is the det-level  binning in Rg
         //I create an extra negative bin where I will put the untagged jets, ie those jets that do not have a SD prong
         Double_t xbinsb[bin_det_rg];
	 xbinsb[0]=-0.05;
	 xbinsb[1]=0.;
	 xbinsb[2]=0.02;
         xbinsb[3]=0.04;
	 xbinsb[4]=0.06;
 	 xbinsb[5]=0.1;
	 xbinsb[6]=0.2;
         xbinsb[7]=0.35;
         //the untagged bin is put somwhere else
        if(flag==9){
         xbinsb[0]=0.;
         xbinsb[1]=0.02;
         xbinsb[2]=0.04;
         xbinsb[3]=0.06;
         xbinsb[4]=0.1;
         xbinsb[5]=0.2;
         xbinsb[6]=0.35;
         xbinsb[7]=0.55;
         }



	//a slight change of binning is considered
         if(flag==3){	
         xbinsb[0]=-0.05;
         xbinsb[1]=0.;
         xbinsb[2]=0.018;
         xbinsb[3]=0.038;
         xbinsb[4]=0.058;
         xbinsb[5]=0.08;
         xbinsb[6]=0.2;
         xbinsb[7]=0.35;
         }



         //xbinsc is the true-level binning in Rg, the range is expanded
         //also the true pT expands from 0 to 160 GeV
	 Double_t xbinsc[bin_true_rg];
         xbinsc[0]=-0.05;
	 xbinsc[1]=0.;
	 xbinsc[2]=0.02;
         xbinsc[3]=0.04;
	 xbinsc[4]=0.06;
 	 xbinsc[5]=0.1;
	 xbinsc[6]=0.2;
         xbinsc[7]=0.35;
	 xbinsc[8]=0.6;

 
	 
   //the raw correlation
   TH2D *h2raw(0);
   h2raw=new TH2D("raw","raw",bin_det_rg-1,xbinsb,bin_det_pt-1,xbins);
   //detector measure level
   TH2D *h2smeared(0);
   h2smeared=new TH2D("smeared","smeared",bin_det_rg-1,xbinsb,bin_det_pt-1,xbins);

   //true correlations 
    TH2D *h2true(0);
    h2true=new TH2D("true","true",bin_true_rg-1,xbinsc,bin_true_pt-1,0,160);
    //fully-efficient true correlation (no det-level cuts)
    TH2D *h2fulleff(0);
    h2fulleff=new TH2D("truef","truef",bin_true_rg-1,xbinsc,bin_true_pt-1,0,160);
   

  
  TH2D *hcovariance(0);
  hcovariance=new TH2D("covariance","covariance",10,0.,1.,10,0,1.);
  TH2D *effnum=(TH2D*)h2fulleff->Clone("effnum");
  TH2D *effdenom=(TH2D*)h2fulleff->Clone("effdenom");
 
   effnum->Sumw2();
   effdenom->Sumw2();
   h2smeared->Sumw2();
   h2true->Sumw2();
   h2raw->Sumw2();
   h2fulleff->Sumw2();
  
   Float_t ptJet,ptJetMatch,zg,EPlane,zgMatch,rg,rgMatch,ktg,ktgMatch,ng,ngMatch,LeadingTrackPt,LeadingTrackPtDet,LeadingTrackPtMatch,ngDet,subjet1,subjet2;

   //******************Fill the data input*******************//
  Int_t nEv=0;; 
  
  nEv=data->GetEntries(); 
  cout<<"entries"<<nEv<<endl;
  data->SetBranchAddress("ptJet", &ptJet); 
  data->SetBranchAddress("zg", &zg); 
  data->SetBranchAddress("rg", &rg); 
  data->SetBranchAddress("ktg", &ktg);
  data->SetBranchAddress("ng", &ng);  
 
   for(int iEntry=0; iEntry< nEv; iEntry++){
   data->GetEntry(iEntry); 
  
   if(ptJet>ptmax_det || ptJet<ptmin_det) continue;
   if(rg>rgmax_det) continue;
   if(flag!=9) if(zg<zcut) rg=-0.025;
   if(flag==9) if(zg<zcut) rg=0.45;
   h2raw->Fill(rg,ptJet);

   }
 
 

   RooUnfoldResponse response;
  
   response.Setup(h2smeared,h2true);
  

 

   //************"Fill the response"**********************//

  ifstream infile3;
  infile3.open(cFiles2.Data());
  char filename3[300];

    while(infile3>>filename3){
    cout<<"one file at a time"<<endl;
    int pthardbin=0;
       double scalefactor=1;
 
      TFile *input2=TFile::Open(filename3);
       if (!input2 || input2->IsZombie()) { cout <<"The file could not be opened!"; delete input2; continue;}

      TTree *mc=(TTree*)input2->Get("JetSubstructure_Jet_AKTChargedR040_tracks_pT0150_E_scheme_TCRawTree_PythiaDef_NoSub_Incl"); 
       if(!mc) continue;
      Int_t nEv=mc->GetEntries(); 
      cout<<nEv<<endl;   
      mc->SetBranchAddress("ptJet", &ptJet); 
      mc->SetBranchAddress("ptJetMatch", &ptJetMatch);
      mc->SetBranchAddress("zg", &zg);
      mc->SetBranchAddress("zgMatch", &zgMatch);
      mc->SetBranchAddress("rg", &rg); 
      mc->SetBranchAddress("rgMatch", &rgMatch);
      mc->SetBranchAddress("ktg", &ktg); 
      mc->SetBranchAddress("ktgMatch", &ktgMatch);
      mc->SetBranchAddress("ng", &ng);
      mc->SetBranchAddress("ngMatch", &ngMatch);
      
      Int_t countm=0;
      Double_t reweight=1;
      Double_t myw=1;
      for(int iEntry=0; iEntry< nEv; iEntry++){
       mc->GetEntry(iEntry);
    
       
      
       
       if(ptJetMatch>ptmax_true ) continue;
       if(rgMatch>rgmax_true) continue;
       if(zgMatch<zcut) rgMatch=-0.025;
        h2fulleff->Fill(rgMatch,ptJetMatch,scalefactor);  
       
  

       if(ptJet>ptmax_det || ptJet<ptmin_det) continue;
          if(rg>rgmax_det) continue;
          if(flag!=9) if(zg<zcut) rg=-0.025;
          if(flag==9) if(zg<zcut) rg=0.45; 

        h2smeared->Fill(rg,ptJet,scalefactor);
        h2true->Fill(rgMatch,ptJetMatch,scalefactor);

       response.Fill(rg,ptJet,rgMatch,ptJetMatch,scalefactor);


      }
    
 
  
    delete mc;
  
    input2->Close();
  
    delete input2;


    }

    
 

     
       TH1F *htrueptd=(TH1F*) h2fulleff->ProjectionX("trueptd",1,-1);
       TH1F *htruept=(TH1F*) h2fulleff->ProjectionY( "truept",1,-1); 
 
     


  
 TH2D* hfold=(TH2D*)h2raw->Clone("hfold");
 hfold->Sumw2();

 //////////Fill the kinematic efficiencies////////////////////////////////////
 TH1D * effok=(TH1D *)h2true->ProjectionX("effok",2,2);
 TH1D * effok1=(TH1D *)h2fulleff->ProjectionX("effok2",2,2);
 effok->Divide(effok1);
 effok->SetName("correff20-40");
 
  TH1D * effok3=(TH1D *)h2true->ProjectionX("effok3",3,3);
  TH1D * effok4=(TH1D *)h2fulleff->ProjectionX("effok4",3,3);
  effok3->Divide(effok4);
  effok3->SetName("correff40-60"); 

  TH1D * effok5=(TH1D *)h2true->ProjectionX("effok5",5,6);
  TH1D * effok6=(TH1D *)h2fulleff->ProjectionX("effok6",5,6);
 effok5->Divide(effok6);
 effok5->SetName("correff80-120"); 

   TH1D * effok7=(TH1D *)h2true->ProjectionX("effok7",4,4);
  TH1D * effok8=(TH1D *)h2fulleff->ProjectionX("effok8",4,4);
 effok7->Divide(effok8);
 effok7->SetName("correff60-80"); 
 
    TH1D * effok9=(TH1D *)h2true->ProjectionX("effok7",6,6);
  TH1D * effok10=(TH1D *)h2fulleff->ProjectionX("effok8",6,6);
 effok9->Divide(effok10);
 effok9->SetName("correff100-120"); 

 
  TFile *fout=new TFile (Form("/data_CMS/cms/lcunquei/UnfoldingSetup/OutputNov2021pp/UnfoldRgzg0.4_total_%d.root",flag),"RECREATE");
  fout->cd();
 effok->Write(); 
 effok3->Write();
 effok5->Write();
 effok7->Write();
 effok9->Write();



  h2raw->SetName("raw");
  h2raw->Write();
  h2smeared->SetName("smeared");
  h2smeared->Write();
  h2fulleff->Write();
  htrueptd->Write();
  TH2D *h2input;
  if(flag!=20) h2input=(TH2D*)h2raw->Clone("h2input");
  if(flag==20) h2input=(TH2D*) h2smeared->Clone("h2input");
  h2input->Write(); 

     for(int jar=1;jar<16;jar++){
      Int_t iter=jar;
      cout<<"iteration"<<iter<<endl;
      cout<<"==============Unfold h1====================="<<endl;

      
      RooUnfoldBayes   unfold(&response, h2input, iter);    // OR
      TH2D* hunf= (TH2D*) unfold.Hreco(errorTreatment);
      //FOLD BACK
      TH2D* hfold= (TH2D*)response.ApplyToTruth(hunf,"");
     


     

  
          TH2D *htempUnf=(TH2D*)hunf->Clone("htempUnf");          
	  htempUnf->SetName(Form("Bayesian_Unfoldediter%d.root",iter));
           
	    TH2D *htempFold=(TH2D*)hfold->Clone("htempFold");          
	  htempFold->SetName(Form("Bayesian_Foldediter%d.root",iter));        


           






     
      	  htempUnf->Write();
	  htempFold->Write();
	  


	  
     


     
     }}

#ifndef __CINT__
int main () { RooUnfoldExample_ppRg(); return 0; }  // Main program when run stand-alone
#endif
