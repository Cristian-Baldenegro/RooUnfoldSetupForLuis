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
#include "sstream"
//#include "RooUnfoldTestHarness2D.h"
#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;


void RooUnfoldExample_GJ(std::string file_mc="/data_CMS/cms/lcunquei/UnfoldingSetup/OutputNov2021pp/mergedMC_split03.root", std::string file_data="/data_CMS/cms/lcunquei/UnfoldingSetup/OutputNov2021pp/mergedData.root", std::string date = "split03",int flag=21)
{



  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
  TFile* f_mc=new TFile(file_mc.c_str());
  TFile* f_data=new TFile(file_data.c_str());
  TH2D *h2raw=(TH2D*)f_data->Get("raw");
  TH2D *h2smeared=(TH2D*)f_mc->Get("smeared");
  TH2D *h2true=(TH2D*)f_mc->Get("true");
  TH2D *h2fulleff=(TH2D*)f_mc->Get("truef");
  RooUnfoldResponse* response = (RooUnfoldResponse*)f_mc->Get("smeared_true");


     TH2D *effnum=(TH2D*)h2fulleff->Clone("effnum");
     TH2D *effdenom=(TH2D*)h2fulleff->Clone("effdenom");
 


  
  TH2D* hfold=(TH2D*)h2raw->Clone("hfold");
  hfold->Sumw2();

  std::stringstream ss;
  ss << "Unfold_nom_" << date << flag<< ".root";
  TFile* fout = new TFile(ss.str().c_str(), "recreate");
  ss.str("");

  h2raw->Write();
  h2smeared->Write();
  h2true->Write();
  h2fulleff->Write();




  

   TH2D *h2input;
   if(flag<20) h2input=(TH2D*)h2raw->Clone("h2input");
   if(flag>=20) h2input=(TH2D*) h2smeared->Clone("h2input");
  //h2input->Write(); 

  //do the unfolding for N iterations

     for(int jar=1;jar<16;jar++){
      Int_t iter=jar;
      cout<<"iteration"<<iter<<endl;
      cout<<"==============Unfold h1====================="<<endl;

      
      RooUnfoldBayes   unfold(response, h2input, iter);    // OR
      TH2D* hunf= (TH2D*) unfold.Hreco(errorTreatment);
      //FOLD BACK
      TH2D* hfold= (TH2D*)response->ApplyToTruth(hunf,"");
     


          TH2D *htempUnf=(TH2D*)hunf->Clone("htempUnf");          
     	  htempUnf->SetName(Form("Bayesian_Unfoldediter%d.root",iter));
           
     	  TH2D *htempFold=(TH2D*)hfold->Clone("htempFold");          
     	  htempFold->SetName(Form("Bayesian_Foldediter%d.root",iter));        

      	  htempUnf->Write();
     	  htempFold->Write();}
          fout->Close();
	     
}

#ifndef __CINT__
int main () { RooUnfoldExample_GJ(); return 0; }  // Main program when run stand-alone
#endif