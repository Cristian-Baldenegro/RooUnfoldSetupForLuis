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


void RooUnfoldExample_matrixinversion(std::string file_mc="/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/ResponseIsRead.root", std::string file_data="/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/ResponseIsRead.root", std::string date = "ConditionNumber",int flag=21)
{

  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
  TFile* f_mc=new TFile(file_mc.c_str());
  TFile* f_data=new TFile(file_data.c_str());
  TH2D *h2raw=(TH2D*)f_mc->Get("smeared");
  TH2D *h2rawTrueBinning=(TH2D*)f_data->Get("rawTrueBinning");
  TH2D *h2smeared=(TH2D*)f_mc->Get("smeared");
  TH2D *h2purity=(TH2D*)f_mc->Get("purity");
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
  h2raw->Divide(h2purity);
  h2raw->Write();
  h2smeared->Write();
  h2true->Write();
  h2fulleff->Write();
 

   TH2D *h2input;
   if(flag<20) h2input=(TH2D*)h2raw->Clone("h2input");
   if(flag>=20) h2input=(TH2D*) h2smeared->Clone("h2input");
  //h2input->Write(); 
//   h2input->Divide(h2purity);
   h2input->Write();
  //do the unfolding for N iterations



      RooUnfoldInvert  unfold(response, h2input);

      TDecompSVD *svd= new TDecompSVD (response->Mresponse()); //this is the singular value decomposition (SVD) matrix
      auto singular_values = svd->GetSig(); //this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
//      svd->Print();
      double singular_value_min;
      for (int i = 0; i < 10 ; i++)
      {
        cout << singular_values[i] << endl;
       if ( singular_values[i] > pow(10,-15)  ) singular_value_min = singular_values[i]; // the pow(10,-15) > requirement is to suppress singular values from empty bins 
      }


      cout << "condition number: "  << singular_values[0]/singular_value_min << endl;

//     RooUnfoldSvd  unfold(response, h2input, 1); 
//      RooUnfoldBayes  unfold(response, h2input, 1);    // OR
      TH2D* hunf= (TH2D*) unfold.Hreco();

          TH2D *htempUnf=(TH2D*)hunf->Clone("htempUnf");          
     	  htempUnf->SetName("MatrixInversionUnfolding.root");
           
      	  htempUnf->Write();
//     	  htempFold->Write();

          fout->Close();
	     
}

#ifndef __CINT__
int main () { RooUnfoldExample_matrixinversion(); return 0; }  // Main program when run stand-alone
#endif
