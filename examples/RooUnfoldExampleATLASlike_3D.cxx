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
#include "TH3D.h"
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

/*
 *
 *
ReadResponseATLASlike_JECdown.root
ReadResponseATLASlike_JECup.root
ReadResponseATLASlike_JERdown.root
ReadResponseATLASlike_JERup.root
ReadResponseATLASlike_PUdown.root
ReadResponseATLASlike_PUup.root
ReadResponseATLASlike_herwig7.root
ReadResponseATLASlike_nom.root
ReadResponseATLASlike_trackIneff.root
 *
 *
 * */

void RooUnfoldExampleATLASlike_3D(std::string file_mc="/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/ReadResponseATLASlike_ak4_v2.root", std::string file_data="/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/ReadDataATLASlike2018_test_B_0_ak4.root", std::string date = "ATLAScutsNom",int flag=19)
{
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
  TFile* f_mc=new TFile(file_mc.c_str());
  TFile* f_data=new TFile(file_data.c_str());
  TH3D *h2raw=(TH3D*)f_data->Get("raw");
  h2raw->SetName("raw");
  TH3D *h2recoall=(TH3D*)f_mc->Get("recoall");
  h2recoall->SetName("recoall");
  TH3D *h2rawTrueBinning=(TH3D*)f_data->Get("rawTrueBinning");
  TH3D *h2smeared=(TH3D*)f_mc->Get("smeared");
  TH3D *h2purity = (TH3D*) h2smeared->Clone();
  h2purity->SetName("purity");
  h2purity->Sumw2();
//  h2purity->Divide(h2recoall);

  h2purity->Divide(h2purity, h2recoall, 1., 1., "B" ); 

//  TH2D *h2recoall=(TH2D*)f_mc->Get("recoall");
//  TH3D *h2purity=(TH3D*)f_mc->Get("purity");

  TH1D *pTjet_det=(TH1D*)f_data->Get("pTjet_det");
  pTjet_det->SetName("pTjet_raw");

  TH1D *pTjet_smeared=(TH1D*)f_mc->Get("pTjet_det");
  pTjet_smeared->SetName("pTjet_smeared");

  TH1D *pTjet_true=(TH1D*)f_mc->Get("pTjet_true");
  pTjet_true->SetName("pTjet_true");


  TH1D *nJets=(TH1D*)f_mc->Get("nJets");
  TH1D *nSplittings=(TH1D*)f_mc->Get("nSplittings");

  TH3D *htruef=(TH3D*)f_mc->Get("truef");
  TH3D *htrue=(TH3D*)f_mc->Get("true");
  htruef->SetName("truefDummy");
  htrue->SetName("trueDummy");
//  htruef->GetYaxis()->SetRangeUser(0.084,6.5);
//  htrue->GetYaxis()->SetRangeUser(0.084,6.5);

//  htruef->GetYaxis()->SetRangeUser(-1,6);
//  htrue->GetYaxis()->SetRangeUser(-1,6);

  htruef->GetZaxis()->SetRange(2,2);
  htrue->GetZaxis()->SetRange(2,2);

  TH2D *efficiency = (TH2D*) htrue->Project3D("yx");
  TH2D *efficiencyDenom = (TH2D*) htruef->Project3D("yx");
  efficiency->Sumw2();
//  efficiency->Divide(efficiencyDenom);
  efficiency->Divide(efficiency, efficiencyDenom, 1., 1., "B" );
//  efficiency->GetYaxis()->SetRangeUser(0.084,6.5);
//  efficiency->GetYaxis()->SetRangeUser(-1,6);  
  efficiency->SetName("efficiency");
  TH3D *h2true=(TH3D*)f_mc->Get("true");
  TH3D *h2fulleff=(TH3D*)f_mc->Get("truef");
  RooUnfoldResponse* response = (RooUnfoldResponse*)f_mc->Get("smeared_true");
  RooUnfoldResponse* pT_response = (RooUnfoldResponse*)f_mc->Get("pTjet_det_pTjet_true");

//  TH2D *effnum=(TH2D*)h2fulleff->Clone("effnum");
//  TH2D *effdenom=(TH2D*)h2fulleff->Clone("effdenom");
  
  TH3D* hfold=(TH3D*)h2raw->Clone("hfold");
  hfold->Sumw2();

  std::stringstream ss;
  ss << "Unfold_nom_" << date << flag<< ".root";
  TFile* fout = new TFile(ss.str().c_str(), "recreate");
  ss.str("");
//  h2raw->Divide(h2purity);
//  if (flag != 20) h2raw->Divide(h2purity);
  h2raw->Write();
  h2smeared->Write();
  h2recoall->Write();
  pTjet_smeared->Write();
  pTjet_det->Write();
  h2true->Write();
  pTjet_true->Write();
  h2fulleff->Write();
  nJets->Write();
  nSplittings->Write();

   TH3D *h2input; TH1D *h2input_pT;
   if(flag<20) h2input=(TH3D*)h2raw->Clone("h2input");
   if(flag>=20) h2input=(TH3D*) h2smeared->Clone("h2input");

   if(flag<20) h2input_pT=(TH1D*) pTjet_det->Clone("h2input_pT");
   if(flag>=20) h2input_pT=(TH1D*) pTjet_smeared->Clone("h2input_pT");

   h2input->Write(); h2input_pT->Write();
   h2input->Multiply(h2purity);
   //do the unfolding for N iterations

//   TMatrixT<double> covariance();

     for(int jar=0;jar<16;jar++)
     {
      Int_t iter=jar;
      cout<<"iteration"<<iter<<endl;
      cout<<"==============Unfold h1====================="<<endl;



      
      RooUnfoldBayes   unfold(response, h2input, iter);
      RooUnfoldBayes   unfold_pT(pT_response, h2input_pT, iter);

//      unfold.IncludeSystematics(2);
//      
//      unfold_pT.IncludeSystematics(2);

      auto covariance = unfold.Ereco();



      TH3D* hunf= (TH3D*) unfold.Hreco(errorTreatment);
      TH1D* hunf_pT = (TH1D*) unfold_pT.Hreco(errorTreatment);
      //FOLD BACK
      TH3D* hfold= (TH3D*)response->ApplyToTruth(hunf,"");
      TH1D* hfold_pT = (TH1D*)pT_response->ApplyToTruth(hunf_pT,"");

     
          TH3D *htempUnf=(TH3D*)hunf->Clone("htempUnf");
          TH1D *htempUnf_pT=(TH1D*)hunf_pT->Clone("htempUnfPt");
     	  htempUnf->SetName(Form("Bayesian_Unfoldediter%d.root",iter));
          htempUnf_pT->SetName(Form("Bayesian_UnfoldedPtiter%d.root",iter));
//          if (jar == 4) covariance.SetName(Form("covariance%d",iter ));
//          htempUnf->Divide(efficiency);
           
     	  TH3D *htempFold=(TH3D*)hfold->Clone("htempFold");          
          TH1D *htempFold_pT=(TH1D*)hfold_pT->Clone("htempFoldPt");
     	  htempFold->SetName(Form("Bayesian_Foldediter%d.root",iter));       
          htempFold_pT->SetName(Form("Bayesian_FoldedPtiter%d.root",iter)); 
          covariance.Write();
      	  htempUnf->Write();
          htempUnf_pT->Write();
          htempFold->Divide(h2purity);
     	  htempFold->Write();
          htempFold_pT->Write();
     }
          h2purity->Write();
          efficiency->Write();
          fout->Close();	     
}

#ifndef __CINT__
int main () { RooUnfoldExampleATLASlike_3D(); return 0; }  // Main program when run stand-alone
#endif
