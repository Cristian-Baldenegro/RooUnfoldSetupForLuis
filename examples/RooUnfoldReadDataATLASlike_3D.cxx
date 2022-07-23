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
#include "LundPlane_LLR/AnalysisFW/interface/QCDEvent.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

//==============================================================================
// Global definitions
//==============================================================================

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


void RooUnfoldReadDataATLASlike_3D(std::string cFiles2 = "",std::string tag = "", std::string date = "", int fileindex=0,  int year = 2018, std::string era = "B",  std::string algo = "")
{


  //errors from the diagonal elements of the covariance matrix
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
 
   int bin_det_lnkt;
   int bin_true_lnkt;
   int bin_det_lntheta;
   int bin_true_lntheta;

   Double_t lnthetamin_det;
   Double_t lnthetamin_true;

   Double_t lnthetamax_det;
   Double_t lnthetamax_true;

   Double_t lnktmax_det;
   Double_t lnktmax_true;

   Double_t lnktmin_det;
   Double_t lnktmin_true;

   if (algo == "ak4")
   {
   bin_det_lnkt=20;
   bin_true_lnkt=20;
   bin_det_lntheta=14;
   bin_true_lntheta=16;

   lnthetamin_det=0.0;
   lnthetamin_true=0;

   lnthetamax_det=2.6;
   lnthetamax_true=3;

   lnktmax_det=0.2;
   lnktmax_true=0.3;

   lnktmin_det=log(0.4);
   lnktmin_true=-6.;
   }

   if (algo == "ak8")
   {
   bin_det_lnkt=12;
   bin_true_lnkt=15;
   bin_det_lntheta=16;
   bin_true_lntheta=18;

   lnthetamin_det=0.0;
   lnthetamin_true=0;

   lnthetamax_det=2.6;
   lnthetamax_true=3;

   lnktmax_det=0.2;
   lnktmax_true=0.3;

   lnktmin_det=log(0.4)+0.5*2.;
   lnktmin_true=-6.;
   }

   /***************************************************/

   //det binning in xj
     Double_t ybins[bin_det_lnkt];
   if (algo == "ak8")
   {
         ybins[0]=log(0.4)+0.5*2.;
         ybins[1]=log(0.4)+0.5*3.;
         ybins[2]=log(0.4)+0.5*4.;
         ybins[3]=log(0.4)+0.5*5.;
         ybins[4]=log(0.4)+0.5*6.;
         ybins[5]=log(0.4)+0.5*7.;
         ybins[6]=log(0.4)+0.5*8.;
         ybins[7]=log(0.4)+0.5*9.;
         ybins[8]=log(0.4)+0.5*10.;
         ybins[9]=log(0.4)+0.5*11.;
         ybins[10]=log(0.4)+0.5*12.;
         ybins[11]=log(0.4)+0.5*15.;

    }

   if (algo == "ak4")
   {
         ybins[0]=0.693147;
         ybins[1]=0.693147+0.27726;
         ybins[2]=0.693147+0.27726*2;
         ybins[3]=0.693147+0.27726*3;
         ybins[4]=0.693147+0.27726*4;
         ybins[5]=0.693147+0.27726*5;
         ybins[6]=0.693147+0.27726*6;
         ybins[7]=0.693147+0.27726*7;
         ybins[8]=0.693147+0.27726*8;
         ybins[9]=0.693147+0.27726*9;
         ybins[10]=0.693147+0.27726*10;
         ybins[11]=0.693147+0.27726*11;
         ybins[12]=0.693147+0.27726*12;
         ybins[13]=0.693147+0.27726*13;
         ybins[14]=0.693147+0.27726*14;
         ybins[15]=0.693147+0.27726*15;
         ybins[16]=0.693147+0.27726*16;
         ybins[17]=0.693147+0.27726*17;
         ybins[18]=0.693147+0.27726*18;
         ybins[19]=0.693147+0.27726*19; //5.822436317
//         ybins[19]=; //5.822436317
    }


          //det binning in Rg
         Double_t xbins[bin_det_lntheta];
      if (algo == "ak4")
      {

         xbins[0]=0;
         xbins[1]=.33333333;
         xbins[2]=.66666667;
         xbins[3]=1.;
         xbins[4]=1.33333333;
         xbins[5]=1.66666667;
         xbins[6]=2.0;
         xbins[7]=2.33333333;
         xbins[8]=2.66666667;
         xbins[9]=3.0;
         xbins[10]=3.33333333;
         xbins[11]=3.66666667;
         xbins[12]=4.0;
         xbins[13]=4.33333333;
      }

      if (algo == "ak8")
      {

         xbins[0]=0;
         xbins[1]=.33333333;
         xbins[2]=.66666667;
         xbins[3]=1.;
         xbins[4]=1.33333333;
         xbins[5]=1.66666667;
         xbins[6]=2.0;
         xbins[7]=2.33333333;
         xbins[8]=2.66666667;
         xbins[9]=3.0;
         xbins[10]=3.33333333;
         xbins[11]=3.66666667;
         xbins[12]=4.0;
         xbins[13]=4.33333333;
         xbins[14]=4.66666667;
         xbins[15]=5.0;
      }

         Double_t pT_det[2];
         pT_det[0]=450.0;
         pT_det[1]=4000.;


   JetCorrectorParameters *ResJetPar ;
   JetCorrectorParameters *L3JetPar ;
   JetCorrectorParameters *L2JetPar ;
   JetCorrectorParameters *L1JetPar;
   JetCorrectionUncertainty *jecUnc ;


if (year == 2016 && era == "BCD" && algo == "ak4")
{
cout << year << era << endl;

   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_L2L3Residual_AK4PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_L3Absolute_AK4PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_L2Relative_AK4PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_L1FastJet_AK4PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");

}

if (year == 2016 && era == "EF" && algo == "ak4")
{
cout << year << era << endl;

   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunEF_V7_DATA_L2L3Residual_AK4PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunEF_V7_DATA_L3Absolute_AK4PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunEF_V7_DATA_L2Relative_AK4PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunEF_V7_DATA_L1FastJet_AK4PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");

}


if (year == 2016 && era == "FGH" && algo == "ak4")
{
cout << year << era << endl;

   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_data/Summer19UL16_RunFGH_V7_DATA_L2L3Residual_AK4PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_data/Summer19UL16_RunFGH_V7_DATA_L3Absolute_AK4PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_data/Summer19UL16_RunFGH_V7_DATA_L2Relative_AK4PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_data/Summer19UL16_RunFGH_V7_DATA_L1FastJet_AK4PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");

}


if (year == 2017 && era == "B" && algo == "ak4")
{
cout << year << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunB_V6_DATA_L2L3Residual_AK4PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunB_V6_DATA_L3Absolute_AK4PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunB_V6_DATA_L2Relative_AK4PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunB_V6_DATA_L1FastJet_AK4PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}

if (year == 2017 && era == "C" && algo == "ak4")
{
cout << year << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunC_V6_DATA_L2L3Residual_AK4PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunC_V6_DATA_L3Absolute_AK4PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunC_V6_DATA_L2Relative_AK4PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunC_V6_DATA_L1FastJet_AK4PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}

if (year == 2017 && era == "D" && algo == "ak4")
{
cout << year << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunD_V6_DATA_L2L3Residual_AK4PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunD_V6_DATA_L3Absolute_AK4PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunD_V6_DATA_L2Relative_AK4PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunD_V6_DATA_L1FastJet_AK4PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}

if (year == 2017 && era == "E" && algo == "ak4")
{
cout << year << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunE_V6_DATA_L2L3Residual_AK4PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunE_V6_DATA_L3Absolute_AK4PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunE_V6_DATA_L2Relative_AK4PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunE_V6_DATA_L1FastJet_AK4PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}

if (year == 2017 && era == "F" && algo == "ak4")
{
cout << year << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunF_V6_DATA_L2L3Residual_AK4PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunF_V6_DATA_L3Absolute_AK4PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunF_V6_DATA_L2Relative_AK4PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunF_V6_DATA_L1FastJet_AK4PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}


if (year == 2018 && era == "A" && algo == "ak4")
{
cout << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunA_V5_DATA_L2L3Residual_AK4PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunA_V5_DATA_L3Absolute_AK4PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunA_V5_DATA_L2Relative_AK4PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunA_V5_DATA_L1FastJet_AK4PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}


if ( year == 2018 && era == "B" && algo == "ak4")
{

cout << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunB_V5_DATA_L2L3Residual_AK4PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunB_V5_DATA_L3Absolute_AK4PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunB_V5_DATA_L2Relative_AK4PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunB_V5_DATA_L1FastJet_AK4PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}


if (year == 2018 && era == "C" && algo == "ak4" )
{
cout << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunC_V5_DATA_L2L3Residual_AK4PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunC_V5_DATA_L3Absolute_AK4PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunC_V5_DATA_L2Relative_AK4PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunC_V5_DATA_L1FastJet_AK4PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}

if ( year == 2018 && era == "D" && algo == "ak4")
{
cout << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunD_V5_DATA_L2L3Residual_AK4PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunD_V5_DATA_L3Absolute_AK4PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunD_V5_DATA_L2Relative_AK4PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunD_V5_DATA_L1FastJet_AK4PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}

if (year == 2016 && era == "BCD" && algo == "ak8")
{
cout << year << era << endl;

   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_L2L3Residual_AK8PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_L3Absolute_AK8PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_L2Relative_AK8PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_L1FastJet_AK8PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");

}

if (year == 2016 && era == "EF" && algo == "ak8")
{
cout << year << era << endl;

   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunEF_V7_DATA_L2L3Residual_AK8PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunEF_V7_DATA_L3Absolute_AK8PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunEF_V7_DATA_L2Relative_AK8PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunEF_V7_DATA_L1FastJet_AK8PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");

}


if (year == 2016 && era == "FGH" && algo == "ak8")
{
cout << year << era << endl;

   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_data/Summer19UL16_RunFGH_V7_DATA_L2L3Residual_AK8PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_data/Summer19UL16_RunFGH_V7_DATA_L3Absolute_AK8PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_data/Summer19UL16_RunFGH_V7_DATA_L2Relative_AK8PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_data/Summer19UL16_RunFGH_V7_DATA_L1FastJet_AK8PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");

}


if (year == 2017 && era == "B" && algo == "ak8")
{
cout << year << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunB_V6_DATA_L2L3Residual_AK8PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunB_V6_DATA_L3Absolute_AK8PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunB_V6_DATA_L2Relative_AK8PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunB_V6_DATA_L1FastJet_AK8PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}

if (year == 2017 && era == "C" && algo == "ak8")
{
cout << year << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunC_V6_DATA_L2L3Residual_AK8PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunC_V6_DATA_L3Absolute_AK8PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunC_V6_DATA_L2Relative_AK8PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunC_V6_DATA_L1FastJet_AK8PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}

if (year == 2017 && era == "D" && algo == "ak8")
{
cout << year << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunD_V6_DATA_L2L3Residual_AK8PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunD_V6_DATA_L3Absolute_AK8PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunD_V6_DATA_L2Relative_AK8PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunD_V6_DATA_L1FastJet_AK8PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}

if (year == 2017 && era == "E" && algo == "ak8")
{
cout << year << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunE_V6_DATA_L2L3Residual_AK8PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunE_V6_DATA_L3Absolute_AK8PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunE_V6_DATA_L2Relative_AK8PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunE_V6_DATA_L1FastJet_AK8PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}

if (year == 2017 && era == "F" && algo == "ak8")
{
cout << year << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunF_V6_DATA_L2L3Residual_AK8PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunF_V6_DATA_L3Absolute_AK8PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunF_V6_DATA_L2Relative_AK8PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_DATA/Summer19UL17_RunF_V6_DATA_L1FastJet_AK8PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}


if (year == 2018 && era == "A" && algo == "ak8")
{
cout << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunA_V5_DATA_L2L3Residual_AK8PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunA_V5_DATA_L3Absolute_AK8PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunA_V5_DATA_L2Relative_AK8PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunA_V5_DATA_L1FastJet_AK8PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}


if ( year == 2018 && era == "B" && algo == "ak8")
{

cout << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunB_V5_DATA_L2L3Residual_AK8PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunB_V5_DATA_L3Absolute_AK8PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunB_V5_DATA_L2Relative_AK8PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunB_V5_DATA_L1FastJet_AK8PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}


if (year == 2018 && era == "C" && algo == "ak8" )
{
cout << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunC_V5_DATA_L2L3Residual_AK8PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunC_V5_DATA_L3Absolute_AK8PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunC_V5_DATA_L2Relative_AK8PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunC_V5_DATA_L1FastJet_AK8PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}

if ( year == 2018 && era == "D" && algo == "ak8")
{
cout << era << endl;
  ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunD_V5_DATA_L2L3Residual_AK8PFPuppi.txt");
  L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunD_V5_DATA_L3Absolute_AK8PFPuppi.txt");
  L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunD_V5_DATA_L2Relative_AK8PFPuppi.txt");
  L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_DATA/Summer19UL18_RunD_V5_DATA_L1FastJet_AK8PFPuppi.txt");
  jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_data/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFPuppi.txt");
}



   vector<JetCorrectorParameters> vPar;

   vPar.push_back(*L1JetPar);
   vPar.push_back(*L2JetPar);
   vPar.push_back(*L3JetPar);
   vPar.push_back(*ResJetPar);
//
//
   FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(vPar);
      
  FactorizedJetCorrector *JetCorrectorRecoLoop = new FactorizedJetCorrector(vPar);

   TFile *fhotMap; TH2F *hotMap;


   if (year == 2016)
   {
   fhotMap = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/hotcoldmaps/hotjets-16.root");
   hotMap = (TH2F*)fhotMap->Get("h2hotfilter");
   }


   if (year == 2017)
   {
   fhotMap = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/hotcoldmaps/hotjets-17.root");
   hotMap = (TH2F*)fhotMap->Get("h2hotfilter");
   }



   if (year == 2018)
   {
   fhotMap = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/hotcoldmaps/hotjets-18.root");
   hotMap = (TH2F*)fhotMap->Get("h2hotfilter");
   }




   TH1D *jet1Pt(0), *jet2Pt(0), *jetPtAsymmetry(0), *DeltaPhi(0), *jet1Rapidity(0), *jet2Rapidity(0), *DeltaRapidity(0), *nJets(0), *nSplittings(0);

   nJets = new TH1D("nJets", "nJets",11,-0.5,10.5);
   nSplittings = new TH1D("nSplittings", "nSplittings",16,-0.5,15.5);
   TH1D* jet1Pt_fired = new TH1D("jet1Pt_fired", "jet1Pt",100,0,1000);
   jet1Pt = new TH1D("jet1Pt", "jet1Pt",98,550,3000);
   jet2Pt = new TH1D("jet2Pt", "jet2Pt",98,550,3000);
   jetPtAsymmetry = new TH1D("jetPtAsymmetry", "jetPtAsymmetry",20,0, 1);
   DeltaPhi = new TH1D("DeltaPhi", "DeltaPhi",32,0, M_PI);
   jet1Rapidity =  new TH1D("jet1Rapidity", "jet1Rapidity",20,-2.5, 2.5);
   jet2Rapidity =  new TH1D("jet2Rapidity", "jet2Rapidity",20,-2.5, 2.5);
   DeltaRapidity =  new TH1D("DeltaRapidity", "DeltaRapidity",20,0, 4.0);

   int nbins = 14;
   int nbins_dR = 14;
   float logkTmax = 6.;
   float nbinsTrue = 10;
   float logkTmaxTrue = 10.0;
   float logkTmin = log(0.4);
   float logRmax = 4.33333333;
   float logRmaxTrue = 6.0 ;


   if (algo == "ak8")
   {
    logkTmin = log(0.8); logRmax = 5.000; logkTmax = 6.5;
   }

      double R = 0.4;

      if (algo == "ak8") R = 0.8;

   //the raw correlation
   TH3D *h2raw(0), *h2rawTrueBinning(0);
   h2raw=new TH3D("raw","raw",bin_det_lntheta-1, xbins,bin_det_lnkt-1,ybins, 1, pT_det );
   h2rawTrueBinning=new TH3D("rawTrueLevelBinning","rawTrueLevelBinning",bin_det_lntheta-1, xbins,bin_det_lnkt-1,ybins,  1, pT_det );
   h2raw->Sumw2();
   h2raw->Sumw2();

   TH1D *pfRho_density(0);
   pfRho_density = new TH1D("pfRho_density", "pfRho_density", 60, 0, 60);

   TH1D * nVtxGood = new TH1D("nvtxGood", "nvtxGood" ,50, 0,50);

   TH1D *pTjet_det(0);
   pTjet_det = new TH1D("pTjet_det", "pTjet_det", 1, pT_det);



    ifstream infile;
    infile.open(cFiles2);
    char filename[300];
    int count=0; int n_splittings = 0;
    while(infile>>filename)
    {

      TFile *input2=TFile::Open(filename);
      cout << filename << endl;
      count=count+1;
      cout<<"reading data file "<<count<<endl;
      if(count==1000) break;
//
//      delete input2;      
      TDirectory *tdir;
      if (algo == "ak4") tdir = input2->GetDirectory("ak4");
      if (algo == "ak8") tdir = input2->GetDirectory("ak8");

      TTree *T = (TTree*)tdir->Get("ProcessedTree");

      QCDEvent        *events;
      TBranch        *b_events;   //!

      // Set object pointer
      events = 0;
      // Set branch addresses and branch pointers
      if (!T) return;     
      T->SetMakeClass(1);     

      Int_t nEvmc=T->GetEntries();
      T->SetBranchAddress("events", &events, &b_events);



   for(int iEntry=0; iEntry< nEvmc; iEntry++)
   {
    T->GetEntry(iEntry);

//      if (1 == 1) continue;


   if ( events->nPFJetsCHS() < 2) continue;
    auto eventHdr = events->evtHdr();
//    auto otherEventInfo = events->evt();
    double pfRho = eventHdr.pfRho();
    auto jet1 = events->pfjetchs(0);
    auto jet2 = events->pfjetchs(1);
    bool HLT_450PFjet = false;
//    for (int i = 0; i < ; i++  )
//    {
//    cout << events->fired(0) << endl;
      if (events->fired(0) == 8 && year == 2016 ) HLT_450PFjet = true;      
      if (algo == "ak4" && events->fired(0) == 9 && ( year == 2018 || year == 2017) ) HLT_450PFjet = true;
//    }


//      if (1 == 1) continue;


    if ( HLT_450PFjet == false ) continue;
    jet1Pt_fired->Fill(jet1.pt());
//    if (events->fired(0) == 0 ) cout << jet1.pt() << endl;


    JetCorrector->setJetEta(jet1.eta());
    JetCorrector->setJetPt(jet1.pt()/jet1.cor());
    JetCorrector->setRho(pfRho);
    JetCorrector->setJetA(jet1.area());

     double newCorrection = JetCorrector->getCorrection();

     jecUnc->setJetEta(jet1.eta());
     jecUnc->setJetPt( jet1.pt()/jet1.cor()*newCorrection);
     double unc = jecUnc->getUncertainty(true); 

     double jet1_pT_corrected = jet1.pt()/jet1.cor()*newCorrection;


     JetCorrector->setJetEta(jet2.eta());
     JetCorrector->setJetPt(jet2.pt()/jet2.cor());
     JetCorrector->setJetA(jet2.area());
     JetCorrector->setRho(pfRho);

      newCorrection = JetCorrector->getCorrection();

     jecUnc->setJetEta(jet2.eta());
     jecUnc->setJetPt( jet2.pt()/jet2.cor()*newCorrection);
     unc = jecUnc->getUncertainty(true);

     double jet2_pT_corrected = jet2.pt()/jet2.cor()*newCorrection;

    pfRho_density->Fill(eventHdr.pfRho());
    nVtxGood->Fill(eventHdr.nVtxGood());
    if (!( jet1_pT_corrected > 675.0 && jet2_pT_corrected > 450 && fabs(jet1.eta()) < 2.1 && fabs(jet2.eta()) < 2.1 && jet1_pT_corrected < 4000 && jet2_pT_corrected < 4000 && jet1.looseID() == true && jet2.looseID() == true && jet2_pT_corrected/jet1_pT_corrected > 2./3.)) continue;
//    if (jet1.looseID() == false || jet2.looseID() == false ) continue;
//    if (!(jet1_pT_corrected > 550 && jet2_pT_corrected > 550 && fabs(jet1.y()) < 2 && fabs(jet2.y()) <2 && jet1_pT_corrected < 4000 && jet2_pT_corrected < 4000 )) continue;
/*    jet1Pt->Fill(jet1_pT_corrected ); //, *jet2Pt(0), *jetPtAsymmetry, *DeltaPhi(0), *jet1Rapidity(0), *jet2Rapidity(0), *DeltaRapidity(0);
    jet2Pt->Fill(jet2_pT_corrected );
    double delta_phi = 0.;
    if (fabs(jet1.phi()-jet2.phi()) < M_PI) delta_phi = fabs(jet1.phi()-jet2.phi());
    else delta_phi = 2*M_PI-fabs(jet1.phi()-jet2.phi());
    DeltaPhi->Fill(delta_phi);
    jet1Rapidity->Fill(jet1.y());
    jet2Rapidity->Fill(jet2.y());
    DeltaRapidity->Fill(fabs(jet1.y()-jet2.y()));
    pfRho_density->Fill(eventHdr.pfRho());
    double pTasymmetry = fabs(jet1_pT_corrected-jet2_pT_corrected )/(jet1_pT_corrected+jet2_pT_corrected);
//    cout << jet1.cor() << endl;
    jetPtAsymmetry->Fill(pTasymmetry);*/



   for (int k = 0; k < 2; k++) // det-level jets
   {
     auto jet_reco = events->pfjetchs(k);

//     if ( hotMap->GetBinContent(hotMap->FindBin(jet_reco.eta(),jet_reco.phi())) > 5 ) cout << hotMap->GetBinContent(hotMap->FindBin(jet_reco.eta(),jet_reco.phi())) << endl; 

     if ( hotMap->GetBinContent(hotMap->FindBin(jet_reco.eta(),jet_reco.phi())) > 5 ) continue;

     auto z_reco = jet_reco.z();
     auto theta_reco = jet_reco.theta();
     auto kT_reco = jet_reco.kT();
     auto eta2_reco = jet_reco.eta2_splitting();
     auto phi2_reco = jet_reco.phi2_splitting();

      JetCorrectorRecoLoop->setJetEta(jet_reco.eta());
      JetCorrectorRecoLoop->setJetPt(jet_reco.pt()/jet_reco.cor());
      JetCorrectorRecoLoop->setRho(pfRho );
      JetCorrectorRecoLoop->setJetA(jet_reco.area());

      newCorrection = JetCorrectorRecoLoop->getCorrection();
//       cout << jet_reco.cor()/newCorrection << endl;
      jecUnc->setJetEta(jet_reco.eta() );
      jecUnc->setJetPt( jet_reco.pt()/jet_reco.cor()*newCorrection);
      unc = jecUnc->getUncertainty(true);

      double jetPtCorrected = jet_reco.pt()/jet_reco.cor()*newCorrection; //rescaling pT with jet energy resolution, new jet energy scale
      if  (!(   jetPtCorrected > 450.0 && fabs(jet_reco.eta()) < 2.1 && jet_reco.looseID() == true && jetPtCorrected < 4000  ) ) continue;
//      if  (!(   z_reco.size() > 0 && jetPtCorrected > 550 && fabs(jet_reco.y()) < 1.7 && jet_reco.looseID() == true && jetPtCorrected < 4000 ) ) continue;



       jet1Pt->Fill(jetPtCorrected );
       jet1Rapidity->Fill(jet_reco.eta());

      int first_splitting = 0;
      n_splittings = 0;
      for (unsigned j = 0; j < z_reco.size(); ++j) //loop over reco-level splittings
      {
       
       if (log(1/z_reco.at(j)) < 0.693147 || log(1./z_reco.at(j)) > 0.693147+0.27726*19 ) continue;
       if (log(R/theta_reco.at(j)) < 0. || log(R/theta_reco.at(j)) > logRmax) continue;
       n_splittings++;
       h2raw->Fill(log(R/theta_reco.at(j)) , log(1./z_reco.at(j) ) , jetPtCorrected );
       h2rawTrueBinning->Fill(log(R/theta_reco.at(j)) , log(1./z_reco.at(j) ) , jetPtCorrected );
      }
       nSplittings->Fill(n_splittings);
       nJets->Fill(1);
       pTjet_det->Fill(jetPtCorrected);
     }
     
    }
  }




//    TFile *fout=new TFile ("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/rawData.root","RECREATE");
    TFile *fout=new TFile (Form("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/ReadDataATLASlike%i_%s_%s_%i_%s.root",year, tag.c_str(),era.c_str(), fileindex, algo.c_str() ),"RECREATE");
    fout->cd();
  jet1Pt->Write(); jet1Pt_fired->Write(); jet2Pt->Write(); jetPtAsymmetry->Write(); DeltaPhi->Write(); jet1Rapidity->Write(); jet2Rapidity->Write(); DeltaRapidity->Write();
    nSplittings->Write(); nJets->Write(); pfRho_density->Write();
    h2rawTrueBinning->SetName("rawTrueBinning");
    h2raw->SetName("raw");
    h2rawTrueBinning->Write();
    h2raw->Write();
    pTjet_det->Write(); nVtxGood->Write();
}

#ifndef __CINT__
int main () { RooUnfoldReadDataATLASlike_3D(); return 0; }  // Main program when run stand-alone
#endif
