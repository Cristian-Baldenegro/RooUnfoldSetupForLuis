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
#include <algorithm>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TH3D.h"


#include "TFile.h"
#include "TVectorD.h"


#include "TLorentzVector.h"
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
//#include "/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/interface/QCDEvent.h"

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


void RooUnfoldReadResponse_3D(std::string cFiles2="",std::string tag = "", std::string date = "", int fileindex=0, int year = 2018, std::string mc_name = "", std::string jec_option = "", std::string jer_option ="", std::string pu_option = "", std::string algo = "" )
{

#ifdef __CINT__
  gSystem->Load("libRooUnfold.so");
//  gSystem->Load("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/interface/QCDEvent.h");
#endif
  //errors from the diagonal elements of the covariance matrix
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;

   cout << Form("ReadResponse%s_%i_JEC%s_JER%s_PU%s.root",mc_name.c_str(), year,jec_option.c_str(),jer_option.c_str(), pu_option.c_str() ) << endl;

   double lumiWeight = 1.;

   if (mc_name == "pythia8")
   {
    if (year == 20161) lumiWeight = 4.185143303;
    if (year == 20162) lumiWeight = 5.391541055;
    if (year == 2017) lumiWeight = 13.32431283;
    if (year == 2018) lumiWeight = 19.38108241;
   }

   if (mc_name == "herwig7")
   {
    if (year == 20161) lumiWeight = 5.4147969;
    if (year == 20162) lumiWeight = 4.742935626;
    if (year == 2017) lumiWeight = 10.67497687;
    if (year == 2018) lumiWeight = 17.726844; 
   }

   JetCorrectorParameters *ResJetPar ;
   JetCorrectorParameters *L3JetPar ;
   JetCorrectorParameters *L2JetPar ;
   JetCorrectorParameters *L1JetPar;
   JetCorrectionUncertainty *jecUnc ;
   if (year == 20161 && algo == "ak4")
   {
   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_MC/Summer19UL16APV_V7_MC_L2L3Residual_AK4PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_MC/Summer19UL16APV_V7_MC_L3Absolute_AK4PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_MC/Summer19UL16APV_V7_MC_L2Relative_AK4PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_MC/Summer19UL16APV_V7_MC_L1FastJet_AK4PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_MC/Summer19UL16APV_V7_MC_Uncertainty_AK4PFPuppi.txt");
   }


   if (year == 20162 && algo == "ak4")
   {
   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_MC/Summer19UL16_V7_MC_L2L3Residual_AK4PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_MC/Summer19UL16_V7_MC_L3Absolute_AK4PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_MC/Summer19UL16_V7_MC_L2Relative_AK4PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_MC/Summer19UL16_V7_MC_L1FastJet_AK4PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_MC/Summer19UL16_V7_MC_Uncertainty_AK4PFPuppi.txt");

   }

   if (year == 2017 && algo == "ak4")
   {
   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_MC/Summer19UL17_V6_MC_L2L3Residual_AK4PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_MC/Summer19UL17_V6_MC_L3Absolute_AK4PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_MC/Summer19UL17_V6_MC_L2Relative_AK4PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_MC/Summer19UL17_V6_MC_L1FastJet_AK4PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_MC/Summer19UL17_V6_MC_Uncertainty_AK4PFPuppi.txt");
   }

   if (year == 2018 && algo == "ak4")
   {
   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_MC/Summer19UL18_V5_MC_L2L3Residual_AK4PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_MC/Summer19UL18_V5_MC_L3Absolute_AK4PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_MC/Summer19UL18_V5_MC_L2Relative_AK4PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_MC/Summer19UL18_V5_MC_L1FastJet_AK4PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_MC/Summer19UL18_V5_MC_Uncertainty_AK4PFPuppi.txt");
   }

   if (year == 20161 && algo == "ak8")
   {
   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_MC/Summer19UL16APV_V7_MC_L2L3Residual_AK8PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_MC/Summer19UL16APV_V7_MC_L3Absolute_AK8PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_MC/Summer19UL16APV_V7_MC_L2Relative_AK8PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_MC/Summer19UL16APV_V7_MC_L1FastJet_AK8PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_preVFP_MC/Summer19UL16APV_V7_MC_Uncertainty_AK8PFPuppi.txt");
   }


   if (year == 20162 && algo == "ak8")
   {
   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_MC/Summer19UL16_V7_MC_L2L3Residual_AK8PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_MC/Summer19UL16_V7_MC_L3Absolute_AK8PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_MC/Summer19UL16_V7_MC_L2Relative_AK8PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_MC/Summer19UL16_V7_MC_L1FastJet_AK8PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2016_postVFP_MC/Summer19UL16_V7_MC_Uncertainty_AK8PFPuppi.txt");

   }

   if (year == 2017 && algo == "ak8")
   {
   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_MC/Summer19UL17_V6_MC_L2L3Residual_AK8PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_MC/Summer19UL17_V6_MC_L3Absolute_AK8PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_MC/Summer19UL17_V6_MC_L2Relative_AK8PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_MC/Summer19UL17_V6_MC_L1FastJet_AK8PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2017_MC/Summer19UL17_V6_MC_Uncertainty_AK8PFPuppi.txt");
   }

   if (year == 2018 && algo == "ak8")
   {
   ResJetPar = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_MC/Summer19UL18_V5_MC_L2L3Residual_AK8PFPuppi.txt");
   L3JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_MC/Summer19UL18_V5_MC_L3Absolute_AK8PFPuppi.txt");
   L2JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_MC/Summer19UL18_V5_MC_L2Relative_AK8PFPuppi.txt");
   L1JetPar  = new JetCorrectorParameters("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_MC/Summer19UL18_V5_MC_L1FastJet_AK8PFPuppi.txt");
   jecUnc = new JetCorrectionUncertainty("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JEC2018_MC/Summer19UL18_V5_MC_Uncertainty_AK8PFPuppi.txt");
   }


   vector<JetCorrectorParameters> vPar;

   vPar.push_back(*L1JetPar);
   vPar.push_back(*L2JetPar);
   vPar.push_back(*L3JetPar);
//   vPar.push_back(*ResJetPar);

   FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(vPar);
   FactorizedJetCorrector *JetCorrectorRecoLoop = new FactorizedJetCorrector(vPar);



  JME::JetResolution resolution ;
  JME::JetResolutionScaleFactor resolution_sf ;
  if (year == 20161 && algo == "ak4")
  {
  resolution = JME::JetResolution("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2016_preVFP_MC/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFPuppi.txt");
  resolution_sf = JME::JetResolutionScaleFactor("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2016_preVFP_MC/Summer20UL16APV_JRV3_MC_SF_AK4PFPuppi.txt");

//  resolution = JME::JetResolution("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2016_preVFP_MC/Summer20UL16_JRV3_MC_PtResolution_AK4PFPuppi.txt");
//  resolution_sf = JME::JetResolutionScaleFactor("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2016_preVFP_MC/Summer20UL16_JRV3_MC_SF_AK4PFPuppi.txt");
  }

  if (year == 20162 && algo == "ak4")
  {
  resolution = JME::JetResolution("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2016_preVFP_MC/Summer20UL16_JRV3_MC_PtResolution_AK4PFPuppi.txt");
  resolution_sf = JME::JetResolutionScaleFactor("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2016_preVFP_MC/Summer20UL16_JRV3_MC_SF_AK4PFPuppi.txt");
  }

  if (year == 2017 && algo == "ak4")
  {
  resolution = JME::JetResolution("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2017/Summer19UL17_JRV3_MC_PtResolution_AK4PFPuppi.txt");
  resolution_sf = JME::JetResolutionScaleFactor("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2017/Summer19UL17_JRV3_MC_SF_AK4PFPuppi.txt");
  }

  if (year == 2018 && algo == "ak4")
  {
  resolution = JME::JetResolution("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2018/Summer19UL18_JRV2_MC_PtResolution_AK4PFPuppi.txt");
  resolution_sf = JME::JetResolutionScaleFactor("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2018/Summer19UL18_JRV2_MC_SF_AK4PFPuppi.txt");
  }


  if (year == 20161 && algo == "ak8")
  {
  resolution = JME::JetResolution("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2016_preVFP_MC/Summer20UL16APV_JRV3_MC_PtResolution_AK8PFPuppi.txt");
  resolution_sf = JME::JetResolutionScaleFactor("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2016_preVFP_MC/Summer20UL16APV_JRV3_MC_SF_AK8PFPuppi.txt");

//  resolution = JME::JetResolution("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2016_preVFP_MC/Summer20UL16_JRV3_MC_PtResolution_AK4PFPuppi.txt");
//  resolution_sf = JME::JetResolutionScaleFactor("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2016_preVFP_MC/Summer20UL16_JRV3_MC_SF_AK4PFPuppi.txt");
  }

  if (year == 20162 && algo == "ak8")
  {
  resolution = JME::JetResolution("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2016_preVFP_MC/Summer20UL16_JRV3_MC_PtResolution_AK8PFPuppi.txt");
  resolution_sf = JME::JetResolutionScaleFactor("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2016_preVFP_MC/Summer20UL16_JRV3_MC_SF_AK8PFPuppi.txt");
  }

  if (year == 2017 && algo == "ak8")
  {
  resolution = JME::JetResolution("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2017/Summer19UL17_JRV3_MC_PtResolution_AK8PFPuppi.txt");
  resolution_sf = JME::JetResolutionScaleFactor("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2017/Summer19UL17_JRV3_MC_SF_AK8PFPuppi.txt");
  }

  if (year == 2018 && algo == "ak8")
  {
  resolution = JME::JetResolution("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2018/Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi.txt");
  resolution_sf = JME::JetResolutionScaleFactor("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/JER2018/Summer19UL18_JRV2_MC_SF_AK8PFPuppi.txt");
  }


   TFile *fPU;
//   fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_weights_pythia8.root");


   TH1F *puReweightNominal;


   TFile *fL1PrefiringFile; TH2F *prefireMap; 

   if (year == 20161 || year == 2018 || year == 20162) //2018 is dummy
   {
   fL1PrefiringFile  = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/prefire/L1prefiring_jetpt_2016BtoH.root");
   prefireMap = (TH2F*)fL1PrefiringFile->Get("L1prefiring_jetpt_2016BtoH");
   }

   if (year == 2017) //2018 is dummy
   { 
   fL1PrefiringFile  = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/prefire/L1prefiring_jetpt_2017BtoF.root");
   prefireMap = (TH2F*)fL1PrefiringFile->Get("L1prefiring_jetpt_2017BtoF");
   }


   TFile *fhotMap; TH2F *hotMap;


   if (year == 20161)
   {
   fhotMap = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/hotcoldmaps/hotjets-16.root");
   hotMap = (TH2F*)fhotMap->Get("h2hotfilter");

   if (mc_name == "herwig7") fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_files/pileup_weights_herwig7_2016preVFP.root");
//   if (mc_name == "herwig7") fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_weights_herwig7.root");
//   if (mc_name == "pythia8") fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_weights_pythia8.root");
   if (mc_name == "pythia8") fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_files/pileup_weights_pythia8_2016preVFP.root");

//   fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_files/pileup_weights_pythia8_2018.root");
   if (pu_option == "nom" ) puReweightNominal=(TH1F*)fPU->Get("puReweightNominal");
   if (pu_option == "down" ) puReweightNominal=(TH1F*)fPU->Get("puReweightDown");
   if (pu_option == "up" ) puReweightNominal=(TH1F*)fPU->Get("puReweightUp");
   }

   if (year == 20162)
   {
   fhotMap = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/hotcoldmaps/hotjets-16.root");
   hotMap = (TH2F*)fhotMap->Get("h2hotfilter");

   if (mc_name == "herwig7") fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_files/pileup_weights_herwig7_2016postVFP.root");

   if (mc_name == "pythia8") fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_files/pileup_weights_pythia8_2016postVFP.root");

   if (pu_option == "nom" ) puReweightNominal=(TH1F*)fPU->Get("puReweightNominal");
   if (pu_option == "down" ) puReweightNominal=(TH1F*)fPU->Get("puReweightDown");
   if (pu_option == "up" ) puReweightNominal=(TH1F*)fPU->Get("puReweightUp");
   }



   if (year == 2017)
   {
   fhotMap = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/hotcoldmaps/hotjets-17.root");
   hotMap = (TH2F*)fhotMap->Get("h2hotfilter");

   if (mc_name == "herwig7") fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_files/pileup_weights_herwig7_2017.root");
   if (mc_name == "pythia8") fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_files/pileup_weights_pythia8_2017.root");

//   fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_files/pileup_weights_pythia8_2018.root");
//   puReweightNominal=(TH1F*)fPU->Get("puReweightNominal");
   if (pu_option == "nom" ) puReweightNominal=(TH1F*)fPU->Get("puReweightNominal");
   if (pu_option == "down" ) puReweightNominal=(TH1F*)fPU->Get("puReweightDown");
   if (pu_option == "up" ) puReweightNominal=(TH1F*)fPU->Get("puReweightUp");

   }



   if (year == 2018)
   {
   fhotMap = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/hotcoldmaps/hotjets-18.root");
   hotMap = (TH2F*)fhotMap->Get("h2hotfilter");

   if (mc_name == "herwig7") fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_files/pileup_weights_herwig7_2018.root");
   if (mc_name == "pythia8") fPU = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/examples/pileup_files/pileup_weights_pythia8_2018.root");
//   puReweightNominal=(TH1F*)fPU->Get("puReweightNominal");
   if (pu_option == "nom" ) puReweightNominal=(TH1F*)fPU->Get("puReweightNominal");
   if (pu_option == "down" ) puReweightNominal=(TH1F*)fPU->Get("puReweightDown");
   if (pu_option == "up" ) puReweightNominal=(TH1F*)fPU->Get("puReweightUp");

   }

//

//  cout << sf << endl;

  TRandom3* rand = new TRandom3(0); 




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
   bin_det_lnkt=13;
   bin_true_lnkt=15;
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
   bin_true_lnkt=14;
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
      
/*         ybins[0]=lnktmin_det;
         ybins[1]=-0.6666667;
	 ybins[2]=-0.3333333;
 	 ybins[3]=0.;
         ybins[4]=0.3333333;
         ybins[5]=0.6666667;
         ybins[6]=1.0;
         ybins[7]=1.3333333;
         ybins[8]=1.6666667;
         ybins[9]=2.0;
         ybins[10]=2.3333333;
         ybins[11]=2.6666667;
         ybins[12]=3.0;
         ybins[13]=3.3333333;
         ybins[14]=3.6666667;
         ybins[15]=4.0;
         ybins[16]=4.3333333;
         ybins[17]=6.0;*/

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
//         ybins[12]=log(0.4)+0.5*15.;

    }

   if (algo == "ak4")
   {
         ybins[0]=lnktmin_det;
         ybins[1]=log(0.4)+0.5;
         ybins[2]=log(0.4)+0.5*2.;
         ybins[3]=log(0.4)+0.5*3.;
         ybins[4]=log(0.4)+0.5*4.;
         ybins[5]=log(0.4)+0.5*5.;
         ybins[6]=log(0.4)+0.5*6.;
         ybins[7]=log(0.4)+0.5*7.;
         ybins[8]=log(0.4)+0.5*8.;
         ybins[9]=log(0.4)+0.5*9.;
         ybins[10]=log(0.4)+0.5*10.;
         ybins[11]=log(0.4)+0.5*11.;
         ybins[12]=6;
    }



       //true binning in xj
        Double_t ybins_true[bin_true_lnkt];
/*         ybins_true[0]=-2.0;
         ybins_true[1]=-1.0;
         ybins_true[2]=-0.75;
         ybins_true[3]=-0.5;
         ybins_true[4]=-0.25;
         ybins_true[5]=0.0;
         ybins_true[6]=0.25;
         ybins_true[7]=0.5;
         ybins_true[8]=0.75;
         ybins_true[9]=1.0;
         ybins_true[10]=1.25;
         ybins_true[11]=1.5;
         ybins_true[12]=1.75;
         ybins_true[13]=2.0;
         ybins_true[14]=2.25;
         ybins_true[15]=2.5;
         ybins_true[16]=2.75;
         ybins_true[17]=3.0;
         ybins_true[18]=3.25;
         ybins_true[19]=3.5;
         ybins_true[20]=3.75;
         ybins_true[21]=4.0;
         ybins_true[22]=4.5;
         ybins_true[23]=6.;
         ybins_true[24]=10.; */

      if (algo == "ak4")
      {
         ybins_true[0]=-6.0;
         ybins_true[1]=log(0.4);
         ybins_true[2]=log(0.4)+0.5;
         ybins_true[3]=log(0.4)+0.5*2.;
         ybins_true[4]=log(0.4)+0.5*3.;
         ybins_true[5]=log(0.4)+0.5*4.;
         ybins_true[6]=log(0.4)+0.5*5.;
         ybins_true[7]=log(0.4)+0.5*6.;
         ybins_true[8]=log(0.4)+0.5*7.;
         ybins_true[9]=log(0.4)+0.5*8.;
         ybins_true[10]=log(0.4)+0.5*9.;
         ybins_true[11]=log(0.4)+0.5*10.;
         ybins_true[12]=log(0.4)+0.5*11.;
         ybins_true[13]=6.0;
         ybins_true[14]=10.;
      }


      if (algo == "ak8")
      {
         ybins_true[0]=-6.0;
         ybins_true[1]=log(0.4)+0.5*2.;
         ybins_true[2]=log(0.4)+0.5*3.;
         ybins_true[3]=log(0.4)+0.5*4.;
         ybins_true[4]=log(0.4)+0.5*5.;
         ybins_true[5]=log(0.4)+0.5*6.;
         ybins_true[6]=log(0.4)+0.5*7.;
         ybins_true[7]=log(0.4)+0.5*8.;
         ybins_true[8]=log(0.4)+0.5*9.;
         ybins_true[9]=log(0.4)+0.5*10.;
         ybins_true[10]=log(0.4)+0.5*11.;
         ybins_true[11]=log(0.4)+0.5*12.;
         ybins_true[12]=log(0.4)+0.5*15.;
         ybins_true[13]=10.;
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

//         xbins[10]=3.33333333;
        //true binning in Rg
         Double_t xbins_true[bin_true_lntheta];

      if (algo == "ak4")
      {
         xbins_true[0]=-1;
         xbins_true[1]=0.;
         xbins_true[2]=.33333333;
         xbins_true[3]=.66666667;
         xbins_true[4]=1.;
         xbins_true[5]=1.33333333;
         xbins_true[6]=1.66666667;
         xbins_true[7]=2.0;
         xbins_true[8]=2.33333333;
         xbins_true[9]=2.66666667;
         xbins_true[10]=3.0;
         xbins_true[11]=3.33333333;
         xbins_true[12]=3.66666667;
         xbins_true[13]=4.0;
         xbins_true[14]=4.33333333;
         xbins_true[15]=6.0;
      }

      if (algo == "ak8")
      {
         xbins_true[0]=-1;
         xbins_true[1]=0.;
         xbins_true[2]=.33333333;
         xbins_true[3]=.66666667;
         xbins_true[4]=1.;
         xbins_true[5]=1.33333333;
         xbins_true[6]=1.66666667;
         xbins_true[7]=2.0;
         xbins_true[8]=2.33333333;
         xbins_true[9]=2.66666667;
         xbins_true[10]=3.0;
         xbins_true[11]=3.33333333;
         xbins_true[12]=3.66666667;
         xbins_true[13]=4.0;
         xbins_true[14]=4.33333333;
         xbins_true[15]=4.66666667;
         xbins_true[16]=5.0;
         xbins_true[17]=6.0;
      }
       //  xbins_true[12]=3.;


         int bin_det_pTjet = 2;
         int bin_true_pTjet = 3;

        //true binning in Rg
         Double_t pT_true[bin_true_pTjet];
         pT_true[0]=600.;
         pT_true[1]=700.;
         pT_true[2]=4000.;

         Double_t pT_det[bin_det_pTjet];
         pT_det[0]=700.;
         pT_det[1]=4000.;


   TH1D *jet1Pt(0), *jet2Pt(0), *jetPtAsymmetry(0), *DeltaPhi(0), *jet1Rapidity(0), *jet2Rapidity(0), *DeltaRapidity(0), *nJets(0), *nSplittings(0);


   TH1D * residualJetPtMINIAOD = new TH1D("residualMINIAOD", "residualMINIAOD",50,-0.5,0.5);
   TH1D * residualJetPtNew = new TH1D("residualNew", "residualNew",50,-0.5,0.5); 
   TH1D * residualJetPtSmeared = new TH1D("residualSmeared", "residualSmeared",50,-0.5,0.5);
   TH1D * residualJetPtUp = new TH1D("residualUp", "residualUp",50,-0.5,0.5);
   TH1D * residualJetPtDown = new TH1D("residualDown", "residualDown",50,-0.5,0.5);
   TH1D * prefireWeights = new TH1D("prefireWeights", "prefireWeight" ,20, 0,1);
   TH1D * scaleFactors = new TH1D("scaleFactors", "scaleFactors" ,100, 0.9,1.1);
   TH1D * weights = new TH1D("weights", "weights" ,200, -12,0);
   TH1D * zTrue = new TH1D("zTrue", "zTrue" ,50, 0,0.5);
   TH1D * nVtxGood = new TH1D("nvtxGood", "nvtxGood" ,50, 0,50);

   nVtxGood->Sumw2();

   TH1D * corMINIAOD = new TH1D("corMINIAOD", "corMINIAOD",50,1,1.5);
   TH1D * corNew = new TH1D("corNew", "corNew",50,1,1.5);
   nJets = new TH1D("nJets", "nJets",11,-0.5,10.5);
   nSplittings = new TH1D("nSplittings", "nSplittings",16,-0.5,15.5);
   TH1D * gen_jet1Pt = new TH1D("genJet1Pt", "genJet1Pt",140,500,4000);
   TH1D * ptHat = new TH1D("ptHat", "ptHat",100,0,3000);
   jet1Pt = new TH1D("jet1Pt", "jet1Pt",132,700,4000);
   jet2Pt = new TH1D("jet2Pt", "jet2Pt",132,700,4000);
   jetPtAsymmetry = new TH1D("jetPtAsymmetry", "jetPtAsymmetry",20,0, 1);
   DeltaPhi = new TH1D("DeltaPhi", "DeltaPhi",32,0, M_PI);
   jet1Rapidity =  new TH1D("jet1Rapidity", "jet1Rapidity",20,-1.7, 1.7);
   jet2Rapidity =  new TH1D("jet2Rapidity", "jet2Rapidity",20,-1.7, 1.7);
   DeltaRapidity =  new TH1D("DeltaRapidity", "DeltaRapidity",20,0, 4.0);
   TH1D * residualKt =  new TH1D("residual_kT_perturbative", "residual_kT_perturbative",20,-1., 1.);

          
   int nbins = 14;
   int nbins_dR = 13;
   float logkTmax = 6.;
   float nbinsTrue = 10;
   float logkTmaxTrue = 10.0;
   float logkTmin = log(0.4);
   float logRmax = 4.33333333;
   float logRmaxTrue = 6.0 ;


   if (algo == "ak8")
   {
    logkTmin = lnktmin_det; logRmax = 5.000; logkTmax = 6.5+log(0.4)+2.*0.5;
    nbins_dR = 15; nbins = 14;
   }

   TH2D *h2pThatratio(0);
   h2pThatratio=new TH2D("pThatratio_vs_pTreco","pThatratio_vs_pTreco", 800, 0, 4000, 40, 0, 6);

   TH2D * h2weightVsPthat = new TH2D("weightVsPthat","weightVsPthat", 400, 0, 4000, 100, -14, 0);
   TH2D * h2weightVsPtreco = new TH2D("weightVsPtreco","weightVsPtreco", 400, 0, 4000, 100, -14, 0);
   TH2D * h2genPtvsPthat =  new TH2D("genPtVsPthat","genPtVsPthat", 400, 0, 4000, 400, 0, 4000); 
   TH2D * h2jetPtVsEtaPrefire = new TH2D("jetPtvsEtaPrefire","jetPtVsEtaPrefire", 10, -2., 2., 25, 550, 4000);
   TH2D * h2jetPtVsEta = new TH2D("jetPtvsEta","jetPtVsEta", 10, -2., 2., 25, 550, 4000);
   TH2D * h2prefireVsEta = new TH2D("prefireVsEta","prefireVsEta", 10, -2., 2., 20, 0., 1.);
   TH2D * h2prefireVsPt = new TH2D("prefireVsPt","prefireVsPt", 25, 550., 4000., 20, 0., 1.);
 
   //MC detector correlation
   TH3D *h2smeared(0);
   h2smeared=new TH3D("smeared","smeared",bin_det_lntheta-1, xbins,bin_det_lnkt-1,ybins, bin_det_pTjet-1, pT_det );

   TH3D *h2recoall(0);
   h2recoall=new TH3D("recoall","recoall",bin_det_lntheta-1, xbins,bin_det_lnkt-1,ybins,  bin_det_pTjet-1, pT_det );


   TH3D *h2purity(0);
   h2purity=new TH3D("purity","purity",bin_det_lntheta-1, xbins,bin_det_lnkt-1,ybins,  bin_det_pTjet-1, pT_det);

// matching efficiency

   TH3D *h2efficiency(0);
   h2efficiency=new TH3D("efficiency","efficiency",bin_true_lntheta-1, xbins_true, bin_true_lnkt-1,ybins_true, bin_true_pTjet-1, pT_true);


   //true correlations with det-level cuts
   TH3D *h2true(0);

    h2true=new TH3D("true","true",bin_true_lntheta-1, xbins_true,bin_true_lnkt-1,ybins_true, bin_true_pTjet-1, pT_true);
    TH3D *h2fulleff(0);

     h2fulleff=new TH3D("truef","truef",bin_true_lntheta-1, xbins_true,bin_true_lnkt-1,ybins_true, bin_true_pTjet-1, pT_true);
     TH3D *effnum=(TH3D*)h2fulleff->Clone("effnum");
     TH3D *effdenom=(TH3D*)h2fulleff->Clone("effdenom");

   TH2D *kT_response(0);

   kT_response=new TH2D("kT_response","kT_response",nbins, logkTmin, logkTmax, nbins, logkTmin, logkTmax);
   TH1D *kT_gen_spectrum(0);
   kT_gen_spectrum = new TH1D("kT_gen_spectrum", "kT_gen_spectrum", nbins, logkTmin, logkTmax);

   TH1D *theta_gen_spectrum(0);
   theta_gen_spectrum = new TH1D("theta_gen_spectrum", "theta_gen_spectrum",nbins_dR, 0, logRmax);

   TH2D *theta_response(0);
   theta_response=new TH2D("theta_response","theta_response",nbins_dR, 0, logRmax,nbins_dR, 0, logRmax);

   TH1D *pfRho_density(0);
   pfRho_density = new TH1D("pfRho_density", "pfRho_density", 60, 0, 60);

   TH1D *puTrue(0);
   puTrue = new TH1D("puTrue", "pfRho_density", 99, 0, 99);

   TH1D *pTjet_true(0);
   pTjet_true = new TH1D("pTjet_true", "pTjet_true", bin_true_pTjet-1, pT_true);

   TH1D *pTjet_det(0);
   pTjet_det = new TH1D("pTjet_det", "pTjet_det", bin_det_pTjet-1, pT_det);


 

   effdenom->Sumw2();
   h2smeared->Sumw2();
   h2true->Sumw2();
  
   h2fulleff->Sumw2();
   //setup the binnings of the 3D response
    RooUnfoldResponse response;
//    response.Setup(h2true,h2true);
    response.Setup(h2smeared,h2true);

    RooUnfoldResponse pT_response;
    pT_response.Setup(pTjet_det, pTjet_true);

 
   //******************Fill the MCresponse*******************//


    ifstream infile;
    infile.open(cFiles2);
    char filename[300];
    int count=1;
    int n_splittings = 0, n_jets = 0;
    while(infile>>filename)
    {

      TFile *input2=TFile::Open(filename);
      cout << filename << endl;
      count=count+1;
//      cout<<"reading data file "<<count<<endl;
      if(count==10) break;

      TDirectory *tdir = input2->GetDirectory("ak4");

      TTree *T = (TTree*)tdir->Get("ProcessedTree");


      QCDEvent        *events;
      TBranch        *b_events;   //!

      // Set object pointer
      events = 0;
      // Set branch addresses and branch pointers
      if (!T) return;
      T->SetMakeClass(1);
//      double logkTmin = -1;
 
      Int_t nEvmc=T->GetEntries();
      T->SetBranchAddress("events", &events, &b_events); 
      double factor = 0.;
      if (jec_option == "nom") factor  = 0.;
      if (jec_option == "up") factor  = 1.;
      if (jec_option == "down") factor  = -1.;

      double R = 0.4;

      if (algo == "ak8") R = 0.8;

   for(int iEntry=0; iEntry< nEvmc; iEntry++)
   {
    T->GetEntry(iEntry);

   if ( events->nPFJetsCHS() < 1) continue;
//    if (jet1.looseID() == false ) continue;
    auto eventHdr = events->evtHdr();
    double weight =  eventHdr.weight();
    auto psWeight =  eventHdr.psWeight();
    double pfRho = eventHdr.pfRho();
    double PU = eventHdr.trpu();

    if (eventHdr.isPVgood() == false) continue;

 
    if ( PU <= 70)
    {
     weight = weight;
     weight = weight*puReweightNominal->GetBinContent(PU+1);
    }

    if ( PU > 70 ) continue;

//    cout << "psWeight: " << psWeight.at(4)/weight << endl;
    weight = weight*lumiWeight*psWeight.at(4)/psWeight.at(1);
    puTrue->Fill(eventHdr.trpu(),weight);
//    if (eventHdr.pthat() < 300 ) continue;
//    weights->Fill(log10(weight));
//    ptHat->Fill(eventHdr.pthat(), weight );
    auto jet1 = events->pfjetchs(0);
//       jet2 = events->pfjetchs(1);

      JetCorrector->setJetEta(jet1.eta());
      JetCorrector->setJetPt(jet1.pt()/jet1.cor());
      JetCorrector->setRho(pfRho );
      JetCorrector->setJetA(jet1.area());
      QCDGenJet jet1_gen, jet2_gen;
      QCDPFJet jet2;

    double HEM1516 = rand->Rndm();
    bool HEM1516_era = HEM1516 < 0.65;

      if ( jet1.genidx() != -1 ) jet1_gen = events->genjet(jet1.genidx());


      double newCorrection = JetCorrector->getCorrection();

      jecUnc->setJetEta(jet1.eta());
      jecUnc->setJetPt( jet1.pt()/jet1.cor()*newCorrection);
      double unc = jecUnc->getUncertainty(true);

      newCorrection = newCorrection*(1+factor*unc);
//      cout << newCorrection << endl;

      JME::JetParameters parameters_1;
      parameters_1.setJetPt(jet1.pt()/jet1.cor()*newCorrection);
      parameters_1.setJetEta(jet1.eta());
      parameters_1.setRho(pfRho);

      float sf = 1.;

      if (jer_option == "nom") sf = resolution_sf.getScaleFactor(parameters_1); 
      if (jer_option == "up") sf = resolution_sf.getScaleFactor(parameters_1, Variation::UP);
      if (jer_option == "down") sf = resolution_sf.getScaleFactor(parameters_1, Variation::DOWN);

      float sim_resolution = resolution.getResolution(parameters_1);
      if (jet1.genidx() != -1 ) sf = 1+(sf-1)*(jet1.pt()/jet1.cor()*newCorrection-jet1_gen.pt())/ ( jet1.pt()/jet1.cor()*newCorrection );
      else sf = 1+gRandom->Gaus(0,sim_resolution)*std::sqrt(sf*sf-1);

     double jet1_pT_corrected = jet1.pt()/jet1.cor()*newCorrection*sf;


//    if (jet1.eta() < -1.3 && jet1.phi() > -1.57 && jet1.phi() < -0.87 && jet1.eta() > -2.5 && HEM1516_era ) jet1_pT_corrected = jet1_pT_corrected*0.8;
    if ( !(jet1_pT_corrected > 700 ) ) continue;
    if ( jet1.looseID() == false ) continue;
    if ( jet1_pT_corrected/eventHdr.pthat() > 2. ) continue;
//    if ( jet1_pT_corrected/eventHdr.pthat() > 1.5 ) continue;
    if ( eventHdr.pthat() < 450.0 ) continue;
//      weight = 1.;


      float sfUp = resolution_sf.getScaleFactor(parameters_1, Variation::UP);
      if (jet1.genidx() != -1 )  sfUp = 1+(sfUp-1)*(jet1.pt()/jet1.cor()*newCorrection-jet1_gen.pt())/( jet1.pt()/jet1.cor()*newCorrection );
      else sfUp = 1+gRandom->Gaus(0,sim_resolution)*std::sqrt(sfUp*sfUp-1);

      float sfDown = resolution_sf.getScaleFactor(parameters_1, Variation::DOWN);
      if (jet1.genidx() != -1 ) sfDown = 1+(sfDown-1)*(jet1.pt()/jet1.cor()*newCorrection-jet1_gen.pt())/( jet1.pt()/jet1.cor()*newCorrection ) ;
      else sfDown = 1+gRandom->Gaus(0,sim_resolution)*std::sqrt(sfDown*sfDown-1);


/*     if ( jet1.genidx() == -1 ) continue;
     if (jet1_gen.pt() < 500 || fabs(jet1_gen.eta() > 1.3 )) continue;
     residualJetPtMINIAOD->Fill( (jet1.pt()-jet1_gen.pt() )/jet1_gen.pt() , weight) ;
     residualJetPtNew->Fill( (jet1_pT_corrected/sf-jet1_gen.pt() )/jet1_gen.pt() , weight) ;
     residualJetPtSmeared->Fill( (jet1_pT_corrected -jet1_gen.pt() )/jet1_gen.pt(), weight ) ;
     residualJetPtUp->Fill((jet1_pT_corrected/sf*sfUp -jet1_gen.pt() )/jet1_gen.pt(), weight ) ;
     residualJetPtDown->Fill( (jet1_pT_corrected/sf*sfDown -jet1_gen.pt() )/jet1_gen.pt()  , weight );*/

     double jet2_pT_corrected;





// prefire weight

   double prefireWeight = 1.;
   for (int k = 0; k < events->nPFJetsCHS(); k++) // det-level jets
   {
     bool mismatch = false;
     auto jet_reco = events->pfjetchs(k);


     if ( jet_reco.genidx() == -1 ) continue; //if det-level jet is not matched, skip event


     auto jet_gen = events->genjet(jet_reco.genidx());

      JetCorrectorRecoLoop->setJetEta(jet_reco.eta());
      JetCorrectorRecoLoop->setJetPt(jet_reco.pt()/jet_reco.cor());
      JetCorrectorRecoLoop->setRho(pfRho );
      JetCorrectorRecoLoop->setJetA(jet_reco.area());

      newCorrection = JetCorrectorRecoLoop->getCorrection();

      jecUnc->setJetEta(jet_reco.eta() );
      jecUnc->setJetPt( jet_reco.pt()/jet_reco.cor()*newCorrection);
      unc = jecUnc->getUncertainty(true);


      newCorrection = newCorrection*(1+factor*unc);

      parameters_1.setJetPt(jet_reco.pt()/jet_reco.cor()*newCorrection);
      parameters_1.setJetEta(jet_reco.eta());

//      sf = resolution_sf.getScaleFactor(parameters_1);

      if (jer_option == "nom") sf = resolution_sf.getScaleFactor(parameters_1);               
      if (jer_option == "up") sf = resolution_sf.getScaleFactor(parameters_1, Variation::UP);
      if (jer_option == "down") sf = resolution_sf.getScaleFactor(parameters_1, Variation::DOWN);

      sf = 1+(sf-1)*(jet_reco.pt()/jet_reco.cor()*newCorrection-jet1_gen.pt())/( jet_reco.pt()/jet1.cor()*newCorrection );
      double jetPtCorrected = jet_reco.pt()/jet_reco.cor()*newCorrection*sf; //rescaling pT with jet energy resolution, new jet energy scale

      if  (!( jetPtCorrected > 20 && fabs(jet_reco.eta()) > 2. && fabs(jet_reco.eta()) < 3.0  ) ) continue;
      prefireWeight *= 1-prefireMap->GetBinContent(prefireMap->FindBin(jet_reco.eta(),jetPtCorrected));

   }
//TODO
//     prefireWeight = 1.;
   if (year == 2018) prefireWeight = 0.999999;
//   if (prefireWeight < 0.5) cout << prefireWeight << endl;

     prefireWeights->Fill(prefireWeight,weight);

     weight = prefireWeight*weight;

     if ( events->nPFJetsCHS() >= 2 )
     {
     jet2 = events->pfjetchs(1);
     JetCorrector->setJetEta(jet2.eta());
     JetCorrector->setJetPt(jet2.pt()/jet2.cor());
     JetCorrector->setJetA(jet2.area());
     JetCorrector->setRho(pfRho);

      newCorrection = JetCorrector->getCorrection(); 

      jecUnc->setJetEta(jet1.eta());
      jecUnc->setJetPt( jet1.pt()/jet1.cor()*newCorrection);
      unc = jecUnc->getUncertainty(true);

      newCorrection = newCorrection*(1+factor*unc);

      parameters_1.setJetPt(jet2.pt()/jet2.cor()*newCorrection);
      parameters_1.setJetEta(jet2.eta());
      parameters_1.setRho(pfRho);

     if ( jet2.genidx() != -1 ) jet2_gen = events->genjet(jet2.genidx());



 //     sf = resolution_sf.getScaleFactor(parameters_1);

      if (jer_option == "nom") sf = resolution_sf.getScaleFactor(parameters_1);               
      if (jer_option == "up") sf = resolution_sf.getScaleFactor(parameters_1, Variation::UP);
      if (jer_option == "down") sf = resolution_sf.getScaleFactor(parameters_1, Variation::DOWN);

      sim_resolution = resolution.getResolution(parameters_1);
      if (jet2.genidx() != -1 ) sf = 1+(sf-1)*(jet2.pt()/jet2.cor()*newCorrection-jet2_gen.pt())/ ( jet2.pt()/jet2.cor()*newCorrection );
      else sf = 1+gRandom->Gaus(0,sim_resolution)*std::sqrt(sf*sf-1);

      jet2_pT_corrected = jet2.pt()/jet2.cor()*newCorrection*sf;

//      cout << jet2_pT_corrected << endl;


     if ( jet1_pT_corrected > 700. && jet2_pT_corrected > 700. && fabs(jet1.y()) < 1.7 && fabs(jet2.y()) < 1.7 && jet1_pT_corrected < 4000 && jet2_pT_corrected < 4000 && jet1.looseID() == true && jet2.looseID() == true )
      {
      double delta_phi = 0.;
      if (fabs(jet1.phi()-jet2.phi()) < M_PI) delta_phi = fabs(jet1.phi()-jet2.phi());
      else delta_phi = 2*M_PI-fabs(jet1.phi()-jet2.phi());
      DeltaPhi->Fill(delta_phi, weight);
//      jet2Rapidity->Fill(jet2.rapidity(), weight);
      DeltaRapidity->Fill(fabs(jet1.y()-jet2.y()),weight);
      double pTasymmetry = fabs(jet1_pT_corrected-jet2_pT_corrected)/(jet1_pT_corrected+jet2_pT_corrected );
      jetPtAsymmetry->Fill(pTasymmetry,weight);

      }
      }


//      weight = weight*prefireWeight;

//    if (jet1.looseID() == false || jet2.looseID() == false ) continue;
//    if (!( jet1_pT_corrected > 550 && jet2_pT_corrected > 550 && fabs(jet1.eta()) < 2 && fabs(jet2.eta()) < 2 && jet1_pT_corrected < 4000 && jet2_pT_corrected < 4000 )) continue;
//    if (weight > 1.69824365E-7 ) continue;


//    h2pThatratio->Fill(jet1_pT_corrected, jet1_pT_corrected/eventHdr.pthat());


    h2pThatratio->Fill(jet1_pT_corrected, jet1_pT_corrected/eventHdr.pthat());
    h2genPtvsPthat->Fill(jet1_gen.pt(), eventHdr.pthat());
    h2weightVsPthat->Fill(eventHdr.pthat() , log10(weight));
    ptHat->Fill(eventHdr.pthat(), weight );
//    if (log10(weight) > -5.8 || log10(weight) < -10  ) continue;


//    h2pThatratio->Fill(jet1_pT_corrected, jet1_pT_corrected/eventHdr.pthat(), weight);

//    ptHat->Fill(eventHdr.pthat(), weight );
    nVtxGood->Fill(eventHdr.nVtxGood(),weight);
//    weights->Fill(log10(weight));
/*     residualJetPtMINIAOD->Fill( (jet1.pt()-jet1_gen.pt() )/jet1_gen.pt() , weight) ;
     residualJetPtNew->Fill( (jet1_pT_corrected/sf-jet1_gen.pt() )/jet1_gen.pt() , weight) ;
     residualJetPtSmeared->Fill( (jet1_pT_corrected -jet1_gen.pt() )/jet1_gen.pt(), weight ) ;*/

    pfRho_density->Fill(eventHdr.pfRho(),weight);

//    jet1Pt->Fill(jet1_pT_corrected , weight); //, *jet2Pt(0), *jetPtAsymmetry, *DeltaPhi(0), *jet1Rapidity(0), *jet2Rapidity(0), *DeltaRapidity(0);
//    jet2Pt->Fill(jet2_pT_corrected , weight);
//

//    if (jet1.looseID() == false || jet2.looseID() == false ) continue;



    n_jets = 0;


    double split = rand->Rndm();
    bool splitPass = split < 0.30;
    bool splitDoesNotPass = split > 0.30;

    splitPass = true; splitDoesNotPass = true;
    //weight = 1.;

   for (int k = 0; k < events->nPFJetsCHS(); k++) // det-level jets
   {
     bool mismatch = false;
     auto jet_reco = events->pfjetchs(k);

     auto z_reco = jet_reco.z();
     auto theta_reco = jet_reco.theta();
     auto kT_reco = jet_reco.kT();
     auto eta2_reco = jet_reco.eta2_splitting();
     auto phi2_reco = jet_reco.phi2_splitting();

     if ( jet_reco.genidx() == -1 ) continue; //if det-level jet is not matched, skip event
  
     auto jet_gen = events->genjet(jet_reco.genidx());
     auto z_gen = jet_gen.z();
     auto theta_gen = jet_gen.theta();
     auto kT_gen = jet_gen.kT();
     auto eta2_gen = jet_gen.eta2_splitting();
     auto phi2_gen = jet_gen.phi2_splitting();

//     cout << hotMap->GetBinContent(hotMap->FindBin(jet_reco.eta(),jet_reco.phi())) << endl;

     if ( hotMap->GetBinContent(hotMap->FindBin(jet_reco.eta(),jet_reco.phi())) > 5 ) continue;

      JetCorrectorRecoLoop->setJetEta(jet_reco.eta());
      JetCorrectorRecoLoop->setJetPt(jet_reco.pt()/jet_reco.cor());
      JetCorrectorRecoLoop->setRho(pfRho );
      JetCorrectorRecoLoop->setJetA(jet_reco.area());

      newCorrection = JetCorrectorRecoLoop->getCorrection();

      jecUnc->setJetEta(jet_reco.eta() );
      jecUnc->setJetPt( jet_reco.pt()/jet_reco.cor()*newCorrection);
      unc = jecUnc->getUncertainty(true);


      newCorrection = newCorrection*(1+factor*unc);

      parameters_1.setJetPt(jet_reco.pt()/jet_reco.cor()*newCorrection);
      parameters_1.setJetEta(jet_reco.eta());

//      sf = resolution_sf.getScaleFactor(parameters_1);

      if (jer_option == "nom") sf = resolution_sf.getScaleFactor(parameters_1);               
      if (jer_option == "up") sf = resolution_sf.getScaleFactor(parameters_1, Variation::UP);
      if (jer_option == "down") sf = resolution_sf.getScaleFactor(parameters_1, Variation::DOWN);

      sf = 1+(sf-1)*(jet_reco.pt()/jet_reco.cor()*newCorrection-jet1_gen.pt())/( jet_reco.pt()/jet1.cor()*newCorrection ); 
      double jetPtCorrected = jet_reco.pt()/jet_reco.cor()*newCorrection*sf; //rescaling pT with jet energy resolution, new jet energy scale

      if (log10(weight) > -5) continue;
//      if (jet_reco.eta() < -1.3 && jet_reco.phi() > -1.57 && jet_reco.phi() < -0.87 && jet_reco.eta() > -2.5 && HEM1516_era ) jetPtCorrected = jetPtCorrected*0.8;


//      if  (!(z_gen.size()> 0 && jet_gen.pt() > 500 && fabs(jet_gen.y()) < 2.0 && z_reco.size() > 0 && jetPtCorrected > 550 && fabs(jet_reco.y()) < 2 && jet_reco.looseID() == true && jet_gen.pt() < 4000. && jetPtCorrected < 4000 && fabs(jet_reco.eta()) < 2.0  ) ) continue;

      if  (!( jet_gen.pt() > 600 && jetPtCorrected > 700 && fabs(jet_reco.y()) < 1.7 && jet_reco.looseID() == true && jet_gen.pt() < 4000. && jetPtCorrected < 4000 ) ) continue;
//      if  (!( jet_gen.pt() > 500 && fabs(jet_gen.y()) < 2.0 && jetPtCorrected > 550 && fabs(jet_reco.y()) < 2 && jet_reco.looseID() == true && jet_gen.pt() < 4000. && jetPtCorrected < 4000 && fabs(jet_reco.eta()) < 2.0  ) ) continue;
//      if  (!(   z_reco.size() > 0 && jetPtCorrected > 550 && fabs(jet_reco.y()) < 2 && jet_reco.looseID() == true && jetPtCorrected < 4000 && fabs(jet_reco.eta()) < 2.0  ) ) continue;
//
//      cout << jet_gen.pt() << endl;
      weights->Fill(log10(weight));
      h2weightVsPtreco->Fill(jetPtCorrected, log10(weight));
      jet1Pt->Fill(jetPtCorrected , weight);
      jet1Rapidity->Fill(jet_reco.eta(), weight);
      gen_jet1Pt->Fill(jet_gen.pt(),weight);
//      cout << jetPtCorrected << endl;

      int first_splitting = 0;
      n_jets++;
      n_splittings = 0;
     if (splitDoesNotPass == true) pT_response.Fill(jetPtCorrected , jet_gen.pt(),  weight);
     if (splitPass == true) pTjet_det->Fill(jetPtCorrected, weight);
     if (splitDoesNotPass == true) pTjet_true->Fill(jet_gen.pt(),weight);

     h2jetPtVsEtaPrefire->Fill(jet_reco.eta(), jetPtCorrected, weight*prefireWeight);
     h2jetPtVsEta->Fill(jet_reco.eta(), jetPtCorrected, weight);

     h2prefireVsPt->Fill(jetPtCorrected, prefireWeight, weight);
     h2prefireVsEta->Fill(jet_reco.eta(), prefireWeight, weight);

     scaleFactors->Fill(sf,weight);

     for (unsigned i = 0; i < z_gen.size(); ++i) // loop over gen-level splittings
     {
//     cout  << z_gen.at(i) << endl;
     if (mismatch == true) cout << "=======================================" << endl;
     if (mismatch == true) cout << "jet" << "\t" << jet_gen.pt() << "\t" << jet_gen.phi() << "\t" << jet_gen.eta() << endl;
//      cout << "=================== " << i << endl;
//      cout << "gen-level splitting " << i << endl;
 //     cout << "=================== " << i << endl;

      if (log(kT_gen.at(i)) < lnktmin_true || log(kT_gen.at(i)) > logkTmaxTrue) continue; //gen-level phase-space
      if (log(R/theta_gen.at(i)) < -1. || log(R/theta_gen.at(i)) > logRmaxTrue ) continue; //gen-level phase-space
      if (splitDoesNotPass )  h2fulleff->Fill(log(R/theta_gen.at(i)) , log(kT_gen.at(i)), jet_gen.pt(), weight );
      n_splittings++;
//       histosTH2F["Lund_plane_gen_all"]->Fill(log(0.4/theta_gen.at(i)),log(kT_gen.at(i) ));

      first_splitting = first_splitting+1;
      float dR_window = 0.1;
      float dR_max = dR_window;
      int index_true = -1;
      int index_reco = -1;
      float kT_min = 0.5;
      for (unsigned j = 0; j < z_reco.size(); ++j) //loop over reco-level splittings
      {
//       cout  << z_reco.at(j) << endl;
       if (log(kT_reco.at(j)) < logkTmin || log(kT_reco.at(j)) > logkTmax ) continue;
       if (log(R/theta_reco.at(j)) < 0. || log(R/theta_reco.at(j)) > logRmax ) continue;
      // cout << " first splitting " << first_splitting << " j " << j << endl;
       if (mismatch == true) cout << "reco " << log(R/theta_reco.at(j)) << "\t" << log(kT_reco.at(j)) << "\t" << phi2_reco.at(j) << "\t" << eta2_reco.at(j) << "\t" << endl;      
       if (first_splitting == 1 && splitPass == true )
       {
        h2recoall->Fill(log(R/theta_reco.at(j)),log(kT_reco.at(j) ), jetPtCorrected ,  weight);

       }

       float dphi = fabs(phi2_reco.at(j) - phi2_gen.at(i));
       float deta = eta2_reco.at(j) - eta2_gen.at(i);
       if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
       float dR = std::sqrt(dphi*dphi + deta*deta);
       float kT_ratio = min(kT_reco.at(j),kT_gen.at(i))/max(kT_reco.at(j),kT_gen.at(i));
       float pT_ratio = ( kT_reco.at(j)/sin(theta_reco.at(j)) ) / ( kT_gen.at(i)/sin(theta_gen.at(i)) );
       if (dR < dR_max) //find best pair
       {
	 dR_max = dR;
         index_true = i;
         index_reco = j;
         kT_min = kT_ratio;
       }//endif


      }//end for

       if ( index_reco == -1 || index_true == -1 ) continue;
       int index_true_mc = -1;
       dR_max = dR_window;
       kT_min = 0.5;
       float dR_final, dEta_final, dPhi_final;
       for (unsigned l = 0; l < z_gen.size(); ++l) // loop over gen-level splitings, check if the previous match is true the other way around
       {

         if (log(kT_gen.at(l)) < lnktmin_true || log(kT_gen.at(l)) > logkTmaxTrue ) continue;
         if (log(R/theta_gen.at(l)) < -1. || log(R/theta_gen.at(l)) > logRmaxTrue ) continue;
         float dphi = fabs(phi2_reco.at(index_reco) - phi2_gen.at(l));
         float deta = eta2_reco.at(index_reco) - eta2_gen.at(l);
         if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
         float dR = std::sqrt(dphi*dphi + deta*deta);
         float kT_ratio = min(kT_reco.at(index_reco),kT_gen.at(l))/max(kT_reco.at(index_reco),kT_gen.at(l));
         float pT_ratio = ( kT_reco.at(index_reco)/sin(theta_reco.at(index_reco)) ) / ( kT_gen.at(i)/sin(theta_gen.at(l)) );
         if (mismatch == true)  cout << "gen " << log(R/theta_gen.at(l)) << "\t" << log(kT_gen.at(l)) << "\t" << phi2_gen.at(l) << "\t" << eta2_gen.at(l) << "\t" << endl; 
         if (dR < dR_max ) //find best pair
         {
           dR_max = dR;
           index_true_mc = l;
           kT_min = kT_ratio;
         }//endif
       }//end for
      if ( index_true != index_true_mc ) continue;

      if (splitPass == true  ) h2smeared->Fill(log(R/theta_reco.at(index_reco)) , log(kT_reco.at(index_reco) )  , jetPtCorrected , weight);
      if (splitDoesNotPass == true  )
      {
       h2true->Fill(log(R/theta_gen.at(index_true)) , log(kT_gen.at(index_true)), jet_gen.pt(), weight);
       theta_response->Fill(log(R/theta_gen.at(index_true)), log(R/theta_reco.at(index_reco)), weight );
       kT_response->Fill(log(kT_gen.at(index_true)) , log(kT_reco.at(index_reco) ) , weight );
       kT_gen_spectrum->Fill(log(kT_gen.at(index_true)), weight);
       theta_gen_spectrum->Fill( log(R/theta_gen.at(index_true)) , weight);
       response.Fill(log(R/theta_reco.at(index_reco)) , log(kT_reco.at(index_reco) )  , jetPtCorrected,  log(R/theta_gen.at(index_true)) , log(kT_gen.at(index_true) ) , jet_gen.pt(),  weight);
       if ( log(kT_gen.at(index_true)) > 3 ) residualKt->Fill( (kT_reco.at(index_reco) -kT_gen.at(index_true))/kT_gen.at(index_true)  ,weight);
//       if ( z_gen.size() == z_reco.size() ) mismatch = false;
       if ( log(kT_gen.at(index_true)) > 3 && (kT_reco.at(index_reco)/theta_reco.at(index_reco)) / ( kT_gen.at(index_true)/theta_gen.at(index_true) ) < 0.2 ) {zTrue->Fill(z_gen.at(index_true),weight ) ; /* mismatch = true;*/}
      }
//     histosTH2F["Lund_plane_reco_matched"]->Fill(log(0.4/theta_reco.at(index_reco)),log(kT_reco.at(index_reco) ));
//     histosTH2F["Lund_plane_gen_matched"]->Fill(log(0.4/theta_gen.at(index_true)),log(kT_gen.at(index_true) ));
      }
     }
    nJets->Fill(n_jets,weight);
    nSplittings->Fill(n_splittings,weight);
   }
  }
   
// TFile *fout = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/herwig7_2018.root","RECREATE");
//    TFile *fout=new TFile (Form("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/ResponseHerwig7_%s_%i_%s_pT500toInf_split30.root",tag.c_str(),fileindex,date.c_str()),"RECREATE");
    TFile *fout=new TFile (Form("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/MiniAODV2_FSRdown%s%s_%i_JEC%s_JER%s_PU%s_%s_%i.root", tag.c_str(), mc_name.c_str(), year,jec_option.c_str(),jer_option.c_str(), pu_option.c_str(), algo.c_str(), fileindex ),"RECREATE");
  fout->cd();

  for (int i = 1; i < nbins+1; i++) // normalization for response matrix
  {

   for (int j = 1; j < nbins+1; j++)
   {
    kT_response->SetBinContent(i,j, float(kT_response->GetBinContent(i,j))/float(kT_gen_spectrum->GetBinContent(i)) );
   }

  }

  for (int i = 1; i < nbins_dR+1; i++) // normalization for response matrix
  {
   for (int j = 1; j < nbins_dR+1; j++)
   {
    theta_response->SetBinContent(i,j, float(theta_response->GetBinContent(i,j))/float(theta_gen_spectrum->GetBinContent(i)));
   }
  }
  residualJetPtMINIAOD->Write(); residualJetPtNew->Write(); residualJetPtSmeared->Write();  residualJetPtUp->Write(); residualJetPtDown->Write();
  corMINIAOD->Write(); corNew->Write(); residualKt->Write();
  pTjet_true->Write(); pTjet_det->Write();
  nSplittings->Write();
  weights->Write(); nVtxGood->Write();
  nJets->Write();
  pfRho_density->Write(); puTrue->Write();
  h2jetPtVsEtaPrefire->Divide(h2jetPtVsEta);
  h2jetPtVsEtaPrefire->Write();
  h2prefireVsPt->Write(); h2prefireVsEta->Write();
  h2purity->Divide(h2recoall,h2smeared);
  h2efficiency->Divide(h2true,h2fulleff);
  h2smeared->SetName("smeared");
  h2smeared->Write();
  h2purity->SetName("purity");
  h2purity->Write();
  h2efficiency->SetName("efficiency");
  h2efficiency->Write();
  h2fulleff->Write();
  h2true->Write();
  h2recoall->Write();
  zTrue->Write();
  h2pThatratio->Write(); h2weightVsPthat->Write(); h2genPtvsPthat->Write();  h2weightVsPtreco->Write();
  ptHat->Write(); gen_jet1Pt->Write();   jet1Pt->Write(); jet2Pt->Write(); jetPtAsymmetry->Write(); DeltaPhi->Write(); jet1Rapidity->Write(); jet2Rapidity->Write(); DeltaRapidity->Write();
  theta_response->Write(); kT_response->Write();
  response.Write();
  pT_response.Write();
  prefireWeights->Write();
  scaleFactors->Write();
}

#ifndef __CINT__
int main () { RooUnfoldReadResponse_3D(); return 0; }  // Main program when run stand-alone
#endif
