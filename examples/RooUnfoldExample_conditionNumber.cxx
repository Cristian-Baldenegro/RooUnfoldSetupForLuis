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
 *ReadResponseRun2_JECdown.root
ReadResponseRun2_JECup.root
ReadResponseRun2_JERdown.root
ReadResponseRun2_JERup.root
ReadResponseRun2_PUdown.root
ReadResponseRun2_PUup.root
ReadResponseRun2_TrackIneff.root
ReadResponseRun2_nom.root

ReadResponseRun2V2_HEM1516_ak4.root
ReadResponseRun2V2_JECdown_ak4.root
ReadResponseRun2V2_JECup_ak4.root
ReadResponseRun2V2_JERdown_ak4.root
ReadResponseRun2V2_JERup_ak4.root
ReadResponseRun2V2_PUdown_ak4.root
ReadResponseRun2V2_PUup_ak4.root
ReadResponseRun2V2_herwig7_ak4.root
ReadResponseRun2V2_nom_ak4.root
ReadResponseRun2V2_trackIneff_ak4.root

for 2016:

ReadResponseLumiWeightTrackIneffV2pythia82016_trackIneff_ak4.root



Redone:ReadResponseLumiWeightTrackIneffV2pythia8_Run2_ak4.root
 *
 *
 * */

void RooUnfoldExample_conditionNumber(std::string file_mc="/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/ConditionNumberV2_ak8.root", 
std::string file_data="/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/ReadData2016.root", std::string date = "herwig7V32016_ak4",int flag=19)
{
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
  TFile* f_mc=new TFile(file_mc.c_str());
  TFile* f_data=new TFile(file_data.c_str());

  RooUnfoldResponse* response = (RooUnfoldResponse*)f_mc->Get("true_true");

//   h2input->Multiply(h2purity);
   //do the unfolding for N iterations

//   TMatrixT<double> covariance();

      TDecompSVD *svd= new TDecompSVD (response->Mresponse()); //this is the singular value decomposition (SVD) matrix
      auto singular_values = svd->GetSig(); //this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
      svd->Print();
      double singular_value_min;
      for (int i = 0; i < 50 ; i++)
      {
        cout << singular_values[i] << endl;
       if ( singular_values[i] > pow(10,-15)  ) singular_value_min = singular_values[i]; // the pow(10,-15) > requirement is to suppress singular values from empty bins 
      }


      cout << "condition number: "  << singular_values[0]/singular_value_min << endl;



}

#ifndef __CINT__
int main () { RooUnfoldExample_conditionNumber(); return 0; }  // Main program when run stand-alone
#endif
