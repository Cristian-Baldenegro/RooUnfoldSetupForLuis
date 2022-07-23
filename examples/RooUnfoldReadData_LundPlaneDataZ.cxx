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


void RooUnfoldReadData_LundPlaneDataZ(std::string cFiles2 = "",std::string tag = "", std::string date = "", int fileindex=0)
{


  //errors from the diagonal elements of the covariance matrix
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
 
   const int bin_det_lnkt=6;
   const int bin_true_lnkt=8;
   const int bin_det_lntheta=12;
   const int bin_true_lntheta=12;

   Double_t lnthetamin_det=0.0;
   Double_t lnthetamin_true=0;

   Double_t lnthetamax_det=2.6;
   Double_t lnthetamax_true=3;

   Double_t lnktmax_det=0.2;
   Double_t lnktmax_true=0.3;

   Double_t lnktmin_det=-1.;
   Double_t lnktmin_true=-2.;

   /***************************************************/

   //det binning in xj
       //true binning in xj
        Double_t ybins_true[bin_true_lnkt];
         ybins_true[0]=0.69;
         ybins_true[1]=1.8;
//         ybins_true[1]=lnktmin_det;
         ybins_true[2]=2.08;
         ybins_true[3]=3.0;
         ybins_true[4]=4.0;
         ybins_true[5]=5.13;
         ybins_true[6]=5.41;
         ybins_true[7]=6.0;

          //det binning in Rg
         Double_t xbins[bin_det_lntheta];
         xbins[0]=0;
         xbins[1]=.33;
         xbins[2]=.67;
         xbins[3]=1.;
         xbins[4]=1.33;
         xbins[5]=1.67;
         xbins[6]=2.0;
         xbins[7]=2.33;
         xbins[8]=2.67;
         xbins[9]=3.0;
         xbins[10]=3.33;
         xbins[11]=3.67;




        //true binning in Rg
         Double_t xbins_true[bin_true_lntheta];
         xbins_true[0]=0;
         xbins_true[1]=.33;
         xbins_true[2]=.67;
         xbins_true[3]=1.;
         xbins_true[4]=1.33;
         xbins_true[5]=1.67;
         xbins_true[6]=2.0;
         xbins_true[7]=2.33;
         xbins_true[8]=2.67;
         xbins_true[9]=3.0;
         xbins_true[10]=3.33;
         xbins_true[11]=3.67;
       //  xbins_true[12]=3.;

   //the raw correlation
   TH2D *h2raw(0), *h2rawTrueBinning(0);
   h2raw=new TH2D("raw","raw",bin_true_lntheta-1,xbins_true,bin_true_lnkt-1,ybins_true);
   h2rawTrueBinning=new TH2D("rawTrueLevelBinning","rawTrueLevelBinning",bin_true_lntheta-1,xbins_true,bin_true_lnkt-1,ybins_true);
   h2raw->Sumw2();
   h2raw->Sumw2();


      TFile *input2=TFile::Open("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/DATA.root");     
      TDirectory *tdir = input2->GetDirectory("ak4");
      TTree *T = (TTree*)tdir->Get("ProcessedTree");


      QCDEvent        *events;
      TBranch        *b_events;   //!

      // Set object pointer
      events = 0;
      // Set branch addresses and branch pointers
      if (!T) return;
      T->SetMakeClass(1);
      double logkTmin = -1;

      Int_t nEvmc=T->GetEntries();
      T->SetBranchAddress("events", &events, &b_events);

   for(int iEntry=0; iEntry< nEvmc; iEntry++)
   {
    T->GetEntry(iEntry);

   if ( events->nPFJetsCHS() < 2) continue;


   for (int k = 0; k < events->nPFJetsCHS(); k++) // det-level jets
   {
     bool mismatch = false;
     auto jet_reco = events->pfjetchs(k);

     auto z_reco = jet_reco.z();
     auto theta_reco = jet_reco.theta();
     auto kT_reco = jet_reco.kT();
     auto eta2_reco = jet_reco.eta2_splitting();
     auto phi2_reco = jet_reco.phi2_splitting();

     if  (!(z_reco.size()> 0 && jet_reco.pt() > 500 && fabs(jet_reco.y()) < 2.0 && z_reco.size() > 0 && fabs(jet_reco.eta()) < 2.0 )) continue;

     if  (!(z_reco.size()> 0 && jet_reco.pt() > 500 && fabs(jet_reco.y()) < 2.0 && z_reco.size() > 0 )) continue;
      int first_splitting = 0;

      for (unsigned j = 0; j < z_reco.size(); ++j) //loop over reco-level splittings
      {
       if (log(1./z_reco.at(j)) < 0.69 || log(1./z_reco.at(j)) > 6.0) continue;
       if (log(0.4/theta_reco.at(j)) < 0. || log(0.4/theta_reco.at(j)) > 3.67) continue;
       h2raw->Fill(log(0.4/theta_reco.at(j)) , log(1./z_reco.at(j) )  );
       h2rawTrueBinning->Fill(log(0.4/theta_reco.at(j)) , log(1/z_reco.at(j) )  );
      }

   }
  }


    TFile *fout=new TFile ("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/rawData.root","RECREATE");
    fout->cd();
    h2rawTrueBinning->SetName("rawTrueBinning");
    h2raw->SetName("raw");
    h2rawTrueBinning->Write();
    h2raw->Write();
}

#ifndef __CINT__
int main () { RooUnfoldReadData_LundPlaneDataZ(); return 0; }  // Main program when run stand-alone
#endif
