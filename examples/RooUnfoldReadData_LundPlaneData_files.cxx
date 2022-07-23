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


void RooUnfoldReadData_LundPlaneData_files(std::string cFiles2 = "",std::string tag = "", std::string date = "", int fileindex=0)
{


  //errors from the diagonal elements of the covariance matrix
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
 
   const int bin_det_lnkt=18;
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
     Double_t ybins[bin_det_lnkt];

         ybins[0]=lnktmin_det;
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
         ybins[17]=6.0;

//         ybins[6]=5.;


       //true binning in xj
        Double_t ybins_true[bin_true_lnkt];
         ybins_true[0]=lnktmin_true;
         ybins_true[1]=lnktmin_det;
         ybins_true[2]=0.;
         ybins_true[3]=1.;
         ybins_true[4]=2.;
         ybins_true[5]=3.;
         ybins_true[6]=4.;
         ybins_true[7]=5.;


          //det binning in Rg
         Double_t xbins[bin_det_lntheta];
         xbins[0]=0;
         xbins[1]=.25;
         xbins[2]=.5;
         xbins[3]=.75;
         xbins[4]=1.;
         xbins[5]=1.25;
         xbins[6]=1.5;
         xbins[7]=1.75;
         xbins[8]=2.;
         xbins[9]=2.25;
         xbins[10]=2.5;
         xbins[11]=2.75;




        //true binning in Rg
         Double_t xbins_true[bin_true_lntheta];
         xbins_true[0]=0;
         xbins_true[1]=.25;
         xbins_true[2]=.5;
         xbins_true[3]=.75;
         xbins_true[4]=1.;
         xbins_true[5]=1.25;
         xbins_true[6]=1.5;
         xbins_true[7]=1.75;
         xbins_true[8]=2.;
         xbins_true[9]=2.25;
         xbins_true[10]=2.5;
         xbins_true[11]=2.75;
       //  xbins_true[12]=3.;




   TH1D *jet1Pt(0), *jet2Pt(0), *jetPtAsymmetry(0), *DeltaPhi(0), *jet1Rapidity(0), *jet2Rapidity(0), *DeltaRapidity(0), *nJets(0), *nSplittings(0);

   nJets = new TH1D("nJets", "nJets",11,-0.5,10.5);
   nSplittings = new TH1D("nSplittings", "nSplittings",16,-0.5,15.5);
   jet1Pt = new TH1D("jet1Pt", "jet1Pt",98,550,3000);
   jet2Pt = new TH1D("jet2Pt", "jet2Pt",98,550,3000);
   jetPtAsymmetry = new TH1D("jetPtAsymmetry", "jetPtAsymmetry",20,0, 1);
   DeltaPhi = new TH1D("DeltaPhi", "DeltaPhi",32,0, M_PI);
   jet1Rapidity =  new TH1D("jet1Rapidity", "jet1Rapidity",20,-2.0, 2.0);
   jet2Rapidity =  new TH1D("jet2Rapidity", "jet2Rapidity",20,-2.0, 2.0);
   DeltaRapidity =  new TH1D("DeltaRapidity", "DeltaRapidity",20,0, 4.0);

   int nbins = 10;
   int nbins_dR = 9;
   float logkTmax = 6.0;
   float nbinsTrue = 10;
   float logkTmaxTrue = 10.0;
   float logkTmin = -1;
   float logRmax = 3.0;


   //the raw correlation
   TH2D *h2raw(0), *h2rawTrueBinning(0);
   h2raw=new TH2D("raw","raw",nbins_dR, 0, logRmax,bin_det_lnkt-1,ybins);
   h2rawTrueBinning=new TH2D("rawTrueLevelBinning","rawTrueLevelBinning",nbins_dR, 0, logRmax,bin_det_lnkt-1,ybins);
   h2raw->Sumw2();
   h2raw->Sumw2();

   TH1D *pfRho_density(0);
   pfRho_density = new TH1D("pfRho_density", "pfRho_density", 60, 0, 60);




    ifstream infile;
    infile.open(cFiles2);
    char filename[300];
    int count=0; int n_splittings = 0;
    while(infile>>filename)
    {

//      TFile *input2=TFile::Open("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/DATA.root");    
      TFile *input2=TFile::Open(filename);
      cout << filename << endl;
      count=count+1;
      cout<<"reading data file "<<count<<endl;
      if(count==1001) break;
//      delete input2;      
      TDirectory *tdir = input2->GetDirectory("ak4");
      TTree *T = (TTree*)tdir->Get("ProcessedTree");

//      delete T;
//      delete events;
//      delete b_events;
//      input2->Close();
//      delete input2;

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

   if ( events->nPFJetsCHS() < 2) continue;
    auto eventHdr = events->evtHdr();
    double pfRho = eventHdr.pfRho();
    pfRho_density->Fill(eventHdr.pfRho());



    auto jet1 = events->pfjetchs(0);
    auto jet2 = events->pfjetchs(1);

    if (jet1.looseID() == false || jet2.looseID() == false ) continue;
    if (!(jet1.pt() > 550 && jet2.pt() > 550 && fabs(jet1.y()) < 2 && fabs(jet2.y()) <2 )) continue;
    jet1Pt->Fill(jet1.pt()); //, *jet2Pt(0), *jetPtAsymmetry, *DeltaPhi(0), *jet1Rapidity(0), *jet2Rapidity(0), *DeltaRapidity(0);
    jet2Pt->Fill(jet2.pt());
    double delta_phi = 0.;
    if (fabs(jet1.phi()-jet2.phi()) < M_PI) delta_phi = fabs(jet1.phi()-jet2.phi());
    else delta_phi = 2*M_PI-fabs(jet1.phi()-jet2.phi());
    DeltaPhi->Fill(delta_phi);
    jet1Rapidity->Fill(jet1.y());
    jet2Rapidity->Fill(jet2.y());
    DeltaRapidity->Fill(fabs(jet1.y()-jet2.y()));
    double pTasymmetry = fabs(jet1.pt()-jet2.pt())/(jet1.pt()+jet2.pt());
    jetPtAsymmetry->Fill(pTasymmetry);
   for (int k = 0; k < 2; k++) // det-level jets
   {
     bool mismatch = false;
     auto jet_reco = events->pfjetchs(k);

     auto z_reco = jet_reco.z();
     auto theta_reco = jet_reco.theta();
     auto kT_reco = jet_reco.kT();
     auto eta2_reco = jet_reco.eta2_splitting();
     auto phi2_reco = jet_reco.phi2_splitting();


     if  (!(z_reco.size()> 0 && jet_reco.pt() > 550 && fabs(jet_reco.y()) < 2.0 )) continue;

     if  (!(z_reco.size()> 0 && jet_reco.pt() > 550 && fabs(jet_reco.y()) < 2.0 && z_reco.size() > 0 )) continue;
      int first_splitting = 0;
      n_splittings = 0;
      for (unsigned j = 0; j < z_reco.size(); ++j) //loop over reco-level splittings
      {
       
       if (log(kT_reco.at(j)) < logkTmin || log(kT_reco.at(j)) > logkTmax) continue;
       if (log(0.4/theta_reco.at(j)) < 0. || log(0.4/theta_reco.at(j)) > logRmax) continue;
       n_splittings++;
       h2raw->Fill(log(0.4/theta_reco.at(j)) , log(kT_reco.at(j) )  );
       h2rawTrueBinning->Fill(log(0.4/theta_reco.at(j)) , log(kT_reco.at(j) )  );
      }
       nSplittings->Fill(n_splittings);
       nJets->Fill(1);
     }
     
    }
  }




//    TFile *fout=new TFile ("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/rawData.root","RECREATE");
    TFile *fout=new TFile (Form("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/ReadData_%s_%i_%s.root",tag.c_str(),fileindex,date.c_str()),"RECREATE");
    fout->cd();
  jet1Pt->Write(); jet2Pt->Write(); jetPtAsymmetry->Write(); DeltaPhi->Write(); jet1Rapidity->Write(); jet2Rapidity->Write(); DeltaRapidity->Write();
    nSplittings->Write(); nJets->Write(); pfRho_density->Write();
    h2rawTrueBinning->SetName("rawTrueBinning");
    h2raw->SetName("raw");
    h2rawTrueBinning->Write();
    h2raw->Write();
}

#ifndef __CINT__
int main () { RooUnfoldReadData_LundPlaneData_files(); return 0; }  // Main program when run stand-alone
#endif
