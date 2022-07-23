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


void RooUnfoldReadResponse_LundPlane2(std::string cFiles2="",std::string tag = "", std::string date = "", int fileindex=0)
{

#ifdef __CINT__
  gSystem->Load("libRooUnfold.so");
//  gSystem->Load("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/interface/QCDEvent.h");
#endif
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
     Double_t ybins[bin_det_lnkt];
      
         ybins[0]=lnktmin_det;
         ybins[1]=0;
	 ybins[2]=1;
 	 ybins[3]=2.;
         ybins[4]=3.;
         ybins[5]=4.;
//         ybins[6]=5.;
       
     
       //true binning in xj
        Double_t ybins_true[bin_true_lnkt];
         ybins_true[0]=lnktmin_det;
         ybins_true[1]=-0.5;
//         ybins_true[1]=lnktmin_det;
	 ybins_true[2]=0.;	
         ybins_true[3]=1.;
         ybins_true[4]=2.;
         ybins_true[5]=3.;
         ybins_true[6]=3.5;
         ybins_true[7]=4.;

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



           
 
  
   //MC detector correlation
   TH2D *h2smeared(0);
//   h2smeared=new TH2D("smeared","smeared",bin_det_lntheta-1,xbins,bin_det_lnkt-1,ybins);

   h2smeared=new TH2D("smeared","smeared",bin_true_lntheta-1,xbins_true,bin_true_lnkt-1,ybins_true);

   TH2D *h2recoall(0);
//   h2recoall=new TH2D("recoall","recoall",bin_det_lntheta-1,xbins,bin_det_lnkt-1,ybins);

   h2recoall=new TH2D("recoall","recoall",bin_true_lntheta-1,xbins_true,bin_true_lnkt-1,ybins_true);


   TH2D *h2purity(0);
//   h2purity=new TH2D("purity","purity",bin_det_lntheta-1,xbins,bin_det_lnkt-1,ybins);
   h2purity=new TH2D("purity","purity",bin_true_lntheta-1,xbins_true,bin_true_lnkt-1,ybins_true);


   //true correlations with det-level cuts
   TH2D *h2true(0);
    h2true=new TH2D("true","true",bin_true_lntheta-1,xbins_true,bin_true_lnkt-1,ybins_true);
     
    //fully-efficient true correlation (no det-level cuts)
    TH2D *h2fulleff(0);
    h2fulleff=new TH2D("truef","truef",bin_true_lntheta-1,xbins_true,bin_true_lnkt-1,ybins_true);
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



 
   //******************Fill the MCresponse*******************//

      TFile *input2=TFile::Open("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/DATA_MG8_testJan312022.root");
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



     if ( jet_reco.genidx() == -1 ) continue; //if det-level jet is not matched, skip

     auto jet_gen = events->genjet(jet_reco.genidx());
     auto z_gen = jet_gen.z();
     auto theta_gen = jet_gen.theta();
     auto kT_gen = jet_gen.kT();
     auto eta2_gen = jet_gen.eta2_splitting();
     auto phi2_gen = jet_gen.phi2_splitting();



     if  (!(z_gen.size()> 0 && jet_gen.pt() > 500 && fabs(jet_gen.y()) < 2.0 && z_reco.size() > 0 && fabs(jet_gen.eta()) < 2.0 )) continue;

     if  (!(z_gen.size()> 0 && jet_gen.pt() > 500 && fabs(jet_gen.y()) < 2.0 && z_reco.size() > 0 )) continue;
      int first_splitting = 0;


     for (unsigned i = 0; i < z_gen.size(); ++i) // loop over gen-level splittings
     {
     if (mismatch == true) cout << "=======================================" << endl;
     if (mismatch == true) cout << "jet" << "\t" << jet_gen.pt() << "\t" << jet_gen.phi() << "\t" << jet_gen.eta() << endl;
//      cout << "=================== " << i << endl;
//      cout << "gen-level splitting " << i << endl;
 //     cout << "=================== " << i << endl;

      if (log(kT_gen.at(i)) < -2 || log(kT_gen.at(i)) > 5) continue; //gen-level phase-space
      if (log(0.4/theta_gen.at(i)) < 0. || log(0.4/theta_gen.at(i)) > 3) continue; //gen-level phase-space
      h2fulleff->Fill(log(0.4/theta_gen.at(i)) , log(kT_gen.at(i)) );
//       histosTH2F["Lund_plane_gen_all"]->Fill(log(0.4/theta_gen.at(i)),log(kT_gen.at(i) ));

      first_splitting = first_splitting+1;
      float dR_window = 0.1;
      float dR_max = dR_window;
      int index_true = -1;
      int index_reco = -1;
      float kT_min = 0.5;
      for (unsigned j = 0; j < z_reco.size(); ++j) //loop over reco-level splittings
      {
       if (log(kT_reco.at(j)) < logkTmin || log(kT_reco.at(j)) > 4) continue;
       if (log(0.4/theta_reco.at(j)) < 0. || log(0.4/theta_reco.at(j)) > 2.75) continue;
      // cout << " first splitting " << first_splitting << " j " << j << endl;
       if (mismatch == true) cout << "reco " << log(0.4/theta_reco.at(j)) << "\t" << log(kT_reco.at(j)) << "\t" << phi2_reco.at(j) << "\t" << eta2_reco.at(j) << "\t" << endl;      
       if (first_splitting == 1)
       {
        h2recoall->Fill(log(0.4/theta_reco.at(j)),log(kT_reco.at(j) ));
//        h_purities_2D_denom->Fill(log(0.4/theta_reco.at(j)),log(kT_reco.at(j) ));
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
//         cout << " we have a match with dR " << dR << " and kTratio " << kT_min << endl;
       }//endif


      }//end for

       if ( index_reco == -1 || index_true == -1 ) continue;
       int index_true_mc = -1;
       dR_max = dR_window;
       kT_min = 0.5;
       float dR_final, dEta_final, dPhi_final;
       for (unsigned l = 0; l < z_gen.size(); ++l) // loop over gen-level splitings, check if the previous match is true the other way around
       {
         if (log(kT_gen.at(l)) < -2 || log(kT_gen.at(l)) > 5) continue;
         if (log(0.4/theta_gen.at(l)) < 0. || log(0.4/theta_gen.at(l)) > 3) continue;
         float dphi = fabs(phi2_reco.at(index_reco) - phi2_gen.at(l));
         float deta = eta2_reco.at(index_reco) - eta2_gen.at(l);
         if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
         float dR = std::sqrt(dphi*dphi + deta*deta);
         float kT_ratio = min(kT_reco.at(index_reco),kT_gen.at(l))/max(kT_reco.at(index_reco),kT_gen.at(l));
         float pT_ratio = ( kT_reco.at(index_reco)/sin(theta_reco.at(index_reco)) ) / ( kT_gen.at(i)/sin(theta_gen.at(l)) ); 
         if (dR < dR_max ) //find best pair
         {
           dR_max = dR;
           index_true_mc = l;
           kT_min = kT_ratio;
         }//endif
       }//end for
      if ( index_true != index_true_mc ) continue;

       h2smeared->Fill(log(0.4/theta_reco.at(index_reco)) , log(kT_reco.at(index_reco) )  );
       h2true->Fill(log(0.4/theta_gen.at(index_true)) , log(kT_gen.at(index_true)));
       response.Fill(log(0.4/theta_reco.at(index_reco)) , log(kT_reco.at(index_reco) )  , log(0.4/theta_gen.at(index_true)) , log(kT_gen.at(index_true) ) );

//     histosTH2F["Lund_plane_reco_matched"]->Fill(log(0.4/theta_reco.at(index_reco)),log(kT_reco.at(index_reco) ));
//     histosTH2F["Lund_plane_gen_matched"]->Fill(log(0.4/theta_gen.at(index_true)),log(kT_gen.at(index_true) ));

     }
   }
  }
   
 TFile *fout = new TFile("/grid_mnt/data__data.polcms/cms/cbaldene/gen_level_angularities_debug2/C10620p1/src/LundPlane_LLR/AnalysisFW/python/UnfoldingTest/RooUnfold/ResponseIsRead_test3.root","RECREATE");
  fout->cd();
  h2purity->Divide(h2recoall,h2smeared);
  h2smeared->SetName("smeared");
  h2smeared->Write();
  h2purity->SetName("purity");
  h2purity->Write();
  h2fulleff->Write();
  h2true->Write();
  response.Write();

}

#ifndef __CINT__
int main () { RooUnfoldReadResponse_LundPlane2(); return 0; }  // Main program when run stand-alone
#endif
