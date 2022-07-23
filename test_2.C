#define test_2_cxx
//#include "test_2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>



using namespace std;


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "LundPlane_LLR/AnalysisFW/interface/QCDEvent.h"

class test_2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   QCDEvent        *events;

   // List of branches
   TBranch        *b_events;   //!

   test_2(TTree *tree=0);
   virtual ~test_2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

test_2::test_2(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DATA_MG8_500MeV_0MeV.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("DATA_MG8_500MeV_0MeV.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("DATA_MG8_500MeV_0MeV.root:/ak4");
      dir->GetObject("ProcessedTree",tree);
   }
   Init(tree);
}

test_2::~test_2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test_2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test_2::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void test_2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   events = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("events", &events, &b_events);
   Notify();
}

Bool_t test_2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test_2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test_2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void test_2::Loop()
{
//   In a ROOT session, you can do:
//      root> .L test_2.C
//      root> test_2 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   map<string,TH1F*> histosTH1F;
   map<string,TH2F*> histosTH2F;



   histosTH1F["reco_efficiency_PF"] = new TH1F("reco_efficiency_PF", "reco_efficiency_PF",11,0.,1.1);
   histosTH1F["reco_efficiency_charged"] = new TH1F("reco_efficiency_charged", "reco_efficiency_charged",11,0.,1.1);

   histosTH2F["Lund plane raw splittings"] = new TH2F("Lund plane raw splittings", "Lund plane raw splittings", 40, 0, 6.0, 40, 0.6, 7);
   histosTH2F["Lund plane raw kT vs theta"] = new TH2F("Lund plane raw kT vs theta", "Lund plane raw kT vs theta", 40, 0, 6.0, 40, -3, 5.5);

   histosTH2F["Lund plane raw splittings charged"] = new TH2F("Lund plane raw splittings charged", "Lund plane raw splittings charged", 40, 0, 6.0, 40, 0.693, 7);
   histosTH2F["Lund plane raw kT vs theta charged"] = new TH2F("Lund plane raw kT vs theta charged", "Lund plane raw kT vs theta charged", 40, 0, 6.0, 40, -3, 5.5);

   histosTH2F["Lund plane raw splittings gen"] = new TH2F("Lund plane raw splittings gen", "Lund plane raw splittings gen", 40, 0, 6.0, 40, 0.693, 7);
   histosTH2F["Lund plane raw kT vs theta gen"] = new TH2F("Lund plane raw kT vs theta gen", "Lund plane raw kT vs theta gen", 40, 0, 6.0, 40, -3, 5.5);

   histosTH2F["Lund plane raw splittings charged gen"] = new TH2F("Lund plane raw splittings charged gen", "Lund plane raw splittings charged gen", 40, 0, 6.0, 40, 0.693, 7);
   histosTH2F["Lund plane raw kT vs theta charged gen"] = new TH2F("Lund plane raw kT vs theta charged gen", "Lund plane raw kT vs theta charged gen", 40, 0, 6.0, 40, -3, 5.5);

   histosTH1F["theta_highKt"] =  new TH1F("theta_highKt", "theta_highKt",10,0., 4.0);
   histosTH1F["theta_lowKt"] =  new TH1F("theta_lowKt", "theta_lowKt",20,0., 6.0);

   histosTH1F["kT_alltheta"] =  new TH1F("kT_alltheta", "kT_alltheta",20,-1.0, 6.0);

   int nbins = 16;
   int nbins_dR = 13;
   float logkTmax = 4.33333;
   float logkTmin = -1;
   float logRmax = 2.6;

   histosTH1F["dEta"] =  new TH1F("dEta", "dEta",200,-0.1, 0.1);
   histosTH1F["dPhi"] =  new TH1F("dPhi", "dPhi",100,0, 0.1);

   histosTH1F["dR"] =  new TH1F("dR", "dR",100,0, 0.1);

   histosTH1F["residual_kT"] =  new TH1F("residual_kT", "residual_kT",200,-1., 1.);
   histosTH1F["residual_kT_perturbative"] =  new TH1F("residual_kT_perturbative", "residual_kT_perturbative",200,-1., 1.);
   histosTH1F["residual_kT_nonperturbative"] =  new TH1F("residual_kT_nonperturbative", "residual_kT_nonperturbative",200,-1., 1.);

   histosTH1F["ratio_response_kT"] =  new TH1F("ratio_response_kT", "ratio_response_kT",200,0., 2.);
   histosTH1F["ratio_response_kT_perturbative"] =  new TH1F("ratio_response_kT_perturbative", "ratio_response_kT_perturbative",200,0., 2.);
   histosTH1F["ratio_response_kT_nonperturbative"] =  new TH1F("ratio_response_kT_nonperturbative", "ratio_response_kT_nonperturbative",200,0., 2.);

   histosTH2F["residual_vs_kTreco"] = new TH2F("residual_vs_kTreco", "", nbins, logkTmin, logkTmax, 25, -1, 1);
   histosTH2F["residual_vs_kTtrue"] = new TH2F("residual_vs_kTtrue", "", nbins, logkTmin, logkTmax, 25, -1, 1);

   histosTH2F["residual_vs_thetareco"] = new TH2F("residual_vs_thetareco", "", nbins_dR, 0, logRmax, 25, -1, 1);
   histosTH2F["residual_vs_thetatrue"] = new TH2F("residual_vs_thetatrue", "", nbins_dR, 0, logRmax, 25, -1, 1);

   histosTH1F["residual_theta"] =  new TH1F("residual_theta", "residual_theta",200,-1., 1.);
   histosTH1F["residual_theta_perturbative"] =  new TH1F("residual_theta_perturbative", "residual_theta_perturbative",200,-1., 1.);
   histosTH1F["residual_theta_nonperturbative"] =  new TH1F("residual_theta_nonperturbative", "residual_theta_nonperturbative",200,-1., 1.);

   histosTH1F["theta_highKt_charged"] =  new TH1F("theta_highKt_charged", "theta_highKt_charged",10,0., 4.0);
   histosTH1F["theta_lowKt_charged"] =  new TH1F("theta_lowKt_charged", "theta_lowKt_charged",20,0., 6.0);

   histosTH1F["kT_alltheta_charged"] =  new TH1F("kT_alltheta_charged", "kT_alltheta_charged",nbins_dR, 0, logRmax);

   histosTH2F["response kT neutral+charged"] = new TH2F("response kT neutral+charged", "", nbins, logkTmin, logkTmax, nbins, logkTmin, logkTmax);
   histosTH2F["response theta neutral+charged"] = new TH2F("response theta neutral+charged", "", nbins_dR, 0, logRmax, nbins_dR, 0, logRmax);

   histosTH1F["kT_gen_spectrum"] = new TH1F("kT_gen_spectrum", "", nbins, logkTmin, logkTmax);
   histosTH1F["theta_gen_spectrum"] = new TH1F("theta_gen_spectrum", "", nbins_dR, 0, logRmax);

   histosTH2F["response kT charged"] = new TH2F("response kT charged", "", nbins, logkTmin, logkTmax, nbins, logkTmin, logkTmax);
   histosTH2F["response theta charged"] = new TH2F("response theta charged", "", nbins_dR, 0, logRmax, nbins_dR, 0, logRmax );


   histosTH1F["kT_charged_gen_spectrum"] = new TH1F("kT_charged_gen_spectrum", "", nbins, logkTmin, logkTmax);
   histosTH1F["theta_charged_gen_spectrum"] = new TH1F("theta_charged_gen_spectrum", "", nbins_dR, 0, logRmax);

   histosTH2F["Lund_plane_reco_all"] = new TH2F("Lund_plane_reco_all", "", nbins_dR, 0, logRmax, nbins, logkTmin, logkTmax);
   histosTH2F["Lund_plane_reco_matched"] = new TH2F("Lund_plane_reco_matched", "", nbins_dR, 0, logRmax, nbins, logkTmin, logkTmax);
   histosTH2F["Lund_plane_gen_all"] = new TH2F("Lund_plane_gen_all", "", nbins_dR, 0, logRmax, nbins, logkTmin, logkTmax);
   histosTH2F["Lund_plane_gen_matched"] = new TH2F("Lund_plane_gen_matched", "", nbins_dR, 0, logRmax, nbins, logkTmin, logkTmax);

   histosTH2F["Lund_plane_reco_all_z"] = new TH2F("Lund_plane_reco_all_z", "", nbins_dR, 0, logRmax, nbins, 0.69, 6);
   histosTH2F["Lund_plane_reco_matched_z"] = new TH2F("Lund_plane_reco_matched_z", "", nbins_dR, 0, logRmax, nbins, 0.69, 6);
   histosTH2F["Lund_plane_gen_all_z"] = new TH2F("Lund_plane_gen_all_z", "", nbins_dR, 0, logRmax, nbins, 0.69, 6);
   histosTH2F["Lund_plane_gen_matched_z"] = new TH2F("Lund_plane_gen_matched_z", "", nbins_dR, 0, logRmax, nbins, 0.69, 6);

   histosTH2F["efficiencies_2D"] = new TH2F("efficiencies_2D", "", nbins_dR, 0, logRmax, nbins, logkTmin, logkTmax);
   TH2F *h_purities_2D = new TH2F("purities_2D", "", nbins_dR, 0, logRmax, nbins, logkTmin, logkTmax);
   TH2F *h_purities_2D_denom = new TH2F("purities_2D_denom", "", nbins_dR, 0, logRmax, nbins, logkTmin, logkTmax );

   TFile* output = new TFile("output.root","RECREATE");



   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
  //    Long64_t ientry = LoadTree(jentry);
 //     if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //
   //QCDPFJet = events[jentry].pfjetchs[0]; 
   //cout << events->pfmjj() << endl;
   if ( events->nPFJetsCHS() < 2) continue;
   for (int i = 0; i < events->nPFJetsCHS(); i++)
   {
     auto jet = events->pfjetchs(i);
//   cout << jet1.lha() << endl;
     auto z = jet.z();
     auto theta = jet.theta();
     auto kT = jet.kT();


     auto z_charged = jet.z_charged();
     auto theta_charged = jet.theta_charged();
     auto kT_charged = jet.kT_charged();

     if (z.size()> 0 && jet.pt() > 500 && fabs(jet.y()) < 2.0 && z_charged.size() > 0 && fabs(jet.eta()) < 2.0 )
      {
        for (unsigned i = 0; i < z.size(); ++i)
        {
          histosTH2F["Lund plane raw splittings"]->Fill(log(0.4/theta.at(i)), log(1/z.at(i)) );
          histosTH2F["Lund plane raw kT vs theta"]->Fill(log(0.4/theta.at(i)), log(kT.at(i)) );

          if (log(kT.at(i))>2) histosTH1F["theta_highKt"]->Fill(log(0.4/theta.at(i)));
          if (log(kT.at(i))<2 && log(kT.at(i)) > -1) histosTH1F["theta_lowKt"]->Fill(log(0.4/theta.at(i)));
          histosTH1F["kT_alltheta"]->Fill(log(kT.at(i)));

        }

        for (unsigned i = 0; i < z_charged.size(); ++i)
        {
          histosTH2F["Lund plane raw splittings charged"]->Fill(log(0.4/theta_charged.at(i)), log(1/z_charged.at(i)) );
          histosTH2F["Lund plane raw kT vs theta charged"]->Fill(log(0.4/theta_charged.at(i)), log(kT_charged.at(i)) );
          if (log(kT_charged.at(i))>2) histosTH1F["theta_highKt_charged"]->Fill(log(0.4/theta_charged.at(i)));
          if (log(kT_charged.at(i))<2 && log(kT_charged.at(i)) > -1 ) histosTH1F["theta_lowKt_charged"]->Fill(log(0.4/theta_charged.at(i)));
          histosTH1F["kT_alltheta_charged"]->Fill(log(kT_charged.at(i)));
        }

     }
    }

//   if ( events->nPFJetsCHS < 2) continue;
   for (int i = 0; i < 2; i++)
   {
     auto jet = events->genjet(i);
   //cout << jet.lha() << endl;
     auto z = jet.z();
     auto theta = jet.theta();
     auto kT = jet.kT();

     auto z_charged = jet.z_charged();
     auto theta_charged = jet.theta_charged();
     auto kT_charged = jet.kT_charged();
//    cout << jet.pt() << endl;

     if (z.size()> 0 && jet.pt() > 500  && z_charged.size() > 0 && fabs(jet.y()) < 2.0 )
      {
        for (unsigned i = 0; i < z.size(); ++i)
        {
          histosTH2F["Lund plane raw splittings gen"]->Fill(log(0.4/theta.at(i)), log(1/z.at(i)) );
          histosTH2F["Lund plane raw kT vs theta gen"]->Fill(log(0.4/theta.at(i)), log(kT.at(i)) );
//          cout << eta2.at(i)-jet.eta() << " " << kT.at(i) << endl; 
        }

        for (unsigned i = 0; i < z_charged.size(); ++i)
        {
          histosTH2F["Lund plane raw splittings charged gen"]->Fill(log(0.4/theta_charged.at(i)), log(1/z_charged.at(i)) );
          histosTH2F["Lund plane raw kT vs theta charged gen"]->Fill(log(0.4/theta_charged.at(i)), log(kT_charged.at(i)) );
        }

     }
  }

   for (int k = 0; k < events->nPFJetsCHS(); k++) // det-level jets
   {
     bool mismatch = false;
     auto jet_reco = events->pfjetchs(k);

     auto z_reco = jet_reco.z();
     auto theta_reco = jet_reco.theta();
     auto kT_reco = jet_reco.kT();
     auto eta2_reco = jet_reco.eta2_splitting();
     auto phi2_reco = jet_reco.phi2_splitting();

     auto eta2_charged_reco = jet_reco.eta2_splitting_charged();
     auto phi2_charged_reco = jet_reco.eta2_splitting_charged();

     auto z_charged_reco = jet_reco.z_charged();
     auto theta_charged_reco = jet_reco.theta_charged();
     auto kT_charged_reco = jet_reco.kT_charged();

     if ( jet_reco.genidx() == -1 ) continue; //if det-level jet is not matched, skip
//     if ( jet_reco.pfDeepCSVb() < 0.50 ) continue; 

     auto jet_gen = events->genjet(jet_reco.genidx());
     auto z_gen = jet_gen.z();
     auto theta_gen = jet_gen.theta();
     auto kT_gen = jet_gen.kT();
     auto eta2_gen = jet_gen.eta2_splitting();
     auto phi2_gen = jet_gen.phi2_splitting();


     auto z_charged_gen = jet_gen.z_charged();
     auto theta_charged_gen = jet_gen.theta_charged();
     auto kT_charged_gen = jet_gen.kT_charged();
     auto eta2_charged_gen = jet_gen.eta2_charged();
     auto phi2_charged_gen = jet_gen.phi2_charged();

     if  (!(z_charged_gen.size()> 0 && jet_gen.pt() > 500 && fabs(jet_gen.y()) < 2.0 && z_charged_reco.size() > 0 && fabs(jet_gen.eta()) < 2.0 )) continue;

     if  (!(z_gen.size()> 0 && jet_gen.pt() > 500 && fabs(jet_gen.y()) < 2.0 && z_reco.size() > 0 )) continue;
      int first_splitting = 0;


     for (unsigned i = 0; i < z_gen.size(); ++i) // loop over gen-level splittings
     {
     if (mismatch == true) cout << "=======================================" << endl;
     if (mismatch == true) cout << "jet" << "\t" << jet_gen.pt() << "\t" << jet_gen.phi() << "\t" << jet_gen.eta() << endl;
//      cout << "=================== " << i << endl;
//      cout << "gen-level splitting " << i << endl;
 //     cout << "=================== " << i << endl;

      if (log(kT_gen.at(i)) < -6 || log(kT_gen.at(i)) > 7) continue; //gen-level phase-space
      if (log(0.4/theta_gen.at(i)) < 0. || log(0.4/theta_gen.at(i)) > 3) continue; //gen-level phase-space
 //     if (log(1/z_gen.at(i)) > 6.9) continue; //gen-level phase-space
       histosTH2F["Lund_plane_gen_all"]->Fill(log(0.4/theta_gen.at(i)),log(kT_gen.at(i) ));
       histosTH2F["Lund_plane_gen_all_z"]->Fill(log(0.4/theta_gen.at(i)),log(1./z_gen.at(i) ));

      first_splitting = first_splitting+1;
      float dR_window = 0.1;
      float dR_max = dR_window;
      int index_true = -1;
      int index_reco = -1;
      float kT_min = 0.5;
      for (unsigned j = 0; j < z_reco.size(); ++j) //loop over reco-level splittings
      {
       if (log(kT_reco.at(j)) < logkTmin || log(kT_reco.at(j)) > 4.33333) continue;
       if (log(0.4/theta_reco.at(j)) < 0. || log(0.4/theta_reco.at(j)) > 2.6) continue;
      // cout << " first splitting " << first_splitting << " j " << j << endl;
       if (mismatch == true) cout << "reco " << log(0.4/theta_reco.at(j)) << "\t" << log(kT_reco.at(j)) << "\t" << phi2_reco.at(j) << "\t" << eta2_reco.at(j) << "\t" << endl;      
       if (first_splitting == 1)
       {
        histosTH2F["Lund_plane_reco_all"]->Fill(log(0.4/theta_reco.at(j)),log(kT_reco.at(j) ));
        histosTH2F["Lund_plane_reco_all_z"]->Fill(log(0.4/theta_reco.at(j)),log(1./z_reco.at(j) ));
        h_purities_2D_denom->Fill(log(0.4/theta_reco.at(j)),log(kT_reco.at(j) ));
       }

       float dphi = fabs(phi2_reco.at(j) - phi2_gen.at(i));
       float deta = eta2_reco.at(j) - eta2_gen.at(i);
       if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
       float dR = std::sqrt(dphi*dphi + deta*deta);
       float kT_ratio = min(kT_reco.at(j),kT_gen.at(i))/max(kT_reco.at(j),kT_gen.at(i));
       float pT_ratio = ( kT_reco.at(j)/sin(theta_reco.at(j)) ) / ( kT_gen.at(i)/sin(theta_gen.at(i)) );
//       cout << "det-level splitting " << j << " jetRecoPhi " <<  jet_reco.phi() << " recoSplitPhi " << phi2_reco.at(j) << "jetGenPhi " << jet_gen.phi() << " genSplitPhi " << phi2_gen.at(i)  << " dPhi " << dphi << endl;
//       cout << "det-level splitting " << j << " dR " <<  dR << " dEta " << deta << " dPhi " << dphi  << " kT ratio " << kT_ratio << endl;
//       cout << "det-level splitting " << j << " jetRecoPhi " <<  jet_reco.jetPhiCA_charged() << " recoSplitPhi " << phi2_reco.at(j) << "jetGenPhi " << jet_gen.jetPhiCA_charged() << "genSplitPhi " << phi2_gen.at(i)  << endl;
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
         if (log(kT_gen.at(l)) < -3 || log(kT_gen.at(l)) > 7) continue;
         if (log(0.4/theta_gen.at(l)) < 0. || log(0.4/theta_gen.at(l)) > 3) continue;
//         if (log(1/z_gen.at(l)) > 6.9) continue; //gen-level phase-space
         if (mismatch == true)  cout << "gen " << log(0.4/theta_gen.at(l)) << "\t" << log(kT_gen.at(l)) << "\t" << phi2_gen.at(l) << "\t" << eta2_gen.at(l) << "\t" << endl;
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

     dPhi_final = fabs(phi2_reco.at(index_reco) - phi2_gen.at(index_true));
     dEta_final = eta2_reco.at(index_reco) - eta2_gen.at(index_true);
     if (dPhi_final > TMath::Pi()) dPhi_final = 2.*TMath::Pi() - dPhi_final;
     dR_final = std::sqrt(dEta_final*dEta_final+dPhi_final*dPhi_final);

//     cout << theta_reco.at(index_reco) << " " << kT_reco(index_reco) << " " << theta_gen.at(index_true) << " " << kT_gen(index_true) << endl;
 

     histosTH1F["dR"]->Fill(dR_final);
     histosTH1F["dEta"]->Fill(dEta_final);
     histosTH1F["dPhi"]->Fill(dPhi_final);
     histosTH2F["Lund_plane_reco_matched"]->Fill(log(0.4/theta_reco.at(index_reco)),log(kT_reco.at(index_reco) ));
     histosTH2F["Lund_plane_gen_matched"]->Fill(log(0.4/theta_gen.at(index_true)),log(kT_gen.at(index_true) ));

     histosTH2F["Lund_plane_reco_matched_z"]->Fill(log(0.4/theta_reco.at(index_reco)),log(1./z_reco.at(index_reco) ));
     histosTH2F["Lund_plane_gen_matched_z"]->Fill(log(0.4/theta_gen.at(index_true)),log(1./z_gen.at(index_true) ));

     histosTH1F["residual_kT"]->Fill((kT_reco.at(index_reco)-kT_gen.at(index_true) ) /kT_gen.at(index_true) );
     if (log(kT_gen.at(index_true)) > 2) histosTH1F["residual_kT_perturbative"]->Fill((kT_reco.at(index_reco)-kT_gen.at(index_true) ) /kT_gen.at(index_true) );
     if (log(kT_gen.at(index_true)) < 2) histosTH1F["residual_kT_nonperturbative"]->Fill((kT_reco.at(index_reco)-kT_gen.at(index_true) ) /kT_gen.at(index_true) );

     histosTH1F["ratio_response_kT"]->Fill(kT_reco.at(index_reco) /kT_gen.at(index_true) );
     if (log(kT_gen.at(index_true)) > 2)
     {
      histosTH1F["ratio_response_kT_perturbative"]->Fill(kT_reco.at(index_reco) /kT_gen.at(index_true) );
      histosTH2F["residual_vs_kTreco"]->Fill(log(kT_reco.at(index_reco)), (kT_reco.at(index_reco)-kT_gen.at(index_true) ) /kT_gen.at(index_true));
      histosTH2F["residual_vs_kTtrue"]->Fill(log(kT_gen.at(index_true)), (kT_reco.at(index_reco)-kT_gen.at(index_true) ) /kT_gen.at(index_true));
      histosTH2F["residual_vs_thetareco"]->Fill(log(0.4/theta_reco.at(index_reco)), (theta_reco.at(index_reco)-theta_gen.at(index_true) ) /theta_gen.at(index_true)  );
      histosTH2F["residual_vs_thetatrue"]->Fill(log(0.4/theta_gen.at(index_true)), (theta_reco.at(index_reco)-theta_gen.at(index_true) ) /theta_gen.at(index_true)  );
      //if ( fabs(log(kT_gen.at(index_true)) - 1.55325) < 0.0001 && fabs(log(0.4/theta_gen.at(index_true)) - 0.659149) < 0.0001) mismatch = true;
      //if (log(kT_reco.at(index_reco)) < -1 &&  fabs(log(kT_gen.at(index_true)) - 2.34441) < 0.001 && fabs(log(0.4/theta_gen.at(index_true)) - 1.85783 ) < 0.001 ) mismatch = true;

     }
     if (log(kT_gen.at(index_true)) < 2) histosTH1F["ratio_response_kT_nonperturbative"]->Fill(kT_reco.at(index_reco) /kT_gen.at(index_true) );


     histosTH1F["residual_theta"]->Fill((theta_reco.at(index_reco)-theta_gen.at(index_true) ) /theta_gen.at(index_true) );
     if (log(kT_gen.at(index_true)) > 2) histosTH1F["residual_theta_perturbative"]->Fill((theta_reco.at(index_reco)-theta_gen.at(index_true) ) /theta_gen.at(index_true)  );
     if (log(kT_gen.at(index_true)) < 2) histosTH1F["residual_theta_nonperturbative"]->Fill((theta_reco.at(index_reco)-theta_gen.at(index_true) ) /theta_gen.at(index_true) );


     h_purities_2D->Fill(log(0.4/theta_reco.at(index_reco)),log(kT_reco.at(index_reco) ));

     histosTH2F["response kT neutral+charged"]->Fill(log(kT_gen.at(index_true)), log(kT_reco.at(index_reco)));
     histosTH2F["response theta neutral+charged"]->Fill(log(0.4/theta_gen.at(index_true)), log(0.4/theta_reco.at(index_reco)));
     histosTH1F["theta_gen_spectrum"]->Fill(log(0.4/theta_gen.at(index_true)));
     histosTH1F["kT_gen_spectrum"]->Fill(log(kT_gen.at(index_true)));
     }///end for


/*     if  (!(z_charged_gen.size()> 0 && jet_gen.pt() > 50 && fabs(jet_gen.y()) < 2.0 && z_charged_reco.size() > 0 && fabs(jet_gen.eta() < 2.0 ))) continue;
*/

  }
 
}

  for (int i = 1; i < nbins+1; i++) // normalization for response matrix
  {

   for (int j = 1; j < nbins+1; j++)
   {
//    histosTH2F["response kT charged"]->SetBinContent(i,j, float(histosTH2F["response kT charged"]->GetBinContent(i,j))/float(histosTH1F["kT_charged_gen_spectrum"]->GetBinContent(i)) );
    histosTH2F["response kT neutral+charged"]->SetBinContent(i,j, float(histosTH2F["response kT neutral+charged"]->GetBinContent(i,j))/float(histosTH1F["kT_gen_spectrum"]->GetBinContent(i)) );    
   }

  }

  for (int i = 1; i < nbins_dR+1; i++) // normalization for response matrix
  {

   for (int j = 1; j < nbins_dR+1; j++)
   {
//    histosTH2F["response theta charged"]->SetBinContent(i,j, float((histosTH2F["response theta charged"]->GetBinContent(i,j) ))/float((histosTH1F["theta_charged_gen_spectrum"]->GetBinContent(i))) );
    histosTH2F["response theta neutral+charged"]->SetBinContent(i,j, float(histosTH2F["response theta neutral+charged"]->GetBinContent(i,j))/float(histosTH1F["theta_gen_spectrum"]->GetBinContent(i)));
   }

  }

//  histosTH2F["purities_2D"] = (TH2F*) histosTH2F["Lund_plane_reco_matched"]->Clone();
//  histosTH2F["purities_2D"]->Divide(histosTH2F["Lund_plane_reco_all"]);
//  histosTH2F["purities_2D"]->SetName("purities_2D");
//  h_purities_2D->Divide(h_purities_2D_denom);
/*  histosTH2F["efficiencies_2D"] = (TH2F*) histosTH2F["Lund_plane_gen_matched"]->Clone();
  histosTH2F["efficiencies_2D"]->Divide(histosTH2F["Lund_plane_gen_all"]);
  histosTH2F["efficiencies_2D"]->SetName("efficiencies_2D");*/

//for some reason Divide is not working as it should, maybe due to empty bins?



  output->cd();
  //histosTH1F["jet1Pt"]->Sumw2();

    h_purities_2D->Write();
    h_purities_2D_denom->Write();

   for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
      it->second->Sumw2();

  for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin();
                                  it_histo != histosTH1F.end(); ++it_histo)
     (*it_histo).second->Write();

   for(map<string,TH2F*>::const_iterator it2 = histosTH2F.begin(); it2 != histosTH2F.end(); ++it2)
      it2->second->Sumw2();

  for(map<string,TH2F*>::iterator it_histo2 = histosTH2F.begin();
                                  it_histo2 != histosTH2F.end(); ++it_histo2)
     (*it_histo2).second->Write();

  output->Close();
}

