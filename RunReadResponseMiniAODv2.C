#include "./examples/RooUnfoldReadResponse_3D.cxx"
void RunReadResponseMiniAODv2(int i )
{
  gSystem->Load("libRooUnfold.so");
  gROOT->ProcessLine(".L ./examples/RooUnfoldReadResponse_3D.cxx");
  std::stringstream ss;   std::stringstream ss_herwig7;   std::stringstream ss_trackIneff;
  std::stringstream mctype; int year = 2018;
  std::stringstream yearPrefix;
  yearPrefix << "2018";
  mctype << "pythia8";
//  ss << "inputtest.dat";
//  ss << "./inputmcfiles/pythia8_preVFP_apr2022";
//   ss << "./inputmcfiles/mc2018_herwig7.dat";
//   ss << "./inputmcfiles/mc2018_pythia8.dat";


//   ss << "./inputmcfiles/mc2016preVFP_pythia8_trackIneff.dat";

//   ss << "./inputmcfiles/mc2016preVFP_pythia8_trackIneff_ak8.dat";
//   ss << "./inputmcfiles/mc2018_pythia8_trackIneff_ak8.dat";
// ss << "./inputmcfiles/mc2017_pythia8_ak8.dat";
//  ss << "./inputmcfiles/mc2017_herwig7_ak8.dat";
//  ss << "./inputmcfiles/mc2018_pythia8_ak8.dat";

/*  ss << "./inputmcfiles/mc" ; ss << yearPrefix.str().c_str(); ss << "_pythia8_";
  ss_trackIneff << "./inputmcfiles/mc"; ss_trackIneff << yearPrefix.str().c_str(); ss_trackIneff <<"_pythia8_trackIneff_"; 
  ss_herwig7 << "./inputmcfiles/mc"; ss_herwig7 << yearPrefix.str().c_str(); ss_herwig7 << "_herwig7_";*/

  ss << "./inputmcfiles_miniAODv2/mc2018_pythia8_ak4_miniAODv2.dat"; 
//   ss << "./inputmcfiles/mc2018_herwig7_ak8.dat";

//   ss << "./inputmcfiles/mc2016preVFP_herwig7_ak8.dat";
//      ss << "./inputmcfiles/mc2016postVFP_pythia8_ak8.dat";

//   ss << "./inputmcfiles/mc2016postVFP_herwig7_ak8.dat";
//
//      ss << "./inputmcfiles/mc2016postVFP_pythia8_ak8.dat";
 
// ss << "./inputmcfiles/mc2017_pythia8.dat";
//  ss << "./inputmcfiles/mc2017_herwig7.dat";

//ss << "./inputmcfiles/mc2017_pythia8.dat";
//   ss << "./inputmcfiles/mc2016preVFP_herwig7.dat";
//   ss << "./inputmcfiles/mc2016preVFP_pythia8.dat";

////   ss << "./inputmcfiles/mc2016postVFP_herwig7.dat";
//   ss << "./inputmcfiles/mc2016postVFP_pythia8.dat";
//   

//
//
//
//    ss << "./inputmcfiles/herwig7_preVFP_apr2022";
//  ss << "./inputmcfiles/pythia8_preVFP_chargeddown.dat";
//  ss << "./inputmcfiles/herwig7_preVFP.dat";
//  if (i < 10) {ss << "0"; ss_herwig7 << "0" ; ss_trackIneff << "0"; }
//  if ((i >= 10) && (i < 100)) ss << "0";
//  ss << i; ss_herwig7 << i ; ss_trackIneff << i;
//  cout << ss.str().c_str() << endl;
//std::string mc_name = "", std::string jec_option = "", std::string jer_option ="", std::string pu_option = ""
  RooUnfoldReadResponse_3D(ss.str().c_str(),"","",i, year, mctype.str().c_str(), "nom", "nom", "nom", "ak4");
/*  RooUnfoldReadResponse_3D(ss.str().c_str(),"","",i, year, mctype.str().c_str(), "nom", "nom", "up", "ak4");
  RooUnfoldReadResponse_3D(ss.str().c_str(),"","",i, year, mctype.str().c_str(), "nom", "nom", "down", "ak4");
  RooUnfoldReadResponse_3D(ss.str().c_str(),"","",i, year, mctype.str().c_str(), "nom", "up", "nom", "ak4");
  RooUnfoldReadResponse_3D(ss.str().c_str(),"","",i, year, mctype.str().c_str(), "nom", "nom", "nom", "ak4");
  RooUnfoldReadResponse_3D(ss.str().c_str(),"","",i, year, mctype.str().c_str(), "nom", "down", "nom", "ak4");
  RooUnfoldReadResponse_3D(ss.str().c_str(),"","",i, year, mctype.str().c_str(), "up", "nom", "nom", "ak4");
  RooUnfoldReadResponse_3D(ss.str().c_str(),"","",i, year, mctype.str().c_str(), "down", "nom", "nom", "ak4");*/

//  RooUnfoldReadResponse_3D(ss_herwig7.str().c_str(),"","",i, year, "herwig7", "nom", "nom", "nom", "ak8");
//  RooUnfoldReadResponse_3D(ss_trackIneff.str().c_str(),"TrackIneff","trackIneff",i, year, "pythia8", "nom", "nom", "nom", "ak4");
//  RooUnfoldReadResponse_3D(ss.str().c_str(),"nom","apr22",i, 2017, "pythia8", "nom", "nom", "up", "ak8");
//  RooUnfoldReadResponse_3D(ss.str().c_str(),"nom","apr22",i, 2017, "pythia8", "nom", "nom", "nom", "ak8");
                                                                               // jec, jer, pu, algo 
  ss.str("");
}
