#include "./examples/RooUnfoldSplitResponse_GJ.cxx"
void RunSplitResponse(int i,double frac)
{
  gSystem->Load("libRooUnfold.so");
  gROOT->ProcessLine(".L ./examples/RooUnfoldSplitResponse_GJ.cxx");
  std::stringstream ss;
 
  ss << "./inputmcfiles/mc_";
  if (i < 10) ss << "00";
  if ((i >= 10) && (i < 100)) ss << "0";
  ss << i;
 
  RooUnfoldSplitResponse_GJ(ss.str().c_str(),"nom","nov_21",i,frac);
  ss.str("");  
}
