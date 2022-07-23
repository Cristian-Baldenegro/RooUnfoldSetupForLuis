#include "./examples/RooUnfoldReadData_GJ.cxx"
void RunReadData(int i)
{
  gSystem->Load("libRooUnfold.so");
  gROOT->ProcessLine(".L ./examples/RooUnfoldReadData_GJ.cxx");
  std::stringstream ss;
 
  ss << "./inputdatafiles/data_";
  if (i < 10) ss << "00";
  if ((i >= 10) && (i < 100)) ss << "0";
  ss << i;
 
  RooUnfoldReadData_GJ(ss.str().c_str(),"nom","nov_21",i);
  ss.str("");  
}
