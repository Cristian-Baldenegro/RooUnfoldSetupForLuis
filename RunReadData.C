#include "./examples/RooUnfoldReadData_3D.cxx"
void RunReadData(int i )
{
  gSystem->Load("libRooUnfold.so");
  gROOT->ProcessLine(".L ./examples/RooUnfoldReadData_3D.cxx");
  std::stringstream ss, era; 
//  ss << "./inputdatafiles/Run16Bl.dat";
//
//    ss << "./inputdatafiles/data_";
//    inputdata2016preVFP
//      if (i < 10) ss << "00";
//        if ((i >= 10) && (i < 100)) ss << "0";
//          ss << i;
//

//  ss << "inputtest.dat";
  era << "C";
//  
//  ss << "./inputdatafiles/2016EF_ak8_";
//  ss << "./inputdatafiles/2018D_ak8.dat";

//  ss << "./inputdatafiles/2016BCD_ak8_";
//  ss << "./inputdatafiles/2016FGH_ak8_";
//  ss << "./inputdatafiles/2016FGH";
//  ss << "./inputdatafiles/2016BCD";
// ss << "./inputdatafiles/2017B";  
 ss << "./inputdatafiles/2017C_ak8_";


//  ss << "./inputdatafiles/inputdata2016preVFP";
//   ss << "./inputdatafiles/2018A.dat";
//   ss << "./inputdatafiles/2018B.dat";
//   ss << "./inputdatafiles/2018C.dat";
//   ss << "./inputdatafiles/2018D.dat";

//   ss << "./inputdatafiles/2017B.dat";
//    ss << "./inputdatafiles/2017C.dat";
//   ss << "./inputdatafiles/2017D.dat";
//   ss << "./inputdatafiles/2017F.dat";
//   ss << "./inputdatafiles/2017F.dat";

//  if (i < 10) ss << "0";
//  if ((i >= 10) && (i < 100)) ss << "0";
  ss << i;

  RooUnfoldReadData_3D(ss.str().c_str(),"test","apr22",i, 2017, era.str(), "ak8");
  ss.str("");  
}
