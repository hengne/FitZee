#include "TFile.h"
//#include "drawMee.h"
#include "drawMeeComb.h"

int _nbins = 240;
double _scale = 1.0; //1.04141815; //1.040002;
int _oddeven = 0; // 0: both, 1: Odd, 2: Even


int main(int argc, char* argv[])
{
  if (argc<4)
  {
    std::cout << argv[0] << " <input_tree_file.root> <calib_table_file.dat> <output_file.root> <max_n_events> <method> <overall_scale> <OddEven> <eta_scale_file>"
    << std::endl;
    return 0;
  }
  
  std::string inrootfilename(argv[1]);
  std::string calibtablename(argv[2]);
  std::string outrootfilename(argv[3]);
 
  int maxevt = -1;
  if (argc>4)
  {
    maxevt = atoi(argv[4]);
  }
  
  int method=31;
  if (argc>5)
  {
    method = atoi(argv[5]);
  }
  
  if (argc>6)
  {
    _scale = atof(argv[6]);
  }

  if (argc>7)
  {
    _oddeven = atoi(argv[7]);
  }

  std::string etascalefilename;
  if (argc>8)
  { 
    etascalefilename = std::string(argv[8]);
  }

  if (method==7)
  {
    if (argc<=8 || etascalefilename=="") 
    {
      std::cout << " Error loading Etascale file for method 7, run " << argv[0] << " to get help." << std::endl; 
      return 1;
    }
  }

  std::cout << "Input file: " << inrootfilename << std::endl;
  std::cout << "Calib Table file: " << calibtablename << std::endl;
  std::cout << "ouput file: " << outrootfilename << std::endl;
  std::cout << "max events: " << maxevt << std::endl;
  std::cout << "method: " << method << std::endl;
  std::cout << "scale: " << _scale << std::endl;
  std::cout << "Both/Odd/Even: " << _oddeven << std::endl;
  if (method==7) std::cout << "Eta scale file : " << etascalefilename << std::endl;

  // reading data
  TChain* tree = new TChain("tree", "tree");
  tree->Add(inrootfilename.c_str());
  
  // output root file
  TFile* fout = new TFile(outrootfilename.c_str(), "recreate");

  // calibTable
  std::vector<calibRecord> calibTable = getCalibTableFromFile(calibtablename.c_str());
  

  // reference scales
  std::vector<EnergyScale> scalesref;
  if (method==7)
  {
    // read reference scales
    std::ifstream reffile(etascalefilename.c_str());
    if (reffile.is_open())
    {
      std::string line;
      EnergyScale scale;
      while (getline(reffile,line))
      {
        std::stringstream sline(line);
        sline >> scale.min 
            >> scale.max 
            >> scale.s
            >> scale.serr;
        scalesref.push_back(scale);
      }
      reffile.close();
    }
  }


  // Set the branches for the TChain/TTree
  SetTreeBranch(tree);

  // Get Original Mee Hist
  TH1D* h0 = (TH1D*)getTH1DOriginalMee(tree, "hMeeOrg", _nbins, maxevt, 1.0, _oddeven);
  h0->SetMarkerColor(2);
  h0->SetLineColor(2);
  
  // Get Mee Hist
  TH1D* h5 = (TH1D*)getTH1DMeeWithEtaScale(tree, calibTable, scalesref, "hMeeEtaScale", _nbins, maxevt, method, 1.0);
  h5->SetMarkerColor(4);
  h5->SetLineColor(4);
  
  // Get Mee Hist
  TH1D* h3 = (TH1D*)getTH1DMeeRegV8Elec(tree, "hMeeRegV8Elec", _nbins, maxevt, 1.0, _oddeven);
  h3->SetMarkerColor(6);
  h3->SetLineColor(6);

  // Get Mee Hist
  TH1D* h4 = (TH1D*)getTH1DMeeRegV8ElecWithEtaScale(tree, scalesref, "hMeeRegV8ElecEtaScale", _nbins, maxevt, 1.0);
  h4->SetMarkerColor(8);
  h4->SetLineColor(8);
 
  // Get Mee Hist
  TH1D* h1 = (TH1D*)getTH1DMee(tree, calibTable, "hMee", _nbins, maxevt, 23, 1.0, _oddeven);
  h1->SetMarkerColor(9);
  h1->SetLineColor(9);
 
  // print
  //std::cout << "nEvents: " << h1->GetEntries() << std::endl;
  
  // write the hist to out root file
  fout->cd();
  h0->Write();
  //h01->Write();
  h1->Write();
  //h2->Write();
  h3->Write();
  h4->Write();
  h5->Write();

  // delete the chain no more need it
  tree->Delete();
  
  // close
  fout->Close();
  
  //
  return 0;
}
