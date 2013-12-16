#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<fstream>
#include"TH2D.h"
#include"TFile.h"
#include"TROOT.h"


void printLine(int idx, int ix, int iy, int iz);
void printLine(int idx, int ix, int iy, int iz, std::ofstream& ofile);

std::vector<std::vector<int> > GetAllCellsInEEP();
std::vector<std::vector<int> > GetAllCellsInEEN();
std::vector<std::vector<int> > GetAllCellsInEE();
std::vector<std::vector<int> > GetAllCellsInEB();

int main()
{
  char name[1000];
  
  // EE+
  std::vector<std::vector<int> > cellsp = GetAllCellsInEEP(); 
  for (int idx=0; idx<(int)cellsp.size(); idx++)
  {
    printLine(idx, cellsp.at(idx).at(0), cellsp.at(idx).at(1), cellsp.at(idx).at(2));
  } 
  // to file
  std::ofstream ofilep("emptyCalibTable_IZ1.dat");
  if (ofilep.is_open())
  {
    for (int idx=0; idx<(int)cellsp.size(); idx++)
    {
      printLine(idx, cellsp.at(idx).at(0), cellsp.at(idx).at(1), cellsp.at(idx).at(2), ofilep);
    }
    ofilep.close();
  }
  // to small files
  // each ix is one file
  for (int ix=1; ix<=100; ix++)
  {
    sprintf(name, "emptyCalibTable_IZ1_IX%d.dat", ix);
    std::ofstream A_ofile(name);
    for (int idx=0; idx<(int)cellsp.size(); idx++)
    {
      if (cellsp.at(idx).at(0)==ix)
      {
        printLine(idx, cellsp.at(idx).at(0), cellsp.at(idx).at(1), cellsp.at(idx).at(2), A_ofile);
      }
    }
    A_ofile.close();
  }
  
  // EE-
  std::vector<std::vector<int> > cellsn = GetAllCellsInEEN(); 
  for (int idx=0; idx<(int)cellsn.size(); idx++)
  {
    printLine(idx, cellsn.at(idx).at(0), cellsn.at(idx).at(1), cellsn.at(idx).at(2));
  }  
  // to file
  std::ofstream ofilen("emptyCalibTable_IZ-1.dat");
  if (ofilen.is_open())
  {
    for (int idx=0; idx<(int)cellsn.size(); idx++)
    {
      printLine(idx, cellsn.at(idx).at(0), cellsn.at(idx).at(1), cellsn.at(idx).at(2), ofilen);
    }
    ofilen.close();
  }
  // to small files
  // each ix is one file
  for (int ix=1; ix<=100; ix++)
  {
    sprintf(name, "emptyCalibTable_IZ-1_IX%d.dat", ix);
    std::ofstream A_ofile(name);
    for (int idx=0; idx<(int)cellsn.size(); idx++)
    {
      if (cellsn.at(idx).at(0)==ix)
      {
        printLine(idx, cellsn.at(idx).at(0), cellsn.at(idx).at(1), cellsn.at(idx).at(2), A_ofile);
      }
    }
    A_ofile.close();
  }
 
  // EB


  for (int iy=1; iy<=360; iy++)
  {
    //emptyCalibTable_IZ${IZ}_IY${IY}_IX${IX}.dat
    sprintf(name, "emptyCalibTable_IZ0_IY%d_IX-.dat", iy);
    std::ofstream A_ofile(name);
    int idx=0;
    for (int ix=-85; ix<=-1; ix++)
    {
      printLine(idx, ix, iy, 0, A_ofile);
    }
    A_ofile.close();

    idx=0;
    sprintf(name, "emptyCalibTable_IZ0_IY%d_IX+.dat", iy);
    std::ofstream B_ofile(name);
    for (int ix=1; ix<=85; ix++)
    {
      printLine(idx, ix, iy, 0, B_ofile);
    }
    B_ofile.close();
  }

  // EB/EE+/EE-
  // to file
  std::vector<std::vector<int> > cells0 = GetAllCellsInEB();
  std::ofstream ofile_all("emptyCalibTable_All.dat");
  if (ofile_all.is_open())
  { 
    //EB
    for (int idx=0; idx<(int)cells0.size(); idx++)
    {
      printLine(idx, cells0.at(idx).at(0), cells0.at(idx).at(1), cells0.at(idx).at(2), ofile_all);
    }
    //EE+
    for (int idx=0; idx<(int)cellsp.size(); idx++)
    {
      printLine(idx, cellsp.at(idx).at(0), cellsp.at(idx).at(1), cellsp.at(idx).at(2), ofile_all);
    }
    //EE-
    for (int idx=0; idx<(int)cellsn.size(); idx++)
    {
      printLine(idx, cellsn.at(idx).at(0), cellsn.at(idx).at(1), cellsn.at(idx).at(2), ofile_all);
    }

    ofile_all.close();
  }
 

  // check
  TFile* fout = new TFile("printEmptyTable.root", "recreate");
  //EE
  TH2D* h2d = new TH2D("h2d", "h2d", 600, -10, 110, 600, -10, 110);
  h2d->Sumw2();
  for (int idx=0; idx<(int)cellsp.size(); idx++)
  {
    h2d->Fill(cellsp.at(idx).at(0), cellsp.at(idx).at(1));
  }
  h2d->Write();
  //EB
  TH2D* h2d0 = new TH2D("h2d0", "h2d0", 900, -90, 90, 1850, -5, 365);
  h2d0->Sumw2();
  for (int idx=0; idx<(int)cells0.size(); idx++)
  {
    h2d0->Fill(cells0.at(idx).at(0), cells0.at(idx).at(1));
  }
  h2d0->Write();

  fout->Close();
  
  return 0;
}


void printLine(int idx, int ix, int iy, int iz)
{
  std::cout << idx << " "
  << ix << " "
  << iy << " "
  << iz << " "
  << 1.0 << " "
  << 0.01 << " "
  << 0 << " "
  << 0 << " "
  << std::endl;
}

void printLine(int idx, int ix, int iy, int iz, std::ofstream& ofile)
{
  ofile << idx << " "
  << ix << " "
  << iy << " "
  << iz << " "
  << 1.0 << " "
  << 0.01 << " "
  << 0 << " "
  << 0 << " "
  << std::endl;
}

std::vector<std::vector<int> > GetAllCellsInEEP()
{
  std::vector<std::vector<int> > cells;
  
  int ix=0;
  int iy=0;
  int iz=+1;
  
  // print cells from ix left to right
  for (ix=1; ix<=100; ix++)
  {
    for (iy=1; iy<=100; iy++)
    {
      // a list of exclusions
      
      // outer
      // left top
      if (ix<4&&iy>60) continue;
      if (ix<6&&iy>65) continue;
      if (ix<9&&iy>75) continue;
      if (ix<14&&iy>80) continue;
      if (ix<16&&iy>85) continue;
      if (ix<21&&iy>87) continue;
      if (ix<26&&iy>92) continue;
      if (ix<36&&iy>95) continue;
      if (ix<41&&iy>97) continue;
      // left bottom
      if (ix<4&&iy<41) continue;
      if (ix<6&&iy<36) continue;
      if (ix<9&&iy<26) continue;
      if (ix<14&&iy<21) continue;
      if (ix<16&&iy<16) continue;
      if (ix<21&&iy<14) continue;
      if (ix<26&&iy<9) continue;
      if (ix<36&&iy<6) continue;
      if (ix<41&&iy<4) continue;
      // right top
      if (ix>60&&iy>97) continue;
      if (ix>65&&iy>95) continue;
      if (ix>75&&iy>92) continue;
      if (ix>80&&iy>87) continue;
      if (ix>85&&iy>85) continue;
      if (ix>87&&iy>80) continue;
      if (ix>92&&iy>75) continue;
      if (ix>95&&iy>65) continue;
      if (ix>97&&iy>60) continue;
      // right bottom
      if (ix>60&&iy<4) continue;
      if (ix>65&&iy<6) continue;
      if (ix>75&&iy<9) continue;
      if (ix>80&&iy<14) continue;
      if (ix>85&&iy<16) continue;
      if (ix>87&&iy<21) continue;
      if (ix>92&&iy<26) continue;
      if (ix>95&&iy<36) continue;
      if (ix>97&&iy<41) continue;

      // inner
      if (iy==40&&ix>45&&ix<56) continue;
      if (iy==41&&ix>43&&ix<58) continue;
      if (iy==42&&ix>42&&ix<59) continue;
      if (iy==43&&ix>41&&ix<60) continue;
      if (iy==44&&ix>40&&ix<61) continue;
      if (iy==45&&ix>40&&ix<61) continue;
      if (iy>45&&iy<56&&ix>39&&ix<62) continue;
      if (iy==56&&ix>40&&ix<61) continue;
      if (iy==57&&ix>40&&ix<61) continue;
      if (iy==58&&ix>41&&ix<60) continue;
      if (iy==59&&ix>42&&ix<59) continue;
      if (iy==60&&ix>43&&ix<58) continue;
      if (iy==61&&ix>45&&ix<56) continue;
      
      std::vector<int> acell;
      acell.push_back(ix);
      acell.push_back(iy);
      acell.push_back(iz);
      cells.push_back(acell);
    }
  }

  //return
  return cells;
  
}

std::vector<std::vector<int> > GetAllCellsInEEN()
{
  std::vector<std::vector<int> > cells;
  
  int ix=0;
  int iy=0;
  int iz=-1;
  
  // print cells from ix left to right
  
  for (ix=1; ix<=100; ix++)
  {
    for (iy=1; iy<=100; iy++)
    {
      // a list of exclusions
      
      // outer
      // left top
      if (ix<4&&iy>60) continue;
      if (ix<6&&iy>65) continue;
      if (ix<9&&iy>75) continue;
      if (ix<14&&iy>80) continue;
      if (ix<16&&iy>85) continue;
      if (ix<21&&iy>87) continue;
      if (ix<26&&iy>92) continue;
      if (ix<36&&iy>95) continue;
      if (ix<41&&iy>97) continue;
      // left bottom
      if (ix<4&&iy<41) continue;
      if (ix<6&&iy<36) continue;
      if (ix<9&&iy<26) continue;
      if (ix<14&&iy<21) continue;
      if (ix<16&&iy<16) continue;
      if (ix<21&&iy<14) continue;
      if (ix<26&&iy<9) continue;
      if (ix<36&&iy<6) continue;
      if (ix<41&&iy<4) continue;
      // right top
      if (ix>60&&iy>97) continue;
      if (ix>65&&iy>95) continue;
      if (ix>75&&iy>92) continue;
      if (ix>80&&iy>87) continue;
      if (ix>85&&iy>85) continue;
      if (ix>87&&iy>80) continue;
      if (ix>92&&iy>75) continue;
      if (ix>95&&iy>65) continue;
      if (ix>97&&iy>60) continue;
      // right bottom
      if (ix>60&&iy<4) continue;
      if (ix>65&&iy<6) continue;
      if (ix>75&&iy<9) continue;
      if (ix>80&&iy<14) continue;
      if (ix>85&&iy<16) continue;
      if (ix>87&&iy<21) continue;
      if (ix>92&&iy<26) continue;
      if (ix>95&&iy<36) continue;
      if (ix>97&&iy<41) continue;

      // inner
      if (iy==40&&ix>45&&ix<56) continue;
      if (iy==41&&ix>43&&ix<58) continue;
      if (iy==42&&ix>42&&ix<59) continue;
      if (iy==43&&ix>41&&ix<60) continue;
      if (iy==44&&ix>40&&ix<61) continue;
      if (iy==45&&ix>40&&ix<61) continue;
      if (iy>45&&iy<56&&ix>39&&ix<62) continue;
      if (iy==56&&ix>40&&ix<61) continue;
      if (iy==57&&ix>40&&ix<61) continue;
      if (iy==58&&ix>41&&ix<60) continue;
      if (iy==59&&ix>42&&ix<59) continue;
      if (iy==60&&ix>43&&ix<58) continue;
      if (iy==61&&ix>45&&ix<56) continue;
      
      std::vector<int> acell;
      acell.push_back(ix);
      acell.push_back(iy);
      acell.push_back(iz);
      cells.push_back(acell);
    }
  }

  //return
  return cells;
  
}

std::vector<std::vector<int> > GetAllCellsInEE()
{
  std::vector<std::vector<int> > cells;
  
  int ix=0;
  int iy=0;
  int iz=-1;
  
  for (iz=-1; iz<=+1; iz=iz+2)
  {
    // print cells from ix left to right
    for (ix=1; ix<=100; ix++)
    {
      for (iy=1; iy<=100; iy++)
      {
        // a list of exclusions
      
        // outer
        // left top
        if (ix<4&&iy>60) continue;
        if (ix<6&&iy>65) continue;
        if (ix<9&&iy>75) continue;
        if (ix<14&&iy>80) continue;
        if (ix<16&&iy>85) continue;
        if (ix<21&&iy>87) continue;
        if (ix<26&&iy>92) continue;
        if (ix<36&&iy>95) continue;
        if (ix<41&&iy>97) continue;
        // left bottom
        if (ix<4&&iy<41) continue;
        if (ix<6&&iy<36) continue;
        if (ix<9&&iy<26) continue;
        if (ix<14&&iy<21) continue;
        if (ix<16&&iy<16) continue;
        if (ix<21&&iy<14) continue; 
        if (ix<26&&iy<9) continue;
        if (ix<36&&iy<6) continue;
        if (ix<41&&iy<4) continue;
        // right top
        if (ix>60&&iy>97) continue;
        if (ix>65&&iy>95) continue;
        if (ix>75&&iy>92) continue;
        if (ix>80&&iy>87) continue;
        if (ix>85&&iy>85) continue;
        if (ix>87&&iy>80) continue;
        if (ix>92&&iy>75) continue;
        if (ix>95&&iy>65) continue;
        if (ix>97&&iy>60) continue;
        // right bottom
        if (ix>60&&iy<4) continue;
        if (ix>65&&iy<6) continue;
        if (ix>75&&iy<9) continue;
        if (ix>80&&iy<14) continue;
        if (ix>85&&iy<16) continue;
        if (ix>87&&iy<21) continue;
        if (ix>92&&iy<26) continue;
        if (ix>95&&iy<36) continue;
        if (ix>97&&iy<41) continue;

        // inner
        if (iy==40&&ix>45&&ix<56) continue;
        if (iy==41&&ix>43&&ix<58) continue;
        if (iy==42&&ix>42&&ix<59) continue;
        if (iy==43&&ix>41&&ix<60) continue;
        if (iy==44&&ix>40&&ix<61) continue;
        if (iy==45&&ix>40&&ix<61) continue;
        if (iy>45&&iy<56&&ix>39&&ix<62) continue;
        if (iy==56&&ix>40&&ix<61) continue;
        if (iy==57&&ix>40&&ix<61) continue;
        if (iy==58&&ix>41&&ix<60) continue;
        if (iy==59&&ix>42&&ix<59) continue;
        if (iy==60&&ix>43&&ix<58) continue;
        if (iy==61&&ix>45&&ix<56) continue;
      
        std::vector<int> acell;
        acell.push_back(ix);
        acell.push_back(iy);
        acell.push_back(iz);
        cells.push_back(acell);
      }
    }
  }
  //return
  return cells;
  
}

std::vector<std::vector<int> > GetAllCellsInEB()
{
  std::vector<std::vector<int> > cells;
  
  int ix=0;  // iEta
  int iy=0;  // iPhi
  int iz=0;  // 0 for EB

  // print cells from ix (iEta) left to right
  for (ix=-85; ix<=+85; ix++)
  {
    // ix (iEta) not to be 0, no iEta=0
    if (ix==0) continue;
     
    // loop over iy (iPhi) 
    for (iy=1; iy<=360; iy++)
    {
      std::vector<int> acell;
      acell.push_back(ix);
      acell.push_back(iy);
      acell.push_back(iz);
      cells.push_back(acell);
    }
  }  

  //return
  return cells;
  
}
