#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include "TROOT.h"
#include "TChain.h"
#include "TH1D.h"
#include "TMath.h"


typedef struct
{
  int idx; // index
  int ix;  // cell ix
  int iy;  // cell iy
  int iz;  // cell iz
  double c; // calibration constant
  double cerr; // error of calibratioin constant
  int fixed;  // true: fixed ; false: not fixed.
  int nfits; // n times of fits
} calibRecord;

typedef struct
{
  double min;
  double max;
  double s; // scale
  double serr; // error of scale
} EnergyScale;

// variables for ROOT Tree
Int_t           runNum;
ULong64_t       evtNum;
Int_t           lumBlk;
UInt_t          runTime;
Int_t           nVtx;
Int_t           tnPar;
Int_t           rStat[2];
Double_t        rR9[2];
Double_t        rPx[2];
Double_t        rPy[2];
Double_t        rPz[2];
Double_t        rPt[2];
Double_t        rE[2];
Double_t        rEta[2];
Double_t        rPhi[2];
Double_t        rVtx[2];
Double_t        rVty[2];
Double_t        rVtz[2];
Double_t        rERaw[2];
Double_t        rNBCl[2];
Double_t        rEtaWidth[2];
Double_t        rPhiWidth[2];
Float_t         rCaloE[2];
Float_t         rEcalE[2];
Float_t         rPresE[2];
Int_t           rEB[2];
Float_t         rSeedE[2];
Int_t           rSeedIX[2];
Int_t           rSeedIY[2];
Int_t           rSeedIZ[2];
Int_t           rNHits[2];
Float_t         rHitE[2][200];
Int_t           rHitIX[2][200];
Int_t           rHitIY[2][200];
Int_t           rHitIZ[2][200];
Int_t           rZStat;

// newly added variables
Int_t           rHLTFire;
Int_t           rNPV;
Int_t           rEleID[2];
Float_t         rERegV8Elec[2];
Float_t         rERegV8Phot[2];
Float_t         rESigmaRegV8Elec[2];
Float_t         rESigmaRegV8Phot[2];


void printCalibRecord(calibRecord& record)
{
  std::cout << record.idx << " "
  << record.ix << " "
  << record.iy << " "
  << record.iz << " "
  << record.c << " "
  << record.cerr << " "
  << record.fixed << " "
  << record.nfits << " "
  << std::endl;
}

std::vector<calibRecord> getCalibTableFromFile(const char* filename)
{
  std::vector<calibRecord> mycalibTable;
  
  std::ifstream myfile(filename);
  std::string line;
  if (myfile.is_open())
  {
    mycalibTable.clear();
    while (getline(myfile,line))
    {
      calibRecord calib;
      std::stringstream sline(line);
      sline >> calib.idx
      >> calib.ix
      >> calib.iy
      >> calib.iz
      >> calib.c
      >> calib.cerr
      >> calib.fixed
      >> calib.nfits ;
      mycalibTable.push_back(calib);
    }
    myfile.close();
  }
  return mycalibTable;
}

void printCalibTableToFile(const std::vector<calibRecord>& calibTable, const char* filename = "calibTable_out.dat")
{
  std::ofstream myfile(filename);
  if (myfile.is_open())
  {
    for (int i=0; i<(int)calibTable.size(); i++)
    {
      myfile << calibTable.at(i).idx << " "
      << calibTable.at(i).ix << " "
      << calibTable.at(i).iy<< " "
      << calibTable.at(i).iz<< " "
      << calibTable.at(i).c << " "
      << calibTable.at(i).cerr << " "
      << calibTable.at(i).fixed << " "
      << calibTable.at(i).nfits << " "
      << std::endl;
    }
    myfile.close();
  }
}

void printCalibTableToScreen(const std::vector<calibRecord>& calibTable)
{
  for (int i=0; i<(int)calibTable.size(); i++)
  {
    std::cout << calibTable.at(i).idx << " "
    << calibTable.at(i).ix << " "
    << calibTable.at(i).iy<< " "
    << calibTable.at(i).iz<< " "
    << calibTable.at(i).c << " "
    << calibTable.at(i).cerr << " "
    << calibTable.at(i).fixed << " "
    << calibTable.at(i).nfits << " "
    << std::endl;
  }

}


// set tree branches
void SetTreeBranch(TChain* tree)
{
  // Set branch addresses.
  tree->SetBranchAddress("runNum",&runNum);
  tree->SetBranchAddress("evtNum",&evtNum);
  tree->SetBranchAddress("lumBlk",&lumBlk);
  tree->SetBranchAddress("runTime",&runTime);
  tree->SetBranchAddress("nVtx",&nVtx);
  tree->SetBranchAddress("tnPar",&tnPar);
  tree->SetBranchAddress("rStat",rStat);
  tree->SetBranchAddress("rR9",rR9);
  tree->SetBranchAddress("rPx",rPx);
  tree->SetBranchAddress("rPy",rPy);
  tree->SetBranchAddress("rPz",rPz);
  tree->SetBranchAddress("rPt",rPt);
  tree->SetBranchAddress("rE",rE);
  tree->SetBranchAddress("rEta",rEta);
  tree->SetBranchAddress("rPhi",rPhi);
  tree->SetBranchAddress("rVtx",rVtx);
  tree->SetBranchAddress("rVty",rVty);
  tree->SetBranchAddress("rVtz",rVtz);
  tree->SetBranchAddress("rERaw",rERaw);
  tree->SetBranchAddress("rNBCl",rNBCl);
  tree->SetBranchAddress("rEtaWidth",rEtaWidth);
  tree->SetBranchAddress("rPhiWidth",rPhiWidth);
  tree->SetBranchAddress("rCaloE",rCaloE);
  tree->SetBranchAddress("rEcalE",rEcalE);
  tree->SetBranchAddress("rPresE",rPresE);
  tree->SetBranchAddress("rEB",rEB);
  tree->SetBranchAddress("rSeedE",rSeedE);
  tree->SetBranchAddress("rSeedIX",rSeedIX);
  tree->SetBranchAddress("rSeedIY",rSeedIY);
  tree->SetBranchAddress("rSeedIZ",rSeedIZ);
  tree->SetBranchAddress("rNHits",rNHits);
  tree->SetBranchAddress("rHitE",rHitE);
  tree->SetBranchAddress("rHitIX",rHitIX);
  tree->SetBranchAddress("rHitIY",rHitIY);
  tree->SetBranchAddress("rHitIZ",rHitIZ);
  tree->SetBranchAddress("rZStat",&rZStat);  
  // newly added branch
  tree->SetBranchAddress("rHLTFire",&rHLTFire);
  tree->SetBranchAddress("rNPV",&rNPV);
  tree->SetBranchAddress("rEleID",rEleID);
  tree->SetBranchAddress("rERegV8Elec",rERegV8Elec);
  tree->SetBranchAddress("rERegV8Phot",rERegV8Phot);
  tree->SetBranchAddress("rESigmaRegV8Elec",rESigmaRegV8Elec);
  tree->SetBranchAddress("rESigmaRegV8Phot",rESigmaRegV8Phot);
}



double calcMass(double E1, double Eta1, double Phi1,
                double E2, double Eta2, double Phi2)
{
  // direction
  double Theta1 = 2.0*atan(exp(-Eta1));
  double px1 = E1*sin(Theta1)*cos(Phi1);
  double py1 = E1*sin(Theta1)*sin(Phi1);
  double pz1 = E1*cos(Theta1);
  double Theta2 = 2.0*atan(exp(-Eta2));
  double px2 = E2*sin(Theta2)*cos(Phi2);
  double py2 = E2*sin(Theta2)*sin(Phi2);
  double pz2 = E2*cos(Theta2);
  // mass
  double px = px1+px2;
  double py = py1+py2;
  double pz = pz1+pz2;
  double E = E1+E2;
  double mass = sqrt(E*E-(px*px+py*py+pz*pz));
  return mass;
}

double getCalibConstFromCalibTable(int ix, int iy, int iz, const std::vector<calibRecord> calibTable)
{
  double calibConst(-100);
  // find the index from the calibTable
  for(int i=0; i<(int)calibTable.size(); i++)
  {
    if (ix==calibTable.at(i).ix && 
        iy==calibTable.at(i).iy &&
        iz==calibTable.at(i).iz)
    {
      calibConst = calibTable.at(i).c;
      return calibConst;
    }
  }
  //
  // Note here, if cannot find an entry in the calibTable,
  // the preset number -100 will be returned.
  return calibConst;
}


double calibCell(double E, int ix, int iy, int iz, const std::vector<calibRecord>& calibTable)
{
  double calibConst = getCalibConstFromCalibTable(ix, iy, iz, calibTable);
  // if cannot find the cell return its original energy
  if (calibConst<=0.0) {
    return E;
  }
  // if good, return the calibrated energy
  return calibConst*E;
}

  
double GetCalibConstFromCalibTable(const int& ix, const int& iy, const int& iz, 
      const std::map<int, std::map<int, std::map<int, calibRecord> > >& calibMap)
{
  std::map<int, std::map<int, std::map<int, calibRecord> > >::const_iterator _it_ix;
  std::map<int, std::map<int, calibRecord> >::const_iterator _it_iy;
  std::map<int, calibRecord>::const_iterator _it_iz;
  // check ix
  _it_ix = calibMap.find(ix);
  if (_it_ix == calibMap.end()) return -1000.0 ;
  // check iy
  _it_iy = (_it_ix->second).find(iy);
  if (_it_iy == (_it_ix->second).end()) return -1000.0;
  // check iz
  _it_iz = (_it_iy->second).find(iz);
  if (_it_iz == (_it_iy->second).end()) return -1000.0;
  // 
  return _it_iz->second.c ;

}

TH1D* getTH1DMeeWithEtaScale(TChain* tree,
                 const std::vector<calibRecord>& calibTable,
                 const std::vector<EnergyScale>& EtaScale,
                 const char* histname = "hist",
                 const int nbins=60,
                 const int maxevt=-1,
                 const int method=7,
                 const double scale=1.0)
{

  // new a histogram
  TH1D* hist = new TH1D(histname, histname, nbins, 60,120);
  hist->Sumw2();
  hist->GetXaxis()->SetTitle("M(ee) (GeV)");
  hist->GetYaxis()->SetTitle("nEvents");
  hist->SetMarkerStyle(20);
  hist->SetTitle(histname);

  // re-format calibTable to calibMap for fast access
  std::map<int, std::map<int, std::map<int, calibRecord> > > calibMap;
  // directly refer to the one that is given
  for (int idx=0; idx<(int)calibTable.size(); idx++)
  {
    calibRecord record = calibTable.at(idx);
    calibMap[record.ix][record.iy][record.iz] = record;
  }

  //
  int nevts = (int)tree->GetEntries();

  // fill data vectors
  for (int i=0; i<nevts; i++)
  {
    // get event from tree
    tree->GetEntry(i);

    // loose+fid cut
    if (rEleID[0]<3||rEleID[1]<3) continue;

    // Energy
    double E1 = rERaw[0]+rPresE[0]; //rE[0];
    double E2 = rERaw[1]+rPresE[1]; //rE[1];

    if (method==7) 
    {
      // similar to method2, but on top of energy regression

      // regression energy
      double EReg1 = rERegV8Elec[0];
      double EReg2 = rERegV8Elec[1];
      double RegScale1 = EReg1/E1;
      double RegScale2 = EReg2/E2;

      // recalib E1
      for (int ii=0; ii<rNHits[0]; ii++)
      {
        double CalibC = GetCalibConstFromCalibTable(rHitIX[0][ii], rHitIY[0][ii], rHitIZ[0][ii], calibMap);
        double Eold = rHitE[0][ii];
        double Enew = Eold * CalibC * scale;

        if (CalibC>-900.0)
        {
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E1 = E1 - Eold + Enew;
        }
      }

      // recalib E2
      for (int ii=0; ii<rNHits[1]; ii++)
      {
        double CalibC = GetCalibConstFromCalibTable(rHitIX[1][ii], rHitIY[1][ii], rHitIZ[1][ii], calibMap);
        double Eold = rHitE[1][ii];
        double Enew = Eold * CalibC * scale;

        if (CalibC>-900.0)
        {
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E2 = E2 - Eold + Enew;
        }
      }

      // apply regression scale
      E1 = E1 * RegScale1;
      E2 = E2 * RegScale2;

      // eta scale
      double scale1(1.0), scale2(1.0);
      for (int ibin=0; ibin<(int)EtaScale.size(); ibin++)
      {
        if (rEta[0]>EtaScale.at(ibin).min&&rEta[0]<EtaScale.at(ibin).max) scale1 = EtaScale.at(ibin).s;
        if (rEta[1]>EtaScale.at(ibin).min&&rEta[1]<EtaScale.at(ibin).max) scale2 = EtaScale.at(ibin).s;
      }

      // calculate mass
      double mee = calcMass(E1*scale1, rEta[0], rPhi[0], E2*scale2, rEta[1], rPhi[1]);

      // fill
      hist->Fill(mee);
    }

    // print
    if(i%10000==0)
    {
      std::cout << "Events " << i << std::endl;
    }

    //
    if (maxevt!=-1&&i>maxevt)
    {
      break;
    }
  }

  return hist;
}
//////////


TH1D* getTH1DMee(TChain* tree, 
                 const std::vector<calibRecord>& calibTable, 
                 const char* histname = "hist", 
                 const int nbins=60,
                 const int maxevt=-1,
                 const int method=2,
                 const double scale=1.0,
                 const int oddeven=0)
{
  
  // new a histogram
  TH1D* hist = new TH1D(histname, histname, nbins, 60,120);
  hist->Sumw2();
  hist->GetXaxis()->SetTitle("M(ee) (GeV)");
  hist->GetYaxis()->SetTitle("nEvents");
  hist->SetMarkerStyle(20);
  hist->SetTitle(histname);
  
  // re-format calibTable to calibMap for fast access
  std::map<int, std::map<int, std::map<int, calibRecord> > > calibMap;
  // directly refer to the one that is given
  for (int idx=0; idx<(int)calibTable.size(); idx++)
  {
    calibRecord record = calibTable.at(idx);
    calibMap[record.ix][record.iy][record.iz] = record;
  }
    
  //
  int nevts = (int)tree->GetEntries();
  // fill data vectors
  for (int i=0; i<nevts; i++)
  {
     // even odd event check
    if (oddeven==1 && i%2==0) continue;
    else if (oddeven==2 && i%2==1) continue;

    // get event from tree
    tree->GetEntry(i);

    // loose+fid cut
    if (rEleID[0]<3||rEleID[1]<3) continue;
    
    // Energy
    double E1 = rERaw[0]+rPresE[0]; //rE[0];
    double E2 = rERaw[1]+rPresE[1]; //rE[1];
    
    
    if (method==2) 
    {
      // recalib E1
      for (int ii=0; ii<rNHits[0]; ii++)
      {
        double CalibC = GetCalibConstFromCalibTable(rHitIX[0][ii], rHitIY[0][ii], rHitIZ[0][ii], calibMap);
        double Eold = rHitE[0][ii];
        double Enew = Eold * CalibC * scale;
      
        if (CalibC>-900.0)
        {
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E1 = E1 - Eold + Enew;
        }
      }
    
      // recalib E2
      for (int ii=0; ii<rNHits[1]; ii++)
      {
        double CalibC = GetCalibConstFromCalibTable(rHitIX[1][ii], rHitIY[1][ii], rHitIZ[1][ii], calibMap);
        double Eold = rHitE[1][ii];
        double Enew = Eold * CalibC * scale;
      
        if (CalibC>-900.0)
        {
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E2 = E2 - Eold + Enew;
        }
      }
      
      // calculate mass
      double mee = calcMass(E1, rEta[0], rPhi[0], E2, rEta[1], rPhi[1]);
      
      // fill
      hist->Fill(mee);
    
    }
    if (method==23)
    {
      // similar to method2, but on top of energy regression

      // regression energy
      double EReg1 = rERegV8Elec[0];
      double EReg2 = rERegV8Elec[1];
      double RegScale1 = EReg1/E1;
      double RegScale2 = EReg2/E2;

      // recalib E1
      for (int ii=0; ii<rNHits[0]; ii++)
      {
        double CalibC = GetCalibConstFromCalibTable(rHitIX[0][ii], rHitIY[0][ii], rHitIZ[0][ii], calibMap);
        double Eold = rHitE[0][ii];
        double Enew = Eold * CalibC * scale;

        if (CalibC>-900.0)
        {
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E1 = E1 - Eold + Enew;
        }
      }

      // recalib E2
      for (int ii=0; ii<rNHits[1]; ii++)
      {
        double CalibC = GetCalibConstFromCalibTable(rHitIX[1][ii], rHitIY[1][ii], rHitIZ[1][ii], calibMap);
        double Eold = rHitE[1][ii];
        double Enew = Eold * CalibC * scale;

        if (CalibC>-900.0)
        {
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E2 = E2 - Eold + Enew;
        }
      }

      // apply regression scale
      E1 = E1 * RegScale1;
      E2 = E2 * RegScale2;

      // calculate mass
      double mee = calcMass(E1, rEta[0], rPhi[0], E2, rEta[1], rPhi[1]);

      // fill
      hist->Fill(mee);

    }
    else if (method==3)
    {
      // method 3:
      // only recalculate one of the two
      // even ievt for E1, odd ievt for E2
      // the other take GsfElectron
      // take E1  or E2
      bool takeE1(false);
      bool takeE2(false);
      // even or odd
      if (i%2==0) takeE2 = true;
      else takeE1 = true;
      
      // dE1 = Sum(-Eold + Enew)
      double dE1(0);
      
      // recalib E1
      if (takeE1) 
      for (int ii=0; ii<rNHits[0]; ii++)
      {
        double CalibC = GetCalibConstFromCalibTable(rHitIX[0][ii], rHitIY[0][ii], rHitIZ[0][ii], calibMap);
        double Eold = rHitE[0][ii];
        double Enew = Eold * CalibC * scale;
      
        if (CalibC>-900.0)
        {
          //calculate new Energy as :
          // dEcorr = 0 - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          dE1 = dE1 - Eold + Enew;
        }
      }
    
      // check if take E1
      if (takeE1) E1 = E1 + dE1;
                
                
      // dE2 = Sum(-Eold + Enew)
      double dE2(0);   
    
      // recalib E2
      if (takeE2)
      for (int ii=0; ii<rNHits[1]; ii++)
      {
        double CalibC = GetCalibConstFromCalibTable(rHitIX[1][ii], rHitIY[1][ii], rHitIZ[1][ii], calibMap);
        double Eold = rHitE[1][ii];
        double Enew = Eold * CalibC * scale;
      
        if (CalibC>-900.0)
        {
          //calculate new Energy as :
          // dEcorr = 0 - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          dE2 = dE2 - Eold + Enew;
        }
      }
      
      // check if take E2
      if (takeE2) E2 = E2 + dE2;
      
      // calculate mass
      double mee = calcMass(E1, rEta[0], rPhi[0], E2, rEta[1], rPhi[1]);
      
      // fill
      hist->Fill(mee);
      
    }
    else if (method==31)
    {
      // method 31:
      // one recalculate one GsfElectron
      // repeat for the two
      // i.e. fill two Mee for each event for the two possible combinations.

      
      // dE1 = Sum(-Eold + Enew)
      double dE1(0);
      
      // recalib E1
      for (int ii=0; ii<rNHits[0]; ii++)
      {
        double CalibC = GetCalibConstFromCalibTable(rHitIX[0][ii], rHitIY[0][ii], rHitIZ[0][ii], calibMap);
        double Eold = rHitE[0][ii];
        double Enew = Eold * CalibC * scale;
      
        if (CalibC>-900.0)
        {
          //calculate new Energy as :
          // dEcorr = 0 - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          dE1 = dE1 - Eold + Enew;
        }
      }
    
      // calculate mass for recalib E1 + GsfElectron of E2
      double mee1 = calcMass(E1+dE1, rEta[0], rPhi[0], E2, rEta[1], rPhi[1]);
      
      // fill
      hist->Fill(mee1);
                
                
      // dE2 = Sum(-Eold + Enew)
      double dE2(0);   
    
      // recalib E2
      for (int ii=0; ii<rNHits[1]; ii++)
      {
        double CalibC = GetCalibConstFromCalibTable(rHitIX[1][ii], rHitIY[1][ii], rHitIZ[1][ii], calibMap);
        double Eold = rHitE[1][ii];
        double Enew = Eold * CalibC * scale;
      
        if (CalibC>-900.0)
        {
          //calculate new Energy as :
          // dEcorr = 0 - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          dE2 = dE2 - Eold + Enew;
        }
      }
      
      // calculate mass for recalib E2 + GsfElectron of E1
      double mee2 = calcMass(E1, rEta[0], rPhi[0], E2+dE2, rEta[1], rPhi[1]);
      
      // fill
      hist->Fill(mee2);
      
    }

    
    // print
    if(i%10000==0)
    {
      std::cout << "Events " << i << std::endl;
    }
    
    //
    if (maxevt!=-1&&i>maxevt)
    {
      break;
    }
  }
  
  return hist;
}




TH1D* getTH1DOriginalMee(TChain* tree, 
                         const char* histname = "hist", 
                         const int nbins=60, 
                         const int maxevt=-1,
                         const double scale=1.0,
                         const int oddeven=0)
{
  
  // new a histogram
  TH1D* hist = new TH1D(histname, histname, nbins, 60,120);
  hist->Sumw2();
  hist->GetXaxis()->SetTitle("M(ee) (GeV)");
  hist->GetYaxis()->SetTitle("nEvents");
  hist->SetMarkerStyle(20);
  hist->SetTitle(histname);
  
  int nevts = (int)tree->GetEntries();
  // fill data vectors
  for (int i=0; i<nevts; i++){

    // even odd event check
    if (oddeven==1 && i%2==0) continue;
    else if (oddeven==2 && i%2==1) continue;

    // get event from tree
    tree->GetEntry(i);
    
    // loose+fid cut
    if (rEleID[0]<3||rEleID[1]<3) continue;

    // calculate mass
    double mee = calcMass(rE[0], rEta[0], rPhi[0], rE[1], rEta[1], rPhi[1]);
    
    // apply scale
    mee = mee*scale;
    
    // fill
    hist->Fill(mee);
    
    //
    if (maxevt!=-1&&i>maxevt)
    {
      break;
    }
  }
  
  return hist;
}


TH1D* getTH1DMeeRegV8ElecWithEtaScale(TChain* tree,
                         const std::vector<EnergyScale>& EtaScale,
                         const char* histname = "hist",
                         const int nbins=60,
                         const int maxevt=-1,
                         const double scale=1.0)
{

  // new a histogram
  TH1D* hist = new TH1D(histname, histname, nbins, 60,120);
  hist->Sumw2();
  hist->GetXaxis()->SetTitle("M(ee) (GeV)");
  hist->GetYaxis()->SetTitle("nEvents");
  hist->SetMarkerStyle(20);
  hist->SetTitle(histname);

  int nevts = (int)tree->GetEntries();
  // fill data vectors
  for (int i=0; i<nevts; i++){
    // get event from tree
    tree->GetEntry(i);

    // loose+fid cut
    if (rEleID[0]<3||rEleID[1]<3) continue;

    double scale1(1.0), scale2(1.0);
    for (int ibin=0; ibin<(int)EtaScale.size(); ibin++)
    {
      if (rEta[0]>EtaScale.at(ibin).min&&rEta[0]<EtaScale.at(ibin).max) scale1 = EtaScale.at(ibin).s; 
      if (rEta[1]>EtaScale.at(ibin).min&&rEta[1]<EtaScale.at(ibin).max) scale2 = EtaScale.at(ibin).s; 
    }

    // calculate mass
    double mee = calcMass(rERegV8Elec[0]*scale1, rEta[0], rPhi[0], rERegV8Elec[1]*scale2, rEta[1], rPhi[1]);

    mee = mee*scale;

    // fill
    hist->Fill(mee);

    //
    if (maxevt!=-1&&i>maxevt)
    {
      break;
    }
  }

  return hist;
}
///////////

TH1D* getTH1DMeeRegV8Elec(TChain* tree,
                         const char* histname = "hist",
                         const int nbins=60,
                         const int maxevt=-1,
                         const double scale=1.0,
                         const int oddeven=0)
{

  // new a histogram
  TH1D* hist = new TH1D(histname, histname, nbins, 60,120);
  hist->Sumw2();
  hist->GetXaxis()->SetTitle("M(ee) (GeV)");
  hist->GetYaxis()->SetTitle("nEvents");
  hist->SetMarkerStyle(20);
  hist->SetTitle(histname);

  int nevts = (int)tree->GetEntries();
  // fill data vectors
  for (int i=0; i<nevts; i++){

    // even odd event check
    if (oddeven==1 && i%2==0) continue;
    else if (oddeven==2 && i%2==1) continue;

    // get event from tree
    tree->GetEntry(i);

    // loose+fid cut
    if (rEleID[0]<3||rEleID[1]<3) continue;

    // calculate mass
    double mee = calcMass(rERegV8Elec[0], rEta[0], rPhi[0], rERegV8Elec[1], rEta[1], rPhi[1]);

    // apply scale
    mee = mee*scale;

    // fill
    hist->Fill(mee);

    //
    if (maxevt!=-1&&i>maxevt)
    {
      break;
    }
  }

  return hist;
}

