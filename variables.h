#include "TROOT.h"
#include "TChain.h"
#include "Minuit2/MnUserParameters.h"
#include <string>
#include <vector>
#include <utility>
#include <fstream>


typedef struct
{
  double min;
  double max;
  double s; // scale
  double serr; // error of scale
} EnergyScale;


// define vectors to store data

// store all events
int nEventsAll, nSignalsAll;
std::vector<double> allE1, allEReg1, allEta1, allPhi1;
std::vector<double> allE2, allEReg2, allEta2, allPhi2;
std::vector<int> allnHits1, allnHits2;
std::vector< std::vector<double> > allHitE1, allHitE2;
std::vector< std::vector<int> > allHitIX1, allHitIY1, allHitIZ1;
std::vector< std::vector<int> > allHitIX2, allHitIY2, allHitIZ2;

// store selected events
int nEvents, nSignals;
std::vector<double*> E1, EReg1, Eta1, Phi1;
std::vector<double*> E2, EReg2, Eta2, Phi2;
std::vector<int*> nHits1, nHits2;
std::vector< std::vector<double>* > HitE1, HitE2;
std::vector< std::vector<int>* > HitIX1, HitIY1, HitIZ1;
std::vector< std::vector<int>* > HitIX2, HitIY2, HitIZ2;
std::vector<bool> UseEle1, UseEle2;
std::vector<int> ScaleBin1, ScaleBin2;

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
Float_t rERegV7Elec[2];
Float_t rERegV7Phot[2];
Float_t rESigmaRegV7Elec[2];
Float_t rESigmaRegV7Phot[2];
Float_t rERegV6Elec[2];
Float_t rERegV6Phot[2];
Float_t rESigmaRegV6Elec[2];
Float_t rESigmaRegV6Phot[2];
Float_t rERegV5Elec[2];
Float_t rERegV5Phot[2];
Float_t rESigmaRegV5Elec[2];
Float_t rESigmaRegV5Phot[2];

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

  tree->SetBranchAddress("rERegV7Elec",rERegV7Elec);
  tree->SetBranchAddress("rERegV7Phot",rERegV7Phot);
  tree->SetBranchAddress("rESigmaRegV7Elec",rESigmaRegV7Elec);
  tree->SetBranchAddress("rESigmaRegV7Phot",rESigmaRegV7Phot);

  tree->SetBranchAddress("rERegV6Elec",rERegV6Elec);
  tree->SetBranchAddress("rERegV6Phot",rERegV6Phot);
  tree->SetBranchAddress("rESigmaRegV6Elec",rESigmaRegV6Elec);
  tree->SetBranchAddress("rESigmaRegV6Phot",rESigmaRegV6Phot);

  tree->SetBranchAddress("rERegV5Elec",rERegV5Elec);
  tree->SetBranchAddress("rERegV5Phot",rERegV5Phot);
  tree->SetBranchAddress("rESigmaRegV5Elec",rESigmaRegV5Elec);
  tree->SetBranchAddress("rESigmaRegV5Phot",rESigmaRegV5Phot);

}

// store all events from the Tree to vectors
void FillAllEvents(TChain* tree, const int debug=0, const std::string regVersion="V8Elec", const bool fitscale=false)
{
  nEventsAll = (int)tree->GetEntries();
  nSignalsAll = 0.99999*nEventsAll; // guess a number
  if (debug>1) std::cout << "fill all data vectors: Using Regression "<< regVersion << std::endl;
  // fill data vectors
  for (int i=0; i<nEventsAll; i++){
    // get event from tree
    tree->GetEntry(i);
    
    // first electron
    //allE1.push_back(rE[0]);
    allE1.push_back(rERaw[0]+rPresE[0]); // SC raw as in Regression
    if(regVersion=="V5Elec") allEReg1.push_back(rERegV5Elec[0]);
    else if(regVersion=="V6Elec") allEReg1.push_back(rERegV6Elec[0]);
    else if(regVersion=="V7Elec") allEReg1.push_back(rERegV7Elec[0]);
    else allEReg1.push_back(rERegV8Elec[0]);    
    allEta1.push_back(rEta[0]);
    allPhi1.push_back(rPhi[0]);
    // if fit scale, then we don't need the following information
    if(!fitscale) 
    {
      allnHits1.push_back(rNHits[0]);
      std::vector<double> hitE1;
      std::vector<int> hitIX1, hitIY1, hitIZ1;
      for (int ii=0; ii<rNHits[0]; ii++){
        hitE1.push_back(rHitE[0][ii]);
        hitIX1.push_back(rHitIX[0][ii]);
        hitIY1.push_back(rHitIY[0][ii]);
        hitIZ1.push_back(rHitIZ[0][ii]);
      }
      allHitE1.push_back(hitE1);
      allHitIX1.push_back(hitIX1);
      allHitIY1.push_back(hitIY1);
      allHitIZ1.push_back(hitIZ1);
    }
    // second electron
    //allE2.push_back(rE[1]);
    allE2.push_back(rERaw[1]+rPresE[1]); // SC raw as in Regression
    if(regVersion=="V5Elec") allEReg2.push_back(rERegV5Elec[1]);
    else if(regVersion=="V6Elec") allEReg2.push_back(rERegV6Elec[1]);
    else if(regVersion=="V7Elec") allEReg2.push_back(rERegV7Elec[1]);
    else allEReg2.push_back(rERegV8Elec[1]);
    allEta2.push_back(rEta[1]);
    allPhi2.push_back(rPhi[1]);
    // if fit scale, then we don't need the following information
    if(!fitscale)
    {
      allnHits2.push_back(rNHits[1]);
      std::vector<double> hitE2;
      std::vector<int> hitIX2, hitIY2, hitIZ2;
      for (int ii=0; ii<rNHits[1]; ii++){
        hitE2.push_back(rHitE[1][ii]);
        hitIX2.push_back(rHitIX[1][ii]);
        hitIY2.push_back(rHitIY[1][ii]);
        hitIZ2.push_back(rHitIZ[1][ii]);
      }
      allHitE2.push_back(hitE2);
      allHitIX2.push_back(hitIX2);
      allHitIY2.push_back(hitIY2);
      allHitIZ2.push_back(hitIZ2);
    }
  }
}


// generate an empty calib table
void GenEmptyCalibTable(std::vector<std::vector<int> > cells, const char* filename = "EmptyCalibTable.dat", int start_idx=0)
{
  // init txt file and print
  std::ofstream myfile(filename);
  if (myfile.is_open())
  {
    for (int i=0; i<(int)cells.size(); i++)
    {
      myfile << i << " "
      << cells.at(i).at(0) << " "
      << cells.at(i).at(1) << " "
      << cells.at(i).at(2) << " "
      << 1.0 << " "
      << 0.1 << " "
      << 0 << " "
      << 0 << " "
      << std::endl;
    }
    myfile.close();
  }
}

// this one is store the selection in E1, Eta1, .. E2, Eta2, .. HitE1, ... vectors.
// select events with at least one hit falling in one particular cell
int SelectEventsInOneCell(int ix, int iy, int iz)
{
  // clear previous vectors
  nEvents = 0;
  nSignals = 0;
  E1.clear();
  EReg1.clear();
  Eta1.clear();
  Phi1.clear();
  E2.clear();
  EReg2.clear();
  Eta2.clear();
  Phi2.clear();
  nHits1.clear();
  nHits2.clear();
  HitE1.clear();
  HitE2.clear();
  HitIX1.clear();
  HitIY1.clear();
  HitIZ1.clear();
  HitIX2.clear();
  HitIY2.clear();
  HitIZ2.clear();
  
  // loop over all events and select events
  for (int i=0; i<nEventsAll; i++)
  {
    bool take=false;
    // loop over all hits of electron1
    for (int ii=0; ii<(int)allHitE1.at(i).size(); ii++ )
    {
      if(take==true) break;
      
      if((allHitIX1.at(i).at(ii)==ix) &&
         (allHitIY1.at(i).at(ii)==iy) &&
         (allHitIZ1.at(i).at(ii)==iz) )
      {
        // if there is even one hit match a cell in selection, take the event
        take = true;
        break;
      }
    }
    
    // loop over all hits of electron2
    for (int ii=0; ii<(int)allHitE2.at(i).size(); ii++ )
    {
      if(take==true) break;
      
      if((allHitIX2.at(i).at(ii)==ix) &&
         (allHitIY2.at(i).at(ii)==iy) &&
         (allHitIZ2.at(i).at(ii)==iz) )
      {
        // if there is even one hit match a cell in selection, take the event
        take = true;
        break;
      }
    }
    // if do not decide to take this event, continue
    if (!take)
    {
      continue;
    }
    
    // if not continue above, it is a useful event to use
    E1.push_back(&(allE1.at(i)));
    EReg1.push_back(&(allEReg1.at(i)));
    Eta1.push_back(&(allEta1.at(i)));
    Phi1.push_back(&(allPhi1.at(i)));
    E2.push_back(&(allE2.at(i)));
    EReg2.push_back(&(allEReg2.at(i)));
    Eta2.push_back(&(allEta2.at(i)));
    Phi2.push_back(&(allPhi2.at(i)));
    nHits1.push_back(&(allnHits1.at(i)));
    HitE1.push_back(&(allHitE1.at(i)));
    HitIX1.push_back(&(allHitIX1.at(i)));
    HitIY1.push_back(&(allHitIY1.at(i)));
    HitIZ1.push_back(&(allHitIZ1.at(i)));
    nHits2.push_back(&(allnHits2.at(i)));
    HitE2.push_back(&(allHitE2.at(i)));
    HitIX2.push_back(&(allHitIX2.at(i)));
    HitIY2.push_back(&(allHitIY2.at(i)));
    HitIZ2.push_back(&(allHitIZ2.at(i)));
  }
  
  // nEvents
  nEvents = (int)E1.size();
  nSignals = nEvents; // assume no background (fix me)
  //
  return nEvents;
}



// this one is store the selection in E1, Eta1, .. E2, Eta2, .. HitE1, ... vectors.
// select events with one hit falling in one particular cell
// and the energy of this hit should be greater than a fraction
int SelectEventsInOneCellWithFraction(int ix, int iy, int iz, double fraction=0.5, int doEvenOdd=0)
{
  // clear previous vectors
  nEvents = 0;
  nSignals = 0;
  E1.clear();
  EReg1.clear();
  Eta1.clear();
  Phi1.clear();
  E2.clear();
  EReg2.clear();
  Eta2.clear();
  Phi2.clear();
  nHits1.clear();
  nHits2.clear();
  HitE1.clear();
  HitE2.clear();
  HitIX1.clear();
  HitIY1.clear();
  HitIZ1.clear();
  HitIX2.clear();
  HitIY2.clear();
  HitIZ2.clear();
  
  // loop over all events and select events
  for (int i=0; i<nEventsAll; i++)
  {
    bool take=false;
    
    // even odd event check
    if (doEvenOdd==1 && i%2==0) continue;
    else if (doEvenOdd==2 && i%2==1) continue; 
    
    // loop over all hits of electron1
    for (int ii=0; ii<(int)allHitE1.at(i).size(); ii++ )
    {
      if(take==true) break;
      
      // if there is even one hit match a cell in selection,
      // also this hit energy is bigger than fraction*E(total),
      // take this event
      if((allHitIX1.at(i).at(ii)==ix) &&
         (allHitIY1.at(i).at(ii)==iy) &&
         (allHitIZ1.at(i).at(ii)==iz) &&
         (allHitE1.at(i).at(ii)>fraction*allE1.at(i)) )
      {

        take = true;
        break;
      }
    }
    
    // loop over all hits of electron2
    for (int ii=0; ii<(int)allHitE2.at(i).size(); ii++ )
    {
      if(take==true) break;
      // if there is even one hit match a cell in selection,
      // also this hit energy is bigger than fraction*E(total),
      // take this event
      if((allHitIX2.at(i).at(ii)==ix) &&
         (allHitIY2.at(i).at(ii)==iy) &&
         (allHitIZ2.at(i).at(ii)==iz) &&
         (allHitE2.at(i).at(ii)>fraction*allE2.at(i)) )
      {
        take = true;
        break;
      }
    }
    // if do not decide to take this event, continue
    if (!take)
    {
      continue;
    }
    
    // if not continue above, it is a useful event to use
    E1.push_back(&(allE1.at(i)));
    EReg1.push_back(&(allEReg1.at(i)));
    Eta1.push_back(&(allEta1.at(i)));
    Phi1.push_back(&(allPhi1.at(i)));
    E2.push_back(&(allE2.at(i)));
    EReg2.push_back(&(allEReg2.at(i)));
    Eta2.push_back(&(allEta2.at(i)));
    Phi2.push_back(&(allPhi2.at(i)));
    nHits1.push_back(&(allnHits1.at(i)));
    HitE1.push_back(&(allHitE1.at(i)));
    HitIX1.push_back(&(allHitIX1.at(i)));
    HitIY1.push_back(&(allHitIY1.at(i)));
    HitIZ1.push_back(&(allHitIZ1.at(i)));
    nHits2.push_back(&(allnHits2.at(i)));
    HitE2.push_back(&(allHitE2.at(i)));
    HitIX2.push_back(&(allHitIX2.at(i)));
    HitIY2.push_back(&(allHitIY2.at(i)));
    HitIZ2.push_back(&(allHitIZ2.at(i)));
  }
  
  // nEvents
  nEvents = (int)E1.size();
  nSignals = nEvents; // assume no background (fix me)
  //
  return nEvents;
}

// this one is store the selection in E1, Eta1, .. E2, Eta2, .. HitE1, ... vectors.
// select events with one hit falling in one particular cell
// and the energy of this hit should be greater than a fraction
// with selection of EBEB, EBEE, or EEEE combination
int SelectEventsInOneCellWithFractionEBorEECombine(int ix, int iy, int iz, double fraction=0.5, int doEvenOdd=0, std::string Combine="")
{
  // clear previous vectors
  nEvents = 0;
  nSignals = 0;
  E1.clear();
  EReg1.clear();
  Eta1.clear();
  Phi1.clear();
  E2.clear();
  EReg2.clear();
  Eta2.clear();
  Phi2.clear();
  nHits1.clear();
  nHits2.clear();
  HitE1.clear();
  HitE2.clear();
  HitIX1.clear();
  HitIY1.clear();
  HitIZ1.clear();
  HitIX2.clear();
  HitIY2.clear();
  HitIZ2.clear();

  // loop over all events and select events
  for (int i=0; i<nEventsAll; i++)
  {
    bool take=false;

    // even odd event check
    if (doEvenOdd==1 && i%2==0) continue;
    else if (doEvenOdd==2 && i%2==1) continue;

    // EB or EE combinatioin check
    if ( Combine=="EBEB" ) // both in EB
    {
      if ( !(fabs(allEta1[i])<1.48&&fabs(allEta2[i])<1.48) ) continue;
    }
    if ( Combine=="EBEE" ) // one in EB one in EE
    {
      if ( !( (fabs(allEta1[i])<1.48&&fabs(allEta2[i])>1.48)||(fabs(allEta2[i])<1.48&&fabs(allEta1[i])>1.48) ) ) continue;
    }
    if ( Combine=="EEEE" ) // both in EE
    {
      if ( !(fabs(allEta1[i])>1.48&&fabs(allEta2[i])>1.48) ) continue;
    }
    if ( Combine=="EE" ) // any one of the two in EE 
    {
      if ( !(fabs(allEta1[i])>1.48||fabs(allEta2[i])>1.48) ) continue;
    }
    if ( Combine=="EB" ) // any one of the two in EB
    {
      if ( !(fabs(allEta1[i])<1.48||fabs(allEta2[i])<1.48) ) continue;
    }

    // loop over all hits of electron1
    for (int ii=0; ii<(int)allHitE1.at(i).size(); ii++ )
    {
      if(take==true) break;

      // if there is even one hit match a cell in selection,
      // also this hit energy is bigger than fraction*E(total),
      // take this event
      if((allHitIX1.at(i).at(ii)==ix) &&
         (allHitIY1.at(i).at(ii)==iy) &&
         (allHitIZ1.at(i).at(ii)==iz) &&
         (allHitE1.at(i).at(ii)>fraction*allE1.at(i)) )
      {

        take = true;
        break;
      }
    }

    // loop over all hits of electron2
    for (int ii=0; ii<(int)allHitE2.at(i).size(); ii++ )
    {
      if(take==true) break;
      // if there is even one hit match a cell in selection,
      // also this hit energy is bigger than fraction*E(total),
      // take this event
      if((allHitIX2.at(i).at(ii)==ix) &&
         (allHitIY2.at(i).at(ii)==iy) &&
         (allHitIZ2.at(i).at(ii)==iz) &&
         (allHitE2.at(i).at(ii)>fraction*allE2.at(i)) )
      {
        take = true;
        break;
      }
    }

    // if do not decide to take this event, continue
    if (!take)
    {
      continue;
    }

    // if not continue above, it is a useful event to use
    E1.push_back(&(allE1.at(i)));
    EReg1.push_back(&(allEReg1.at(i)));
    Eta1.push_back(&(allEta1.at(i)));
    Phi1.push_back(&(allPhi1.at(i)));
    E2.push_back(&(allE2.at(i)));
    EReg2.push_back(&(allEReg2.at(i)));
    Eta2.push_back(&(allEta2.at(i)));
    Phi2.push_back(&(allPhi2.at(i)));
    nHits1.push_back(&(allnHits1.at(i)));
    HitE1.push_back(&(allHitE1.at(i)));
    HitIX1.push_back(&(allHitIX1.at(i)));
    HitIY1.push_back(&(allHitIY1.at(i)));
    HitIZ1.push_back(&(allHitIZ1.at(i)));
    nHits2.push_back(&(allnHits2.at(i)));
    HitE2.push_back(&(allHitE2.at(i)));
    HitIX2.push_back(&(allHitIX2.at(i)));
    HitIY2.push_back(&(allHitIY2.at(i)));
    HitIZ2.push_back(&(allHitIZ2.at(i)));
  }

  // nEvents
  nEvents = (int)E1.size();
  nSignals = nEvents; // assume no background (fix me)
  //
  return nEvents;

}

//
// Store the selection in E1, Eta1, .. E2, Eta2, .. HitE1, ... vectors.
// select events according to the following rule:
// - center at ix,iy,iz among a group of 3x3 = 9 cells
// - if one event has one electron falls inside these 9cells and deposits >90% 
//   of its energy inside these 9 cells, take the event.
// - no need to have all the 9 cells be fired by the same electron in the same 
//   event, but >90% electron energy is required.
// - no need to be exactly 9 cells if the center cell is at (or near) an edge.
// - store those cells among these 9 cells that follow some requirements.
// - if not an electron/event among all the events deposits >10% of the electron 
//   energy into this center cell ix, iy, iz, skip this cell and also all the 
//   other 8 cells around it.
// - if not an electron/event among all the events deposits >2% of its energy 
//   inside the 8 cells around the center cell, skip this surrounding cell.
// 
// Has selection of EBEB, EBEE, or EEEE combination
int SelectEventsInOneCellWith3x3Others(int ix, int iy, int iz, 
                   std::vector<int>& ixx, std::vector<int>& iyy, std::vector<int>& izz, 
                   std::string Combine="")
{
  // clear previous vectors
  nEvents = 0;
  nSignals = 0;
  E1.clear();
  EReg1.clear();
  Eta1.clear();
  Phi1.clear();
  E2.clear();
  EReg2.clear();
  Eta2.clear();
  Phi2.clear();
  nHits1.clear();
  nHits2.clear();
  HitE1.clear();
  HitE2.clear();
  HitIX1.clear();
  HitIY1.clear();
  HitIZ1.clear();
  HitIX2.clear();
  HitIY2.clear();
  HitIZ2.clear();

  // define default 3x3 celles taking ix,iy,iz as center
  // positions
  std::vector<int> iixx, iiyy, iizz;
  // largest energy fraction in all events
  std::vector<double> iimaxefrac;
  // push_back the fitted cell as the first item in the vector
  iixx.push_back(ix);
  iiyy.push_back(iy);
  iizz.push_back(iz);
  iimaxefrac.push_back(0.0);
  // push_back the rest 8 cells
  for (int jx=ix-1; jx<=ix+1; jx++)
  {
    for (int jy=iy-1; jy<=iy+1; jy++)
    {
      // skip the fitting cell that has already been booked.
      if (jx==ix&&jy==iy) continue;
      iiyy.push_back(jy);
      iizz.push_back(iz);
      iimaxefrac.push_back(0.0);
      if (iz==0&&jx==0&&ix==1) 
      {
        iixx.push_back(jx-1);
      }
      else if (iz==0&&jx==0&&ix==-1)
      {
        iixx.push_back(jx+1);
      }
      else
      {
        iixx.push_back(jx);
      } 
    }
  }

  // some debug
  std::cout << " Debug:: print 3x3 cells " << std::endl;
  for (int icell=0; icell<(int)iixx.size(); icell++)
  {
    std::cout << "(" << iixx.at(icell) << "," << iiyy.at(icell) << "," << iizz.at(icell) << ")" << std::endl;
  }

  // loop over all events and select events
  for (int i=0; i<nEventsAll; i++)
  {

    // EB or EE combinatioin check
    if ( Combine=="EBEB" ) // both in EB
    {
      if ( !(fabs(allEta1[i])<1.48&&fabs(allEta2[i])<1.48) ) continue;
    }
    if ( Combine=="EBEE" ) // one in EB one in EE
    {
      if ( !( (fabs(allEta1[i])<1.48&&fabs(allEta2[i])>1.48)||(fabs(allEta2[i])<1.48&&fabs(allEta1[i])>1.48) ) ) continue;
    }
    if ( Combine=="EEEE" ) // both in EE
    {
      if ( !(fabs(allEta1[i])>1.48&&fabs(allEta2[i])>1.48) ) continue;
    }
    if ( Combine=="EE" ) // any one of the two in EE
    {
      if ( !(fabs(allEta1[i])>1.48||fabs(allEta2[i])>1.48) ) continue;
    }
    if ( Combine=="EB" ) // any one of the two in EB
    {
      if ( !(fabs(allEta1[i])<1.48||fabs(allEta2[i])<1.48) ) continue;
    }
    // note, if not any of the string above, this event will always pass this check.

    // 
    double E1_9cells(0);
    // loop over all hits of electron1
    for (int ii=0; ii<(int)allHitE1.at(i).size(); ii++ )
    {
      double hitfrac = allHitE1.at(i).at(ii)/allE1.at(i);
      // loop over all 9 cells
      for (int icell=0; icell<(int)iixx.size(); icell++)
      {
        // check matching
        if( (allHitIX1.at(i).at(ii)==iixx.at(icell)) &&
            (allHitIY1.at(i).at(ii)==iiyy.at(icell)) &&
            (allHitIZ1.at(i).at(ii)==iizz.at(icell)) )
        {
          E1_9cells += allHitE1.at(i).at(ii);
          if (hitfrac>iimaxefrac.at(icell))
          {
            iimaxefrac.at(icell) = hitfrac;
          }
        } // if ( (allHitIX1.at(i)...
      } // for (int icell=0; ..
    } // for (int ii=0; ..
   
    //
    double E2_9cells(0);
    // loop over all hits of electron2
    for (int ii=0; ii<(int)allHitE2.at(i).size(); ii++ )
    {
      double hitfrac = allHitE2.at(i).at(ii)/allE2.at(i);
      // loop over all 9 cells
      for (int icell=0; icell<(int)iixx.size(); icell++)
      {
        // check matching
        if( (allHitIX2.at(i).at(ii)==iixx.at(icell)) &&
            (allHitIY2.at(i).at(ii)==iiyy.at(icell)) &&
            (allHitIZ2.at(i).at(ii)==iizz.at(icell)) )
        {
          E2_9cells += allHitE2.at(i).at(ii);
          if (hitfrac>iimaxefrac.at(icell))
          {
            iimaxefrac.at(icell) = hitfrac;
          }
        } // if ( (allHitIX2.at(i)...
      } // for (int icell=0; ..
    } // for (int ii=0; ..

    // if the 9cells energy in both the two electrons are 
    //  less than 90% of the electrons' energy, I skip this 
    //  event.
    if ( E1_9cells/allE1.at(i)<0.6 && 
         E2_9cells/allE2.at(i)<0.6 )
    {
      continue;
    }

    // if not continue above, it is a useful event to use
    E1.push_back(&(allE1.at(i)));
    EReg1.push_back(&(allEReg1.at(i)));
    Eta1.push_back(&(allEta1.at(i)));
    Phi1.push_back(&(allPhi1.at(i)));
    E2.push_back(&(allE2.at(i)));
    EReg2.push_back(&(allEReg2.at(i)));
    Eta2.push_back(&(allEta2.at(i)));
    Phi2.push_back(&(allPhi2.at(i)));
    nHits1.push_back(&(allnHits1.at(i)));
    HitE1.push_back(&(allHitE1.at(i)));
    HitIX1.push_back(&(allHitIX1.at(i)));
    HitIY1.push_back(&(allHitIY1.at(i)));
    HitIZ1.push_back(&(allHitIZ1.at(i)));
    nHits2.push_back(&(allnHits2.at(i)));
    HitE2.push_back(&(allHitE2.at(i)));
    HitIX2.push_back(&(allHitIX2.at(i)));
    HitIY2.push_back(&(allHitIY2.at(i)));
    HitIZ2.push_back(&(allHitIZ2.at(i)));
  }

  // nEvents
  nEvents = (int)E1.size();
  nSignals = nEvents; // assume no background (fix me)

  // making selection of the 3x3 cells according to the following rules:
  // - if not an electron/event among all the events deposits >10% of the electron
  //   energy into this center cell ix, iy, iz, skip this cell and also all the
  //   other 8 cells around it.
  // - if not an electron/event among all the events deposits >2% of its energy
  //   inside the 8 cells around the center cell, skip this surrounding cell.

  // clean up the output cell position vectors
  ixx.clear();
  iyy.clear();
  izz.clear();

  // check fitting cell at first
  if ( iimaxefrac.at(0)<0.1 ) 
  {
    return 0; // take this as a flag to skip this cell
  }
  else
  {
    // same, the first entry is the fitting cell.
    ixx.push_back(iixx.at(0));
    iyy.push_back(iiyy.at(0));
    izz.push_back(iizz.at(0));
  }

  // check the rest 8 cells;
  for (int icell=1; icell<(int)iixx.size(); icell++)
  {
    // if the maximum cell energy fraction among all events are <2%, do not use this cell
    if (iimaxefrac.at(icell)<0.02) continue;
    // if not use this cell
    ixx.push_back(iixx.at(icell));
    iyy.push_back(iiyy.at(icell));
    izz.push_back(iizz.at(icell));
  }

  //
  return nEvents;

}
//

// select events in one Eta bin
int SelectEventsInOneEtaBin(double bin_min, double bin_max, std::string Combine="")
{
  // clear previous vectors
  nEvents = 0;
  nSignals = 0;
  E1.clear();
  EReg1.clear();
  Eta1.clear();
  Phi1.clear();
  E2.clear();
  EReg2.clear();
  Eta2.clear();
  Phi2.clear();

  // loop over all events and select events
  for (int i=0; i<nEventsAll; i++)
  {

    // EB or EE combinatioin check
    if ( Combine=="EBEB" ) // both in EB
    {
      if ( !(fabs(allEta1[i])<1.48&&fabs(allEta2[i])<1.48) ) continue;
    }
    if ( Combine=="EBEE" ) // one in EB one in EE
    {
      if ( !( (fabs(allEta1[i])<1.48&&fabs(allEta2[i])>1.48)||(fabs(allEta2[i])<1.48&&fabs(allEta1[i])>1.48) ) ) continue;
    }
    if ( Combine=="EEEE" ) // both in EE
    {
      if ( !(fabs(allEta1[i])>1.48&&fabs(allEta2[i])>1.48) ) continue;
    }
    if ( Combine=="EE" ) // any one of the two in EE
    {
      if ( !(fabs(allEta1[i])>1.48||fabs(allEta2[i])>1.48) ) continue;
    }
    if ( Combine=="EB" ) // any one of the two in EB
    {
      if ( !(fabs(allEta1[i])<1.48||fabs(allEta2[i])<1.48) ) continue;
    }
    // note, if not any of the string above, this event will always pass this check.

    // check which electron in this eta-ring
    bool takeEle1(false), takeEle2(false);
    if (allEta1.at(i)>bin_min&&allEta1.at(i)<bin_max) takeEle1 = true;
    else takeEle1 = false;
    if (allEta2.at(i)>bin_min&&allEta2.at(i)<bin_max) takeEle2 = true;
    else takeEle2 = false;
 
    // continue if do not take any one
    if ( (!takeEle1) && (!takeEle2) ) continue;

    // if not continue above, it is a useful event to use
    E1.push_back(&(allE1.at(i)));
    EReg1.push_back(&(allEReg1.at(i)));
    Eta1.push_back(&(allEta1.at(i)));
    Phi1.push_back(&(allPhi1.at(i)));
    UseEle1.push_back(takeEle1);
    E2.push_back(&(allE2.at(i)));
    EReg2.push_back(&(allEReg2.at(i)));
    Eta2.push_back(&(allEta2.at(i)));
    Phi2.push_back(&(allPhi2.at(i)));
    UseEle2.push_back(takeEle2);
  }

  // nEvents
  nEvents = (int)E1.size();
  nSignals = nEvents; // assume no background (fix me)
  
  return nEvents;

}


int AddEtaBinNumberToElectrons(const std::vector<double>& EtaBins, const std::string Combine="")
{
  // clear previous vectors
  nEvents = 0;
  nSignals = 0;
  E1.clear();
  EReg1.clear();
  Eta1.clear();
  Phi1.clear();
  E2.clear();
  EReg2.clear();
  Eta2.clear();
  Phi2.clear();
  // clear energy scale binning number vectors
  ScaleBin1.clear();
  ScaleBin2.clear();


  // loop over all events and select events
  for (int i=0; i<nEventsAll; i++)
  {

    // EB or EE combinatioin check
    if ( Combine=="EBEB" ) // both in EB
    {
      if ( !(fabs(allEta1[i])<1.48&&fabs(allEta2[i])<1.48) ) continue;
    }
    if ( Combine=="EBEE" ) // one in EB one in EE
    {
      if ( !( (fabs(allEta1[i])<1.48&&fabs(allEta2[i])>1.48)||(fabs(allEta2[i])<1.48&&fabs(allEta1[i])>1.48) ) ) continue;
    }
    if ( Combine=="EEEE" ) // both in EE
    {
      if ( !(fabs(allEta1[i])>1.48&&fabs(allEta2[i])>1.48) ) continue;
    }
    if ( Combine=="EE" ) // any one of the two in EE
    {
      if ( !(fabs(allEta1[i])>1.48||fabs(allEta2[i])>1.48) ) continue;
    }
    if ( Combine=="EB" ) // any one of the two in EB
    {
      if ( !(fabs(allEta1[i])<1.48||fabs(allEta2[i])<1.48) ) continue;
    }
    // note, if not any of the string above, this event will always pass this check.

    // findout which bins it belongs to
    int bin1(-1), bin2(-1); 
    for (int ibin=0; ibin<(int)EtaBins.size()-1; ibin++)
    {
      double bin_min = EtaBins.at(ibin);
      double bin_max = EtaBins.at(ibin+1);
 
      if (allEta1.at(i)>bin_min&&allEta1.at(i)<bin_max) bin1 = ibin;
      if (allEta2.at(i)>bin_min&&allEta2.at(i)<bin_max) bin2 = ibin;
    }

    // check if both bin1 and bin2 have NOT been given a value, we simply skip this event
    if (bin1<0&&bin2<0) continue;

    // give bin number, -1 means no e-scale is going to be applied.
    ScaleBin1.push_back(bin1);
    ScaleBin2.push_back(bin2);

    // if not continue above, it is a useful event to use
    E1.push_back(&(allE1.at(i)));
    EReg1.push_back(&(allEReg1.at(i)));
    Eta1.push_back(&(allEta1.at(i)));
    Phi1.push_back(&(allPhi1.at(i)));
    E2.push_back(&(allE2.at(i)));
    EReg2.push_back(&(allEReg2.at(i)));
    Eta2.push_back(&(allEta2.at(i)));
    Phi2.push_back(&(allPhi2.at(i)));
  }

  // nEvents
  nEvents = (int)E1.size();
  nSignals = nEvents; // assume no background (fix me)

  return nEvents;

}

int AddEtaBinNumberToElectrons(const std::vector<EnergyScale>& EtaScale, const std::string Combine="")
{
  // clear previous vectors
  nEvents = 0;
  nSignals = 0;
  E1.clear();
  EReg1.clear();
  Eta1.clear();
  Phi1.clear();
  E2.clear();
  EReg2.clear();
  Eta2.clear();
  Phi2.clear();
  // clear energy scale binning number vectors
  ScaleBin1.clear();
  ScaleBin2.clear();


  // loop over all events and select events
  for (int i=0; i<nEventsAll; i++)
  {

    // EB or EE combinatioin check
    if ( Combine=="EBEB" ) // both in EB
    {
      if ( !(fabs(allEta1[i])<1.48&&fabs(allEta2[i])<1.48) ) continue;
    }
    if ( Combine=="EBEE" ) // one in EB one in EE
    {
      if ( !( (fabs(allEta1[i])<1.48&&fabs(allEta2[i])>1.48)||(fabs(allEta2[i])<1.48&&fabs(allEta1[i])>1.48) ) ) continue;
    }
    if ( Combine=="EEEE" ) // both in EE
    {
      if ( !(fabs(allEta1[i])>1.48&&fabs(allEta2[i])>1.48) ) continue;
    }
    if ( Combine=="EE" ) // any one of the two in EE
    {
      if ( !(fabs(allEta1[i])>1.48||fabs(allEta2[i])>1.48) ) continue;
    }
    if ( Combine=="EB" ) // any one of the two in EB
    {
      if ( !(fabs(allEta1[i])<1.48||fabs(allEta2[i])<1.48) ) continue;
    }
    // note, if not any of the string above, this event will always pass this check.

    // findout which bins it belongs to
    int bin1(-1), bin2(-1);
    for (int ibin=0; ibin<(int)EtaScale.size(); ibin++)
    {
      if (allEta1.at(i)>EtaScale.at(ibin).min&&allEta1.at(i)<EtaScale.at(ibin).max) bin1 = ibin;
      if (allEta2.at(i)>EtaScale.at(ibin).min&&allEta2.at(i)<EtaScale.at(ibin).max) bin2 = ibin;
    }


    // check if both bin1 and bin2 have NOT been given a value, we simply skip this event
    if (bin1<0&&bin2<0) continue;

    // give bin number, -1 means no e-scale is going to be applied.
    ScaleBin1.push_back(bin1);
    ScaleBin2.push_back(bin2);

    // if not continue above, it is a useful event to use
    E1.push_back(&(allE1.at(i)));
    EReg1.push_back(&(allEReg1.at(i)));
    Eta1.push_back(&(allEta1.at(i)));
    Phi1.push_back(&(allPhi1.at(i)));
    E2.push_back(&(allE2.at(i)));
    EReg2.push_back(&(allEReg2.at(i)));
    Eta2.push_back(&(allEta2.at(i)));
    Phi2.push_back(&(allPhi2.at(i)));
  }

  // nEvents
  nEvents = (int)E1.size();
  nSignals = nEvents; // assume no background (fix me)

  return nEvents;

}


void ApplyEtaScaleToAllEvents(const std::vector<double>& EtaBins,const std::vector<double>& EtaScales)
{
  // loop over all events 
  for (int i=0; i<nEventsAll; i++)
  {
    // findout which bins it belongs to
    int bin1(-1), bin2(-1);
    for (int ibin=0; ibin<(int)EtaBins.size()-1; ibin++)
    {
      double bin_min = EtaBins.at(ibin);
      double bin_max = EtaBins.at(ibin+1);

      if (allEta1.at(i)>bin_min&&allEta1.at(i)<bin_max) bin1 = ibin;
      if (allEta2.at(i)>bin_min&&allEta2.at(i)<bin_max) bin2 = ibin;
    }

    // apply eta scale
    if (bin1>=0) allEReg1.at(i) *= EtaScales.at(bin1);
    if (bin2>=0) allEReg2.at(i) *= EtaScales.at(bin2);

  }

}

void ApplyEtaScaleToAllEvents(const std::vector<EnergyScale>& EtaScale)
{
  // loop over all events
  for (int i=0; i<nEventsAll; i++)
  {
    // findout which bins it belongs to
    int bin1(-1), bin2(-1);
    for (int ibin=0; ibin<(int)EtaScale.size(); ibin++)
    {
      if (allEta1.at(i)>EtaScale.at(ibin).min&&allEta1.at(i)<EtaScale.at(ibin).max) bin1 = ibin;
      if (allEta2.at(i)>EtaScale.at(ibin).min&&allEta2.at(i)<EtaScale.at(ibin).max) bin2 = ibin;
    }

    // apply eta scale
    if (bin1>=0) allEReg1.at(i) *= EtaScale.at(bin1).s;
    if (bin2>=0) allEReg2.at(i) *= EtaScale.at(bin2).s;

  }

}

//////////


// select groups of events according to a vector of cells,
// such as a supercluster.
int SelectSubsetEvents(std::vector<std::vector<int> > cells)
{
  // cells to be selected
  if ((int)cells.size()==0)
  {
    return 0;
  }
  // clear previous vectors
  nEvents = 0;
  nSignals = 0;
  E1.clear();
  EReg1.clear();
  Eta1.clear();
  Phi1.clear();
  E2.clear();
  EReg2.clear();
  Eta2.clear();
  Phi2.clear();
  nHits1.clear();
  nHits2.clear();
  HitE1.clear();
  HitE2.clear();
  HitIX1.clear();
  HitIY1.clear();
  HitIZ1.clear();
  HitIX2.clear();
  HitIY2.clear();
  HitIZ2.clear();
  
  // loop over all events and select events
  for (int i=0; i<nEventsAll; i++)
  {
    bool take=false;
    // loop over all hits of electron1
    for (int ii=0; ii<(int)allHitE1.at(i).size(); ii++ )
    {
      // if one hit match, then no need to loop over the rest hits
      if (take) break;
      // loop over all cells to be selected
      for (int icell=0; icell<(int)cells.size(); icell++)
      {
        if((allHitIX1.at(i).at(ii)==cells.at(icell).at(0)) &&
           (allHitIY1.at(i).at(ii)==cells.at(icell).at(1)) && 
           (allHitIZ1.at(i).at(ii)==cells.at(icell).at(2)))
        {
          // if there is even one hit match a cell in selection, take the event
          take = true;
          break;
        }
      }
    }
    // loop over all hits of electron2
    for (int ii=0; ii<(int)allHitE2.at(i).size(); ii++ )
    {
      // if one hit match, then no need to loop over the rest hits
      if (take) break;
      // loop over all cells to be selected
      for (int icell=0; icell<(int)cells.size(); icell++)
      {
        if((allHitIX2.at(i).at(ii)==cells.at(icell).at(0)) &&
           (allHitIY2.at(i).at(ii)==cells.at(icell).at(1)) && 
           (allHitIZ2.at(i).at(ii)==cells.at(icell).at(2)))
        {
          // if there is even one hit match a cell in selection, take the event
          take = true;
          break;
        }
      }
    }
    
    // if do not decide to take this event, continue
    if (!take) continue;
    
    // if not continue above, it is a useful event ot use
    E1.push_back(&(allE1.at(i)));
    EReg1.push_back(&(allEReg1.at(i)));
    Eta1.push_back(&(allEta1.at(i)));
    Phi1.push_back(&(allPhi1.at(i)));
    E2.push_back(&(allE2.at(i)));
    EReg2.push_back(&(allEReg2.at(i)));
    Eta2.push_back(&(allEta2.at(i)));
    Phi2.push_back(&(allPhi2.at(i)));
    nHits1.push_back(&(allnHits1.at(i)));
    HitE1.push_back(&(allHitE1.at(i)));
    HitIX1.push_back(&(allHitIX1.at(i)));
    HitIY1.push_back(&(allHitIY1.at(i)));
    HitIZ1.push_back(&(allHitIZ1.at(i)));
    nHits2.push_back(&(allnHits2.at(i)));
    HitE2.push_back(&(allHitE2.at(i)));
    HitIX2.push_back(&(allHitIX2.at(i)));
    HitIY2.push_back(&(allHitIY2.at(i)));
    HitIZ2.push_back(&(allHitIZ2.at(i)));
  }
  
  // nEvents
  nEvents = (int)E1.size();
  nSignals = nEvents; // assume no background (fix me)
  //
  return nEvents;
}


// select groups of events according to a vector of cells,
// such as a supercluster,
// with a choice of EBEB, EBEE, or EEEE combination
int SelectSubsetEventsWithEBOrEECombine(std::vector<std::vector<int> > cells, const std::string Combine="")
{
  // cells to be selected
  if ((int)cells.size()==0)
  {
    return 0;
  } 
  // clear previous vectors
  nEvents = 0;
  nSignals = 0;
  E1.clear();
  EReg1.clear();
  Eta1.clear();
  Phi1.clear();
  E2.clear();
  EReg2.clear();
  Eta2.clear();
  Phi2.clear();
  nHits1.clear();
  nHits2.clear();
  HitE1.clear();
  HitE2.clear();
  HitIX1.clear();
  HitIY1.clear();
  HitIZ1.clear();
  HitIX2.clear();
  HitIY2.clear();
  HitIZ2.clear();

  // loop over all events and select events
  for (int i=0; i<nEventsAll; i++)
  {
    // EB or EE combinatioin check
    if ( Combine=="EBEB" ) // both in EB
    {
      if ( !(fabs(allEta1[i])<1.48&&fabs(allEta2[i])<1.48) ) continue;
    }
    if ( Combine=="EBEE" ) // one in EB one in EE
    {
      if ( !( (fabs(allEta1[i])<1.48&&fabs(allEta2[i])>1.48)||(fabs(allEta2[i])<1.48&&fabs(allEta1[i])>1.48) ) ) continue;
    }
    if ( Combine=="EEEE" ) // both in EE
    {
      if ( !(fabs(allEta1[i])>1.48&&fabs(allEta2[i])>1.48) ) continue;
    }
    if ( Combine=="EE" ) // any one of the two in EE
    {
      if ( !(fabs(allEta1[i])>1.48||fabs(allEta2[i])>1.48) ) continue;
    }
    if ( Combine=="EB" ) // any one of the two in EB
    {
      if ( !(fabs(allEta1[i])<1.48||fabs(allEta2[i])<1.48) ) continue;
    }


    bool take=false;
    // loop over all hits of electron1
    for (int ii=0; ii<(int)allHitE1.at(i).size(); ii++ )
    {
      // if one hit match, then no need to loop over the rest hits
      if (take) break;
      // loop over all cells to be selected
      for (int icell=0; icell<(int)cells.size(); icell++)
      {
        if((allHitIX1.at(i).at(ii)==cells.at(icell).at(0)) &&
           (allHitIY1.at(i).at(ii)==cells.at(icell).at(1)) &&
           (allHitIZ1.at(i).at(ii)==cells.at(icell).at(2)))
        {
          // if there is even one hit match a cell in selection, take the event
          take = true;
          break;
        }
      }
    }
    // loop over all hits of electron2
    for (int ii=0; ii<(int)allHitE2.at(i).size(); ii++ )
    {
      // if one hit match, then no need to loop over the rest hits
      if (take) break;
      // loop over all cells to be selected
      for (int icell=0; icell<(int)cells.size(); icell++)
      {
        if((allHitIX2.at(i).at(ii)==cells.at(icell).at(0)) &&
           (allHitIY2.at(i).at(ii)==cells.at(icell).at(1)) &&
           (allHitIZ2.at(i).at(ii)==cells.at(icell).at(2)))
        {
          // if there is even one hit match a cell in selection, take the event
          take = true;
          break;
        }
      }
    }

    // if do not decide to take this event, continue
    if (!take) continue;

    // if not continue above, it is a useful event ot use
    E1.push_back(&(allE1.at(i)));
    EReg1.push_back(&(allEReg1.at(i)));
    Eta1.push_back(&(allEta1.at(i)));
    Phi1.push_back(&(allPhi1.at(i)));
    E2.push_back(&(allE2.at(i)));
    EReg2.push_back(&(allEReg2.at(i)));
    Eta2.push_back(&(allEta2.at(i)));
    Phi2.push_back(&(allPhi2.at(i)));
    nHits1.push_back(&(allnHits1.at(i)));
    HitE1.push_back(&(allHitE1.at(i)));
    HitIX1.push_back(&(allHitIX1.at(i)));
    HitIY1.push_back(&(allHitIY1.at(i)));
    HitIZ1.push_back(&(allHitIZ1.at(i)));
    nHits2.push_back(&(allnHits2.at(i)));
    HitE2.push_back(&(allHitE2.at(i)));
    HitIX2.push_back(&(allHitIX2.at(i)));
    HitIY2.push_back(&(allHitIY2.at(i)));
    HitIZ2.push_back(&(allHitIZ2.at(i)));

  }

  // nEvents
  nEvents = (int)E1.size();
  nSignals = nEvents; // assume no background (fix me)
  //
  return nEvents;

}


// cells in 1.83<eta<2.17
// a band between 25 to 35 crystals apart from the center.
std::vector<std::vector<int> > GetCellsEta1p83To2p17()
{
  std::vector<std::vector<int> > cells;
  
  // loop over all cells and select those in the region
  for (int iz=-1; iz<=1; iz=iz+2 )
  {
    for (int ix=1; ix<=100; ix++)
    {
      for (int iy=1; iy<=100; iy++)
      {
        // R between 25 to 35 crystals
        // thus R=sqrt(xx*xx+yy*yy)
        // xx = (ix-0.5)-50.0
        // yy = (iy-0.5)-50.0
        double xx = (double(ix)-0.5)-50.0;
        double yy = (double(iy)-0.5)-50.0;
        double R = sqrt(xx*xx+yy*yy);
        if (R<=35&&R>=25)
        {
          std::vector<int> acell;
          acell.push_back(ix);
          acell.push_back(iy);
          acell.push_back(iz);
          cells.push_back(acell);
        }
      }
    }
  }
  return cells;
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

std::vector<std::vector<int> > GetAllCellsFromCalibTable(const char* filename)
{
  std::vector<std::vector<int> > cells;
  std::ifstream myfile(filename);
  std::string line;
  int idx, ix, iy, iz, fixed, nfits;
  double c, cerr;
  if (myfile.is_open())
  {
    while (getline(myfile,line))
    {
      calibRecord calib;
      std::stringstream sline(line);
      sline >> idx
             >> ix
             >> iy
             >> iz
             >> c
             >> cerr
             >> fixed
             >> nfits ;
      std::vector<int> acell;
      acell.push_back(ix);
      acell.push_back(iy);
      acell.push_back(iz);
      cells.push_back(acell);
    }
    myfile.close();
  }
  
  //return
  return cells;
  
}
