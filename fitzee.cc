// Authors: Hengne Li of UVa at CERN in 2013   
 
#include "functions.h"
#include "variables.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnParameterScan.h"
#include "Minuit2/MnContours.h"
#include "TFile.h"
#include "TH1D.h"
#include "TChain.h"
#include "TGraph.h"
#include <stdio.h>
#include <iostream>


using namespace ROOT::Minuit2;

int mode = 0;
// 0: first fit, single cell, one parameter
// 1: second fit, single cell, one free parameter, but all parameters to be used.

// Mee calculation method
int method = 2;
// 0: only recaculate E of the second electron, using a calibConst in a reference calibTable
//       for cells not to be fitted.
// 1: recalculate the E of both the two electrons
// 2: recalculate the E of both the two electrons, using a calibConst in a reference calibTable
//    for cells not to be fitted.

// signal fraction
double signalFraction = 0.82;

// Z mass
double _Zmass = 91.5586; //91.513; // 91.188;

// fit Z mass window High
double _FitWindowHigh = 100.0;
double _FitWindowLow = 80.0;

// Gassian Resolution of the Mee distribution
double gaus_reso = 2.0;

// energy fraction of one crystal in fit
double hitEnergyFraction = 0.3;

int debug = 0;
char rootfile_in[300];
char rootfile_out[300];
char calibtable_ref[300];
char calibtable_in[300];
char calibtable_out[300];
char etascale_ref[300];

// even or odd events
int doEvenOdd = 0; // 0 for not to split, 1 for odd, 2 for even 

// combine, can be EBEB, EBEE, EEEE, EE (any one in EE), EB (any one in EB)
std::string _combine = "EE";

int main(int argc, char* argv[])
{
  if (argc<8)
  {
    std::cout << argv[0] << " <mode> <input_file.root> <output_file.root> \\\n"
              << " <input_calibTable.dat> <out_calibTable.dat> <ref_calibTable.dat>\\\n" 
              << " <signalFraction> <method> <GaussResolution> <hitEnergyFraction> <debug> \\\n"
              << " <doEvenOdd> <etascale_reference_file.dat>"
    << std::endl;
    return 0;
  }
  mode = atoi(argv[1]);
  sprintf(rootfile_in, "%s", argv[2]);
  sprintf(rootfile_out, "%s", argv[3]);
  sprintf(calibtable_in, "%s", argv[4]);
  sprintf(calibtable_out, "%s", argv[5]);
  sprintf(calibtable_ref, "%s", argv[6]);
  
  if (argc>7)
  {
    signalFraction = atof(argv[7]);
  }
  
  if (argc>8)
  {
    method = atoi(argv[8]);
  }
  
  if (argc>9)
  {
    gaus_reso = atof(argv[9]);
  }

  if (argc>10)
  {
    hitEnergyFraction = atof(argv[10]);
  }

  if (argc>11)
  {
    debug = atoi(argv[11]);
  }
  
  if (argc>12)
  {
    doEvenOdd = atoi(argv[12]);
  }

  if (argc>13)
  {
    sprintf(etascale_ref, "%s", argv[13]);
  } 

  // check
  if (mode==31||mode==32)
  {
    if (argc<=13)
    {
      std::cout << "Missing etascale reference file. Please run " << argv[0] << " to print usage information. " << std::endl;
      return 1;
    }
  }
 
  std::cout << "fitzee.exe " << std::endl;
  std::cout << "  mode = " << mode << std::endl;
  std::cout << "  method = " << method << std::endl;
  std::cout << "  rootfile_in = " << rootfile_in << std::endl;
  std::cout << "  rootfile_out = " << rootfile_out << std::endl;
  std::cout << "  calibtable_in = " << calibtable_in << std::endl;
  std::cout << "  calibtable_out = " << calibtable_out << std::endl;
  std::cout << "  calibtable_ref = " << calibtable_ref << std::endl;
  std::cout << "  signalFraction = " << signalFraction << std::endl;
  std::cout << "  GaussianResolution = " << gaus_reso << std::endl;
  std::cout << "  hitEnergyFraction = " << hitEnergyFraction << std::endl;
  std::cout << "  doEvenOdd = " << doEvenOdd << std::endl;
  if (mode==31||mode==32) std::cout << "  EtaScale Reference File = " << etascale_ref << std::endl; 
  std::cout << " Start program. " << std::endl;
  
  // reading data
  TChain* tree = new TChain("tree", "tree");
  tree->Add(rootfile_in);
  
  // output root file
  TFile* fout = new TFile(rootfile_out, "recreate");
  
  // Set the branches for the TChain/TTree
  SetTreeBranch(tree);

  // Fill all events into vectors
  FillAllEvents(tree, 0);
  
  // delete the chain no more need it
  tree->Delete();
  
  // test fit
  
  // get all cells to be fit
  std::vector<std::vector<int> > cells = GetAllCellsFromCalibTable(calibtable_in);

  
  // select events for
  //nEvents = SelectSubsetEvents(cells);
  nEvents = SelectSubsetEventsWithEBOrEECombine(cells, _combine);
  
  std::cout << "nEvents = " << nEvents << std::endl;
 
  // reference scales
  std::vector<double> scalesref;
  if (mode==31)
  {
    // read reference scales
    std::ifstream reffile(etascale_ref);
    if (reffile.is_open())
    {
      std::string line;
      double s, e;
      while (getline(reffile,line))
      {
        std::stringstream sline(line);
        sline >> s >> e;
        scalesref.push_back(s);
      }
      reffile.close();
    }
  }

  // reference scale for mode 32
  std::vector<EnergyScale> EtaScaleRef;
  if (mode==32)
  {
    // read reference scales
    std::ifstream reffile(etascale_ref);
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
        EtaScaleRef.push_back(scale);
      }
      reffile.close();
    }
  }
 
  if (!nEvents)
  {
    std::cout << "No event selected, quit the job." << std::endl;
    return 0;
  }
  
  if (debug>0) std::cout << " Step 1: Initialize function overall " << std::endl;
  
  // define the fitting function
  BWGSLikelihoodFCN fcn(nEvents, nSignals,
                        E1, EReg1, Eta1, Phi1, nHits1, HitE1, HitIX1, HitIY1, HitIZ1,
                        E2, EReg2, Eta2, Phi2, nHits2, HitE2, HitIX2, HitIY2, HitIZ2,
                        debug, method);

  if (debug>0) std::cout << " Step 1: Initialize PDF overall " << std::endl;
  // initialize PDF
  fcn.initBWGSParameters(_FitWindowHigh, // windowHigh
                         _FitWindowLow, // windowLow
                         _Zmass, //voigtMass
                         gaus_reso, //voightResolution
                         2.4952, //voigtWidth
                         nSignals,
                         nEvents);
  
  if (debug>0) std::cout << " Step 1: Initialize Data overall " << std::endl;
  // initialize calibTable
  fcn.initCalibTable(calibtable_in);
  
  if (debug>0) std::cout << " Step 1: Initialize CalibTableRef overall " << std::endl;
  // initialize calibTableRef
  fcn.initCalibTableRef(calibtable_ref);
  
  // define the parameters
  MnUserParameters pars;
  
  //if (debug>0) std::cout << " Step 1: Initialize MnUserParameters overall " << std::endl;
  // initialize the parameters from the calibTable in BWGSLikelihoodFCN
  fcn.initMnUserParametersFromCurrentCalibTable(pars);
  
  //if (debug>0) std::cout << " Step 1: Store original overall Mee hist" << std::endl;
  // store hist using original calibration
  //TH1D* h0 = (TH1D*)fcn.GetMeeTH1DFromOriginalCalibration("hMeeOriginal");
  //fout->cd();
  //h0->Write();
  
  //if (debug>0) std::cout << " Step 1: Store overall Mee hist before fit" << std::endl;
  // store hist before fit
  //TH1D* h1 = (TH1D*)fcn.GetMeeTH1DFromMnUserParameters(pars, "hMeeBefore");
  //fout->cd();
  //h1->Write();

  
  // now start the real fit for mode==0
  if (debug>0) std::cout << " Step 2: Start fit cell-by-cell" << std::endl;
  // calibRecord vector to store the calib records
  std::vector<calibRecord> calibTable;
  
  if (mode==0)
  {
    // loop over all the cells,
    // select only those events in this cell,
    // define only this cell one single parameter to be fitted
    for (int ic=0; ic<(int)cells.size(); ic++)
    {
      // cell ix, iy, iz
      const int ix = cells.at(ic).at(0);
      const int iy = cells.at(ic).at(1);
      const int iz = cells.at(ic).at(2);
      
      // select events for only this cell
      //nEvents = SelectEventsInOneCell(ix, iy, iz);
      //nEvents = SelectEventsInOneCellWithFraction(ix, iy, iz, hitEnergyFraction);
      nEvents = SelectEventsInOneCellWithFractionEBorEECombine(ix, iy, iz, hitEnergyFraction, doEvenOdd, _combine);
      // define the signal fraction
      nSignals = int(signalFraction*(double)nEvents);
      
      std::cout << "ic = " << ic << ", " << "nEvents = " << nEvents << std::endl;
      
      // initialize fcn using this set of events
      fcn.initData(nEvents, nSignals,
                   E1, EReg1, Eta1, Phi1, nHits1, HitE1, HitIX1, HitIY1, HitIZ1,
                   E2, EReg2, Eta2, Phi2, nHits2, HitE2, HitIX2, HitIY2, HitIZ2,
                   debug, method);
      // initialize PDF
      fcn.initBWGSParameters(_FitWindowHigh, // windowHigh
                             _FitWindowLow, // windowLow
                             _Zmass, //voigtMass
                             gaus_reso, //voightResolution
                             2.4952, //voigtWidth
                             nSignals,
                             nEvents);
      // set a calibTable with only one cell
      calibRecord AcalibRecord = {ic, ix, iy, iz, 1.0, 0.1, 0, 0};
      std::vector<calibRecord> AcalibTable;
      AcalibTable.push_back(AcalibRecord);
      
      // initialize calibTable in fcn using the only-one-cell calibTable above
      fcn.initCalibTableFromVector(AcalibTable);
      
      // also initialize calibTableRef in fcn using the only-one-cell calibTable above
      fcn.initCalibTableRefFromVector(AcalibTable);
      
      // also define a one-par MnUserParameters
      MnUserParameters Apars;
      char name[100];
      sprintf(name, "par_ix%d_iy%d_iz%d", ix, iy, iz);
      Apars.Add(name, 1.0, 0.1);
      // define scan
      //MnScan scan(fcn, Apars);
      
      //FunctionMinimum min1 = scan(0);
      
      // define migrad
      MnMigrad migrad(fcn, Apars);
      
      // minimize
      FunctionMinimum min = migrad();
      
      // print min
      std::cout << "minimum of cell (ix, iy, iz) = (" 
              << ix << "," << iy << "," << iz << "): \n"
              << min << std::endl;
      
      // update MnUserParameters
      Apars = migrad.Parameters();
      //Apars = scan.Parameters();
      
      // update calibTable, only has one parameter.
      AcalibTable.at(0).c = Apars.Value(0);
      AcalibTable.at(0).cerr = Apars.Error(0);
      
      // push back the fitted result, again, only has one parameter
      calibTable.push_back(AcalibTable.at(0));
    }
  }
  else if (mode==1)
  {
    if (debug>0) std::cout << " Step 2: Fit cell-by-cell, inside mode 1" << std::endl;
    // but initialize calibTableRef using the full table, for method 2
    fcn.initCalibTableRef(calibtable_ref);
    
    // loop over all the cells,
    // select only those events in this cell,
    // define only this cell one single parameter to be fitted
    for (int ic=0; ic<(int)cells.size(); ic++)
    {
      // cell ix, iy, iz
      const int ix = cells.at(ic).at(0);
      const int iy = cells.at(ic).at(1);
      const int iz = cells.at(ic).at(2);

      // select events for only this cell
      //nEvents = SelectEventsInOneCell(ix, iy, iz);
      //nEvents = SelectEventsInOneCellWithFraction(ix, iy, iz, hitEnergyFraction, doEvenOdd);
      nEvents = SelectEventsInOneCellWithFractionEBorEECombine(ix, iy, iz, hitEnergyFraction, doEvenOdd, _combine);
      // define the signal fraction
      nSignals = int(signalFraction*(double)nEvents);
      
      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, for Cell: " << std::endl;
      std::cout << "ic = " << ic << ", " << "nEvents = " << nEvents << std::endl;

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, initialize data " << std::endl;
      // initialize fcn using this set of events
      fcn.initData(nEvents, nSignals,
                   E1, EReg1, Eta1, Phi1, nHits1, HitE1, HitIX1, HitIY1, HitIZ1,
                   E2, EReg2, Eta2, Phi2, nHits2, HitE2, HitIX2, HitIY2, HitIZ2,
                   debug, method);

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, initialize PDF " << std::endl;
      // initialize PDF
      fcn.initBWGSParameters(_FitWindowHigh, // windowHigh
                             _FitWindowLow, // windowLow
                             _Zmass, //voigtMass
                             gaus_reso, //voightResolution
                             2.4952, //voigtWidth
                             nSignals,
                             nEvents);
      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, set one-cell calibTable" << std::endl;
      // set a calibTable with only one cell
      //double calibC_old = fcn.GetCalibConstFromCalibTableRef(ix, iy, iz);
      //double calibCerr_old = fcn.GetCalibConstErrorFromCalibTableRef(ix, iy, iz);
      //calibRecord AcalibRecord = {ic, ix, iy, iz, calibC_old, calibCerr_old, 0, 0};
      calibRecord AcalibRecord = {ic, ix, iy, iz, 1.0, 0.01, 0, 0};
      std::vector<calibRecord> AcalibTable;
      AcalibTable.push_back(AcalibRecord);
      
      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, initialize one-cell calibTable" << std::endl;
      // initialize calibTable in fcn using the only-one-cell calibTable above
      fcn.initCalibTableFromVector(AcalibTable);
      
      // also define a one-par MnUserParameters
      MnUserParameters Apars;
      char name[100];
      sprintf(name, "par_ix%d_iy%d_iz%d", ix, iy, iz);
      //Apars.Add(name, calibC_old, calibCerr_old);
      Apars.Add(name, 1.0, 0.01);
      
      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, define MnMigrad" << std::endl;
      // define migrad
      MnMigrad migrad(fcn, Apars);
      
      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, do fit" << std::endl;
      // minimize
      FunctionMinimum min = migrad();
      
      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, print output" << std::endl;
      // print min
      std::cout << "minimum of cell (ix, iy, iz) = (" << ix << "," << iy << "," << iz << "): \n"
      << min << std::endl;
      
      // update MnUserParameters
      Apars = migrad.Parameters();
      
      // print scan parameter
      // plot
      //MnPlot plot;
      // scan parameters
      //MnParameterScan parscan(fcn, Apars);
      //std::vector< std::pair<double, double> > points_scan = parscan(0);
      //if (min.IsValid()) plot(points_scan);
      
      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, update calibTable and store output" << std::endl;
      // update calibTable, only has one parameter.
      AcalibTable.at(0).c = Apars.Value(0);
      AcalibTable.at(0).cerr = Apars.Error(0);
      
      // push back the fitted result, again, only has one parameter
      calibTable.push_back(AcalibTable.at(0));
      
    }
  }
  else if (mode==2)
  {
    // mode 2 is similar to mode 1, but assuming model Enew = Eold * (Ci + Fi * <dC>),
    // where Fi = (Etotal-Eold) / Eold , and <dC> is the mean delta CalibC bias averaging
    //   all crystals.
    // still fit one calibC one time, but par[0] = <dC> and only par[1] is the CabliC to be fitted.
    if (debug>0) std::cout << " Step 2: Fit cell-by-cell, inside mode 2" << std::endl;
    // but initialize calibTableRef using the full table, for method 2
    fcn.initCalibTableRef(calibtable_ref);

    // loop over all the cells,
    // select only those events in this cell,
    // define only this cell one single parameter to be fitted
    for (int ic=0; ic<(int)cells.size(); ic++)
    {
      // cell ix, iy, iz
      const int ix = cells.at(ic).at(0);
      const int iy = cells.at(ic).at(1);
      const int iz = cells.at(ic).at(2);

      // select events for only this cell
      //nEvents = SelectEventsInOneCell(ix, iy, iz);
      //nEvents = SelectEventsInOneCellWithFraction(ix, iy, iz, hitEnergyFraction, doEvenOdd);
      nEvents = SelectEventsInOneCellWithFractionEBorEECombine(ix, iy, iz, hitEnergyFraction, doEvenOdd, _combine);
      // define the signal fraction
      nSignals = int(signalFraction*(double)nEvents);

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, for Cell: " << std::endl;
      std::cout << "ic = " << ic << ", " << "nEvents = " << nEvents << std::endl;

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, initialize data " << std::endl;
      // initialize fcn using this set of events
      fcn.initData(nEvents, nSignals,
                   E1, EReg1, Eta1, Phi1, nHits1, HitE1, HitIX1, HitIY1, HitIZ1,
                   E2, EReg2, Eta2, Phi2, nHits2, HitE2, HitIX2, HitIY2, HitIZ2,
                   debug, method);

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, initialize PDF " << std::endl;
      // initialize PDF
      fcn.initBWGSParameters(_FitWindowHigh, // windowHigh
                             _FitWindowLow, // windowLow
                             _Zmass, //voigtMass
                             gaus_reso, //voightResolution
                             2.4952, //voigtWidth
                             nSignals,
                             nEvents);
      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, set one-cell calibTable" << std::endl;
      // set a calibTable with only one cell
      calibRecord AcalibRecord = {ic, ix, iy, iz, 1.0, 0.01, 0, 0};
      std::vector<calibRecord> AcalibTable;
      AcalibTable.push_back(AcalibRecord);

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, initialize one-cell calibTable" << std::endl;
      // initialize calibTable in fcn using the only-one-cell calibTable above
      fcn.initCalibTableFromVector(AcalibTable);

      // also define a one-par MnUserParameters
      MnUserParameters Apars;
      char name[100];
      sprintf(name, "mdC_ix%d_iy%d_iz%d", ix, iy, iz);
      Apars.Add(name, 0.0, 0.001); // mean delta C:  <dC>
      sprintf(name, "par_ix%d_iy%d_iz%d", ix, iy, iz);
      Apars.Add(name, 1.0, 0.01); //  calibC


      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, define MnMigrad" << std::endl;
      // define migrad
      MnMigrad migrad(fcn, Apars);

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, do fit" << std::endl;
      // minimize
      FunctionMinimum min = migrad();

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, print output" << std::endl;
      // print min
      std::cout << "minimum of cell (ix, iy, iz) = (" << ix << "," << iy << "," << iz << "): \n"
      << min << std::endl;

      // fit a second time with the bias <dC> fixed
      migrad.Fix((unsigned int)0);
      min = migrad();

      // hesse
      MnHesse hesse;
      hesse(fcn, min);

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, 2nd time, print output" << std::endl;
      // print min
      std::cout << "minimum of cell (ix, iy, iz) = (" << ix << "," << iy << "," << iz << "): \n"
      << min << std::endl;

      // update MnUserParameters
      //Apars = migrad.Parameters();
      Apars = min.UserParameters();

      // print scan parameter
      //if (min.IsValid())
      //{
      //  // plot
      //  MnPlot plot;
      //  // scan parameters
      //  MnParameterScan parscan(fcn, Apars);
      //  std::vector< std::pair<double, double> > points_scan = parscan(0);
      //  plot(points_scan);
      //}

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, update calibTable and store output" << std::endl;
      // update calibTable, par[0] is <dC> , par[1] is Ci.
      AcalibTable.at(0).c = Apars.Value(1);
      AcalibTable.at(0).cerr = Apars.Error(1);

      // push back the fitted result, again, only has one parameter
      calibTable.push_back(AcalibTable.at(0));

    }
  }
  else if (mode==3 || mode==31 ||  mode==32)
  {
    // mode 3 is similar to mode 1, but trying to fit several crystalls at once,
    //  automatically find 3x3 cells around it
    // mode 31 is to apply pre-defined eta-scales before fits
    // mode 32 is the same as mode 31 just use a different eta-scale format

    if (debug>0) std::cout << " Step 2: Fit cell-by-cell, inside mode 3 or 31" << std::endl;
    // but initialize calibTableRef using the full table, for method 2
    fcn.initCalibTableRef(calibtable_ref);

    // apply eta-scale
    if (mode==31)
    {
      // Eta bins
      std::vector<double> EtaBins;
      for (int ibin=0; ibin<=50; ibin++)
      {
        EtaBins.push_back(-2.5 + 0.1*ibin);
      }

      if (debug>0) std::cout << " Step 2: Apply Eta-scale to all events, inside mode 31" << std::endl;
      // apply ref eta-scale
      ApplyEtaScaleToAllEvents(EtaBins, scalesref);
    }
    else if (mode==32)
    {
      if (debug>0) std::cout << " Step 2: Apply Eta-scale to all events, inside mode 32" << std::endl;
      // apply ref eta-scale
      ApplyEtaScaleToAllEvents(EtaScaleRef);
    }

    // loop over all the cells,
    // select only those events in this cell,
    // define only this cell one single parameter to be fitted
    for (int ic=0; ic<(int)cells.size(); ic++)
    {
      // cell ix, iy, iz
      const int ix = cells.at(ic).at(0);
      const int iy = cells.at(ic).at(1);
      const int iz = cells.at(ic).at(2);

      // fitting cell above and 8 cells around it, in total at most 3x3=9 cells
      std::vector<int> ixx, iyy, izz;

      // select events for only this cell
      nEvents = SelectEventsInOneCellWith3x3Others(ix, iy, iz, ixx, iyy, izz, _combine);

      // if nEvents==0, it means no electrons in these events give its energy of 10% to this fitting cells, 
      //   therefore, I skip this cell here, set it to be a default value 1.0 0.1, and continue the next cell
      if (nEvents==0) 
      {
        std::cout << "Skip cell (ix, iy, iz) = (" << ix << "," << iy << "," << iz << "): No sufficient events. \n" << std::endl;
        calibRecord AcalibRecord = {ic, ix, iy, iz, 1.0, 0.1, 0, 0};
        calibTable.push_back(AcalibRecord);
        continue;
      }
      // 
      // some debug
      std::cout << " Debug:: print 3x3 cells After Selection" << std::endl;
      for (int icell=0; icell<(int)ixx.size(); icell++)
      {
        std::cout << "(" << ixx.at(icell) << "," << iyy.at(icell) << "," << izz.at(icell) << ")" << std::endl;
      }

      // define the signal fraction
      nSignals = int(signalFraction*(double)nEvents);

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, for Cell: " << std::endl;
      std::cout << "ic = " << ic << ", " << "nEvents = " << nEvents << std::endl;

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, initialize data " << std::endl;
      // initialize fcn using this set of events
      fcn.initData(nEvents, nSignals,
                   E1, EReg1, Eta1, Phi1, nHits1, HitE1, HitIX1, HitIY1, HitIZ1,
                   E2, EReg2, Eta2, Phi2, nHits2, HitE2, HitIX2, HitIY2, HitIZ2,
                   debug, method);

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, initialize PDF " << std::endl;
      // initialize PDF
      fcn.initBWGSParameters(_FitWindowHigh, // windowHigh
                             _FitWindowLow, // windowLow
                             _Zmass, //voigtMass
                             gaus_reso, //voightResolution
                             2.4952, //voigtWidth
                             nSignals,
                             nEvents);

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, set one-cell calibTable" << std::endl;
      // set a calibTable with this cell at the center of the 3x3 cells
      //  the fitting cell is the first in the ixx, iyy, izz vectors
      std::vector<calibRecord> AcalibTable;
      for (int icell=0; icell<(int)ixx.size(); icell++)
      {
         calibRecord AcalibRecord = {ic, ixx.at(icell), iyy.at(icell), izz.at(icell), 1.0, 0.005, 0, 0};
         AcalibTable.push_back(AcalibRecord);
      }

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, initialize one-cell calibTable" << std::endl;
      // initialize calibTable in fcn using the only-one-cell calibTable above
      fcn.initCalibTableFromVector(AcalibTable);

      // also define MnUserParameters
      MnUserParameters Apars;
      char name[100];
      // the fitting cell is the first in the ixx, iyy, izz vectors
      for (int icell=0; icell<(int)ixx.size(); icell++)
      {
        sprintf(name, "par_ix%d_iy%d_iz%d", ixx.at(icell), iyy.at(icell), izz.at(icell));
        Apars.Add(name, 1.0, 0.005); //  calibC
      }

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, define MnMigrad" << std::endl;
      // define migrad
      MnMigrad migrad(fcn, Apars);

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, do fit" << std::endl;
      // minimize
      FunctionMinimum min = migrad();

      // hesse
      MnHesse hesse;
      hesse(fcn, min);

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, print output" << std::endl;
      // print min
      std::cout << "minimum of cell (ix, iy, iz) = (" << ix << "," << iy << "," << iz << "): \n"
      << min << std::endl;

      // fit a second time with the other calib consts fixed
      for (int icell=1; icell<(int)ixx.size(); icell++)
      {
        migrad.Fix((unsigned int)(icell));
      }
      min = migrad();
      hesse(fcn, min);

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, 2nd time, print output" << std::endl;
      // print min
      std::cout << "minimum of cell (ix, iy, iz) = (" << ix << "," << iy << "," << iz << "): \n"
      << min << std::endl;

      // update MnUserParameters
      Apars = min.UserParameters();

      //
      if (debug>0)
      {
        // print scan parameter
        // plot
        MnPlot plot;
        // scan parameters
        MnParameterScan parscan(fcn, Apars);
        if (fabs(Apars.Value(0))<3&&fabs(Apars.Error(0))<1.0)
        {
          std::cout << "Scan of Parameter " << Apars.GetName(0) << " : " << std::endl;
          std::vector< std::pair<double, double> > points_scan = parscan(0);
          plot(points_scan);
        }
        else
        {
          std::cout << "Scan of Parameter " << Apars.GetName(0) << " : Not pass result quality check.. " << std::endl;
        }
      }

      if (debug>0) std::cout << " Step 2: Fit cell-by-cell, update calibTable and store output" << std::endl;
      // update calibTable, we only take this single cell.
      AcalibTable.at(0).c = Apars.Value(0);
      AcalibTable.at(0).cerr = Apars.Error(0);

      // push back the fitted result, again, only has one parameter
      calibTable.push_back(AcalibTable.at(0));

    }
  }
//////// 
 
  if (debug>0) std::cout << " Step 3: Check the results" << std::endl;
  
  if (debug>0) std::cout << " Step 3: Check the results, Print parameters." << std::endl;
  // print parameters
  for (int ic=0; ic<(int)calibTable.size(); ic++)
  {
    printCalibRecord(calibTable.at(ic));
  }
  
  // draw Mee plot again with calibrated constants

  //if (debug>0) std::cout << " Step 3: Check the results, select events." << std::endl;
  //// select events
  //nEvents = SelectSubsetEvents(cells);
  
  //if (debug>0) std::cout << " Step 3: Check the results, init data" << std::endl;
  //// init data
  //fcn.initData(nEvents, nSignals,
  //             E1, EReg1, Eta1, Phi1, nHits1, HitE1, HitIX1, HitIY1, HitIZ1,
  //             E2, EReg2, Eta2, Phi2, nHits2, HitE2, HitIX2, HitIY2, HitIZ2,
  //             debug, method);
  
  if (debug>0) std::cout << " Step 3: Check the results, init calibTable" << std::endl;
  // initialize calibTable in fcn from the fitted calibTable
  fcn.initCalibTableFromVector(calibTable);
  
  //if (debug>0) std::cout << " Step 3: Check the results, init calibTableRef" << std::endl;
  //// initialize calibTableRef
  //fcn.initCalibTableRef(calibtable_ref);
  
  if (debug>0) std::cout << " Step 3: Check the results, print calibTable" << std::endl;
  // print calibTable
  fcn.printCalibTable(calibtable_out);
  
  //// initialize the parameters from from the calibTable in BWGSLikelihoodFCN
  //fcn.initMnUserParametersFromCurrentCalibTable(pars);
  
  //if (debug>0) std::cout << " Step 3: Check the results, store Mee hist after fit" << std::endl;
  //// store hist after fit
  //TH1D* h2 = (TH1D*)fcn.GetMeeTH1DFromMnUserParameters(pars, "hMeeAfter");
  //fout->cd();
  //h2->Write();
  
  
  // close fout
  fout->Close();

  if (debug>1) std::cout << "End program:: return 0" << std::endl;
  return 0;
}


