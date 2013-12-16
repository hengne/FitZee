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

int mode = 6;
// 6: fit energy scale

// Mee calculation method
int method = 6;
// 6: fit energy scale at electron level

// Energy regression Version
std::string RegVersion="V6Elec";

// signal fraction
double signalFraction = 0.99;

// Z mass
double _Zmass = 91.5586; //91.513; // 91.188;

// fit Z mass window High
double _FitWindowHigh = 100.0;
double _FitWindowLow = 80.0;

// Gassian Resolution of the Mee distribution
double gaus_reso = 1.0;

int debug = 0;
char rootfile_in[1000];
char rootfile_out[1000];

char parfile_ref[1000];

// combine, can be EBEB, EBEE, EEEE, EE (any one in EE), EB (any one in EB)
std::string _combine = "EE";

// fitScale : fit calibration or energy scale, calibration is per-cell, energy scale is per-electron
bool fitscale = true;

int main(int argc, char* argv[])
{
  if (argc<4)
  {
    std::cout << argv[0] << " <mode> <input_file.root> <output_file.root> \\\n"
              << "            <signalFraction> <method> <GaussResolution> <RegVersion> <debug> <parameter_file_reference.dat>\n"
    << std::endl;
    return 0;
  }
  mode = atoi(argv[1]);
  sprintf(rootfile_in, "%s", argv[2]);
  sprintf(rootfile_out, "%s", argv[3]);

  if (argc>4)
  {
    signalFraction = atof(argv[4]);
  }
  
  if (argc>5)
  {
    method = atoi(argv[5]);
  }
  
  if (argc>6)
  {
    gaus_reso = atof(argv[6]);
  }

  if (argc>7)
  {
    RegVersion = std::string(argv[7]);
  }

  if (argc>8)
  {
    debug = atoi(argv[8]);
  }

  if (argc>9)
  {
    sprintf(parfile_ref, "%s", argv[9]);
  }
 
  // check 
  if (mode==73||mode==75||mode==77)
  {
    if (argc<=9) 
    {
      std::cout << "Missing parameters reference file. Please run " << argv[0] << " to print usage information. " << std::endl;  
      return 1;
    } 
  }
 
  std::cout << argv[0] << std::endl;
  std::cout << "  mode = " << mode << std::endl;
  std::cout << "  method = " << method << std::endl;
  std::cout << "  rootfile_in = " << rootfile_in << std::endl;
  std::cout << "  rootfile_out = " << rootfile_out << std::endl;
  std::cout << "  signalFraction = " << signalFraction << std::endl;
  std::cout << "  GaussianResolution = " << gaus_reso << std::endl;
  std::cout << "  RegVersion = " << RegVersion << std::endl;
  if (mode==73||mode==75||mode==77) std::cout << "  Parameter Reference File = " << parfile_ref << std::endl;
  
  std::cout << " Start program. " << std::endl;
  
  // reading data
  TChain* tree = new TChain("tree", "tree");
  tree->Add(rootfile_in);
  
  // output root file
  TFile* fout = new TFile(rootfile_out, "recreate");
  
  // Set the branches for the TChain/TTree
  SetTreeBranch(tree);

  // all events
  nEvents = tree->GetEntries();
  nSignals = nEvents;

  // Fill all events into vectors
  FillAllEvents(tree, 2, RegVersion, fitscale);
  if (debug>0) std::cout << " Step 1: fill all events: " << nEvents << std::endl;
  
  // delete the chain no more need it
  tree->Delete();
  
  // define the fitting function
  BWGSLikelihoodFCN fcn;

  fcn.setdebug(debug);

  if (debug>0) std::cout << " Step 1: Initialize PDF overall " << std::endl;
  // initialize PDF
  fcn.initBWGSParameters(_FitWindowHigh, // windowHigh
                         _FitWindowLow, // windowLow
                         _Zmass, //voigtMass
                         gaus_reso, //voightResolution
                         2.4952, //voigtWidth
                         nSignals,
                         nEvents);
  
  // define the energy scales as parameters
  MnUserParameters scales;
  
  // reference scales 
  std::vector<double> scalesref;

  if (mode==73) 
  {
    // read reference scales
    std::ifstream reffile(parfile_ref);
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
 
  // reference scale for mode 75 77
  std::vector<EnergyScale> EtaScaleRef;
  if (mode==75||mode==77)
  {
    // read reference scales
    std::ifstream reffile(parfile_ref);
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

  if (mode==61)
  {
    // mode 61 is to do eta-scale fit

    if (debug>0) std::cout << " Step 2: do fit mode==61 " << std::endl;

    // Eta bins
    std::vector<double> EtaBins;
    for (int ibin=0; ibin<51; ibin++)
    { 
      EtaBins.push_back(-2.5 + 0.1*ibin);
    }  

    // loop over Bins and do fits
    for (int ibin=0; ibin<(int)EtaBins.size()-1; ibin++)
    {
      double bin_min = EtaBins.at(ibin);
      double bin_max = EtaBins.at(ibin+1);

      nEvents = SelectEventsInOneEtaBin(bin_min, bin_max, _combine);

      // define the signal fraction
      nSignals = int(signalFraction*(double)nEvents);

      // initialize fcn using this set of events
      fcn.initDataScale(nEvents, nSignals,
                   E1, EReg1, Eta1, Phi1, UseEle1,
                   E2, EReg2, Eta2, Phi2, UseEle2,
                   debug, method);

      // initialize PDF
      fcn.initBWGSParameters(_FitWindowHigh, // windowHigh
                             _FitWindowLow, // windowLow
                             _Zmass, //voigtMass
                             gaus_reso, //voightResolution
                             2.4952, //voigtWidth
                             nSignals,
                             nEvents);

      // also define a one-par MnUserParameters
      MnUserParameters Apars;
      char name[100];
      sprintf(name, "par_eta_%d", ibin);
      Apars.Add(name, 1.0, 0.001); //  calibC
      // give the same initial par value to the container of all the scales
      scales.Add(name, 1.0, 0.001); //

      // define migrad
      MnMigrad migrad(fcn, Apars);

      // minimize
      FunctionMinimum min = migrad();

      // print min
      std::cout << "minimum of E-scale of EtaBin " << ibin << " (eta_min, eta_max) = ("<< bin_min << ", " << bin_max << "): \n"
      << min << std::endl;

      // hesse
      MnHesse hesse;
      hesse(fcn, min);

      // get parmeters
      Apars = min.UserParameters();
      scales.SetValue(ibin, Apars.Value(0));
      scales.SetError(ibin, Apars.Error(0));

      std::cout << "Bin " << ibin << " : " << Apars.Value(0) << " " << Apars.Error(0) << std::endl;

    } // for (int ibin=0;...

    std::cout << "Fitting Results: " << std::endl;
    for (int ibin=0; ibin<(int)EtaBins.size()-1; ibin++)
    {
      std::cout << "EtaBin[" << ibin << "] (" << EtaBins.at(ibin) << "<Eta<" << EtaBins.at(ibin+1) << ") : " 
       << scales.Value(ibin) << " +/- " << scales.Error(ibin) << std::endl; 
    }

    fout->cd();

    TH1D* hist = new TH1D("hEtaScale", "hEtaScale", 50, -2.5, 2.5);
    hist->Sumw2();
    hist->SetMarkerStyle(20);

    for (int ibin=0; ibin<(int)EtaBins.size()-1; ibin++)
    {
      hist->SetBinContent(ibin+1, scales.Value(ibin));
      hist->SetBinError(ibin+1, scales.Error(ibin));
    }

    hist->Write();
  } // mode==61
  else if (mode==71 || mode==72)
  {
    // mode 71 is to do eta-scale fit, instead of fit one bin after another, fit all of them together

    if (debug>0) std::cout << " Step 2: do fit mode==71 " << std::endl;

    // Eta bins
    std::vector<double> EtaBins;
    for (int ibin=0; ibin<=50; ibin++)
    {
      EtaBins.push_back(-2.5 + 0.1*ibin);
    }

    nEvents = AddEtaBinNumberToElectrons(EtaBins, _combine);
 
    // define the signal fraction
    nSignals = int(signalFraction*(double)nEvents);

    if (debug>0) std::cout << " Step 2: init data and pars : nEvents = " << nEvents << std::endl;

    // initialize fcn using this set of events
    fcn.initDataScale(nEvents, nSignals,
                   E1, EReg1, Eta1, Phi1, ScaleBin1,
                   E2, EReg2, Eta2, Phi2, ScaleBin2,
                   debug, method);

    // initialize PDF
    fcn.initBWGSParameters(_FitWindowHigh, // windowHigh
                             _FitWindowLow, // windowLow
                             _Zmass, //voigtMass
                             gaus_reso, //voightResolution
                             2.4952, //voigtWidth
                             nSignals,
                             nEvents);    

    // also define MnUserParameters
    MnUserParameters Apars;

    char name[100];
    for (int ibin=0; ibin<(int)EtaBins.size()-1; ibin++)
    {
      sprintf(name, "par_eta_%d", ibin);
      Apars.Add(name, 1.0, 0.0001); //  calibC
    }

    if (debug>0) std::cout << " Step 2: define migrad " << std::endl;
    // define migrad
    MnMigrad migrad(fcn, Apars);

    if (debug>0) std::cout << " Step 2: do fit " << std::endl;
    // minimize
    FunctionMinimum min = migrad();

    // print min
    std::cout << "minimum of E-scale of EtaBin : \n"
    << min << std::endl;

    if (debug>0) std::cout << " Step 2: hesse estimation of uncertainties " << std::endl;
    // hesse
    MnHesse hesse;
    hesse(fcn, min);

    if (debug>0) std::cout << " Step 2: get parameters " << std::endl;
    // get parmeters
    Apars = min.UserParameters();


    // if mode==72, fit again each parameters one by one with all the other parameters fixed.
    if (mode==72)
    {

      if (debug>0) std::cout << " Step 2-1: Refit with fixed other parameters. " << std::endl;

      // fix all parameters
      for (int ip=0; ip<(int)migrad.Params().size(); ip++)
      {
        migrad.Fix(ip); 
      }
      
      // loop all the parameters one by one while fixing all the rests
      for (int ip=0; ip<(int)migrad.Params().size(); ip++)
      {
        migrad.Release(ip);
        min = migrad();
        Apars.SetValue(ip, migrad.Value(ip));
        Apars.SetError(ip, migrad.Error(ip));
        migrad.Fix(ip);
      }      
    
    }

    // pass the "Apars" to "scales" 
    scales = Apars;

    // 
    if (debug>0)
    {
      // print scan parameter
      // plot
      MnPlot plot;
      // scan parameters
      MnParameterScan parscan(fcn, Apars);
      for (int ip=0; ip<(int)Apars.Params().size(); ip++)
      {
        if (fabs(Apars.Value(ip))<3&&fabs(Apars.Error(ip))<1.0)
        {
          std::cout << "Scan of Parameter " << Apars.GetName(ip) << " : " << std::endl;
          std::vector< std::pair<double, double> > points_scan = parscan(ip);
          plot(points_scan);
        }
        else
        {
          std::cout << "Scan of Parameter " << Apars.GetName(ip) << " : Not pass result quality check.. " << std::endl;
        }
      }
    }


    std::cout << "Fitting Results: " << std::endl;
    for (int ibin=0; ibin<(int)EtaBins.size()-1; ibin++)
    {
      std::cout << "EtaBin[" << ibin << "] (" << EtaBins.at(ibin) << "<Eta<" << EtaBins.at(ibin+1) << ") : "
         << Apars.Value(ibin) << " +/- " << Apars.Error(ibin) << std::endl;
    }

    fout->cd();

    TH1D* hist = new TH1D("hEtaScale", "hEtaScale", 50, -2.5, 2.5);
    hist->Sumw2();
    hist->SetMarkerStyle(20);

    for (int ibin=0; ibin<(int)EtaBins.size()-1; ibin++)
    {
      hist->SetBinContent(ibin+1, Apars.Value(ibin));
      hist->SetBinError(ibin+1, Apars.Error(ibin));
    }

    hist->Write();
  } // mode==71 or 72
  else if (mode==73)
  {
    // mode 73 is similar to mode 71, is to do eta-scale fit, instead of fit one bin after another, fit all of them together,
    // but different from mode 71 is that mode 73 applies a pre-defined eta-scale before the fits. 

    if (debug>0) std::cout << " Step 2: do fit mode==73 " << std::endl;

    // Eta bins
    std::vector<double> EtaBins;
    for (int ibin=0; ibin<=50; ibin++)
    {
      EtaBins.push_back(-2.5 + 0.1*ibin);
    }

    // apply ref eta-scale
    ApplyEtaScaleToAllEvents(EtaBins, scalesref);

    // give eta-bin numbers
    nEvents = AddEtaBinNumberToElectrons(EtaBins, _combine);

    // define the signal fraction
    nSignals = int(signalFraction*(double)nEvents);

    if (debug>0) std::cout << " Step 2: init data and pars : nEvents = " << nEvents << std::endl;

    // initialize fcn using this set of events
    fcn.initDataScale(nEvents, nSignals,
                   E1, EReg1, Eta1, Phi1, ScaleBin1,
                   E2, EReg2, Eta2, Phi2, ScaleBin2,
                   debug, method);

    // initialize PDF
    fcn.initBWGSParameters(_FitWindowHigh, // windowHigh
                             _FitWindowLow, // windowLow
                             _Zmass, //voigtMass
                             gaus_reso, //voightResolution
                             2.4952, //voigtWidth
                             nSignals,
                             nEvents);

    // also define MnUserParameters
    MnUserParameters Apars;

    char name[100];
    for (int ibin=0; ibin<(int)EtaBins.size()-1; ibin++)
    {
      sprintf(name, "par_eta_%d", ibin);
      Apars.Add(name, 1.0, 0.0001); //  calibC
    }

    if (debug>0) std::cout << " Step 2: define migrad " << std::endl;
    // define migrad
    MnMigrad migrad(fcn, Apars);

    if (debug>0) std::cout << " Step 2: do fit " << std::endl;
    // minimize
    FunctionMinimum min = migrad();

    // print min
    std::cout << "minimum of E-scale of EtaBin : \n"
    << min << std::endl;

    if (debug>0) std::cout << " Step 2: hesse estimation of uncertainties " << std::endl;
    // hesse
    MnHesse hesse;
    hesse(fcn, min);

    if (debug>0) std::cout << " Step 2: get parameters " << std::endl;
    // get parmeters
    Apars = min.UserParameters();

    // pass the "Apars" to "scales"
    scales = Apars;

    //
    if (debug>0)
    {
      // print scan parameter
      // plot
      MnPlot plot;
      // scan parameters
      MnParameterScan parscan(fcn, Apars);
      for (int ip=0; ip<(int)Apars.Params().size(); ip++)
      {
        if (fabs(Apars.Value(ip))<3&&fabs(Apars.Error(ip))<1.0)
        {
          std::cout << "Scan of Parameter " << Apars.GetName(ip) << " : " << std::endl;
          std::vector< std::pair<double, double> > points_scan = parscan(ip);
          plot(points_scan);
        }
        else
        {
          std::cout << "Scan of Parameter " << Apars.GetName(ip) << " : Not pass result quality check.. " << std::endl;
        }
      }
    }


    std::cout << "Fitting Results: " << std::endl;
    for (int ibin=0; ibin<(int)EtaBins.size()-1; ibin++)
    {
      std::cout << "EtaBin[" << ibin << "] (" << EtaBins.at(ibin) << "<Eta<" << EtaBins.at(ibin+1) << ") : "
         << Apars.Value(ibin) << " +/- " << Apars.Error(ibin) << std::endl;
    }

    fout->cd();

    TH1D* hist = new TH1D("hEtaScale", "hEtaScale", 50, -2.5, 2.5);
    hist->Sumw2();
    hist->SetMarkerStyle(20);

    for (int ibin=0; ibin<(int)EtaBins.size()-1; ibin++)
    {
      hist->SetBinContent(ibin+1, Apars.Value(ibin));
      hist->SetBinError(ibin+1, Apars.Error(ibin));
    }

    hist->Write();
  } // mode==73
  else if (mode==74 || mode==75)
  {
    // mode 74 is same to mode 71, just with a more flexable binning
    // mode 75 same as mode 74, just with a Etascale applied in advance
    if (debug>0) std::cout << " Step 2: do fit mode==73 " << std::endl;

    // Eta bins
    std::vector<EnergyScale> EtaScale;
    // 20 EE- bins
    for (int ibin=0; ibin<20; ibin++)
    {
      EnergyScale Ascale = {-2.5+0.05*ibin, -2.5+0.05*(ibin+1), 1.0, 0.001};
      EtaScale.push_back(Ascale);
    }
    // 10 EB bins
    for (int ibin=0; ibin<10; ibin++)
    {
      EnergyScale Ascale = {-1.5+0.3*ibin, -1.5+0.3*(ibin+1), 1.0, 0.001};
      EtaScale.push_back(Ascale);
    }
    // 20 EE+ bins
    for (int ibin=0; ibin<20; ibin++)
    {
      EnergyScale Ascale = {1.5+0.05*ibin, 1.5+0.05*(ibin+1), 1.0, 0.001};
      EtaScale.push_back(Ascale);
    }
   
    if (mode==75)
    {
      // apply ref eta-scale
      ApplyEtaScaleToAllEvents(EtaScaleRef);
    }

    // give eta-bin numbers
    nEvents = AddEtaBinNumberToElectrons(EtaScale, _combine);

    // define the signal fraction
    nSignals = int(signalFraction*(double)nEvents);

    if (debug>0) std::cout << " Step 2: init data and pars : nEvents = " << nEvents << std::endl;

    // initialize fcn using this set of events
    fcn.initDataScale(nEvents, nSignals,
                   E1, EReg1, Eta1, Phi1, ScaleBin1,
                   E2, EReg2, Eta2, Phi2, ScaleBin2,
                   debug, method);

    // initialize PDF
    fcn.initBWGSParameters(_FitWindowHigh, // windowHigh
                             _FitWindowLow, // windowLow
                             _Zmass, //voigtMass
                             gaus_reso, //voightResolution
                             2.4952, //voigtWidth
                             nSignals,
                             nEvents);

    // also define MnUserParameters
    MnUserParameters Apars;

    char name[100];
    for (int ibin=0; ibin<(int)EtaScale.size(); ibin++)
    {
      sprintf(name, "par_eta_%d", ibin);
      Apars.Add(name, 1.0, 0.0001); //  calibC
    }

    if (debug>0) std::cout << " Step 2: define migrad " << std::endl;
    // define migrad
    MnMigrad migrad(fcn, Apars);

    if (debug>0) std::cout << " Step 2: do fit " << std::endl;
    // minimize
    FunctionMinimum min = migrad();

    // print min
    std::cout << "minimum of E-scale of EtaBin : \n"
    << min << std::endl;

    if (debug>0) std::cout << " Step 2: hesse estimation of uncertainties " << std::endl;
    // hesse
    //MnHesse hesse;
    //hesse(fcn, min);

    if (debug>0) std::cout << " Step 2: get parameters " << std::endl;
    // get parmeters
    Apars = min.UserParameters();

    // pass the "Apars" to "scales"
    scales = Apars;

    //
    if (debug>0)
    {
      // print scan parameter
      // plot
      MnPlot plot;
      // scan parameters
      MnParameterScan parscan(fcn, Apars);
      for (int ip=0; ip<(int)Apars.Params().size(); ip++)
      {
        if (fabs(Apars.Value(ip))<3&&fabs(Apars.Error(ip))<1.0)
        {
          std::cout << "Scan of Parameter " << Apars.GetName(ip) << " : " << std::endl;
          std::vector< std::pair<double, double> > points_scan = parscan(ip);
          plot(points_scan);
        }
        else
        {
          std::cout << "Scan of Parameter " << Apars.GetName(ip) << " : Not pass result quality check.. " << std::endl;
        }
      }
    }


    std::cout << "Fitting Results: " << std::endl;
    for (int ibin=0; ibin<(int)EtaScale.size(); ibin++)
    {
      std::cout << "EtaBin[" << ibin << "] (" << EtaScale.at(ibin).min << "<Eta<" << EtaScale.at(ibin).max << ") : "
         << Apars.Value(ibin) << " +/- " << Apars.Error(ibin) << std::endl;
    }

    fout->cd();

    double histBins[1000];
    for (int ibin=0; ibin<(int)EtaScale.size(); ibin++)
    {
      histBins[ibin] = EtaScale.at(ibin).min;
      if (ibin==(int)EtaScale.size()-1)
      {
        histBins[ibin+1] = EtaScale.at(ibin).max;
      }  
    }
    

    TH1D* hist = new TH1D("hEtaScale", "hEtaScale", (int)EtaScale.size(), histBins);
    hist->Sumw2();
    hist->SetMarkerStyle(20);

    for (int ibin=0; ibin<(int)Apars.Params().size(); ibin++)
    {
      hist->SetBinContent(ibin+1, Apars.Value(ibin));
      hist->SetBinError(ibin+1, Apars.Error(ibin));
    }

    hist->Write();
  } // mode==74 || mode==75
  else if (mode==76 || mode==77)
  {
    // mode 76 is same to mode 74 having a more flexable binning compared to 71, but 76 has more finner binning than 74
    // mode 77 same as mode 76, just with a Etascale applied in advance
    if (debug>0) std::cout << " Step 2: do fit mode==76 " << std::endl;

    // Eta bins
    std::vector<EnergyScale> EtaScale;
    // 40 EE- bins
    for (int ibin=0; ibin<40; ibin++)
    {
      EnergyScale Ascale = {-2.5+0.025*ibin, -2.5+0.025*(ibin+1), 1.0, 0.001};
      EtaScale.push_back(Ascale);
    }
    // 10 EB bins
    for (int ibin=0; ibin<10; ibin++)
    {
      EnergyScale Ascale = {-1.5+0.3*ibin, -1.5+0.3*(ibin+1), 1.0, 0.001};
      EtaScale.push_back(Ascale);
    }
    // 40 EE+ bins
    for (int ibin=0; ibin<40; ibin++)
    {
      EnergyScale Ascale = {1.5+0.025*ibin, 1.5+0.025*(ibin+1), 1.0, 0.001};
      EtaScale.push_back(Ascale);
    }

    if (mode==77)
    {
      // apply ref eta-scale
      ApplyEtaScaleToAllEvents(EtaScaleRef);
    }

    // give eta-bin numbers
    nEvents = AddEtaBinNumberToElectrons(EtaScale, _combine);

    // define the signal fraction
    nSignals = int(signalFraction*(double)nEvents);

    if (debug>0) std::cout << " Step 2: init data and pars : nEvents = " << nEvents << std::endl;

    // initialize fcn using this set of events
    fcn.initDataScale(nEvents, nSignals,
                   E1, EReg1, Eta1, Phi1, ScaleBin1,
                   E2, EReg2, Eta2, Phi2, ScaleBin2,
                   debug, method);

    // initialize PDF
    fcn.initBWGSParameters(_FitWindowHigh, // windowHigh
                             _FitWindowLow, // windowLow
                             _Zmass, //voigtMass
                             gaus_reso, //voightResolution
                             2.4952, //voigtWidth
                             nSignals,
                             nEvents);

    // also define MnUserParameters
    MnUserParameters Apars;

    char name[100];
    for (int ibin=0; ibin<(int)EtaScale.size(); ibin++)
    {
      sprintf(name, "par_eta_%d", ibin);
      Apars.Add(name, 1.0, 0.0001); //  calibC
    }

    if (debug>0) std::cout << " Step 2: define migrad " << std::endl;
    // define migrad
    MnMigrad migrad(fcn, Apars);

    if (debug>0) std::cout << " Step 2: do fit " << std::endl;
    // minimize
    FunctionMinimum min = migrad();

    // print min
    std::cout << "minimum of E-scale of EtaBin : \n"
    << min << std::endl;

    if (debug>0) std::cout << " Step 2: hesse estimation of uncertainties " << std::endl;
    // hesse
    //MnHesse hesse;
    //hesse(fcn, min);

    if (debug>0) std::cout << " Step 2: get parameters " << std::endl;
    // get parmeters
    Apars = min.UserParameters();

    // pass the "Apars" to "scales"
    scales = Apars;

    //
    if (debug>0)
    {
      // print scan parameter
      // plot
      MnPlot plot;
      // scan parameters
      MnParameterScan parscan(fcn, Apars);
      for (int ip=0; ip<(int)Apars.Params().size(); ip++)
      {
        if (fabs(Apars.Value(ip))<3&&fabs(Apars.Error(ip))<1.0)
        {
          std::cout << "Scan of Parameter " << Apars.GetName(ip) << " : " << std::endl;
          std::vector< std::pair<double, double> > points_scan = parscan(ip);
          plot(points_scan);
        }
        else
        {
          std::cout << "Scan of Parameter " << Apars.GetName(ip) << " : Not pass result quality check.. " << std::endl;
        }
      }
    }


    std::cout << "Fitting Results: " << std::endl;
    for (int ibin=0; ibin<(int)EtaScale.size(); ibin++)
    {
      std::cout << "EtaBin[" << ibin << "] (" << EtaScale.at(ibin).min << "<Eta<" << EtaScale.at(ibin).max << ") : "
         << Apars.Value(ibin) << " +/- " << Apars.Error(ibin) << std::endl;
    }

    fout->cd();

    double histBins[1000];
    for (int ibin=0; ibin<(int)EtaScale.size(); ibin++)
    {
      histBins[ibin] = EtaScale.at(ibin).min;
      if (ibin==(int)EtaScale.size()-1)
      {
        histBins[ibin+1] = EtaScale.at(ibin).max;
      }
    }


    TH1D* hist = new TH1D("hEtaScale", "hEtaScale", (int)EtaScale.size(), histBins);
    hist->Sumw2();
    hist->SetMarkerStyle(20);

    for (int ibin=0; ibin<(int)Apars.Params().size(); ibin++)
    {
      hist->SetBinContent(ibin+1, Apars.Value(ibin));
      hist->SetBinError(ibin+1, Apars.Error(ibin));
    }

    hist->Write();
  } // mode==76 77


    
  // close fout
  fout->Close();

  return 0;
}


