// Authors: Hengne Li of UVa at CERN in 2013

#include "calibRecord.h"
#include "voigt.h"
#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/MnUserParameters.h"
#include "TComplex.h"
#include "TH1D.h"
#include "math.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <time.h>


namespace ROOT {
  
  namespace Minuit2 {

// likelihood function of Breit-Wigner x Gaussian function
class BWGSLikelihoodFCN : public FCNBase {

public:

  BWGSLikelihoodFCN(void) {};

  BWGSLikelihoodFCN(int nEvents, int nSignals,
                    const std::vector<double*>& E1,
                    const std::vector<double*>& EReg1,
                    const std::vector<double*>& Eta1,
                    const std::vector<double*>& Phi1,
                    const std::vector<int*>& nHits1,
                    const std::vector< std::vector<double>* >& HitE1,
                    const std::vector< std::vector<int>* >& HitIX1,
                    const std::vector< std::vector<int>* >& HitIY1,
                    const std::vector< std::vector<int>* >& HitIZ1,
                    const std::vector<double*>& E2,
                    const std::vector<double*>& EReg2,
                    const std::vector<double*>& Eta2,
                    const std::vector<double*>& Phi2,
                    const std::vector<int*>& nHits2,
                    const std::vector< std::vector<double>* >& HitE2,
                    const std::vector< std::vector<int>* >& HitIX2,
                    const std::vector< std::vector<int>* >& HitIY2,
                    const std::vector< std::vector<int>* >& HitIZ2,
                    const int debug=0,
                    const int method=0)
  {
    // debug
    _debug = debug; 
 
    // initialize parameters for the Likelihood function
    initBWGSParameters(120.0, // windowHigh
                       60.0, // windowLow
                       91.188, //voigtMass
                       2.2, //voightResolution
                       2.4952, //voigtWidth
                       nSignals,
                       nEvents);
    
    // initialize data
    initData(nEvents, nSignals,
             E1, EReg1, Eta1, Phi1, nHits1, HitE1, HitIX1, HitIY1, HitIZ1,
             E2, EReg2, Eta2, Phi2, nHits2, HitE2, HitIX2, HitIY2, HitIZ2,
             debug, method);
    
    // initialize calibration constant table
    initCalibTable();
    
  }
  
  ~BWGSLikelihoodFCN() {}
  
  void setdebug(int debug) { _debug = debug; }
 
  double operator()(const std::vector<double>& par) const
  {
    // recalculate masses
    recalcM(_method, par);
   
    // calculate negative log likelihood
    double nll(0.0);
    double mass(0.0);
    double probval(0.0);
    double weight(1.0);
    for (int i=0; i<_nEvents; i++)
    {
      mass = _MassRecalc.at(i);
      
      if (_method==5) 
      {
        weight = _MassRecalcWeight.at(i);
      }

      if (_debug>2)
      {
        std::cout << "mass = " << mass << std::endl;
        if (_method==5) 
        {
          std::cout << "weight = " << weight << std::endl;
        }
      }
      
      if (mass < _windowLow) mass = _windowLow;
      if (mass > _windowHigh) mass = _windowHigh;
      //if (mass<_windowLow||mass>_windowHigh) continue;
      

      probval = (1.0-_voigtFraction)/(_windowHigh-_windowLow)
                    + _voigtFraction/_voigtNorm*(_voigt.val(mass)); 
      if (_method==5) 
      {
        nll+=-weight*log(probval);
      }
      else
      {
        nll+=-log(probval);
      }
    }
    
    if (_debug>2) std::cout << "LogL = " << nll << std::endl;
    return nll;
    
  }
  
  double Up() const {return 0.5;}
 
  void initBWGSParameters(double windowHigh, double windowLow, double voigtMass,
                          double voigtResolution, double voigtWidth, int nSignals,
                          int nEvents)
  {
    _windowHigh = windowHigh;
    _windowLow = windowLow;
    _voigtMass = voigtMass;
    _voigtResolution = voigtResolution;
    _voigtWidth = voigtWidth;
    _nSignals = nSignals;
    _nEvents = nEvents;
    _voigtFraction = (double)_nSignals/(double)_nEvents;
    _voigt.setParams(_voigtMass, _voigtResolution, _voigtWidth);
    //_voigtNorm = GetVoigtNorm(_windowLow, _windowHigh, _voigtMass, _voigtResolution, _voigtWidth);
    _voigtNorm = _voigt.integral(_windowLow, _windowHigh);
  }
 
  void initDataScale(int nEvents, int nSignals,
                const std::vector<double*>& E1,
                const std::vector<double*>& EReg1,
                const std::vector<double*>& Eta1,
                const std::vector<double*>& Phi1,
                const std::vector<bool>& UseEle1,
                const std::vector<double*>& E2,
                const std::vector<double*>& EReg2,
                const std::vector<double*>& Eta2,
                const std::vector<double*>& Phi2,
                const std::vector<bool>& UseEle2,
                const int debug=0,
                const int method=0)
  {
    // debug
    _debug = debug;

    // method to calculate Mee
    _method = method;

    // nEvents and nSignals
    _nEvents = nEvents;
    _nSignals = nSignals;

    // initialize data
    _E1 = E1;
    _EReg1 = EReg1;
    _Eta1 = Eta1;
    _Phi1 = Phi1;
    _UseEle1 = UseEle1;
    _E2 = E2;
    _EReg2 = EReg2;
    _Eta2 = Eta2;
    _Phi2 = Phi2;
    _UseEle2 = UseEle2;

    // check nEvents correct or not
    if (_nEvents != (int)E1.size()) {
      std::cout << "BWGSLikelihoodFCN:: "
      << "The nEvents given and the nEvents from the data vectors (e.g. E1[nEvents]) are not the same. \n"
      << "    Use the nEvents in data vectors as the nEvents."
      << std::endl;
      _nEvents = (int)E1.size();
    }

    // clear calculated data vectors and re-calculate again right after
    _Mass.clear();
    _MassRecalc.clear();
    _MassRecalcWeight.clear();

    // calculate data
    for (int i=0; i<(int)E1.size(); i++){
      double mass = CalcMass(*(E1.at(i)), *(Eta1.at(i)), *(Phi1.at(i)), *(E2.at(i)), *(Eta2.at(i)), *(Phi2.at(i)));
      _Mass.push_back(mass);
      _MassRecalc.push_back(mass);
      _MassRecalcWeight.push_back(1.0);
    }
  }

  void initDataScale(int nEvents, int nSignals,
                const std::vector<double*>& E1,
                const std::vector<double*>& EReg1,
                const std::vector<double*>& Eta1,
                const std::vector<double*>& Phi1,
                const std::vector<int>& ScaleBin1,
                const std::vector<double*>& E2,
                const std::vector<double*>& EReg2,
                const std::vector<double*>& Eta2,
                const std::vector<double*>& Phi2,
                const std::vector<int>& ScaleBin2,
                const int debug=0,
                const int method=0)
  {
    // debug
    _debug = debug;

    // method to calculate Mee
    _method = method;

    // nEvents and nSignals
    _nEvents = nEvents;
    _nSignals = nSignals;

    // initialize data
    _E1 = E1;
    _EReg1 = EReg1;
    _Eta1 = Eta1;
    _Phi1 = Phi1;
    _ScaleBin1 = ScaleBin1;
    _E2 = E2;
    _EReg2 = EReg2;
    _Eta2 = Eta2;
    _Phi2 = Phi2;
    _ScaleBin2 = ScaleBin2;

    // check nEvents correct or not
    if (_nEvents != (int)E1.size()) {
      std::cout << "BWGSLikelihoodFCN:: "
      << "The nEvents given and the nEvents from the data vectors (e.g. E1[nEvents]) are not the same. \n"
      << "    Use the nEvents in data vectors as the nEvents."
      << std::endl;
      _nEvents = (int)E1.size();
    }

    // clear calculated data vectors and re-calculate again right after
    _Mass.clear();
    _MassRecalc.clear();
    _MassRecalcWeight.clear();     

    // clear calculated data vectors and re-calculate again right after
    _Mass.clear();
    _MassRecalc.clear();
    _MassRecalcWeight.clear();

    // calculate data
    for (int i=0; i<(int)E1.size(); i++){
      double mass = CalcMass(*(E1.at(i)), *(Eta1.at(i)), *(Phi1.at(i)), *(E2.at(i)), *(Eta2.at(i)), *(Phi2.at(i)));
      _Mass.push_back(mass);
      _MassRecalc.push_back(mass);
      _MassRecalcWeight.push_back(1.0);
    }

  }

  void initDataSeed(int nEvents, int nSignals,
                const std::vector<double*>& E1,
                const std::vector<double*>& EReg1,
                const std::vector<double*>& Eta1,
                const std::vector<double*>& Phi1,
                const std::vector<bool>& UseEle1,
                const std::vector<double*>& RawEEcal1,
                const std::vector<double*>& E2,
                const std::vector<double*>& EReg2,
                const std::vector<double*>& Eta2,
                const std::vector<double*>& Phi2,
                const std::vector<bool>& UseEle2,
                const std::vector<double*>& RawEEcal2,
                const int debug=0,
                const int method=0)
  {
    // debug
    _debug = debug;

    // method to calculate Mee
    _method = method;

    // nEvents and nSignals
    _nEvents = nEvents;
    _nSignals = nSignals;

    // initialize data
    _E1 = E1;
    _EReg1 = EReg1;
    _Eta1 = Eta1;
    _Phi1 = Phi1;
    _UseEle1 = UseEle1;
    _RawEEcal1 = RawEEcal1;
    _E2 = E2;
    _EReg2 = EReg2;
    _Eta2 = Eta2;
    _Phi2 = Phi2;
    _UseEle2 = UseEle2;
    _RawEEcal2 = RawEEcal2;

    // check nEvents correct or not
    if (_nEvents != (int)E1.size()) {
      std::cout << "BWGSLikelihoodFCN:: "
      << "The nEvents given and the nEvents from the data vectors (e.g. E1[nEvents]) are not the same. \n"
      << "    Use the nEvents in data vectors as the nEvents."
      << std::endl;
      _nEvents = (int)E1.size();
    }

    // clear calculated data vectors and re-calculate again right after
    _Mass.clear();
    _MassRecalc.clear();
    _MassRecalcWeight.clear();

    // clear calculated data vectors and re-calculate again right after
    _Mass.clear();
    _MassRecalc.clear();
    _MassRecalcWeight.clear();

    // calculate data
    for (int i=0; i<(int)E1.size(); i++){
      double mass = CalcMass(*(E1.at(i)), *(Eta1.at(i)), *(Phi1.at(i)), *(E2.at(i)), *(Eta2.at(i)), *(Phi2.at(i)));
      _Mass.push_back(mass);
      _MassRecalc.push_back(mass);
      _MassRecalcWeight.push_back(1.0);
    }

  }

 
  void initData(int nEvents, int nSignals,
                const std::vector<double*>& E1,
                const std::vector<double*>& EReg1,
                const std::vector<double*>& Eta1,
                const std::vector<double*>& Phi1,
                const std::vector<int*>& nHits1,
                const std::vector< std::vector<double>* >& HitE1,
                const std::vector< std::vector<int>* >& HitIX1,
                const std::vector< std::vector<int>* >& HitIY1,
                const std::vector< std::vector<int>* >& HitIZ1,
                const std::vector<double*>& E2,
                const std::vector<double*>& EReg2,
                const std::vector<double*>& Eta2,
                const std::vector<double*>& Phi2,
                const std::vector<int*>& nHits2,
                const std::vector< std::vector<double>* >& HitE2,
                const std::vector< std::vector<int>* >& HitIX2,
                const std::vector< std::vector<int>* >& HitIY2,
                const std::vector< std::vector<int>* >& HitIZ2,
                const int debug=0,
                const int method=0)
  {
    // debug
    _debug = debug;
    
    // method to calculate Mee
    _method = method;
    
    // nEvents and nSignals
    _nEvents = nEvents;
    _nSignals = nSignals;
    
    // initialize data
    _E1 = E1;
    _EReg1 = EReg1;
    _Eta1 = Eta1;
    _Phi1 = Phi1;
    _nHits1 = nHits1;
    _HitE1 = HitE1;
    _HitIX1 = HitIX1;
    _HitIY1 = HitIY1;
    _HitIZ1 = HitIZ1;
    _E2 = E2;
    _EReg2 = EReg2;
    _Eta2 = Eta2;
    _Phi2 = Phi2;
    _nHits2 = nHits2;
    _HitE2 = HitE2;
    _HitIX2 = HitIX2;
    _HitIY2 = HitIY2;
    _HitIZ2 = HitIZ2;
    
    // check nEvents correct or not
    if (_nEvents != (int)E1.size()) {
      std::cout << "BWGSLikelihoodFCN:: "
      << "The nEvents given and the nEvents from the data vectors (e.g. E1[nEvents]) are not the same. \n"
      << "    Use the nEvents in data vectors as the nEvents."
      << std::endl;
      _nEvents = (int)E1.size();
    }
    
    // clear calculated data vectors and re-calculate again right after
    _Mass.clear();
    _MassRecalc.clear();
    _MassRecalcWeight.clear();
    
    // calculate data
    for (int i=0; i<(int)E1.size(); i++){
      double mass = CalcMass(*(E1.at(i)), *(Eta1.at(i)), *(Phi1.at(i)), *(E2.at(i)), *(Eta2.at(i)), *(Phi2.at(i)));
      _Mass.push_back(mass);
      _MassRecalc.push_back(mass);
      _MassRecalcWeight.push_back(1.0);
    }
  }
  
  void initCalibTable(const char* filename = "calibTable_in.dat")
  {
    std::ifstream myfile(filename);
    std::string line;
    if (myfile.is_open())
    {
      _calibTable.clear();
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
        _calibTable.push_back(calib);
      }
      myfile.close();
    }
    
    
  }
  
  void initCalibTableRef(const char* filename = "calibTableRef_in.dat")
  {
    std::ifstream myfile(filename);
    std::string line;
    if (myfile.is_open())
    {
      _calibTableRef.clear();
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

        _calibTableRef[calib.ix][calib.iy][calib.iz] = calib;
        
      }
      myfile.close();
    }
    
    
  }
  
  void initCalibTableFromVector(std::vector<calibRecord>& calibTable)
  {
    // clear its own calibTable
    _calibTable.clear();
    // directly refer to the one that is given
    _calibTable = calibTable;
  }

  
  void initCalibTableRefFromVector(std::vector<calibRecord>& calibTableRef)
  {
    // clear its own calibTableRef
    _calibTableRef.clear();
    // directly refer to the one that is given
    for (int idx=0; idx<(int)calibTableRef.size(); idx++)
    {
      calibRecord record = calibTableRef.at(idx);
      _calibTableRef[record.ix][record.iy][record.iz] = record;
    }
  }
  
  std::vector<calibRecord> getCalibTable()
  {
    return _calibTable;
  }
  
  std::vector<calibRecord> getCalibTableRef()
  {
    std::map<int, std::map<int, std::map<int, calibRecord> > >::iterator _it_ix;
    std::map<int, std::map<int, calibRecord> >::iterator _it_iy;
    std::map<int, calibRecord>::iterator _it_iz;
    std::vector<calibRecord> aTable;
    for (_it_ix=_calibTableRef.begin(); _it_ix!=_calibTableRef.end(); ++_it_ix)
    {
      for (_it_iy = _it_ix->second.begin(); _it_iy != _it_ix->second.end(); ++_it_iy)
      {
        for (_it_iz = _it_iy->second.begin(); _it_iz != _it_iy->second.end(); ++_it_iz)
        {
          aTable.push_back(_it_iz->second);
        }
      }    
    }
    return aTable;
  }
  
  void printCalibTable(const char* filename = "calibTable_out.dat")
  {
    std::ofstream myfile(filename);
    if (myfile.is_open())
    {
      for (int i=0; i<(int)_calibTable.size(); i++)
      {
        myfile << _calibTable.at(i).idx << " "
               << _calibTable.at(i).ix << " "
               << _calibTable.at(i).iy<< " "
               << _calibTable.at(i).iz<< " "
               << _calibTable.at(i).c << " "
               << _calibTable.at(i).cerr << " "
               << _calibTable.at(i).fixed << " "
               << _calibTable.at(i).nfits << " "
               << std::endl;
      }
      myfile.close();
    }
  }
  
  void initMnUserParametersFromCurrentCalibTable(MnUserParameters& Pars)
  {
    // size of the reference MnUserParameters
    int SizeOfMnPars = (int)Pars.Parameters().size();
    // size of the current calibTable
    int SizeOfCalibTable = (int)_calibTable.size();
    if (_debug>2) std::cout << "SizeOfMnPars: " << SizeOfMnPars 
              << "; SizeOfCalibTable: " << SizeOfCalibTable << std::endl;
    // add parameters to MnUserParameters if it doesn't have enough parameters
    if (SizeOfMnPars<SizeOfCalibTable)
    {
      for (int i=SizeOfMnPars; i<SizeOfCalibTable; i++)
      {
        char name[200];
        sprintf(name, "par%d_ix%d_iy%d_iz%d", i, _calibTable.at(i).ix, _calibTable.at(i).iy, _calibTable.at(i).iz);
        Pars.Add(std::string(name), 1.0, 0.1);
        if (_debug>2) std::cout << "Init Par: " << name << std::endl;
      }
    }
    // and fix parameters that more than that in calibTables
    else if (SizeOfMnPars>SizeOfCalibTable)
    {
      for (int i=SizeOfCalibTable; i<SizeOfMnPars; i++)
      {
        Pars.Fix(i);
      }
    }
    // give all the values to MnUserParameters
    for (int i=0; i<SizeOfCalibTable; i++)
    {
      Pars.SetValue(i, _calibTable.at(i).c);
      Pars.SetError(i, _calibTable.at(i).cerr);
      if(_calibTable.at(i).fixed==1) {
        Pars.Fix(i);
      }
    }
    
  }
  
  void updateCalibTableFromMnUserParameters(MnUserParameters& Pars)
  {
    // size of the reference MnUserParameters
    int SizeOfMnPars = (int)Pars.Parameters().size();
    // size of the current calibTable
    int SizeOfCalibTable = (int)_calibTable.size();

    // give all the values to MnUserParameters
    for (int i=0; i<SizeOfMnPars; i++)
    {
      _calibTable.at(i).c = Pars.Value(i);
      _calibTable.at(i).cerr = Pars.Error(i);
      if(Pars.Parameter(i).IsFixed())
      {
        _calibTable.at(i).fixed = 1;
      }
      else
      {
        _calibTable.at(i).fixed = 0;
      }
    }
    
  }
  
  void recalcM(const int method, const std::vector<double>& par) const 
  {
  
    // event print inteval
    int n_interval=1;
    // time
    double time0 = (double)clock();
    double time1 = (double)clock();
    
    std::map<int, calibRecord>::iterator itiz;
    // Reproduce mass from raw data, applying calibration (if requested)
    for (Int_t i=0; i<_nEvents; i++) {
      // print
      if(_debug>3 && i%n_interval==0) 
      {
        std::cout << "recalcM: Events " << i << " dtime = " << clock()-time1<< std::endl;
         time1 = clock();
      }
      // Energy
      double E1 = *(_E1.at(i));
      double E2 = *(_E2.at(i));
      if (method==0) 
      {
        //// Note here, we temporarily only want to calibrate E2, all E1 are good GsfElectrons,
        //// we simply take it's energy for E1.
        //// recalib E2
        //for (int ii=0; ii<*(_nHits2.at(i)); ii++)
        //{
        //  double Eold = _HitE2.at(i)->at(ii);
        //  double Enew = calibCell(Eold, _HitIX2.at(i)->at(ii),
        //                           _HitIY2.at(i)->at(ii), 
        //                           _HitIZ2.at(i)->at(ii), par, method);
        //  // Attention:
        //  // Below is a protection assuming all cell energy are positive.
        //  // This is not true in some cases.
        //  // One reason is calibCell return a value of -100.0 if I don't want to fit this cell.
        //  // i.e. not in the calibTable.
        //  // And I only want to calib those that in the calibTable.
        //  // Do nothing if not in the calibTable.
        //  if (Enew>0.0){
        //    //calculate new Energy as :
        //    // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
        //    E2 = E2 - Eold + Enew;
        //  }
        //}
        
        
        if(_debug>3 && i%n_interval==0) 
        {
          std::cout << "recalcM: Using method 0 " << " dtime = " << clock()-time1 <<  std::endl;
           time1 = clock();
        }
        if(_debug>3 && i%n_interval==0) 
        {
          std::cout << "recalcM: Recalib E2 " << " dtime = " << clock()-time1 << std::endl;
          time1 = clock();
        }
        // recalib E2
        for (int ii=0; ii<*(_nHits2.at(i)); ii++)
        {          
          // index in _calibTable
          int index = GetCalibTableIndex(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTable " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // calib const
          double calibC = GetCalibConstFromCalibTableRef(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTableRef " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // coninue the loop if cannot find the cell in the _calibTableRef
          if (calibC<-90.0) continue;
          // if found in _calibTable, use the one in MINUIT par
          if (index>-90)
          {
            calibC = par.at(index); 
          }
          // Enew = Eold * calibC
          double Eold = _HitE2.at(i)->at(ii);
          double Enew = Eold * calibC;
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E2 = E2 - Eold + Enew;
          
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Applied calib. " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
        }
        
      }
      else if (method==1)
      {
        // This method, we calibrate both E1 and E2
        // recalib E1
        for (int ii=0; ii<*(_nHits1.at(i)); ii++){
          int index = GetCalibTableIndex(_HitIX1.at(i)->at(ii), _HitIY1.at(i)->at(ii), _HitIZ1.at(i)->at(ii));
          // coninue the loop if cannot find the cell in the _calibTable
          if (index==-100) continue;
          // Enew = Eold * calibC
          double Eold = _HitE1.at(i)->at(ii);
          double Enew = Eold * (par.at(index));
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E1 = E1 - Eold + Enew;
        }
        // recalib E2
        for (int ii=0; ii<*(_nHits2.at(i)); ii++){
          int index = GetCalibTableIndex(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          // coninue the loop if cannot find the cell in the _calibTable
          if (index==-100) continue;
          // Enew = Eold * calibC
          double Eold = _HitE2.at(i)->at(ii);
          double Enew = Eold * (par.at(index));
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E2 = E2 - Eold + Enew;
        }
      }
      else if (method==2)
      {
        if(_debug>3 && i%n_interval==0) 
        {
          std::cout << "recalcM: Using method 2 " << " dtime = " << clock()-time1 <<  std::endl;
           time1 = clock();
        }
        // This method, we also calibrate both E1 and E2
        // But for the cells not to be calibrated (not free parameters),
        //   we use calibration constants from a separate file.
        if(_debug>3 && i%n_interval==0) 
        {
          std::cout << "recalcM: Recalib E1 " << " dtime = " << clock()-time1 << std::endl;
           time1 = clock();
        }
        // recalib E1
        for (int ii=0; ii<*(_nHits1.at(i)); ii++)
        {
          // index in _calibTable
          int index = GetCalibTableIndex(_HitIX1.at(i)->at(ii), _HitIY1.at(i)->at(ii), _HitIZ1.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTable " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // calib const
          double calibC = GetCalibConstFromCalibTableRef(_HitIX1.at(i)->at(ii), _HitIY1.at(i)->at(ii), _HitIZ1.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTableRef " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // coninue the loop if cannot find the cell in the _calibTableRef
          if (calibC<-90.0) continue;
          // if found in _calibTable, use the one in MINUIT par
          if (index>-90)
          {
            calibC = par.at(index); 
          }
          // Enew = Eold * calibC
          double Eold = _HitE1.at(i)->at(ii);
          double Enew = Eold * calibC;
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E1 = E1 - Eold + Enew;
          
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Applied calib. " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
        }
                
        if(_debug>3 && i%n_interval==0) 
        {
          std::cout << "recalcM: Recalib E2 " << " dtime = " << clock()-time1 << std::endl;
          time1 = clock();
        }
        // recalib E2
        for (int ii=0; ii<*(_nHits2.at(i)); ii++)
        {          
          // index in _calibTable
          int index = GetCalibTableIndex(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTable " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // calib const
          double calibC = GetCalibConstFromCalibTableRef(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTableRef " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // coninue the loop if cannot find the cell in the _calibTableRef
          if (calibC<-90.0) continue;
          // if found in _calibTable, use the one in MINUIT par
          if (index>-90)
          {
            calibC = par.at(index); 
          }
          // Enew = Eold * calibC
          double Eold = _HitE2.at(i)->at(ii);
          double Enew = Eold * calibC;
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E2 = E2 - Eold + Enew;
          
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Applied calib. " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
        }
      }
      else if (method==21)
      {
        if(_debug>3 && i%n_interval==0)
        {
          std::cout << "recalcM: Using method 21 " << " dtime = " << clock()-time1 <<  std::endl;
           time1 = clock();
        }
        // This method similar to method2, we calibrate both E1 and E2.
        // For the cells not to be calibrated (not free parameters),
        // we use calibration constants from a separate file.
        // But more than that, a new feature is introduced to remove the correlation bias.
        // each time Enew = Eold * (Ci + Fi*<dC>)
        // Ci is the unbiased calibration constant,
        //   Fi = 1/fi - 1 = (Eold - Ehit)/Ehit,
        //   <dC> is the mean bias on average of all the other crystals.
        // In par, par[0] will be the <dC>, and starting form (including) par[1], are the
        //   calibration constants, par[1] stores the Ci was stored by par[0] previously.

        // also added is the method in method 23, that it is using energy regression.


        // regression energy
        double EReg1 = *(_EReg1.at(i));
        double EReg2 = *(_EReg2.at(i));
        double RegScale1 = EReg1/E1;
        double RegScale2 = EReg2/E2;

        if(_debug>3 && i%n_interval==0)
        {
          std::cout << "recalcM: Recalib E1 " << " dtime = " << clock()-time1 << std::endl;
           time1 = clock();
        }
        // dE1 = Sum(-Eold + Enew)
        double dE1(0);
        // recalib E1
        for (int ii=0; ii<*(_nHits1.at(i)); ii++)
        {
          // index in _calibTable
          int index = GetCalibTableIndex(_HitIX1.at(i)->at(ii), _HitIY1.at(i)->at(ii), _HitIZ1.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Got index _calibTable " << " index = " << index << ", dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // calib const
          double calibC = GetCalibConstFromCalibTableRef(_HitIX1.at(i)->at(ii), _HitIY1.at(i)->at(ii), _HitIZ1.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Got CalibC from _calibTableRef " << "calibC = " << calibC << ", dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // coninue the loop if cannot find the cell in the _calibTableRef
          if (calibC<-90.0) continue;
          // by default, if the hit is not the one to be fitted,
          //  calculate Enew = Eold * calibC
          double Eold = _HitE1.at(i)->at(ii);
          double Enew = Eold * calibC;

          // if found in _calibTable, use the one in MINUIT par
          //  and calculate Enew = Eold * (CalibC + Fi*<dC>)
          if (index>-90)
          {
            calibC = par.at(index+1);
            // <dC> = par.at(0)
            // Fi = (Etotal/Eold -1)
            //Etotal = (*(_E1.at(i)))
            Enew = Eold * ( calibC + ( E1/Eold-1.0 ) * par.at(0) );
          }

          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          //E1 = E1 - Eold + Enew;
          dE1 = dE1 - Eold + Enew;

          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Applied calib. " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
        }

        if(_debug>3 && i%n_interval==0)
        {
          std::cout << "recalcM: Recalib E2 " << " dtime = " << clock()-time1 << std::endl;
          time1 = clock();
        }

        // dE2 = Sum(-Eold + Enew)
        double dE2(0);
        // recalib E2
        for (int ii=0; ii<*(_nHits2.at(i)); ii++)
        {
          // index in _calibTable
          int index = GetCalibTableIndex(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Got index _calibTable " << " index = " << index << ", dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // calib const
          double calibC = GetCalibConstFromCalibTableRef(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Got CalibC from _calibTableRef " << "calibC = " << calibC << ", dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // coninue the loop if cannot find the cell in the _calibTableRef
          if (calibC<-90.0) continue;


          // by default, if the hit is not the one to be fitted,
          //  calculate Enew = Eold * calibC
          double Eold = _HitE2.at(i)->at(ii);
          double Enew = Eold * calibC;

          // if found in _calibTable, use the one in MINUIT par
          //  and calculate Enew = Eold * (CalibC + Fi*<dC>)
          if (index>-90)
          {
            calibC = par.at(index+1);
            // <dC> = par.at(0)
            // Fi = (Etotal/Eold -1)
            //Etotal = (*(_E2.at(i)))
            Enew = Eold * ( calibC + ( E2/Eold-1.0 ) * par.at(0) );
          }

          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          //E2 = E2 - Eold + Enew;
          dE2 = dE2 - Eold + Enew;

          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Applied calib. " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
        }

        // recalculate E1 and E2
        E1 = E1 + dE1;
        E2 = E2 + dE2;

        // re-apply regression energy scale
        E1 = E1*RegScale1;
        E2 = E2*RegScale2;

      }
      else if (method==23)
      {
        if(_debug>3 && i%n_interval==0)
        {
          std::cout << "recalcM: Using method 23 " << " dtime = " << clock()-time1 <<  std::endl;
           time1 = clock();
        }
        // Method 23 is the same as method2, calibrate both E1 and E2
        // For the cells not to be calibrated (not free parameters),
        //   we use calibration constants from a separate file.
        // but the method 23 use the Energy regression for the electron energy

        // regression energy
        double EReg1 = *(_EReg1.at(i));
        double EReg2 = *(_EReg2.at(i));
        double RegScale1 = EReg1/E1;
        double RegScale2 = EReg2/E2;

        if(_debug>3 && i%n_interval==0) 
        {
          std::cout << "recalcM: Recalib E1 " << " dtime = " << clock()-time1 << std::endl;
           time1 = clock();
        }
        // recalib E1
        for (int ii=0; ii<*(_nHits1.at(i)); ii++)
        {
          // index in _calibTable
          int index = GetCalibTableIndex(_HitIX1.at(i)->at(ii), _HitIY1.at(i)->at(ii), _HitIZ1.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Got index _calibTable " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // calib const
          double calibC = GetCalibConstFromCalibTableRef(_HitIX1.at(i)->at(ii), _HitIY1.at(i)->at(ii), _HitIZ1.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Got index _calibTableRef " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // if found in _calibTable, use the one in MINUIT par
          if (index>-90)
          {
            calibC = par.at(index);
          }
          // if cannot be found in the _calibTable and also cannot be found in the _calibTableRef I skip this hit.
          else if (calibC<-90.0)
          {
            continue;
          }
          // Enew = Eold * calibC
          double Eold = _HitE1.at(i)->at(ii);
          double Enew = Eold * calibC;
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E1 = E1 - Eold + Enew;

          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Applied calib. " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
        }

        if(_debug>3 && i%n_interval==0)
        {
          std::cout << "recalcM: Recalib E2 " << " dtime = " << clock()-time1 << std::endl;
          time1 = clock();
        }
        // recalib E2
        for (int ii=0; ii<*(_nHits2.at(i)); ii++)
        {
          // index in _calibTable
          int index = GetCalibTableIndex(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Got index _calibTable " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // calib const
          double calibC = GetCalibConstFromCalibTableRef(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Got index _calibTableRef " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // if found in _calibTable, use the one in MINUIT par
          if (index>-90)
          {
            calibC = par.at(index);
          }
          // if cannot be found in the _calibTable and also cannot be found in the _calibTableRef I skip this hit.
          else if (calibC<-90.0)
          {
            continue;
          }
          // Enew = Eold * calibC
          double Eold = _HitE2.at(i)->at(ii);
          double Enew = Eold * calibC;
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          E2 = E2 - Eold + Enew;

          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Applied calib. " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
        }

        // re-apply regression energy scale
        E1 = E1*RegScale1;
        E2 = E2*RegScale2;

      }
////////////////////
      else if (method==3) 
      {
        // This method, we also calibrate both E1 or E2, may or may not both
        // If an electron has the cells to be calibrated, it recalculate its energy,
        // if not, we use the GsfElectron to provide its energy.
        // In case of recalculation, and the calibration constants are not in the _calibTable,
        // we use calibration constants from a separate file as stored in _calibTableRef.
        
        if(_debug>3 && i%n_interval==0) 
        {
          std::cout << "recalcM: Using method 3 " << " dtime = " << clock()-time1 <<  std::endl;
           time1 = clock();
        }
        
        // take E1  or not
        bool takeE1(false);
        // dE1 = Sum(-Eold + Enew)
        double dE1(0);
        
        if(_debug>3 && i%n_interval==0) 
        {
          std::cout << "recalcM: Recalib E1 " << " dtime = " << clock()-time1 << std::endl;
           time1 = clock();
        }
        // recalib E1
        for (int ii=0; ii<*(_nHits1.at(i)); ii++)
        {
          // index in _calibTable
          int index = GetCalibTableIndex(_HitIX1.at(i)->at(ii), _HitIY1.at(i)->at(ii), _HitIZ1.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTable " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // calib const
          double calibC = GetCalibConstFromCalibTableRef(_HitIX1.at(i)->at(ii), _HitIY1.at(i)->at(ii), _HitIZ1.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTableRef " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // coninue the loop if cannot find the cell in the _calibTableRef
          if (calibC<-90.0) continue;
          // if found in _calibTable, use the one in MINUIT par
          if (index>-90)
          {
            calibC = par.at(index);
            takeE1 = true; 
          }
          // Enew = Eold * calibC
          double Eold = _HitE1.at(i)->at(ii);
          double Enew = Eold * calibC;
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          //E1 = E1 - Eold + Enew;
          dE1 = dE1 - Eold + Enew;
          
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Applied calib. " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
        }
        
        // check if take E1
        if (takeE1) E1 = E1 + dE1;
                
                
        // take E2  or not
        bool takeE2(false);
        // dE2 = Sum(-Eold + Enew)
        double dE2(0);
        
        if(_debug>3 && i%n_interval==0) 
        {
          std::cout << "recalcM: Recalib E2 " << " dtime = " << clock()-time1 << std::endl;
          time1 = clock();
        }
        // recalib E2
        for (int ii=0; ii<*(_nHits2.at(i)); ii++)
        {          
          // index in _calibTable
          int index = GetCalibTableIndex(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTable " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // calib const
          double calibC = GetCalibConstFromCalibTableRef(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTableRef " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // coninue the loop if cannot find the cell in the _calibTableRef
          if (calibC<-90.0) continue;
          // if found in _calibTable, use the one in MINUIT par
          if (index>-90)
          {
            calibC = par.at(index); 
            takeE2 = true;
          }
          // Enew = Eold * calibC
          double Eold = _HitE2.at(i)->at(ii);
          double Enew = Eold * calibC;
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          //E2 = E2 - Eold + Enew;
          dE2 = dE2 - Eold + Enew;
          
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Applied calib. " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          
          
        }
        
        // check if take E2
        if (takeE2) E2 = E2 + dE2;
        
      }
      else if (method==5)
      {
        if(_debug>3 && i%n_interval==0) 
        {
          std::cout << "recalcM: Using method 5 " << " dtime = " << clock()-time1 <<  std::endl;
           time1 = clock();
        }
        // This method, same as method 2, just add a weight

        if(_debug>3 && i%n_interval==0) 
        {
          std::cout << "recalcM: Recalib E1 " << " dtime = " << clock()-time1 << std::endl;
           time1 = clock();
        }
        
        // recalib E1
        // E1 = E1 + dE1
        // dE1 = Sum(-Eold + Enew)
        double dE1(0);
        // Sum of energies of cells to be fitted
        double ESumOfCalibCells1(0); 
        for (int ii=0; ii<*(_nHits1.at(i)); ii++)
        {
          // old energy of the hit
          double Eold = _HitE1.at(i)->at(ii);
          // index in _calibTable
          int index = GetCalibTableIndex(_HitIX1.at(i)->at(ii), _HitIY1.at(i)->at(ii), _HitIZ1.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTable " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // calib const
          double calibC = GetCalibConstFromCalibTableRef(_HitIX1.at(i)->at(ii), _HitIY1.at(i)->at(ii), _HitIZ1.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTableRef " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // coninue the loop if cannot find the cell in the _calibTableRef
          if (calibC<-90.0) continue;
          // if found in _calibTable, use the one in MINUIT par
          if (index>-90)
          {
            calibC = par.at(index); 
            ESumOfCalibCells1 += Eold; 
          }
          // Enew = Eold * calibC
          double Enew = Eold * calibC;
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          dE1 = dE1 - Eold + Enew;
          
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Applied calib. " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
        }
      
                
        if(_debug>3 && i%n_interval==0) 
        {
          std::cout << "recalcM: Recalib E2 " << " dtime = " << clock()-time1 << std::endl;
          time1 = clock();
        }
        
        // recalib E2
        // E2 = E2 + dE2
        // dE2 = Sum(-Eold + Enew)
        double dE2(0);
        // Sum of energies of cells to be fitted
        double ESumOfCalibCells2(0); 
        
        for (int ii=0; ii<*(_nHits2.at(i)); ii++)
        {          
          // old energy of the hit
          double Eold = _HitE2.at(i)->at(ii);
          // index in _calibTable
          int index = GetCalibTableIndex(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTable " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // calib const
          double calibC = GetCalibConstFromCalibTableRef(_HitIX2.at(i)->at(ii), _HitIY2.at(i)->at(ii), _HitIZ2.at(i)->at(ii));
          if(_debug>3 && i%n_interval==0) 
          {
            std::cout << "recalcM: Got index _calibTableRef " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
          // coninue the loop if cannot find the cell in the _calibTableRef
          if (calibC<-90.0) continue;
          // if found in _calibTable, use the one in MINUIT par
          if (index>-90)
          {
            calibC = par.at(index);
            ESumOfCalibCells2 += Eold; 
          }
          // Enew = Eold * calibC
          double Enew = Eold * calibC;
          //calculate new Energy as :
          // Ecorr = E - Sum(Eraw_i) + Sum(Eraw_i * C_i)
          dE2 = dE2 - Eold + Enew;
          
          if(_debug>3 && i%n_interval==0)
          {
            std::cout << "recalcM: Applied calib. " << " dtime = " << clock()-time1 << std::endl;
            time1 = clock();
          }
        }
        
        // calculate weight of the event (while the mass will be calculated later independent of option 5)
        // define:
        // Energy fraction = [ E1(hits)/E1(total) + E2(hits)/E2(total) ] / 2.0 
        //   note: divide by 2 is to get it in range 0 to 1
        // weight = 1/fraction
        //   note, because we are minimizing -ln(L), a large fraction should be turn 
        //    into a small weight of make a stronger decision in the minimization.
        _MassRecalcWeight.at(i) = 2.0/(ESumOfCalibCells1/E1 + ESumOfCalibCells2/E2);
        // note above the E1 and E2 are the old energy before recalibration.
        
        // re-calculate E1 
        E1 = E1 + dE1;
        
        // re-calculate E2 
        E2 = E2 + dE2;
        
      }
      else if (method==6)
      {
        // Method 6 is not to recalib, but to determine the energy scale at the electron level,
        //  i.e. not at the cell level.

        // regression energy
        double EReg1 = *(_EReg1.at(i));
        double EReg2 = *(_EReg2.at(i));

        // apply energy scale
        if (_UseEle1.at(i)) EReg1 *= par.at(0);
        if (_UseEle2.at(i)) EReg2 *= par.at(0);

        // energy after energy scale
        E1 = EReg1;
        E2 = EReg2;

      }
      else if (method==7)
      {
        // Method 7 is not to recalib, but to determine the energy scale at the electron level,
        //  i.e. not at the cell level.
        // different from method==6, method 7 tries to fit eta scales of serveral parameters at once.

        // regression energy
        double EReg1 = *(_EReg1.at(i));
        double EReg2 = *(_EReg2.at(i));

        // apply energy scale
        if(_ScaleBin1.at(i)>=0) EReg1 *= par.at(_ScaleBin1.at(i));
        if(_ScaleBin2.at(i)>=0) EReg2 *= par.at(_ScaleBin2.at(i));

        //std::cout << "method 7: _ScaleBin1.at(i) = " << _ScaleBin1.at(i) 
        //             << "; _ScaleBin2.at(i) = " << _ScaleBin2.at(i) << std::endl;
        // energy after energy scale
        E1 = EReg1;
        E2 = EReg2;

      }
      else if (method==8)
      {
        // Method 8 is to fit the IC by applying the IC of the Seed crystal to the rawEnergy of the whole SC.
        //   Each time fit only one crystal in the following steps:
        //   1.) select events with electrons having there Seed hits in the crystal of the IC to be fitted.
        //   2.) the minimization varies the IC and recalculates the electron's energy by applying the IC of
        //        the cyrstal being fitted to the raw energy of the whole SC (many crystals/hits).

        // regression energy
        double EReg1 = *(_EReg1.at(i));
        double EReg2 = *(_EReg2.at(i));
        // raw ecal energy
        double RawEEcal1 = *(_RawEEcal1.at(i));
        double RawEEcal2 = *(_RawEEcal2.at(i));
        // 
        // apply energy scale
        if (_UseEle1.at(i)) EReg1 = (EReg1/E1)*(RawEEcal1*par.at(0)+(E1-RawEEcal1));
        if (_UseEle2.at(i)) EReg2 = (EReg2/E2)*(RawEEcal2*par.at(0)+(E2-RawEEcal2));

        // energy after energy scale
        E1 = EReg1;
        E2 = EReg2;

      }

      if(_debug>3 && i%n_interval==0) 
      {
        std::cout << "recalcM: Recalculate Mee" << " dtime = " << clock()-time1 << std::endl;
        time1 = clock();
      }
      // calculate mass
      _MassRecalc.at(i) = CalcMass(E1, *(_Eta1.at(i)), *(_Phi1.at(i)), E2, *(_Eta2.at(i)), *(_Phi2.at(i)));
      
      if(_debug>3 && i%n_interval==0) 
      {
        std::cout << "recalcM: Recalculate Mee done." << " dtime = " << clock()-time1 << std::endl;
        time1 = clock();
      }
    }
    
    // if method==5, weighted -ln(L) = Sum{ - w_i * ln(P_i) }, we need to modify the weight 
    // to avoid the weight being too large.
    // modify the weight to be w'_i = N * w_i / Sum(w_i)
    if (method==5)
    {
       // get sum of the weights
       double sum_weight(0);
       for (int i=0; i<(int)_MassRecalcWeight.size(); i++)
       {
         sum_weight += _MassRecalcWeight.at(i);
       }
       
       // recalculate weight to be w'_i = N * w_i / Sum(w_i)
       for (int i=0;  i<(int)_MassRecalcWeight.size(); i++)
       {
         _MassRecalcWeight.at(i) = (int)_MassRecalcWeight.size() *
                                  _MassRecalcWeight.at(i) / sum_weight;
       }
    }
  }
  
  double calibCell(double E, int ix, int iy, int iz, const std::vector<double>& par, int method=0) const
  {
    int index = GetCalibTableIndex(ix, iy, iz);
    // abort if cannot find the cell in the _calibTable
    if (index==-100) {
      if (_debug>2)
        std::cout << "WARNING:: Cannot find this cell (ix, iy, iz) = ("
                << ix << "," << iy << "," << iz 
                << ") in the _calibTable." << std::endl;
      // abort();
      return -100.0; // if so, do not calib, and return -100
    }
    // if good, return the calibrated energy
    return E*par.at(index);
  }
  
  double GetCalibConstFromCalibTable(int ix, int iy, int iz) const
  {
    int index = GetCalibTableIndex(ix, iy, iz);
    double CalibC(1.0);
    if (index>-100) 
    { 
      CalibC = _calibTable.at(index).c;
    }
    return CalibC;
  }
  
  double GetCalibConstErrorFromCalibTable(int ix, int iy, int iz) const
  {
    int index = GetCalibTableIndex(ix, iy, iz);
    double CalibCErr(0.1);
    if (index>-100) 
    { 
      CalibCErr = _calibTable.at(index).cerr;
    }
    return CalibCErr;
  }

  double GetCalibConstFromCalibTableRef(const int& ix, const int& iy, const int& iz) const
  {
    std::map<int, std::map<int, std::map<int, calibRecord> > >::const_iterator _it_ix;
    std::map<int, std::map<int, calibRecord> >::const_iterator _it_iy;
    std::map<int, calibRecord>::const_iterator _it_iz;
    // check ix
    _it_ix = _calibTableRef.find(ix);
    if (_it_ix == _calibTableRef.end()) return -1000.0 ;
    // check iy
    _it_iy = (_it_ix->second).find(iy);
    if (_it_iy == (_it_ix->second).end()) return -1000.0;
    // check iz
    _it_iz = (_it_iy->second).find(iz);
    if (_it_iz == (_it_iy->second).end()) return -1000.0;
    // 
    return _it_iz->second.c ;

  }
  
  double GetCalibConstErrorFromCalibTableRef(const int& ix, const int& iy, const int& iz) const
  {
    std::map<int, std::map<int, std::map<int, calibRecord> > >::const_iterator _it_ix;
    std::map<int, std::map<int, calibRecord> >::const_iterator _it_iy;
    std::map<int, calibRecord>::const_iterator _it_iz;
    // check ix
    _it_ix = _calibTableRef.find(ix);
    if (_it_ix == _calibTableRef.end()) return -1000.0 ;
    // check iy
    _it_iy = _it_ix->second.find(iy);
    if (_it_iy == _it_ix->second.end()) return -1000.0;
    // check iz
    _it_iz = _it_iy->second.find(iz);
    if (_it_iz == _it_iy->second.end()) return -1000.0;
    // 
    return _it_iz->second.cerr ;

  }  
  
  int GetCalibTableIndex(int ix, int iy, int iz) const
  {
    // find the index from the _calibTable
    for(int i=0; i<(int)_calibTable.size(); i++)
    {
      if (ix==_calibTable.at(i).ix && 
          iy==_calibTable.at(i).iy &&
          iz==_calibTable.at(i).iz)
      {
        return i;
      }
    }
    // if not find anything in the _calibTable, return -100, means cannot find
    return -100;
  }
  
  double GetVoigtNorm(double windowLow, double windowHigh, double voigtMass,
                      double voigtResolution, double voigtWidth)
  {
    double norm(0);
    double mass(0);
    const int NSteps = 1000000;
    _voigt.setParams(voigtMass, voigtResolution, voigtWidth);

    for (int i=0; i<NSteps; i++)
    {
      mass = windowLow + ((double)i+0.5)*(windowHigh-windowLow)/(double)NSteps;
      norm += (windowHigh-windowLow)/(double)NSteps*(_voigt.val(mass));
    }
    
    std::cout << "GetVoigtNorm:: Norm = " << norm << std::endl;
    return norm;
  }
  
  TComplex ComplexErrFunc(const TComplex& Z) const
  {
    // Return CERNlib complex error function
    //
    // This code is translated from the fortran version in the CERN mathlib.
    // (see ftp://asisftp.cern.ch/cernlib/share/pro/src/mathlib/gen/c/cwerf64.F)
    
    TComplex ZH,S,T,V;
    static TComplex R[38];
    
    static const Double_t Z1= 1, HF= Z1/2, Z10= 10;
    static const Double_t C1= 74/Z10, C2= 83/Z10, C3= Z10/32, C4 = 16/Z10;
    static const Double_t C = 1.12837916709551257, P = pow(2*C4,33);
    static const TComplex zero(0);
    
    Double_t X(Z.Re());
    Double_t Y(Z.Im());
    Double_t XA(fabs(X));
    Double_t YA(fabs(Y));
    int N;
    if((YA < C1) && (XA < C2))
    {
      ZH= TComplex(YA+C4,XA);
      R[37]=zero;
      N= 36;
      while(N > 0) {
        TComplex dummy=TComplex::Conjugate(R[N+1]);
        dummy*=N;
        T=ZH+dummy;
        R[N--]=(T*HF)/T.Rho2();
      }
      Double_t XL=P;
      S=zero;
      N= 33;
      while(N > 0) {
        XL=C3*XL;
        S=R[N--]*(S+XL);
      }
      V=S*C;
    }
    else
    {
      ZH=TComplex(YA,XA);
      R[1]=zero;
      N= 9;
      while(N > 0) {
        TComplex dummy=TComplex::Conjugate(R[1]);
        dummy*=N;
        T=ZH+dummy;
        R[1]=(T*HF)/T.Rho2();
        N--;
      }
      V=R[1]*C;
    }
    if(YA==0) V=TComplex(exp(-(XA*XA)),V.Im());
    
    if(Y < 0)
    {
      TComplex tmp(XA,YA);
      tmp= -tmp*tmp;
      TComplex dummy=TComplex::Exp(tmp);
      dummy*=2;
      V=dummy-V;
      if(X > 0) V= TComplex::Conjugate(V);
    }
    else
    {
      if(X < 0) V= TComplex::Conjugate(V);
    }
    return V;
    
  }
  
  double Voigtian(double x, double mean, double sigma, double width) const
  {
    // Evaluate Voigtian
    double s = (sigma>0.0) ? sigma : -sigma;
    double w = (width>0.0) ? width : -width;
    double coef = -0.5/(s*s);
    double arg = x - mean;
    
    // Return constant for zero width and sigma
    if (s==0.0&&w==0.0) return 1.0;
    
    // Gauss for zero width
    if (w==0.0) return exp(coef*arg*arg);
    
    // Actual Voigtian function for non-trivial width and sigma
    double c = 1.0/(sqrt(2.0)*s);
    double a = 0.5*c*w;
    double u = c*arg;
    TComplex z(u,a);
    TComplex v(0.0);
    
    v = ComplexErrFunc(z);
    
    return c/sqrt(M_PI)*v.Re();
  }
  
  double CalcMass(double E1, double Eta1, double Phi1,
                  double E2, double Eta2, double Phi2) const 
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

  
  TH1D* GetMeeTH1DFromOriginalCalibration(const char* histname = "hist")
  {
    // new a histogram
    TH1D* hist = new TH1D(histname, histname, 60, 60,120);
    hist->Sumw2();
    hist->GetXaxis()->SetTitle("M(ee) (GeV)");
    hist->GetYaxis()->SetTitle("nEvents");
    hist->SetMarkerStyle(20);
    hist->SetTitle(histname);
    
    // fill the hist
    for (int i=0; i<(int)_Mass.size(); i++)
    {
      hist->Fill(_Mass.at(i));
    }
    
    return hist;
  }
  
  TH1D* GetMeeTH1DFromMnUserParameters(MnUserParameters& Pars, const char* histname = "hist")
  {
    // recalculate Mass
    recalcM(_method, Pars.Params());
    
    // new a histogram
    TH1D* hist = new TH1D(histname, histname, 60, 60,120);
    hist->Sumw2();
    hist->GetXaxis()->SetTitle("M(ee) (GeV)");
    hist->GetYaxis()->SetTitle("nEvents");
    hist->SetMarkerStyle(20);
    hist->SetTitle(histname);
    
    // fill the hist
    for (int i=0; i<(int)_MassRecalc.size(); i++)
    {
      hist->Fill(_MassRecalc.at(i));
    }
    
    return hist;
  }
  
private:
  
  // method
  int _method;
  
  // debug
  int _debug;
  
  // parameters of BWGS function
  double _windowHigh;
  double _windowLow;
  double _voigtMass;
  double _voigtResolution;
  double _voigtWidth;
  double _voigtFraction;
  double _voigtNorm;

  // voigtian function
  voigtian _voigt; 

  // data
  int _nEvents;
  int _nSignals;
  std::vector<double*> _E1;
  std::vector<double*> _EReg1;
  std::vector<double*> _Eta1;
  std::vector<double*> _Phi1;
  std::vector<bool>_UseEle1;
  std::vector<int*>_nHits1;
  std::vector< std::vector<double>* > _HitE1;
  std::vector< std::vector<int>* > _HitIX1;
  std::vector< std::vector<int>* > _HitIY1;
  std::vector< std::vector<int>* > _HitIZ1;
  std::vector<double*> _E2;
  std::vector<double*> _EReg2;
  std::vector<double*> _Eta2;
  std::vector<double*> _Phi2;
  std::vector<bool>_UseEle2;
  std::vector<int*>_nHits2;
  std::vector< std::vector<double>* > _HitE2;
  std::vector< std::vector<int>* > _HitIX2;
  std::vector< std::vector<int>* > _HitIY2;
  std::vector< std::vector<int>* > _HitIZ2;
 
  std::vector<int> _ScaleBin1;
  std::vector<int> _ScaleBin2;
  
  std::vector<double*> _RawEEcal1; 
  std::vector<double*> _RawEEcal2;
 
  std::vector<double> _Mass;
  mutable std::vector<double> _MassRecalc;
  mutable std::vector<double> _MassRecalcWeight;
  
  // calibration constants
  std::vector<calibRecord> _calibTable;
  // reference calibration constants for _method 2
  std::map<int, std::map<int, std::map<int, calibRecord> > > _calibTableRef;

};
    
  } // namespace Minuit2
  
} // namespace ROOT


