#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>

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

calibRecord CalibRecord(int idx, int ix, int iy, int iz, double c, double cerr, int fixed, int nfits)
{
  calibRecord record = {idx, ix, iy, iz, c, cerr, fixed, nfits};
  return record;
}

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
