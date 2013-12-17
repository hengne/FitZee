#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
  if (argc<4)
  {
    std::cout << argv[0] << "<calibTable_in.dat> <calibTable_out.dat> <scale>"
              << std::endl;
    return 0;
  }


  const char* filename = argv[1];
  const char* filename_out = argv[2];
  const double scale = atof(argv[3]);
  
  std::cout << "Read file : " << filename << std::endl;
  std::cout << "Write modifled file to : " << filename_out << std::endl;
  std::cout << "Apply scale : " << scale << std::endl;

  int idx, ix, iy, iz, fixed, nfits;
  double c, cerr;

  std::ifstream myfile;
  myfile.open(filename);

  std::ofstream myfile_out;
  myfile_out.open(filename_out);

  std::string line;
  
  if (myfile.is_open()&&myfile_out.is_open())
  {
    int index=0;
    while (getline(myfile,line))
    {
      std::stringstream sline(line);
      sline >> idx
      >> ix
      >> iy
      >> iz
      >> c
      >> cerr
      >> fixed
      >> nfits ;
     
      myfile_out << index << " " 
                   << ix << " "
                   << iy << " "
                   << iz << " "
                   << c*scale << " "
                   << cerr*scale << " "
                   << fixed << " "
                   << nfits << " "
                   << std::endl;
      index++;
    }
    myfile.close();
    myfile_out.close();
  }
}
