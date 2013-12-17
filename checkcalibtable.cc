#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h> 


int main(int argc, char* argv[])
{
  if (argc<6)
  {
    std::cout << "> check <file_to_be_checked.txt> <c_min> <c_max> <cerr_min> <cerr_max> <file_modifled.txt>"
              << std::endl;
    return 0;
  }


  const char* filename = argv[1];
  const double c_min = atof(argv[2]);
  const double c_max = atof(argv[3]);
  const double cerr_min = atof(argv[4]);
  const double cerr_max = atof(argv[5]);
  const char* filename_out;
  
  if (argc>6) 
  {
    filename_out = argv[6];  
  }

  std::cout << "Read file : " << filename << std::endl;
  std::cout << " calib. consts range: " << c_min << " to " << c_max << std::endl;
  std::cout << " calib. consts error range: " << cerr_min << " to " << cerr_max << std::endl;
  if (argc>6) 
  {
    std::cout << "Write modifled file to : " << filename_out << std::endl;
  }

  int idx, ix, iy, iz, fixed, nfits;
  double c, cerr;
  std::string cstr, cerrstr;

  std::ifstream myfile;
  myfile.open(filename);

  std::string line;

  std::ofstream myfile_out;
  if (argc>6) 
  {
    myfile_out.open(filename_out);
  }
  
  std::cout << "The following lines are not satisfy the cuts:" << std::endl;
  if (myfile.is_open()&&
     ( (argc>6&&myfile_out.is_open())||(argc==6) ) )
  {
    while (getline(myfile,line))
    {
      std::stringstream sline(line);
      sline >> idx
      >> ix
      >> iy
      >> iz
      >> cstr
      >> cerrstr
      >> fixed
      >> nfits ;
      
      c = atof(cstr.c_str());
      cerr = atof(cerrstr.c_str());

      if (cstr=="nan"||cerrstr=="nan" || 
           c>c_max||c<c_min||
          cerr>cerr_max||cerr<cerr_min)
      {
        std::cout << line << std::endl;
        if (argc>6) 
        {
          myfile_out << idx << " " 
                   << ix << " "
                   << iy << " "
                   << iz << " "
                   << "1.0" << " "
                   << "0.1" << " "
                   << fixed << " "
                   << nfits << " "
                   << std::endl;
        }
      }
      else 
      {
        if (argc>6) 
        {
          myfile_out << line << std::endl; 
        }
      }
    }
    myfile.close();
    myfile_out.close();
  }

  return 0;
}
