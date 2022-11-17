#include <string>
#include <vector>
#include <fstream>
#include "header.H"
#include <array>

void outputFile(std::string outputName, std::vector<std::array<double, 3>>& u, double x0, double dx)
{
  std::ofstream output(outputName);
  
  for (int i=1; i<u.size(); i++)
    {
      double x = x0 + (i-1) * dx;
      output<<x<<" ";
      for (int j=0; j<u[i].size(); j++)
	{
	  output<<u[i][j]<<" ";
	}
      output<<std::endl;
    } 
}
