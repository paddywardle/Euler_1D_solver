#include <vector>
#include <iostream>
#include <array>
#include "header.H"

void initial_conditions(std::vector<std::array<double, 3>>& u, double x0, double dx)
{ 
  for (int i=0; i<u.size(); i++)
    {
      double x = x0 + i * dx;
      
      if (x <= 0.5)
	{
	  // density
	  u[i][0] = 1.0;
	  // velocity
	  u[i][1] = 0.0;
	  // pressure
	  u[i][2] = 1.0;
	}
      else
	{
	  // density
	  u[i][0] = 0.125;
	  // velocity
	  u[i][1] = 0.0;
	  // pressure
	  u[i][2] = 0.1;
	}
    }
}
