#include <vector>
#include <array>
#include <cmath>
#include "header.H"

std::vector<std::array<double, 3>> prim_to_cons(std::vector<std::array<double, 3>> u_prim, double gamma, int nCells)
{

  std::vector<std::array<double, 3>> u_cons;
  u_cons.resize(nCells+2);
  
  for (int i=0; i<u_prim.size(); i++)
    {
      u_cons[i][0] = u_prim[i][0];
      u_cons[i][1] = u_prim[i][0] * u_prim[i][1];
      u_cons[i][2] = u_prim[i][2]/(gamma-1.0) + 0.5 * u_prim[i][0] * pow(u_prim[i][1], 2.0);
    }

  return u_cons;
}

std::vector<std::array<double, 3>> cons_to_prim(std::vector<std::array<double, 3>> u_cons, double gamma, int nCells)
{
  std::vector<std::array<double, 3>> u_prim;
  u_prim.resize(nCells+2);

  for (int i=0; i<u_cons.size(); i++)
    {
      u_prim[i][0] = u_cons[i][0];
      u_prim[i][1] = u_cons[i][1] / u_cons[i][0];
      u_prim[i][2] = (u_cons[i][2] - 0.5 * u_prim[i][0]*pow(u_prim[i][1],2.0)) * (gamma - 1.0);
    }

  return u_prim;
}
