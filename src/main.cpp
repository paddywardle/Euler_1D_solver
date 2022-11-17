#include <vector>
#include <array>
#include <iostream>
#include "header.H"
int main()
{
  double x0 = 0.0;
  double x1 = 1.0;
  double tStart = 0.0;
  double tStop = 0.25;
  double gamma = 1.4;
  double C = 0.8;
  int nCells = 1000;
  double dx = (x1 - x0)/nCells;
  

  std::vector<std::array<double, 3>> u;
  u.resize(nCells+2);
  std::vector<std::array<double, 3>> uPlus1;
  uPlus1.resize(nCells+2);

  initial_conditions(u, x0, dx);
  initial_conditions(uPlus1, x0, dx);

  u = prim_to_cons(u, gamma, nCells);
  uPlus1 = prim_to_cons(uPlus1, gamma, nCells);

  solvers(u, uPlus1, tStart, tStop, dx, nCells, C, gamma);

  u = cons_to_prim(u, gamma, nCells);

  outputFile("../data/test_LF.dat", u, x0, dx);
}
