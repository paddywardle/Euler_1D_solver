#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "header.H"

double flux_fn(double u)
{
  return 0.5 * pow(u, 2.0);
}

double flux_fn_rho(double rho, double v)
{
  return rho * v;
}

double flux_fn_mom(double rho, double v, double p)
{
  return rho * pow(v, 2.0) + p;
}

double flux_fn_E(double E, double v, double p)
{
  return (E + p) * v;
} 

double lax_friedrich_flux(std::array<double, 3> u_i, std::array<double, 3>u_iPlus1, std::array<double, 3>u_prim_i, std::array<double, 3> u_prim_iPlus1, int col, double dt, double dx)
{
  double fhalf;
  
  if (col == 0)
    {
      fhalf = 0.5 * (dx / dt) * (u_i[0] - u_iPlus1[0]) + 0.5 * (flux_fn_rho(u_iPlus1[0], u_prim_iPlus1[1]) + flux_fn_rho(u_i[0], u_prim_i[2]));
    }
  else if (col == 1)
    {
      fhalf = 0.5 * (dx / dt) * (u_i[1] - u_iPlus1[1]) + 0.5 * (flux_fn_mom(u_iPlus1[0], u_prim_iPlus1[1], u_prim_iPlus1[2]) + flux_fn_mom(u_i[0], u_prim_i[1], u_prim_i[2]));
    }
  else if (col == 2)
    {
      fhalf = 0.5 * (dx / dt) * (u_i[2] - u_iPlus1[2]) + 0.5 * (flux_fn_E(u_iPlus1[2], u_prim_iPlus1[1], u_prim_iPlus1[2]) + flux_fn_E(u_i[2], u_prim_i[1], u_prim_i[2]));
    }

  return fhalf;
}

double richtmyer_flux(std::array<double, 3> u_i, std::array<double, 3>u_iPlus1, std::array<double, 3>u_prim_i, std::array<double, 3> u_prim_iPlus1, int col, double dt, double dx)
{
  double fhalf;

  double uhalf_rho = 0.5 * (u_iPlus1[0] + u_i[0]) - 0.5 * (dt/dx) * (flux_fn_rho(u_iPlus1[0], u_prim_iPlus1[1]) - flux_fn_rho(u_i[0], u_prim_i[2]));
  double uhalf_mom = 0.5 * (u_prim_iPlus1[1] + u_prim_i[1]) - 0.5 * (dt/dx) * (flux_fn_mom(u_iPlus1[0], u_prim_iPlus1[1], u_prim_iPlus1[2]) - flux_fn_mom(u_i[0], u_prim_i[1], u_prim_i[2]));
  double uhalf_E = 0.5 * (u_prim_iPlus1[2] + u_prim_i[2]) - 0.5 * (dt/dx) * (flux_fn_E(u_iPlus1[2], u_prim_iPlus1[1], u_prim_iPlus1[2]) + flux_fn_E(u_i[2], u_prim_i[1], u_prim_i[2])); 
  
  if (col == 0)
    {
      double uhalf_v = uhalf_mom / uhalf_rho;
      fhalf = flux_fn_rho(uhalf_rho, uhalf_v);
    }
  else if (col == 1)
    {
      double uhalf_v = uhalf_mom / uhalf_rho;
      double uhalf_p = uhalf_E - 0.5 * uhalf_rho * pow(uhalf_v, 2.0) * (1.4-1.0);
      fhalf = flux_fn_mom(uhalf_rho, uhalf_v, uhalf_p);
    }
  else if (col == 2)
    {
      double uhalf_v = uhalf_mom / uhalf_rho;
      double uhalf_p = uhalf_E - 0.5 * uhalf_rho * pow(uhalf_v, 2.0) * (1.4-1.0);
      fhalf = flux_fn_E(uhalf_E, uhalf_v, uhalf_p);
    }

  return fhalf;
}

double FORCE_flux(std::array<double, 3> u_i, std::array<double, 3> u_iPlus1, std::array<double, 3>u_prim_i, std::array<double, 3> u_prim_iPlus1, int col, double dt, double dx)
{
  double fhalf = 0.5 * (lax_friedrich_flux(u_i, u_iPlus1, u_prim_i, u_prim_iPlus1, col, dt, dx) + richtmyer_flux(u_i, u_iPlus1, u_prim_i, u_prim_iPlus1, col, dt, dx));

  return fhalf;
}

double calculate_timestep(std::vector<std::array<double, 3>> u, int nCells, double C, double dx, double gamma, int counter)
{
  std::vector<std::array<double, 3>> u_prim;
  u_prim.resize(nCells+2);

  u_prim = cons_to_prim(u, gamma, nCells);
  
  std::vector<double> wave_speed;
  wave_speed.resize(nCells+2);
  
  for (int i=0; i<u.size(); i++)
    {
      double cs = sqrt((gamma*u_prim[i][2])/u_prim[i][0]);
      wave_speed[i] = abs(u_prim[i][1]) + cs;
    }
  
  double a_max = *std::max_element(wave_speed.begin(), wave_speed.end());
  
  double dt = C * (dx / a_max);
  
  return dt;
}

void solvers(std::vector<std::array<double, 3>>& u, std::vector<std::array<double, 3>>& uPlus1, double tStart, double tStop, double dx, int nCells, double C, double gamma)
{
  double t = tStart;
  std::vector<std::array<double, 3>> flux;
  flux.resize(nCells+1);

  std::vector<std::array<double, 3>> u_prim;
  u_prim.resize(nCells+2);
  
  int counter = 0;
  do {
    double dt = calculate_timestep(u, nCells, C, dx, gamma, counter);

    u_prim = cons_to_prim(u, gamma, nCells);

    counter+=1;

    t += dt;

    // add transmissive boundary
    for (int i=0; i<u[0].size(); i++)
      {
	u[0][i] = u[1][i];
	u[nCells+1][i] = u[nCells][i];
      }

    for (int i=0; i<flux.size(); i++)
      {
	int counter = 0;
	for (int j=0; j<flux[i].size(); j++)
	  {
	    flux[i][j] = FORCE_flux(u[i], u[i+1], u_prim[i], u_prim[i+1], j, dt, dx);
	    counter+=1;
	    //std::cout<<counter<<std::endl;
	  }
      }
    
    for (int i = 1; i<nCells+1; i++)
      {
	for (int j=0; j<u[i].size(); j++)
	  {
	    uPlus1[i][j] = u[i][j] - (dt/dx) * (flux[i][j] - flux[i-1][j]);
	  }
      }

    u = uPlus1;
  } while (t < tStop);
  std::cout<<"FUCK OFF YOU PRICK"<<std::endl;
}
