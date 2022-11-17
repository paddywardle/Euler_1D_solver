#include <vector>
#include <string>
#include <array>

void initial_conditions(std::vector<std::array<double, 3>>&, double, double);

std::vector<std::array<double, 3>> prim_to_cons(std::vector<std::array<double, 3>>, double, int);

std::vector<std::array<double, 3>> cons_to_prim(std::vector<std::array<double, 3>>, double, int);

void solvers(std::vector<std::array<double, 3>>&, std::vector<std::array<double, 3>>&, double, double, double, int, double, double);

void outputFile(std::string, std::vector<std::array<double, 3>>&, double, double);

/*
void solver(std::vector<double>&, std::vector<double>&, double, double, double, double, int, double);

void outputFile(std::string, std::vector<double>&, int, double, double);

class FORCE{
public:

  FORCE(int, double, double, double, double, double, double, std::vector<std::vector<double>, std::vector<double>, std::vector<double>>, std::vector<std::vector<double>, std::vector<double>, std::vector<double>>);
  
  void initial_conditions(std::vector<std::array<double, 10>>& u, double x0, double dx);
  
  void solvers(std::vector<std::vector<double>, std::vector<double>, std::vector<double>>&, std::vector<std::vector<double>, std::vector<double>, std::vector<double>>&, double, double, double, double, int, double);

  void outputFile(std::string, std::vector<std::vector<double>, std::vector<double>, std::vector<double>>&, int, double, double);

private:

  int nCells;

  double x0;

  double x1;

  double tStart;

  double tStop;

  double a;

  double C;

  std::vector<double> u;

  std::vector<double> uPlus1;
}
*/
