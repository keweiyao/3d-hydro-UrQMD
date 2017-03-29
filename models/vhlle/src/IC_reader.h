#ifndef IC_READER
#define IC_READER
#include <vector>

class IC_reader{
 public:
  IC_reader(char* filename, EoS *_eos, int IC_Nxy, int IC_Neta, double ic_dxy, double ic_deta);
  void setIC(Fluid *f, double tau);
 private:
  EoS *eos;
  const int IC_NX_, IC_NY_, IC_NZ_; // IC grid size
  const double IC_dxy_, IC_dz_, IC_xmin_, IC_ymin_, IC_zmin_;   // IC grid step
  std::vector<double> t00vec_;
  std::vector<double> t01vec_;
  std::vector<double> t02vec_;
  std::vector<double> t03vec_;

  void readFile(char* filename, EoS *eos);
  void readFile_hdf5(char* filename, EoS *eos);
  int findIndex(int ix, int iy, int ieta);
  double interpolate(double x, double y, double eta, int row);
  void setQ(double tau, double x, double y, double eta, double *Q);
};
#endif
