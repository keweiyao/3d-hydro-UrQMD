#include <vector>

// initial state from Monte Carlo Glauber
// plus pre-equilibrium streaming
class MCGstreamIC {
 public:
  MCGstreamIC(double normalization, char* filename, EoS *_eos, int event_number = 0);
  void setIC(Fluid *f, double tau);
 private:
  EoS *eos;
  double e_norm_;
  int IC_NX_, IC_NY_, IC_NZ_; // IC grid size
  int limitindex_;
  double IC_dxy_, IC_deta_; // IC grid step
  std::vector<double> t00vec_;
  std::vector<double> t01vec_;
  std::vector<double> t02vec_;
  std::vector<double> t03vec_;
  std::ofstream ftest_initial;
  void readFile(char* filename, EoS *eos);
  void readFile_hdf5(char* filename, EoS *eos, int event_number);
  int findIndex(int ix, int iy, int ieta);
  double interpolate(double x, double y, double eta,
		     std::vector<double> *tvec);
  void setQ(double tau, double x, double y, double eta, double *Q);
};
