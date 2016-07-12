#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "fld.h"
#include "eos.h"
#include "inc.h"
#include "rmn.h"

#include "mcgstream.h"

using namespace std;

MCGstreamIC::MCGstreamIC(double normalization, char* filename, EoS *_eos, int event_number) {
  eos=_eos;
  e_norm_ = normalization;
  readFile_hdf5(filename, eos, event_number);
}

int MCGstreamIC::findIndex(int ix, int iy, int ieta) 
{
     return (IC_NX_*IC_NY_)*ieta + ix*IC_NY_+iy;
}

double MCGstreamIC::interpolate(double x, double y, double eta,
				std::vector<double> *tvec) {
  int IC_x_half = (IC_NX_%2==0) ? IC_NX_/2 : (IC_NX_-1)/2 ;
  int IC_y_half = (IC_NY_%2==0) ? IC_NY_/2 : (IC_NY_-1)/2 ;
  int IC_eta_half = (IC_NZ_%2==0) ? IC_NZ_/2 : (IC_NZ_-1)/2 ;
 //   eta = eta - 0.465; 
  const double xmin = - IC_x_half * IC_dxy_;
  const double ymin = - IC_y_half * IC_dxy_;
  const double etamin = - IC_eta_half * IC_deta_;
  const int ix = (int)((x - xmin) / IC_dxy_);
  const int iy = (int)((y - ymin) / IC_dxy_);
  const int ieta = (int)((eta - etamin) / IC_deta_);

  if (ix < 0 || ix > IC_NX_ - 1
      || iy < 0 || iy > IC_NY_ - 1
      || ieta < 0 || ieta > IC_NZ_ - 1 ) {
    //    std::cout << ix<<"/"<<IC_NX_ << " " << iy<<"/"<<IC_NY_ << " " << ieta<<"/"<<IC_NZ_  << std::endl;
    return 0.0;
  }
  const double sx = std::abs(x - xmin - ix * IC_dxy_);
  const double sy = std::abs(y - ymin - iy * IC_dxy_);
  const double seta = std::abs(eta - etamin - ieta * IC_deta_);
  double wx[2] = {(1. - sx / IC_dxy_), sx / IC_dxy_};
  double wy[2] = {(1. - sy / IC_dxy_), sy / IC_dxy_};
  double weta[2] = {(1. - seta / IC_deta_), seta / IC_deta_};
  double result = 0.0;
  for (int jx = 0; jx < 2; jx++) {
    for (int jy = 0; jy < 2; jy++) {
      for (int jeta = 0; jeta < 2; jeta++) {
	int i = findIndex(ix + jx, iy + jy, ieta + jeta);
	if (i <= limitindex_) {
	  result += wx[jx] * wy[jy] * weta[jeta] * (*tvec)[i];
	} 
	else {
	  //	  std::cout << "large i=" << i << ", set by zero" << std::endl;
	  return 0.0;
	}
      }
    }
  }
  return result;
}

void MCGstreamIC::setQ(double tau0, double x, double y, double eta, double *Q) 
{
  double t00, t01, t02, t03;
  t00 = interpolate(x, y, eta, &t00vec_);
////YX  t01 = interpolate(x, y, eta, &t01vec_);
////YX  t02 = interpolate(x, y, eta, &t02vec_);
////YX  t03 = interpolate(x, y, eta, &t03vec_);
  Q[0] = t00;
////  Q[1] = t01;
////  Q[2] = t02;
////  Q[3] = t03 * tau0;
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
  Q[5] = 0.0;
  Q[6] = 0.0;
}

void MCGstreamIC::setIC(Fluid *f, double tau) {
  double e = 0., p = 0., nb = 0., nq = 0., ns = 0.;
  double vx = 0., vy = 0., vz = 0.;
  double Q[7];
  Cell *c;
  double avv_num = 0., avv_den = 0.;
  double Etotal = 0.0;
  for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++)
      for (int iz = 0; iz < f->getNZ(); iz++) {
        c = f->getCell(ix, iy, iz);
        double x = f->getX(ix);
        double y = f->getY(iy);
        double eta = f->getZ(iz);
	setQ(tau, x, y, eta, Q);
	transformPV(eos, Q, e, p, nb, nq, ns, vx, vy, vz);
	avv_num += sqrt(vx * vx + vy * vy) * e;
	avv_den += e;
	c->setPrimVar(eos, tau, e, nb, nq, 0., vx, vy, vz);
	double _p = eos->p(e, nb, nq, 0.);
	const double gamma2 = 1.0 / (1.0 - vx * vx - vy * vy - vz * vz);
	Etotal +=
	  ((e + _p) * gamma2 * (cosh(eta) + vz * sinh(eta)) - _p * cosh(eta));
	c->saveQprev();
	if (e > 0.) c->setAllM(1.);

      }
  std::cout << "average initial flow = " << avv_num / avv_den << std::endl;
////YX  std::cout << "total energy = "
////YX	    << Etotal *f->getDx() * f->getDy() * f->getDz() * tau
////YX	    << std::endl;

}

void MCGstreamIC::readFile(char* filename, EoS *eos) {
  std::ifstream inputfile(filename);
  std::ofstream ftest_initial("test_initial.dat", ios::app);
  if (!inputfile) {
    std::cout << "Error! File " << filename << " not found." << std::endl;
    exit(1);
  }
  int line = 0, errors = 0;
  int ix, iy, ieta;
  double t00, t01, t02, t03, check;
  double sd, ed, nb;
  while (!inputfile.eof()) {
    inputfile >> ix >> iy >> ieta >> sd;
    nb = 0.0;
    IC_NX_=101;
    IC_NY_=101;
    IC_NZ_=5;
    IC_dxy_ = 0.1;
    IC_deta_ = 0.2;
    sd *= e_norm_;
    if(sd < 1.0e-8) {ed = 0.0;}
    else {eos->gete(sd, ed, nb);}
   
    t00=ed;
    t01=0.;
    t02=0.;
    t03=0.;
    check=0;
    //ftest_initial << ix<<'\t'  << iy<<'\t'  <<ieta <<'\t' << ed <<'\t' <<sd <<'\t'<<nb<< std::endl; 
    if ((int)check == 0) {
      t00vec_.push_back(t00);
      t01vec_.push_back(t01);
      t02vec_.push_back(t02);
      t03vec_.push_back(t03);

    } else {
      std::cout << "Warning: Line " << line << " in file " << filename << std::endl;
      std::cout << "Nonzero check value: " << check << std::endl;
      errors++;
      if (errors > 100) {
	std::cout << "Too many failed checks, exiting." << std::endl;
	exit(1);
      }
    }
    line++;
    if( (line % 1000) == 0)
    {std::cout << "reading initial lines" << line<<std::endl;}
  }
  std::cout <<"read initial successful :)"<< std::endl;
}

#include "H5Cpp.h"
using namespace H5;
void MCGstreamIC::readFile_hdf5(char* filename, EoS *eos, int event_number) 
{
  limitindex_ = -1;
  const H5std_string FILE_NAME(filename);
  char buffer [50];
  sprintf(buffer, "event_3d_%d", event_number);
  std::cout << "event_number" << event_number << " " << buffer << std::endl ;  
  const H5std_string DATASET_NAME(buffer);
    
  H5File file(FILE_NAME, H5F_ACC_RDONLY);
  DataSet dataset = file.openDataSet(DATASET_NAME);

  DataSpace data_space = dataset.getSpace();
  int const rank = data_space.getSimpleExtentNdims();
  hsize_t dims_out[rank];
  data_space.getSimpleExtentDims(dims_out, NULL);

  const int NX = dims_out[0];
  const int NY = dims_out[1];
  const int NZ = 1;
  
  double data_out[100][100][1] = {0.0};
  hsize_t dimsm[rank];
  dimsm[0] = NX;
  dimsm[1] = NY;
  dimsm[2] = NZ;
  IC_NX_ = dims_out[0];
  IC_NY_ = dims_out[1];
  IC_NZ_ = dims_out[2];
  IC_dxy_ = 0.2;
  IC_deta_ = 0.46875;
	
  DataSpace mem_space(rank,dimsm);

  hsize_t offset[3]; offset[0] = 0; offset[1] = 0; offset[2] = 0;
  hsize_t count[3]; count[0] = NX; count[1] = NY; count[2] = NZ;
  hsize_t offset_scan[3]; offset_scan[0] = 0; offset_scan[1] = 0; offset_scan[2] = 0;
  
  double t00, t01, t02, t03, check;
  double sd, ed, nb;
  double s_tot = 0;
  for(int ieta = 0; ieta<dims_out[2];ieta++)
  {
    if(ieta<dims_out[2])
      {
	offset_scan[2] = ieta;
	data_space.selectHyperslab(H5S_SELECT_SET, count, offset_scan);
	mem_space.selectHyperslab(H5S_SELECT_SET, count, offset);
	dataset.read(data_out, PredType::NATIVE_DOUBLE, mem_space, data_space);
      }
    for (int ix = 0; ix<NX; ix++)
    {
      for (int iy=0; iy<NY; iy++)
      {
	s_tot += data_out[ix][iy][0];
	sd = e_norm_*data_out[ix][iy][0];
	nb = 0.0;
	if(sd < 1.0e-8) {ed = 0.0;}
	else {eos->gete(sd, ed, nb);}
   
	t00=ed;
	t01=0.;
	t02=0.;
	t03=0.;
	
	t00vec_.push_back(t00);
	t01vec_.push_back(t01);
	t02vec_.push_back(t02);
	t03vec_.push_back(t03);
	limitindex_++;
      }
    }
    //std::cout << "eta slice " << ieta << "done; IC cells # = " << limitindex_+1 <<std::endl;
  }
  std::cout <<"IC read finished, s_tot = " << s_tot << std::endl;
}

