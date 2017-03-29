#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "fld.h"
#include "eos.h"
#include "inc.h"
#include "rmn.h"

#include "IC_reader.h"

using namespace std;

IC_reader::IC_reader(char* filename, EoS *_eos, int ic_nxy, int ic_neta, double ic_dxy, double ic_deta)
: IC_NX_(ic_nxy),
  IC_NY_(ic_nxy),
  IC_NZ_(ic_neta),
  IC_dxy_(ic_dxy),
  IC_dz_(ic_deta),
  IC_xmin_(-0.5*(ic_nxy - 1)*ic_dxy),
  IC_ymin_(-0.5*(ic_nxy - 1)*ic_dxy),
  IC_zmin_(-0.5*(ic_neta)*ic_deta)
{
  eos = _eos;
  readFile_hdf5(filename, eos);
}

int IC_reader::findIndex(int ix, int iy, int ieta) 
{
     return (IC_NX_*IC_NY_)*ieta + ix*IC_NY_ + iy;
}

double IC_reader::interpolate(double x, double y, double eta, int row) {
  int ix = std::floor((x - IC_xmin_) / IC_dxy_);
  int iy = std::floor((y - IC_ymin_) / IC_dxy_);
  int iz = std::floor((eta - IC_zmin_) / IC_dz_);

  if (ix < 0 || ix >= IC_NX_ - 1
      || iy < 0 || iy >= IC_NY_ - 1
      || iz < 0 || iz >= IC_NZ_ - 1 ) 
  { return 0.0; }
  
  double sx = x - IC_xmin_ - ix * IC_dxy_;
  double sy = y - IC_ymin_ - iy * IC_dxy_;
  double sz = eta - IC_zmin_ - iz * IC_dz_;

  double wx[2] =   { (1. - sx/IC_dxy_),    sx/IC_dxy_  };
  double wy[2] =   { (1. - sy/IC_dxy_),    sy/IC_dxy_  };
  double wz[2] =   { (1. - sz/IC_dz_),     sz/IC_dz_   };
  double result = 0.0;
  for (int jx = 0; jx < 2; jx++) {
    for (int jy = 0; jy < 2; jy++) {
      for (int jz = 0; jz < 2; jz++) {
	int i = findIndex(ix + jx, iy + jy, iz + jz);
        if (row == 0)
	  result += wx[jx] * wy[jy] * wz[jz] * t00vec_[i];
	if (row == 1)
          result += wx[jx] * wy[jy] * wz[jz] * t01vec_[i];
	if (row == 2)
          result += wx[jx] * wy[jy] * wz[jz] * t02vec_[i];
	if (row == 3)
          result += wx[jx] * wy[jy] * wz[jz] * t03vec_[i];
      }
    }
  }
  return result;
}

void IC_reader::setQ(double tau0, double x, double y, double eta, double *Q) 
{
  double t00, t01, t02, t03;
  t00 = interpolate(x, y, eta, 0);
  t01 = interpolate(x, y, eta, 1);
  t02 = interpolate(x, y, eta, 2);
  t03 = interpolate(x, y, eta, 3);
  Q[0] = t00;
  Q[1] = t01;
  Q[2] = t02;
  Q[3] = t03 * tau0;
  Q[4] = 0.0;
  Q[5] = 0.0;
  Q[6] = 0.0;
}

void IC_reader::setIC(Fluid *f, double tau) {
  double e = 0., p = 0., nb = 0., nq = 0., ns = 0.;
  double vx = 0., vy = 0., vz = 0.;
  double Q[7];
  Cell *c;
  ofstream file("check.txt");
  double avv_num = 0., avv_den = 0.;
  double Etotal = 0.0;
  for (int ix = 0; ix < f->getNX(); ix++){
    for (int iy = 0; iy < f->getNY(); iy++){
      for (int iz = 0; iz < f->getNZ(); iz++){
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
	if (iz == 20){
		double temp, h1, h2, h3, h4;
		eos->eos(e, nb, nq, 0.0, temp, h1, h2, h3, h4);
		file << temp << " ";
	}
	const double gamma2 = 1.0 / (1.0 - vx * vx - vy * vy - vz * vz);
	Etotal +=
	  ((e + _p) * gamma2 * (cosh(eta) + vz * sinh(eta)) - _p * cosh(eta));
	c->saveQprev();
	if (e > 0.) c->setAllM(1.);
      }
    }
  }
  std::cout << "average initial flow = " << avv_num / avv_den << std::endl;
  std::cout << "total energy = "
	    << Etotal *f->getDx() * f->getDy() * f->getDz() * tau
	    << std::endl;
}

void IC_reader::readFile(char* filename, EoS *eos) {
  std::ifstream inputfile(filename);
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
void IC_reader::readFile_hdf5(char* filename, EoS *eos) 
{
  const H5std_string FILE_NAME(filename);
  char buffer [50];
  sprintf(buffer, "event_0");
  std::cout << "dataset name = " << buffer << std::endl ;  
  const H5std_string DATASET_NAME(buffer);
    
  H5File file(FILE_NAME, H5F_ACC_RDONLY);
  DataSet dataset = file.openDataSet(DATASET_NAME);

  DataSpace data_space = dataset.getSpace();
  int const rank = data_space.getSimpleExtentNdims();
  hsize_t dims_out[rank];
  data_space.getSimpleExtentDims(dims_out, NULL);

  if (dims_out[0] != IC_NX_ || dims_out[1] != IC_NY_ || dims_out[2] != IC_NZ_)
    {
      std::cout << "Wrong IC dims" << std::endl;
    }
  
  double data_out[IC_NX_][IC_NY_][1] = {0.0};
  hsize_t dimsm[rank];
  dimsm[0] = IC_NX_;
  dimsm[1] = IC_NY_;
  dimsm[2] = 1;
	
  DataSpace mem_space(rank,dimsm);

  hsize_t offset[3]; offset[0] = 0; offset[1] = 0; offset[2] = 0;
  hsize_t count[3]; count[0] = IC_NX_; count[1] = IC_NY_; count[2] = 1;
  hsize_t offset_scan[3]; offset_scan[0] = 0; offset_scan[1] = 0; offset_scan[2] = 0;
  
  double t00, t01, t02, t03;
  double sd, ed, nb;
  for(int ieta = 0; ieta<dims_out[2]; ieta++)
  {
    if(ieta<dims_out[2])
      {
	offset_scan[2] = ieta;
	data_space.selectHyperslab(H5S_SELECT_SET, count, offset_scan);
	mem_space.selectHyperslab(H5S_SELECT_SET, count, offset);
	dataset.read(data_out, PredType::NATIVE_DOUBLE, mem_space, data_space);
      }
    
    for (int ix = 0; ix<dims_out[0]; ix++)
    {
      for (int iy=0; iy<dims_out[1]; iy++)
      {
	sd = data_out[ix][iy][0];
	nb = 0.0;
	if(sd < 1.0e-8) {ed = 0.0;}
	else {eos->gete(sd, ed, nb);}
	//	std::cout << "check "<< ed << std::endl;
	t00=ed;
	t01=0.;
	t02=0.;
	t03=0.;
	
	t00vec_.push_back(t00);
	t01vec_.push_back(t01);
	t02vec_.push_back(t02);
	t03vec_.push_back(t03);
      }
    }
  }
}

