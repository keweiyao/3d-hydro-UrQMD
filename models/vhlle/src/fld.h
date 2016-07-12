#include <iosfwd>
#include "cll.h"

class EoS ;
class TransportCoeff ;
class Cornelius ;

class Fluid
{
private:
    EoS *eos;
    TransportCoeff *trcoeff ;
    Cornelius *cornelius ;
    Cell *cell ;
    Cell *cell0 ;
    int nx, ny, nz ;
    double minx, maxx, miny, maxy, minz, maxz ;
    double dx, dy, dz, dt ;
    double tcrit ;
  	double ecrit ;
  	double vEff, EtotSurf ;
	std::ofstream foutkw, foutkw_dim, foutxvisc, foutyvisc, foutdiagvisc,
	foutx, fouty, foutdiag, foutz, fout_aniz, fout2d, ffreeze, fvhlle_hydro, fstep, fnormal;
  	int compress2dOut ;
public:
	Fluid(EoS *_eos, TransportCoeff *_trcoeff, int _nx, int _ny, int _nz,
   	double _minx, double _maxx, double _miny, double _maxy, double _minz, double _maxz, double dt,  double eCrit);
	~Fluid();
	void initOutput(char *dir, int maxstep, double tau0, int cmpr2dOut) ;
	int output_xy_spacing ;
	int output_eta_points ;
	int output_tau_spacing ;
	int output_eta_spacing ;

	int output_nt, output_nx, output_ny ;

	inline int getNX () {return nx ; }
	inline int getNY () {return ny ; }
	inline int getNZ () {return nz ; }
	inline double getDx () {return dx ; }
	inline double getDy () {return dy ; }
	inline double getDz () {return dz ; }
	inline double getX(int ix) { return minx + ix*dx ; }
	inline double getY(int iy) { return miny + iy*dy ; }
	inline double getZ(int iz) { return minz + iz*dz ; }

	void getCMFvariables(Cell *c, double tau, double &e, double &nb, double &nq, double &ns, double &vx, double &vy, double &Y) ;
	void getCMFvariables_new(Cell *c, double tau, double &e, double &nb, double &nq, double &ns, double &vx, double &vy, double &vz) ;

  	inline Cell* getCell(int ix, int iy, int iz) 
	{ 
		ix = ix>0 ? ix : 0; ix = ix<nx ? ix : nx-1 ;
    	iy = iy>0 ? iy : 0; iy = iy<ny ? iy : ny-1 ;
    	iz = iz>0 ? iz : 0; iz = iz<nz ? iz : nz-1 ;
		return &cell[ix+nx*iy+nx*ny*iz] ; 
	}

	void correctImagCells(void) ; // only ideal hydro part
	void correctImagCellsFull(void) ; // correct ideal+visc
	void updateM(double tau, double dt) ;

	void outputPDirections(double tau) ;
  	void outputGnuplot(double tau) ;
	void outputSurface(double tau) ;
  	void outputCorona(double tau) ;
};
