#ifndef CDATA_H
#define CDATA_H

#include <vector>
#include <string>
#include <map>
#include <complex>
#include <iostream>
#include "fwd.h"

using namespace std;

struct particle_data
{
  // independent info
  int type;
  double px, py, pz, p0;

  // derived info
  double eta, pt, phi, pcos, psin;
};

struct sample
{
    vector<particle_data> in_list, fin_list;
    int id;
};

class cdata
{
    private:
        const VarMap vm;
	void read(void);
	void calc(void);
	vector<sample> event;
    public:
	void write2file(void);        
	cdata(const VarMap& vm_);
	~cdata(){};
};

#endif
