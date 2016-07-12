#include "cdata.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include "H5Cpp.h"

using namespace std;

typedef struct p_output
{
  double pt;
  double eta;
  double psin;
  double pcos;
  int type;
} p_output;

cdata::cdata(const VarMap & vm_)
:
  vm(vm_)
{
  read();
  cout << "A total of " << event.size() << " samples." << endl;
}

void cdata::calc(void)
{
  
}

void cdata::write2file(void)
{
  particle_data p;
  //define h5 particle type
  const H5std_string mem1("pt"), mem2("eta"), mem3("psin"), mem4("pcos"), mem5("type");
  H5::CompType h5particle( sizeof(p_output) );
  h5particle.insertMember(mem1, HOFFSET(p_output, pt), H5::PredType::NATIVE_DOUBLE);
  h5particle.insertMember(mem2, HOFFSET(p_output, eta), H5::PredType::NATIVE_DOUBLE);
  h5particle.insertMember(mem3, HOFFSET(p_output, psin), H5::PredType::NATIVE_DOUBLE);
  h5particle.insertMember(mem4, HOFFSET(p_output, pcos), H5::PredType::NATIVE_DOUBLE);
  h5particle.insertMember(mem5, HOFFSET(p_output, type), H5::PredType::NATIVE_INT);
  
  //create hdf5 file
  const H5std_string file_name = "data.hdf5";
  H5::H5File* file = new H5::H5File(file_name, H5F_ACC_TRUNC);
  cout << "Output file name = " << file_name << endl;
  
  //loop over samples
  for (int is=0; is < event.size(); is++)
  {
    //define sample name
    char buff[50];
    sprintf(buff, "sample-%d", is);
    const H5std_string dataset_name(buff);
    //define sample dataspace
    hsize_t dim[] = {event[is].fin_list.size()};
    const int RANK = 1;
    cout << "sample multiplicity = " << dim[0] << endl;
    H5::DataSpace space(RANK, dim);
    //create sample dataset
    H5::DataSet* dataset =  new H5::DataSet(file->createDataSet(dataset_name, h5particle, space));
    //load sample into a contagious buffer
    p_output * samplebuff = new p_output[dim[0]]; 
    for(int ip=0; ip<dim[0]; ip++)
    {
      p = event[is].fin_list[ip];
      samplebuff[ip].pt = p.pt;
      samplebuff[ip].eta = p.eta;
      samplebuff[ip].psin = p.psin;
      samplebuff[ip].pcos = p.pcos;
      samplebuff[ip].type = p.type;
    }
    //write this sample
    dataset->write(samplebuff, h5particle);
    delete dataset;
    delete[] samplebuff;
  }
  delete file;
}

void cdata::read(void)
{
    event.clear();
    cout << "read event" << endl;
    // read particle list
    fs::ifstream f_readin (vm["inputfile"].as<fs::path>());
    string line, entry;
    double pabs;
    vector<string> spl;
    const string subsample[] = {"before urqmd", "after urqmd"};
    int sample_index, sub_index;
    while(!(f_readin.eof()))
    {
      // 0 1 2 3 4  5  6  7 8  9   10  11  12
      // x y z t px py pz E id mid chg bar str

        getline(f_readin, line);
        if(line.length() == 0)
        {
            continue;
        }
	istringstream sline(line);
	spl.clear();
	do
	{
	  sline >> entry;
	  spl.push_back(entry);
	}while(sline);
        //cout << "l size = " <<spl.size() << endl;
	if(spl[0] == "#event")
	{
	  sample sample1;
	  event.push_back(sample1);
	  sub_index = -1;
	  sample_index = stoi(spl[1]);
	  event.back().id = sample_index;
	  cout << "read sample # " << event.back().id << endl;
	  
	}
	if(spl[0] == "#cols:")
	{
	  sub_index ++;
	  cout << " --> read sub sample: " << subsample[sub_index] << endl;
	}
        if(spl.size() == 14)
	{
	  if (stoi(spl[10]) == 0) {continue;}
	  particle_data p1;
	  p1.px = stof(spl[4]);
	  p1.py = stof(spl[5]);
	  p1.pz = stof(spl[6]);
	  p1.p0 = stof(spl[7]);
	  p1.type = stoi(spl[8]);
	  p1.pt = sqrt(p1.px*p1.px + p1.py*p1.py);
	  p1.psin = p1.py/p1.pt;
	  p1.pcos = p1.px/p1.pt;
	  pabs = sqrt(p1.pt*p1.pt + p1.pz*p1.pz);
	  p1.eta = 0.5*log((pabs + p1.pz)/(pabs - p1.pz));
	  
	  if (sub_index == 0) event[sample_index].in_list.push_back(p1);
	  else if (sub_index == 1) event[sample_index].fin_list.push_back(p1);
	  else cout << "sub event index wrong ..." << endl;
	}
    }
    f_readin.close();
}



