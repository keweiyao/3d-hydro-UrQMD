#include "cdata.h"
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>
#include "fwd.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void echo_configure(const VarMap& vm)
{
	cout << "Echo configure:" << endl;
  for (const auto& it : vm) {
    cout << it.first.c_str() << " ";
    auto& value = it.second.value();
    if (auto v = boost::any_cast<int>(&value))
      cout << *v << endl;
    else if (auto v = boost::any_cast<bool>(&value))
      cout << *v << endl;
    else if (auto v = boost::any_cast<double>(&value))
      cout << *v << endl;
    else if (auto v = boost::any_cast<fs::path>(&value))
      cout << *v << endl;
    else
      cout << "error" << endl;
  }
}

int main(int argc, char * argv[])
{
  fs::path fpath=argv[1];

  po::options_description config_file_opts{};
  config_file_opts.add_options()
    ("inputfile", po::value<fs::path>(), "URQMD input file")
    ("outputfile", po::value<std::string>(), "analyse output file")
    ("mult", po::value<bool>(), "switch for multiplicity calc")
    ("m_eta_L", po::value<double>(), "lower eta cut for multiplicity")
    ("m_eta_H", po::value<double>(), "higher eta cut for multiplicity")
    ("m_pt_L", po::value<double>(), "lower pt cut for multiplicity")
    ("m_pt_H", po::value<double>(), "higher pt cut for multiplicity")
    ("dndy", po::value<bool>(), "switch for dNdy calc")
    ("dn_eta_L", po::value<double>(), "lower eta cut for dNdy")
    ("dn_eta_H", po::value<double>(), "higher eta cut for dNdy")
    ("dn_eta_bins", po::value<int>(), "eta bins for dNdy")
    ("dn_pt_L", po::value<double>(), "lower pt cut for dNdy")
    ("dn_pt_H", po::value<double>(), "higher pt cut for dNdy")
    ("qns", po::value<bool>(), "switch for Q1, Q2, Q3, Q4, Q6 calc")
    ("q_eta_L", po::value<double>(), "lower eta cut for qn")
    ("q_eta_H", po::value<double>(), "higher eta cut for qn")
    ("q_eta_bins", po::value<int>(), "eta bins for qn")
    ("q_pt_L", po::value<double>(), "lower pt cut for qn")
    ("q_pt_H", po::value<double>(), "higher pt cut for qn");
  VarMap var_map{};
  fs::ifstream conf_file{fpath};
  po::store(po::parse_config_file(conf_file, config_file_opts), var_map);
  echo_configure(var_map);
 

  cdata ds(var_map);
  ds.write2file();
  return 0;
}
