#ifndef FWL
#define FWL
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

namespace boost{
	namespace filesystem {}
	namespace program_options{
		class variables_map;
	}
}
namespace fs = boost::filesystem;
namespace po = boost::program_options;
using VarMap = po::variables_map;
#endif
