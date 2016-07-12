// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "output.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options/variables_map.hpp>

#include "event.h"

// Compile HDF5 output if available.
#ifdef TRENTO_HDF5

#include "H5Cpp.h"

namespace trento {

namespace {

/// Simple functor to write many events to an HDF5 file.
class HDF5Writer {
 public:
  /// Prepare an HDF5 file for writing.
  HDF5Writer(const fs::path& filename);

  /// Write an event.
  void operator()(int num, double impact_param, const Event& event) const;

 private:
  /// Internal storage of the file object.
  H5::H5File file_;
};

// Map types to HDF5 predefined datatypes.
// Must create an explicit template specialization for each used type.
using H5Type = H5::PredType;
template <typename T> const H5Type& h5_type();
template <> const H5Type& h5_type<double>() { return H5Type::NATIVE_DOUBLE; }
template <> const H5Type& h5_type<int>()    { return H5Type::NATIVE_INT; }

// Add a simple scalar attribute to an HDF5 dataset.
template <typename T>
void h5_add_scalar_attr(
    const H5::DataSet& dataset, const std::string& name, const T& value) {
  const auto& datatype = h5_type<T>();
  auto attr = dataset.createAttribute(name, datatype, H5::DataSpace{});
  attr.write(datatype, &value);
}

HDF5Writer::HDF5Writer(const fs::path& filename)
    : file_(filename.string(), H5F_ACC_TRUNC)
{}

void HDF5Writer::operator()(
    int num, double impact_param, const Event& event) const {
  // Prepare arguments for new HDF5 dataset.
    if (!event.is3d())
    {
        // The dataset name is a prefix plus the event number.
        const std::string name{"event_" + std::to_string(num)};

        // Cache a reference to the event grid -- will need it several times.
        const auto& grid = event.reduced_thickness_grid();

        // Define HDF5 datatype and dataspace to match the grid.
        const auto& datatype = h5_type<Event::Grid::element>();
        constexpr auto num_dims = Event::Grid::dimensionality;
        hsize_t shape[num_dims];
        std::copy(grid.shape(), grid.shape() + num_dims, shape);
        H5::DataSpace dataspace{num_dims, shape};

        // Set dataset storage properties.
        H5::DSetCreatPropList proplist{};
        // Set chunk size to the entire grid.  For typical grid sizes (~100x100), this
        // works out to ~80 KiB, which is pretty optimal.  Anyway, it makes logical
        // sense to chunk this way, since chunks must be read contiguously and there's
        // no reason to read a partial grid.
        proplist.setChunk(num_dims, shape);
        // Set gzip compression level.  4 is the default in h5py.
        proplist.setDeflate(4);

        // Create the new dataset and write the grid.
        auto dataset = file_.createDataSet(name, datatype, dataspace, proplist);
        dataset.write(grid.data(), datatype);

        // Write event attributes.
        h5_add_scalar_attr(dataset, "b", impact_param);
        h5_add_scalar_attr(dataset, "npart", event.npart());
        h5_add_scalar_attr(dataset, "mult", event.multiplicity());
        for (const auto& ecc : event.eccentricity())
            h5_add_scalar_attr(dataset, "e" + std::to_string(ecc.first), ecc.second);
    }
    else if (event.out3d())
    {
        std::cout << "3d output" << std::endl;
        // The dataset name is a prefix plus the event number.
        const std::string name{"event_3d_" + std::to_string(num)};
        
        // Cache a reference to the event grid -- will need it several times.
        const auto& grid_3d = event.reduced_thickness_grid_3d();
        
        // Define HDF5 datatype and dataspace to match the grid.
        const auto& datatype_3d = h5_type<Event::Grid_3d::element>();
        constexpr auto num_dims_3d = Event::Grid_3d::dimensionality;
        hsize_t shape_3d[num_dims_3d];
        std::copy(grid_3d.shape(), grid_3d.shape() + num_dims_3d, shape_3d);
        H5::DataSpace dataspace_3d{num_dims_3d, shape_3d};
        
        const unsigned long chunk_num_dims_3d = num_dims_3d;
        hsize_t chunk_shape_3d[num_dims_3d];
        chunk_shape_3d[0] = shape_3d[0];
        chunk_shape_3d[1] = shape_3d[1];
        chunk_shape_3d[2] = shape_3d[2];
        
        // Set dataset storage properties.
        H5::DSetCreatPropList proplist_3d{};
        // Set chunk size to the entire grid.  For typical grid sizes (~100x100), this
        // works out to ~80 KiB, which is pretty optimal.  Anyway, it makes logical
        // sense to chunk this way, since chunks must be read contiguously and there's
        // no reason to read a partial grid.
        proplist_3d.setChunk(num_dims_3d, shape_3d);
        // Set gzip compression level.  4 is the default in h5py.
        proplist_3d.setDeflate(4);
        
        // Create the new dataset and write the grid.
        auto dataset_3d = file_.createDataSet(name, datatype_3d, dataspace_3d, proplist_3d);
        dataset_3d.write(grid_3d.data(), datatype_3d);
        
        // Write event attributes.
        h5_add_scalar_attr(dataset_3d, "b", impact_param);
        h5_add_scalar_attr(dataset_3d, "npart", event.npart());
        h5_add_scalar_attr(dataset_3d, "mult", event.multiplicity());
        for (const auto& ecc : event.eccentricity())
            h5_add_scalar_attr(dataset_3d, "e" + std::to_string(ecc.first), ecc.second);

    }
    else
    {
        std::cout << "reduced 3d output" << std::endl;
        // The dataset name is a prefix plus the event number.
        const std::string name{"event_3d_reduced_" + std::to_string(num)};
        
        // Cache a reference to the event grid -- will need it several times.
        const auto& grid_1d = event.reduced_rapidity_profile();
        
        // Define HDF5 datatype and dataspace to match the grid.
        const auto& datatype_1d = h5_type<Event::Grid_1d::element>();
        constexpr auto num_dims_1d = Event::Grid_1d::dimensionality;
        hsize_t shape_1d[num_dims_1d];
        std::copy(grid_1d.shape(), grid_1d.shape() + num_dims_1d, shape_1d);
        H5::DataSpace dataspace_1d{num_dims_1d, shape_1d};
        
        const unsigned long chunk_num_dims_1d = num_dims_1d;
        hsize_t chunk_shape_1d[num_dims_1d];
        chunk_shape_1d[0] = shape_1d[0];
        
        // Set dataset storage properties.
        H5::DSetCreatPropList proplist_1d{};
        // Set chunk size to the entire grid.  For typical grid sizes (~100x100), this
        // works out to ~80 KiB, which is pretty optimal.  Anyway, it makes logical
        // sense to chunk this way, since chunks must be read contiguously and there's
        // no reason to read a partial grid.
        proplist_1d.setChunk(num_dims_1d, shape_1d);
        // Set gzip compression level.  4 is the default in h5py.
        proplist_1d.setDeflate(4);
        
        // Create the new dataset and write the grid.
        auto dataset_1d = file_.createDataSet(name, datatype_1d, dataspace_1d, proplist_1d);
        dataset_1d.write(grid_1d.data(), datatype_1d);
        
        // Write event attributes.
        h5_add_scalar_attr(dataset_1d, "b", impact_param);
        h5_add_scalar_attr(dataset_1d, "npart", event.npart());
        h5_add_scalar_attr(dataset_1d, "mult", event.multiplicity());
        for (const auto& ecc : event.eccentricity())
            h5_add_scalar_attr(dataset_1d, "e" + std::to_string(ecc.first), ecc.second);
        
        
/////////////
        // The dataset name is a prefix plus the event number.
        const std::string nameEP{"EP_3d_" + std::to_string(num)};
        
        // Cache a reference to the event grid -- will need it several times.
        const auto& grid_EP = event.get_event_plane();
        
        // Define HDF5 datatype and dataspace to match the grid.
        const auto& datatype_EP = h5_type<Event::Grid::element>();
        constexpr auto num_dims_EP = Event::Grid::dimensionality;
        hsize_t shape_EP[num_dims_EP];
        std::copy(grid_EP.shape(), grid_EP.shape() + num_dims_EP, shape_EP);
        H5::DataSpace dataspace_EP{num_dims_EP, shape_EP};
        
        const unsigned long chunk_num_dims_EP = num_dims_EP;
        hsize_t chunk_shape_EP[num_dims_EP];
        chunk_shape_EP[0] = shape_EP[0];
        chunk_shape_EP[1] = shape_EP[1];
        
        // Set dataset storage properties.
        H5::DSetCreatPropList proplist_EP{};
        // Set chunk size to the entire grid.  For typical grid sizes (~100x100), this
        // works out to ~80 KiB, which is pretty optimal.  Anyway, it makes logical
        // sense to chunk this way, since chunks must be read contiguously and there's
        // no reason to read a partial grid.
        proplist_EP.setChunk(num_dims_EP, shape_EP);
        // Set gzip compression level.  4 is the default in h5py.
        proplist_EP.setDeflate(4);
        
        // Create the new dataset and write the grid.
        auto dataset_EP = file_.createDataSet(nameEP, datatype_EP, dataspace_EP, proplist_EP);
        dataset_EP.write(grid_EP.data(), datatype_EP);

    }
}

}  // unnamed namespace

}  // namespace trento

#endif  // TRENTO_HDF5

namespace trento {

namespace {

// These output functions are invoked by the Output class.

void write_stream(std::ostream& os, int width,
    int num, double impact_param, const Event& event) {
  using std::fixed;
  using std::setprecision;
  using std::setw;
  using std::scientific;

  // Write a nicely-formatted line of event properties.
  os << setprecision(10)
     << setw(width)            << num
     << setw(15) << fixed      << impact_param
     << setw(5)                << event.npart()
     << setw(18) << scientific << event.multiplicity()
     << fixed;

  for (const auto& ecc : event.eccentricity())
    os << setw(14)             << ecc.second;

  os << '\n';
}

void write_text_file(const fs::path& output_dir, int width,
    int num, double impact_param, const Event& event) {
  // Open a numbered file in the output directory.
  // Pad the filename with zeros.
  std::ostringstream padded_fname{};
  padded_fname << std::setw(width) << std::setfill('0') << num << ".dat";
  fs::ofstream ofs{output_dir / padded_fname.str()};

  // Write a commented header of event properties as key = value pairs.
  ofs << std::setprecision(10)
      << "# event "   << num                  << '\n'
      << "# b     = " << impact_param         << '\n'
      << "# npart = " << event.npart()        << '\n'
      << "# mult  = " << event.multiplicity() << '\n';

  for (const auto& ecc : event.eccentricity())
    ofs << "# e" << ecc.first << "    = " << ecc.second << '\n';

  // Write IC profile as a block grid.  Use C++ default float format (not
  // fixed-width) so that trailing zeros are omitted.  This significantly
  // increases output speed and saves disk space since many grid elements are
  // zero.
  for (const auto& row : event.reduced_thickness_grid()) {
    auto&& iter = row.begin();
    // Write all row elements except the last with a space delimiter afterwards.
    do {
      ofs << *iter << ' ';
    } while (++iter != --row.end());

    // Write the last element and a linebreak.
    ofs << *iter << '\n';
  }
}

// Determine if a filename is an HDF5 file based on the extension.
bool is_hdf5(const fs::path& path) {
  if (!path.has_extension())
    return false;

  auto hdf5_exts = {".hdf5", ".hdf", ".hd5", ".h5"};
  auto result = std::find(hdf5_exts.begin(), hdf5_exts.end(), path.extension());

  return result != hdf5_exts.end();
}

}  // unnamed namespace

Output::Output(const VarMap& var_map) {
  // Determine the required width (padding) of the event number.  For example if
  // there are 10 events, the numbers are 0-9 and so no padding is necessary.
  // However given 11 events, the numbers are 00-10 with padded 00, 01, ...
  auto nevents = var_map["number-events"].as<int>();
  auto width = static_cast<int>(std::ceil(std::log10(nevents)));

  // Write to stdout unless the quiet option was specified.
  if (!var_map["quiet"].as<bool>()) {
    writers_.emplace_back(
      [width](int num, double impact_param, const Event& event) {
        write_stream(std::cout, width, num, impact_param, event);
      }
    );
  }

  // Possibly write to text or HDF5 files.
  if (var_map.count("output")) {
    const auto& output_path = var_map["output"].as<fs::path>();
    if (is_hdf5(output_path)) {
#ifdef TRENTO_HDF5
      if (fs::exists(output_path))
        throw std::runtime_error{"file '" + output_path.string() +
                                 "' exists, will not overwrite"};
      writers_.emplace_back(HDF5Writer{output_path});
#else
      throw std::runtime_error{"HDF5 output was not compiled"};
#endif  // TRENTO_HDF5
    } else {
      // Text files are all written into a single directory.  Require the
      // directory to be initially empty to avoid accidental overwriting and/or
      // spewing files into an already-used location.  If the directory does not
      // exist, create it.
      if (fs::exists(output_path)) {
        if (!fs::is_empty(output_path)) {
          throw std::runtime_error{"output directory '" + output_path.string() +
                                   "' must be empty"};
        }
      } else {
        fs::create_directories(output_path);
      }
      writers_.emplace_back(
        [output_path, width](int num, double impact_param, const Event& event) {
          write_text_file(output_path, width, num, impact_param, event);
        }
      );
    }
  }
}

}  // namespace trento
