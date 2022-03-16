// Blacklight simulation reader header

#ifndef SIMULATION_READER_H_
#define SIMULATION_READER_H_

// C++ headers
#include <fstream>  // ifstream
#include <string>   // string

// Blacklight headers
#include "../blacklight.hpp"                 // enums
#include "../input_reader/input_reader.hpp"  // InputReader
#include "../utils/array.hpp"                // Array

//--------------------------------------------------------------------------------------------------

// Simulation reader
struct SimulationReader
{
  // Constructors and destructor
  SimulationReader(const InputReader *p_input_reader_);
  SimulationReader(const SimulationReader &source) = delete;
  SimulationReader &operator=(const SimulationReader &source) = delete;
  ~SimulationReader();

  // Pointers to other objects
  const InputReader *p_input_reader;

  // Input data - general
  ModelType model_type;

  // Input data - simulation parameters
  SimulationFormat simulation_format;
  std::string simulation_file;
  bool simulation_multiple;
  int simulation_start;
  int simulation_end;
  Coordinates simulation_coord;
  double simulation_a;
  std::string simulation_kappa_name;

  // Input data - slow-light parameters
  bool slow_light_on;
  int slow_chunk_size;
  double slow_t_start;
  double slow_dt;

  // Input data - plasma parameters
  PlasmaModel plasma_model;

  // Flags for tracking function calls
  bool first_time = true;
  bool first_time_root_object_header = true;

  // Metadata
  std::ifstream data_stream;
  unsigned long int root_object_header_address;
  unsigned long int root_btree_address;
  unsigned long int root_name_heap_address;
  unsigned long int root_data_segment_address;
  std::string metric;
  double metric_a, metric_h;
  std::string *dataset_names;
  int num_dataset_names = 0;
  std::string *variable_names;
  int num_variable_names = 0;
  Array<int> num_variables;
  int ind_hydro;
  int ind_bb;
  int ind_rho, ind_pgas, ind_kappa;
  int ind_uu1, ind_uu2, ind_uu3;
  int ind_bb1, ind_bb2, ind_bb3;
  double adiabatic_gamma;
  int num_arrays;
  int latest_file_number;
  const double extrapolation_tolerance = 1.0;

  // Data
  int n_3_root;
  Array<int> levels;
  Array<int> locations;
  Array<double> x1f, x2f, x3f;
  Array<double> x1v, x2v, x3v;
  Array<double> x2v_alt;
  double *time;
  Array<float> *prim;
  Array<float> prim_transpose;

  // External function
  double Read(int snapshot);

  // Internal functions - simulation_reader.cpp
  std::string FormatFilename(int file_number);
  void VerifyVariablesAthena();
  void VerifyVariablesHarm();

  // Internal functions - simulation_geometry.cpp
  void ConvertCoordinates();
  void ConvertPrimitives(Array<float> &primitives);

  // Internal functions - hdf5_format_structure.cpp
  void ReadHDF5Superblock();
  void ReadHDF5RootGroupSymbolTableEntry();
  unsigned long int ReadHDF5Heap(unsigned long int heap_address);
  void ReadHDF5RootObjectHeader();
  void ReadHDF5FloatAttribute(const char *attribute_name, float *p_val);

  // Internal functions - hdf5_format_metadata.cpp
  unsigned long int ReadHDF5DatasetHeaderAddress(const char *name, unsigned long int btree_address,
      unsigned long int data_segment_address);
  void ReadHDF5DataObjectHeader(unsigned long int data_object_header_address,
      unsigned char **p_datatype_raw, unsigned char **p_dataspace_raw, unsigned char **p_data_raw);
  static void ReadHDF5DataspaceDims(const unsigned char *dataspace_raw, unsigned long int **p_dims,
      int *p_num_dims);

  // Internal functions - hdf5_format_arrays.cpp
  void ReadHDF5StringArray(const char *name, bool allocate, std::string **p_string_array,
      int *p_array_length);
  void ReadHDF5IntArray(const char *name, Array<int> &int_array);
  void ReadHDF5FloatArray(const char *name, Array<float> &float_array);
  void ReadHDF5FloatArray(const char *name, Array<double> &double_array);
  void ReadHDF5DoubleArray(const char *name, Array<double> &double_array);
  static void SetHDF5StringArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, bool allocate,
      std::string **string_array, int *p_array_length);
  static void SetHDF5IntArray(const unsigned char *datatype_raw, const unsigned char *dataspace_raw,
      const unsigned char *data_raw, Array<int> &int_array);
  static void SetHDF5FloatArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, Array<float> &float_array);
  static void SetHDF5FloatArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw,
      Array<double> &double_array);
  static void SetHDF5DoubleArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw,
      Array<double> &double_array);
};

#endif
