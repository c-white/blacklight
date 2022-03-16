// Blacklight simulation reader

// C++ headers
#include <cstdio>    // snprintf
#include <fstream>   // ifstream
#include <ios>       // ios_base
#include <optional>  // optional
#include <sstream>   // ostringstream
#include <string>    // stoi, string

// Library headers
#include <omp.h>  // omp_get_wtime

// Blacklight headers
#include "simulation_reader.hpp"
#include "../blacklight.hpp"                 // enums
#include "../input_reader/input_reader.hpp"  // InputReader
#include "../utils/array.hpp"                // Array
#include "../utils/exceptions.hpp"           // BlacklightException, BlacklightWarning

//--------------------------------------------------------------------------------------------------

// Simulation reader constructor
// Inputs:
//   p_input_reader_: pointer to object containing input parameters
// Notes:
//   File is not opened for writing until Read() function is called because the file name might be
//       reformatted after this constructor is called.
SimulationReader::SimulationReader(const InputReader *p_input_reader_)
  : p_input_reader(p_input_reader_)
{
  // Copy general input data
  model_type = p_input_reader->model_type.value();

  // Copy simulation parameters
  if (model_type == ModelType::simulation)
  {
    simulation_format = p_input_reader->simulation_format.value();
    simulation_file = p_input_reader->simulation_file.value();
    simulation_multiple = p_input_reader->simulation_multiple.value();
    if (simulation_multiple)
    {
      simulation_start = p_input_reader->simulation_start.value();
      if (simulation_start < 0)
        throw BlacklightException("Must have nonnegative index simulation_start.");
      simulation_end = p_input_reader->simulation_end.value();
      if (simulation_end < simulation_start)
        throw
            BlacklightException("Must have simulation_end at least as large as simulation_start.");
    }
    simulation_coord = p_input_reader->simulation_coord.value();
    simulation_a = p_input_reader->simulation_a.value();
  }

  // Copy slow-light parameters
  if (model_type == ModelType::simulation)
  {
    slow_light_on = p_input_reader->slow_light_on.value();
    if (slow_light_on)
    {
      if (not simulation_multiple)
        throw BlacklightException("Must enable simulation_multiple to use slow light.");
      slow_chunk_size = p_input_reader->slow_chunk_size.value();
      if (slow_chunk_size < 2)
        throw BlacklightException("Must have slow_chunk_size be at least 2.");
      if (slow_chunk_size > simulation_end - simulation_start + 1)
        throw BlacklightException("Not enough simulation files for given slow_chunk_size.");
      slow_t_start = p_input_reader->slow_t_start.value();
      slow_dt = p_input_reader->slow_dt.value();
      if (slow_dt <= 0.0)
        throw BlacklightException("Must have positive time interval slow_dt.");
    }
  }
  if (model_type != ModelType::simulation and p_input_reader->slow_light_on.has_value()
      and p_input_reader->slow_light_on.value())
    throw BlacklightException("Can only use slow light with simulation data.");

  // Copy plasma parameters
  if (model_type == ModelType::simulation)
  {
    plasma_model = p_input_reader->plasma_model.value();
    if (plasma_model == PlasmaModel::code_kappa)
      simulation_kappa_name = p_input_reader->simulation_kappa_name.value();
  }

  // Determine how many files will be held in memory simultaneously
  num_arrays = 0;
  if (model_type == ModelType::simulation)
    num_arrays = slow_light_on ? slow_chunk_size : 1;

  // Allocate array of time values
  if (num_arrays > 0)
    time = new double[num_arrays];

  // Allocate arrays of Arrays of cell variables
  if (num_arrays > 0)
    prim = new Array<float>[num_arrays];
}

//--------------------------------------------------------------------------------------------------

// Simulation reader destructor
SimulationReader::~SimulationReader()
{
  if (num_dataset_names > 0)
    delete[] dataset_names;
  if (num_variable_names > 0)
    delete[] variable_names;
  if (num_arrays > 0)
  {
    for (int n = 0; n < num_arrays; n++)
      prim[n].Deallocate();
    delete[] time;
    delete[] prim;
  }
}

//--------------------------------------------------------------------------------------------------

// Simulation reader read and initialize function
// Inputs:
//   snapshot: index (starting at 0) of which snapshot is about to be prepared
// Outputs:
//   returned value: execution time in seconds
// Notes:
//   Does nothing if model does not need to be read from file.
//   The output file offset is always equal to snapshot; the input file offset is equal to snapshot
//       if slow_light_on == false.
//   Opens and closes stream for reading.
//   Initializes all member objects.
//   Implements a subset of the HDF5 standard:
//       portal.hdfgroup.org/display/HDF5/File+Format+Specification
double SimulationReader::Read(int snapshot)
{
  // Only proceed if needed
  if (model_type != ModelType::simulation)
    return 0.0;
  double time_start = omp_get_wtime();

  // Prepare default number of files to read
  int num_read = 1;

  // Determine which files to read with slow light
  if (slow_light_on)
  {
    // Calculate time at camera for current snapshot
    double snapshot_time = slow_t_start + slow_dt * snapshot;

    // Initialize most recent file time and number
    double latest_time = snapshot_time - 2.0 * extrapolation_tolerance;
    if (not first_time)
      latest_time = time[0];
    int latest_file_number_old = -1;
    if (first_time)
      latest_file_number = simulation_start + slow_chunk_size - 2;
    else
      latest_file_number_old = latest_file_number;

    // Go through files until sufficiently late time is found
    while (latest_time < snapshot_time and latest_file_number < simulation_end)
    {
      // Determine file name
      latest_file_number++;
      std::string simulation_file_formatted = FormatFilename(latest_file_number);

      // Open input file
      data_stream =
          std::ifstream(simulation_file_formatted, std::ios_base::in | std::ios_base::binary);
      if (not data_stream.is_open())
        throw BlacklightException("Could not open file for reading.");

      // Read basic data about file
      ReadHDF5Superblock();
      root_data_segment_address = ReadHDF5Heap(root_name_heap_address);
      ReadHDF5RootObjectHeader();

      // Read time
      if (simulation_format == SimulationFormat::athena)
      {
        float time_temp;
        ReadHDF5FloatAttribute("Time", &time_temp);
        latest_time = time_temp;
      }
      else if (simulation_format == SimulationFormat::harm)
      {
        Array<double> time_temp(1);
        ReadHDF5DoubleArray("t", time_temp);
        latest_time = time_temp(0);
      }
    }

    // Check range of files covers desired time
    if (latest_time < snapshot_time - extrapolation_tolerance)
    {
      std::ostringstream message;
      message << "Snapshot " << snapshot << " at time " << snapshot_time;
      message << " would require significant extrapolation beyond file " << simulation_end << ".";
      throw BlacklightException(message.str().c_str());
    }
    else if (latest_time < snapshot_time)
    {
      std::ostringstream message;
      message << "Snapshot " << snapshot << " at time " << snapshot_time;
      message << " requires moderate extrapolation.";
      BlacklightWarning(message.str().c_str());
    }

    // Account for no new data needed
    if (latest_file_number == latest_file_number_old)
      num_read = 0;

    // Shift existing data
    else if (latest_file_number - slow_chunk_size + 1 <= latest_file_number_old)
    {
      num_read = latest_file_number - latest_file_number_old;
      for (int n = slow_chunk_size - 1; n >= num_read; n--)
      {
        prim[n].Swap(prim[n-num_read]);
        time[n] = time[n-num_read];
      }
    }

    // Account for all data being replaced
    else
      num_read = slow_chunk_size;
  }

  // Determine which file to read in the case of multiple files without slow light
  else if (simulation_multiple)
    latest_file_number = simulation_start + snapshot;

  // Determine which file to read in the case of a single file
  else
    latest_file_number = -1;

  // Read new files
  for (int n = 0; n < num_read; n++)
  {
    // Determine file name
    std::string simulation_file_formatted = simulation_file;
    if (latest_file_number >= 0)
      simulation_file_formatted = FormatFilename(latest_file_number - n);

    // Open input file
    data_stream =
        std::ifstream(simulation_file_formatted, std::ios_base::in | std::ios_base::binary);
    if (not data_stream.is_open())
      throw BlacklightException("Could not open file for reading.");

    // Read basic data about file
    ReadHDF5Superblock();
    root_data_segment_address = ReadHDF5Heap(root_name_heap_address);
    ReadHDF5RootObjectHeader();

    // Read time
    if (simulation_format == SimulationFormat::athena)
    {
      float time_temp;
      ReadHDF5FloatAttribute("Time", &time_temp);
      time[n] = time_temp;
    }
    else if (simulation_format == SimulationFormat::harm)
    {
      Array<double> time_temp(1);
      ReadHDF5DoubleArray("t", time_temp);
      time[n] = time_temp(0);
    }

    // Read metric
    if (first_time and simulation_format == SimulationFormat::harm)
    {
      std::string *p_temp_metric;
      int temp_count;
      ReadHDF5StringArray("header/metric", true, &p_temp_metric, &temp_count);
      metric = *p_temp_metric;
      delete[] p_temp_metric;
      if (simulation_coord == Coordinates::sph_ks)
      {
        if (metric != "MKS")
        {
          std::ostringstream message;
          message << "Given metric mks does not match file value of " << metric;
          message << "; ignoring the latter.";
          BlacklightWarning(message.str().c_str());
        }
        Array<double> a_temp, h_temp;
        ReadHDF5DoubleArray("header/geom/mks/a", a_temp);
        ReadHDF5DoubleArray("header/geom/mks/hslope", h_temp);
        metric_a = a_temp(0);
        if (metric_a != simulation_a)
        {
          std::ostringstream message;
          message << "Given spin of " << simulation_a << " does not match file value of ";
          message << metric_a << "; ignoring the latter.";
          BlacklightWarning(message.str().c_str());
        }
        metric_h = h_temp(0);
      }
      else
        throw BlacklightException("Invalid simulation_coord for Harm format.");
    }

    // Read block layout
    if (first_time)
    {
      if (simulation_format == SimulationFormat::athena)
      {
        ReadHDF5IntArray("Levels", levels);
        ReadHDF5IntArray("LogicalLocations", locations);
      }
      else if (simulation_format == SimulationFormat::harm)
      {
        levels.Allocate(1);
        levels(0) = 0;
        locations.Allocate(1, 3);
        locations(0,0) = 0;
        locations(0,1) = 0;
        locations(0,2) = 0;
      }
    }

    // Read coordinates
    if (first_time)
    {
      if (simulation_format == SimulationFormat::athena)
      {
        ReadHDF5FloatArray("x1f", x1f);
        ReadHDF5FloatArray("x2f", x2f);
        ReadHDF5FloatArray("x3f", x3f);
        ReadHDF5FloatArray("x1v", x1v);
        ReadHDF5FloatArray("x2v", x2v);
        ReadHDF5FloatArray("x3v", x3v);
      }
      else if (simulation_format == SimulationFormat::harm)
      {
        Array<int> num_cells;
        Array<double> x_start, dx;
        ReadHDF5IntArray("header/n1", num_cells);
        ReadHDF5DoubleArray("header/geom/startx1", x_start);
        ReadHDF5DoubleArray("header/geom/dx1", dx);
        x1f.Allocate(1, num_cells(0) + 1);
        x1v.Allocate(1, num_cells(0));
        for (int i = 0; i < num_cells(0); i++)
        {
          x1f(i+1) = x_start(0) + (i + 1) * dx(0);
          x1v(i) = 0.5 * (x1f(i) + x1f(i+1));
        }
        ReadHDF5IntArray("header/n2", num_cells);
        ReadHDF5DoubleArray("header/geom/startx2", x_start);
        ReadHDF5DoubleArray("header/geom/dx2", dx);
        x2f.Allocate(1, num_cells(0) + 1);
        x2v.Allocate(1, num_cells(0));
        for (int j = 0; j < num_cells(0); j++)
        {
          x2f(j+1) = x_start(0) + (j + 1) * dx(0);
          x2v(j) = 0.5 * (x2f(j) + x2f(j+1));
        }
        ReadHDF5IntArray("header/n3", num_cells);
        ReadHDF5DoubleArray("header/geom/startx3", x_start);
        ReadHDF5DoubleArray("header/geom/dx3", dx);
        x3f.Allocate(1, num_cells(0) + 1);
        x3v.Allocate(1, num_cells(0));
        for (int k = 0; k < num_cells(0); k++)
        {
          x3f(k+1) = x_start(0) + (k + 1) * dx(0);
          x3v(k) = 0.5 * (x3f(k) + x3f(k+1));
        }
        ConvertCoordinates();
      }
    }

    // Read cell data
    if (simulation_format == SimulationFormat::athena)
    {
      VerifyVariablesAthena();
      if (first_time)
      {
        int n5 = num_variables(ind_hydro) + num_variables(ind_bb);
        int n4 = levels.n1;
        int n3 = x3v.n1;
        int n2 = x2v.n1;
        int n1 = x1v.n1;
        prim[n].Allocate(n5, n4, n3, n2, n1);
      }
      Array<float> hydro(prim[n]);
      hydro.Slice(5, 0, num_variables(ind_hydro) - 1);
      ReadHDF5FloatArray("prim", hydro);
      Array<float> bb(prim[n]);
      bb.Slice(5, num_variables(ind_hydro), num_variables(ind_hydro) + num_variables(ind_bb) - 1);
      ReadHDF5FloatArray("B", bb);
    }
    else if (simulation_format == SimulationFormat::harm)
    {
      VerifyVariablesHarm();
      if (first_time)
      {
        int n5 = num_variables(0);
        int n4 = levels.n1;
        int n3 = x3v.n1;
        int n2 = x2v.n1;
        int n1 = x1v.n1;
        prim[n].Allocate(n5, n4, n3, n2, n1);
        prim_transpose.Allocate(n1, n2, n3, n5);
      }
      ReadHDF5FloatArray("prims", prim_transpose);
      for (int n_variable = 0; n_variable < num_variables(0); n_variable++)
        for (int k = 0; k < x3v.n1; k++)
          for (int j = 0; j < x2v.n1; j++)
            for (int i = 0; i < x1v.n1; i++)
              prim[n](n_variable,0,k,j,i) = prim_transpose(i,j,k,n_variable);
      for (int k = 0; k < x3v.n1; k++)
        for (int j = 0; j < x2v.n1; j++)
          for (int i = 0; i < x1v.n1; i++)
            prim[n](ind_pgas,0,k,j,i) *= static_cast<float>(adiabatic_gamma - 1.0);
      ConvertPrimitives(prim[n]);
    }

    // Close input file
    data_stream.close();

    // Update first time flag
    first_time = false;
  }

  // Calculate elapsed time
  return omp_get_wtime() - time_start;
}

//--------------------------------------------------------------------------------------------------

// Function to construct filename formatted with file number
// Inputs:
//   file_number: number of simulation file to construct
// Outputs:
//   returned_value: string containing formatted filename
std::string SimulationReader::FormatFilename(int file_number)
{
  // Locate braces
  std::string::size_type simulation_pos_open = simulation_file.find_first_of('{');
  if (simulation_pos_open == std::string::npos)
    throw BlacklightException("Invalid simulation_file for multiple runs.");
  std::string::size_type simulation_pos_close
      = simulation_file.find_first_of('}', simulation_pos_open);
  if (simulation_pos_close == std::string::npos)
    throw BlacklightException("Invalid simulation_file for multiple runs.");

  // Parse integer format string
  if (simulation_file[simulation_pos_close-1] != 'd')
    throw BlacklightException("Invalid simulation_file for multiple runs.");
  int simulation_field_length = 0;
  if (simulation_pos_close - simulation_pos_open > 2)
    simulation_field_length = std::stoi(simulation_file.substr(simulation_pos_open + 1,
        simulation_pos_close - simulation_pos_open - 2));
  int file_number_length = std::snprintf(nullptr, 0, "%d", file_number);
  if (file_number_length < 0)
    throw BlacklightException("Could not format file name.");
  int num_zeros = 0;
  if (file_number_length < simulation_field_length)
    num_zeros = simulation_field_length - file_number_length;

  // Create filename
  std::ostringstream simulation_filename;
  simulation_filename << simulation_file.substr(0, simulation_pos_open);
  for (int n = 0; n < num_zeros; n++)
    simulation_filename << "0";
  simulation_filename << file_number;
  simulation_filename << simulation_file.substr(simulation_pos_close + 1);
  std::string simulation_file_formatted = simulation_filename.str();
  return simulation_file_formatted;
}

//--------------------------------------------------------------------------------------------------

// Function to check that needed Athena++ variables are located as expected
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Sets indices locating specific variables among primitives.
//   Assumes metadata set.
void SimulationReader::VerifyVariablesAthena()
{
  // Check that array of all primitives is present
  int prim_offset = 0;
  for (ind_hydro = 0; ind_hydro < num_dataset_names; ind_hydro++)
    if (dataset_names[ind_hydro] == "prim")
      break;
    else
      prim_offset += num_variables(ind_hydro);
  if (ind_hydro == num_dataset_names)
    throw BlacklightException("Unable to locate array \"prim\" in data file.");

  // Check that all necessary primitives are present
  for (ind_rho = prim_offset; ind_rho < prim_offset + num_variables(ind_hydro); ind_rho++)
    if (variable_names[ind_rho] == "rho")
      break;
  if (ind_rho == num_variables(ind_hydro))
    throw BlacklightException("Unable to locate \"rho\" slice of \"prim\" in data file.");
  for (ind_pgas = prim_offset; ind_pgas < prim_offset + num_variables(ind_hydro); ind_pgas++)
    if (variable_names[ind_pgas] == "press")
      break;
  if (ind_pgas == num_variables(ind_hydro))
    throw BlacklightException("Unable to locate \"press\" slice of \"prim\" in data file.");
  if (plasma_model == PlasmaModel::code_kappa)
  {
    for (ind_kappa = prim_offset; ind_kappa < prim_offset + num_variables(ind_hydro); ind_kappa++)
      if (variable_names[ind_kappa] == simulation_kappa_name)
        break;
    if (ind_kappa == num_variables(ind_hydro))
      throw
          BlacklightException("Unable to locate electron entropy slice of \"prim\" in data file.");
  }
  for (ind_uu1 = prim_offset; ind_uu1 < prim_offset + num_variables(ind_hydro); ind_uu1++)
    if (variable_names[ind_uu1] == "vel1")
      break;
  if (ind_uu1 == num_variables(ind_hydro))
    throw BlacklightException("Unable to locate \"vel1\" slice of \"prim\" in data file.");
  for (ind_uu2 = prim_offset; ind_uu2 < prim_offset + num_variables(ind_hydro); ind_uu2++)
    if (variable_names[ind_uu2] == "vel2")
      break;
  if (ind_uu2 == num_variables(ind_hydro))
    throw BlacklightException("Unable to locate \"vel2\" slice of \"prim\" in data file.");
  for (ind_uu3 = prim_offset; ind_uu3 < prim_offset + num_variables(ind_hydro); ind_uu3++)
    if (variable_names[ind_uu3] == "vel3")
      break;
  if (ind_uu3 == num_variables(ind_hydro))
    throw BlacklightException("Unable to locate \"vel3\" slice of \"prim\" in data file.");

  // Check that array of all magnetic field components is present
  int bb_offset = 0;
  for (ind_bb = 0; ind_bb < num_dataset_names; ind_bb++)
    if (dataset_names[ind_bb] == "B")
      break;
    else
      bb_offset += num_variables(ind_bb);
  if (ind_bb == num_dataset_names)
    throw BlacklightException("Unable to locate array \"B\" in data file.");

  // Check that all necessary magnetic field components are present
  for (ind_bb1 = bb_offset; ind_bb1 < bb_offset + num_variables(ind_bb); ind_bb1++)
    if (variable_names[ind_bb1] == "Bcc1")
      break;
  if (ind_bb1 == bb_offset + num_variables(ind_bb))
    throw BlacklightException("Unable to locate \"Bcc1\" slice of \"prim\" in data file.");
  for (ind_bb2 = bb_offset; ind_bb2 < bb_offset + num_variables(ind_bb); ind_bb2++)
    if (variable_names[ind_bb2] == "Bcc2")
      break;
  if (ind_bb2 == bb_offset + num_variables(ind_bb))
    throw BlacklightException("Unable to locate \"Bcc2\" slice of \"prim\" in data file.");
  for (ind_bb3 = bb_offset; ind_bb3 < bb_offset + num_variables(ind_bb); ind_bb3++)
    if (variable_names[ind_bb3] == "Bcc3")
      break;
  if (ind_bb3 == bb_offset + num_variables(ind_bb))
    throw BlacklightException("Unable to locate \"Bcc3\" slice of \"prim\" in data file.");
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to check that needed Harm variables are located as expected
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Sets indices locating specific variables among primitives.
//   Sets adiabatic_gamma.
//   Assumes metadata set.
void SimulationReader::VerifyVariablesHarm()
{
  // Read number of primitives
  ReadHDF5IntArray("header/n_prim", num_variables);

  // Read names of primitives
  ReadHDF5StringArray("header/prim_names", num_variable_names == 0, &variable_names,
      &num_variable_names);
  if (num_variables(0) != num_variable_names)
    throw BlacklightException("Inconsistency in number of primitive variables.");

  // Check that all necessary primitives are present
  for (ind_rho = 0; ind_rho < num_variable_names; ind_rho++)
    if (variable_names[ind_rho] == "RHO")
      break;
  if (ind_rho == num_variable_names)
    throw BlacklightException("Unable to locate \"RHO\" slice of \"prims\" in data file.");

  for (ind_pgas = 0; ind_pgas < num_variable_names; ind_pgas++)
    if (variable_names[ind_pgas] == "UU")
      break;
  if (ind_pgas == num_variable_names)
    throw BlacklightException("Unable to locate \"UU\" slice of \"prims\" in data file.");
  if (plasma_model == PlasmaModel::code_kappa)
  {
    for (ind_kappa = 0; ind_kappa < num_variable_names; ind_kappa++)
      if (variable_names[ind_kappa] == simulation_kappa_name)
        break;
    if (ind_kappa == num_variable_names)
      throw
          BlacklightException("Unable to locate electron entropy slice of \"prims\" in data file.");
  }
  for (ind_uu1 = 0; ind_uu1 < num_variable_names; ind_uu1++)
    if (variable_names[ind_uu1] == "U1")
      break;
  if (ind_uu1 == num_variable_names)
    throw BlacklightException("Unable to locate \"U1\" slice of \"prims\" in data file.");
  for (ind_uu2 = 0; ind_uu2 < num_variable_names; ind_uu2++)
    if (variable_names[ind_uu2] == "U2")
      break;
  if (ind_uu2 == num_variable_names)
    throw BlacklightException("Unable to locate \"U2\" slice of \"prims\" in data file.");
  for (ind_uu3 = 0; ind_uu3 < num_variable_names; ind_uu3++)
    if (variable_names[ind_uu3] == "U3")
      break;
  if (ind_uu3 == num_variable_names)
    throw BlacklightException("Unable to locate \"U3\" slice of \"prims\" in data file.");
  for (ind_bb1 = 0; ind_bb1 < num_variable_names; ind_bb1++)
    if (variable_names[ind_bb1] == "B1")
      break;
  if (ind_bb1 == num_variable_names)
    throw BlacklightException("Unable to locate \"B1\" slice of \"prims\" in data file.");
  for (ind_bb2 = 0; ind_bb2 < num_variable_names; ind_bb2++)
    if (variable_names[ind_bb2] == "B2")
      break;
  if (ind_bb2 == num_variable_names)
    throw BlacklightException("Unable to locate \"B2\" slice of \"prims\" in data file.");
  for (ind_bb3 = 0; ind_bb3 < num_variable_names; ind_bb3++)
    if (variable_names[ind_bb3] == "B3")
      break;
  if (ind_bb3 == num_variable_names)
    throw BlacklightException("Unable to locate \"B3\" slice of \"prims\" in data file.");

  // Read adiabatic index
  Array<double> gamma;
  ReadHDF5DoubleArray("header/gam", gamma);
  adiabatic_gamma = gamma(0);
  return;
}
