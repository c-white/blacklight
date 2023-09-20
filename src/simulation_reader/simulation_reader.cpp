// Blacklight simulation reader

// C++ headers
#include <algorithm>  // remove
#include <cmath>      // abs
#include <cstdint>    // int32_t
#include <cstdio>     // snprintf
#include <cstring>    // strncmp, strtok
#include <fstream>    // ifstream
#include <ios>        // ios_base
#include <iosfwd>     // streampos
#include <optional>   // optional
#include <sstream>    // ostringstream
#include <string>     // getline, stod, stoi, string

// Library headers
#include <omp.h>  // pragmas, omp_get_wtime

// Blacklight headers
#include "simulation_reader.hpp"
#include "../blacklight.hpp"                 // Math, enums
#include "../input_reader/input_reader.hpp"  // InputReader
#include "../utils/array.hpp"                // Array
#include "../utils/exceptions.hpp"           // BlacklightException, BlacklightWarning
#include "../utils/file_io.hpp"              // ReadBinary

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
    simulation_m_msun = p_input_reader->simulation_m_msun.value();
    simulation_rho_cgs = p_input_reader->simulation_rho_cgs.value();
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
    plasma_mu = p_input_reader->plasma_mu.value();
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
  if (model_type == ModelType::simulation and simulation_format == SimulationFormat::athenak
      and athenak_variable_size == 8)
    delete[] athenak_cell_data_double;
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
      if (simulation_format == SimulationFormat::athena
          or simulation_format == SimulationFormat::iharm3d)
      {
        ReadHDF5Superblock();
        root_data_segment_address = ReadHDF5Heap(root_name_heap_address);
        ReadHDF5RootObjectHeader();
      }
      else if (simulation_format == SimulationFormat::athenak)
        ReadAthenaKHeader();

      // Read time
      if (simulation_format == SimulationFormat::athena)
      {
        float time_temp;
        ReadHDF5FloatAttribute("Time", &time_temp);
        latest_time = time_temp;
      }
      else if (simulation_format == SimulationFormat::athenak)
        latest_time = athenak_time;
      else if (simulation_format == SimulationFormat::iharm3d)
      {
        Array<double> time_temp(1);
        ReadHDF5DoubleArray("t", time_temp);
        latest_time = time_temp(0);
      }
      else if (simulation_format == SimulationFormat::harm3d)
        data_stream >> latest_time;
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
    if (simulation_format == SimulationFormat::athena
        or simulation_format == SimulationFormat::iharm3d)
    {
      ReadHDF5Superblock();
      root_data_segment_address = ReadHDF5Heap(root_name_heap_address);
      ReadHDF5RootObjectHeader();
    }
    else if (simulation_format == SimulationFormat::athenak)
    {
      ReadAthenaKHeader();
      if (first_time)
        ReadAthenaKInputs();
    }

    // Read time
    if (simulation_format == SimulationFormat::athena)
    {
      float time_temp;
      ReadHDF5FloatAttribute("Time", &time_temp);
      time[n] = time_temp;
    }
    else if (simulation_format == SimulationFormat::athenak)
      time[n] = athenak_time;
    else if (simulation_format == SimulationFormat::iharm3d)
    {
      Array<double> time_temp(1);
      ReadHDF5DoubleArray("t", time_temp);
      time[n] = time_temp(0);
    }
    else if (simulation_format == SimulationFormat::harm3d)
      data_stream >> time[n];

    // Read metric
    if (first_time and simulation_format == SimulationFormat::iharm3d)
    {
      std::string *p_temp_metric;
      int temp_count;
      ReadHDF5StringArray("header/metric", true, &p_temp_metric, &temp_count);
      metric = *p_temp_metric;
      delete[] p_temp_metric;
      if (simulation_coord == Coordinates::sks or simulation_coord == Coordinates::fmks)
      {
        std::string metric_lower = metric;
        std::transform(metric_lower.begin(), metric_lower.end(), metric_lower.begin(),
                       [](unsigned char c){ return std::tolower(c); });

        if (metric != "MKS" and metric != "MMKS" and metric != "FMKS")
        {
          std::ostringstream message;
          message << "Given metric mks does not match file value of " << metric;
          message << "; ignoring the latter.";
          BlacklightWarning(message.str().c_str());
        }
        Array<double> a_temp, h_temp;
        ReadHDF5DoubleArray(("header/geom/" + metric_lower + "/a").c_str(), a_temp);
        ReadHDF5DoubleArray(("header/geom/" + metric_lower + "/hslope").c_str(), h_temp);
        metric_a = a_temp(0);
        if (metric_a != simulation_a)
        {
          std::ostringstream message;
          message << "Given spin of " << simulation_a << " does not match file value of ";
          message << metric_a << "; ignoring the latter.";
          BlacklightWarning(message.str().c_str());
        }
        metric_h = h_temp(0);
        if (metric == "MMKS" or metric == "FMKS") {
          Array<double> poly_xt_temp, poly_alpha_temp, mks_smooth_temp, rin_temp;
          ReadHDF5DoubleArray(("header/geom/" + metric_lower + "/poly_xt").c_str(), poly_xt_temp);
          ReadHDF5DoubleArray(("header/geom/" + metric_lower + "/poly_alpha").c_str(), poly_alpha_temp);
          ReadHDF5DoubleArray(("header/geom/" + metric_lower + "/mks_smooth").c_str(), mks_smooth_temp);
          try
          {
            ReadHDF5DoubleArray(("header/geom/" + metric_lower + "/r_in").c_str(), rin_temp);
          }
          catch (...)
          {
            try
            {
              ReadHDF5DoubleArray(("header/geom/" + metric_lower + "/Rin").c_str(), rin_temp);
            }
            catch (...) {
              throw BlacklightException("Unable to identify r_in parameter for iharm3d-format file.");
            }
          }
          metric_rin = rin_temp(0);
          metric_poly_xt = poly_xt_temp(0);
          metric_poly_alpha = poly_alpha_temp(0);
          metric_mks_smooth = mks_smooth_temp(0);
          metric_derived_poly_norm = 1.0 / (1.0 + 1.0/(metric_poly_alpha+1.0)/std::pow(metric_poly_xt, metric_poly_alpha));
          metric_derived_poly_norm *= 0.5 * Math::pi;
          //native_x1in = std::log(metric_rin);
          //native_deltax1 =
          //native_deltax2 = 1.0 / n2;
        }
      }
      else
        throw BlacklightException("Invalid simulation_coord for Harm format.");
    }

    // Read AthenaK data
    if (simulation_format == SimulationFormat::athenak)
    {
      // Count blocks
      if (first_time)
      {
        data_stream.seekg(athenak_data_offset);
        int32_t block_indices[6];
        data_stream.read(reinterpret_cast<char *>(block_indices), 24);
        athenak_block_nx = block_indices[1] - block_indices[0] + 1;
        athenak_block_ny = block_indices[3] - block_indices[2] + 1;
        athenak_block_nz = block_indices[5] - block_indices[4] + 1;
        athenak_cells_per_block = athenak_block_nz * athenak_block_ny * athenak_block_nx;
        athenak_block_size_bytes = 24 + 16 + 6 * athenak_location_size
            + num_variable_names * athenak_cells_per_block * athenak_variable_size;
        athenak_num_blocks = 0;
        data_stream.seekg(athenak_data_offset);
        while (not data_stream.eof())
        {
          data_stream.ignore(athenak_block_size_bytes);
          athenak_num_blocks++;
        }
        athenak_num_blocks--;
      }

      // Allocate arrays
      if (first_time)
      {
        levels.Allocate(athenak_num_blocks);
        locations.Allocate(athenak_num_blocks, 3);
        x1f.Allocate(athenak_num_blocks, athenak_block_nx + 1);
        x2f.Allocate(athenak_num_blocks, athenak_block_ny + 1);
        x3f.Allocate(athenak_num_blocks, athenak_block_nz + 1);
        x1v.Allocate(athenak_num_blocks, athenak_block_nx);
        x2v.Allocate(athenak_num_blocks, athenak_block_ny);
        x3v.Allocate(athenak_num_blocks, athenak_block_nz);
        int n5 = plasma_model == PlasmaModel::code_kappa ? 9 : 8;
        for (int nn = 0; nn < num_read; nn++)
          prim[nn].Allocate(n5, athenak_num_blocks, athenak_block_nz, athenak_block_ny,
              athenak_block_nx);
      }

      // Go through blocks
      data_stream.seekg(athenak_data_offset);
      for (int block = 0; block < athenak_num_blocks; block++)
      {
        // Skip block indices
        data_stream.ignore(24);

        // Read block layout
        if (first_time)
        {
          data_stream.read(reinterpret_cast<char *>(&locations(block,0)), 12);
          data_stream.read(reinterpret_cast<char *>(&levels(block)), 4);
        }
        else
          data_stream.ignore(16);

        // Read coordinates
        if (first_time)
        {
          double face_coordinates[6];
          if (athenak_location_size == 4)
          {
            float face_coordinates_single[6];
            data_stream.read(reinterpret_cast<char *>(face_coordinates_single), 24);
            for (int ind = 0; ind < 6; ind++)
              face_coordinates[ind] = face_coordinates_single[ind];
          }
          else if (athenak_location_size == 8)
            data_stream.read(reinterpret_cast<char *>(face_coordinates), 48);
          x1f(block,0) = face_coordinates[0];
          x1f(block,athenak_block_nx) = face_coordinates[1];
          double dx = (face_coordinates[1] - face_coordinates[0]) / athenak_block_nx;
          for (int i = 1; i < athenak_block_nx; i++)
            x1f(block,i) = face_coordinates[0] + i * dx;
          for (int i = 0; i < athenak_block_nx; i++)
            x1v(block,i) = 0.5 * (x1f(block,i) + x1f(block,i+1));
          x2f(block,0) = face_coordinates[2];
          x2f(block,athenak_block_ny) = face_coordinates[3];
          double dy = (face_coordinates[3] - face_coordinates[2]) / athenak_block_ny;
          for (int j = 1; j < athenak_block_ny; j++)
            x2f(block,j) = face_coordinates[2] + j * dy;
          for (int j = 0; j < athenak_block_nx; j++)
            x2v(block,j) = 0.5 * (x2f(block,j) + x2f(block,j+1));
          x3f(block,0) = face_coordinates[4];
          x3f(block,athenak_block_nz) = face_coordinates[5];
          double dz = (face_coordinates[5] - face_coordinates[4]) / athenak_block_nz;
          for (int k = 1; k < athenak_block_nz; k++)
            x3f(block,k) = face_coordinates[4] + k * dz;
          for (int k = 0; k < athenak_block_nx; k++)
            x3v(block,k) = 0.5 * (x3f(block,k) + x3f(block,k+1));
        }
        else
          data_stream.ignore(6 * athenak_location_size);

        // Read cell data
        if (first_time)
          VerifyVariablesAthenaK();
        if (first_time and athenak_variable_size == 8)
          athenak_cell_data_double = new double[athenak_cells_per_block];
        std::streampos cell_data_begin = data_stream.tellg();
        int athenak_inds[8] = {athenak_ind_rho, athenak_ind_uu1, athenak_ind_uu2, athenak_ind_uu3,
            athenak_ind_pgas, athenak_ind_bb1, athenak_ind_bb2, athenak_ind_bb3};
        int prim_inds[8] =
            {ind_rho, ind_uu1, ind_uu2, ind_uu3, ind_pgas, ind_bb1, ind_bb2, ind_bb3};
        for (int ind_ind = 0; ind_ind < 8; ind_ind++)
        {
          data_stream.seekg(cell_data_begin);
          int offset = athenak_inds[ind_ind] * athenak_cells_per_block * athenak_variable_size;
          data_stream.seekg(offset, std::ios_base::cur);
          if (athenak_variable_size == 4)
            data_stream.read(reinterpret_cast<char *>(&prim[n](prim_inds[ind_ind],block,0,0,0)),
                athenak_cells_per_block * 4);
          else if (athenak_variable_size == 8)
          {
            data_stream.read(reinterpret_cast<char *>(athenak_cell_data_double),
                athenak_cells_per_block * 8);
            for (int k = 0, ind = 0; k < athenak_block_nz; k++)
              for (int j = 0; j < athenak_block_ny; j++)
                for (int i = 0; i < athenak_block_nx; i++, ind++)
                  prim[n](prim_inds[ind_ind],block,k,j,i) =
                      static_cast<float>(athenak_cell_data_double[ind]);
          }
        }
        if (plasma_model == PlasmaModel::code_kappa)
        {
          data_stream.seekg(cell_data_begin);
          int offset = athenak_ind_kappa * athenak_cells_per_block * athenak_variable_size;
          data_stream.seekg(offset, std::ios_base::cur);
          if (athenak_variable_size == 4)
            data_stream.read(reinterpret_cast<char *>(&prim[n](ind_kappa,block,0,0,0)),
                athenak_cells_per_block * 4);
          else if (athenak_variable_size == 8)
          {
            data_stream.read(reinterpret_cast<char *>(athenak_cell_data_double),
                athenak_cells_per_block * 8);
            for (int k = 0, ind = 0; k < athenak_block_nz; k++)
              for (int j = 0; j < athenak_block_ny; j++)
                for (int i = 0; i < athenak_block_nx; i++, ind++)
                  prim[n](ind_kappa,block,k,j,i) =
                      static_cast<float>(athenak_cell_data_double[ind]);
          }
        }
      }

      // Convert internal energy to pressure
      #pragma omp parallel for schedule(static) collapse(3)
      for (int block = 0; block < athenak_num_blocks; block++)
        for (int k = 0; k < athenak_block_nz; k++)
          for (int j = 0; j < athenak_block_ny; j++)
            for (int i = 0; i < athenak_block_nx; i++)
              prim[n](ind_pgas,block,k,j,i) *= static_cast<float>(adiabatic_gamma - 1.0);
    }

    // Read block layout
    if (first_time)
    {
      if (simulation_format == SimulationFormat::athena)
      {
        ReadHDF5IntArray("Levels", levels);
        ReadHDF5IntArray("LogicalLocations", locations);
      }
      else if (simulation_format == SimulationFormat::iharm3d
          or simulation_format == SimulationFormat::harm3d)
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
      else if (simulation_format == SimulationFormat::iharm3d)
      {
        // TODO check that this makes sense.
        Array<int> num_cells;
        Array<double> x_start, dx;
        ReadHDF5IntArray("header/n1", num_cells);
        ReadHDF5DoubleArray("header/geom/startx1", x_start);
        ReadHDF5DoubleArray("header/geom/dx1", dx);
        x1f.Allocate(1, num_cells(0) + 1);
        x1v.Allocate(1, num_cells(0));
        x1f(0,0) = x_start(0);
        for (int i = 0; i < num_cells(0); i++)
        {
          // TODO: note the legacy code here assumed x1f started with 0, which is not always true!
          x1f(0,i+1) = x_start(0) + (i + 1) * dx(0);
          x1v(0,i) = 0.5 * (x1f(0,i) + x1f(0,i+1));
        }
        ReadHDF5IntArray("header/n2", num_cells);
        ReadHDF5DoubleArray("header/geom/startx2", x_start);
        ReadHDF5DoubleArray("header/geom/dx2", dx);
        x2f.Allocate(1, num_cells(0) + 1);
        x2v.Allocate(1, num_cells(0));
        for (int j = 0; j < num_cells(0); j++)
        {
          x2f(0,j+1) = x_start(0) + (j + 1) * dx(0);
          x2v(0,j) = 0.5 * (x2f(0,j) + x2f(0,j+1));
        }
        ReadHDF5IntArray("header/n3", num_cells);
        ReadHDF5DoubleArray("header/geom/startx3", x_start);
        ReadHDF5DoubleArray("header/geom/dx3", dx);
        x3f.Allocate(1, num_cells(0) + 1);
        x3v.Allocate(1, num_cells(0));
        for (int k = 0; k < num_cells(0); k++)
        {
          x3f(0,k+1) = x_start(0) + (k + 1) * dx(0);
          x3v(0,k) = 0.5 * (x3f(0,k) + x3f(0,k+1));
        }
        ConvertCoordinates();
      }
      else if (simulation_format == SimulationFormat::harm3d)
      {
        int num_cells_1, num_cells_2, num_cells_3;
        data_stream >> num_cells_1 >> num_cells_2 >> num_cells_3;
        double x1_start, x2_start, x3_start;
        data_stream >> x1_start >> x2_start >> x3_start;
        double dx1, dx2, dx3;
        data_stream >> dx1 >> dx2 >> dx3;
        x1f.Allocate(1, num_cells_1 + 1);
        x1v.Allocate(1, num_cells_1);
        for (int i = 0; i < num_cells_1; i++)
        {
          x1f(0,i+1) = x1_start + (i + 1) * dx1;
          x1v(0,i) = 0.5 * (x1f(0,i) + x1f(0,i+1));
        }
        x2f.Allocate(1, num_cells_2 + 1);
        x2v.Allocate(1, num_cells_2);
        for (int j = 0; j < num_cells_2; j++)
        {
          x2f(0,j+1) = x2_start + (j + 1) * dx2;
          x2v(0,j) = 0.5 * (x2f(0,j) + x2f(0,j+1));
        }
        x3f.Allocate(1, num_cells_3 + 1);
        x3v.Allocate(1, num_cells_3);
        for (int k = 0; k < num_cells_3; k++)
        {
          x3f(0,k+1) = x3_start + (k + 1) * dx3;
          x3v(0,k) = 0.5 * (x3f(0,k) + x3f(0,k+1));
        }
        data_stream >> metric_a;
        if (metric_a != simulation_a)
        {
          std::ostringstream message;
          message << "Given spin of " << simulation_a << " does not match file value of ";
          message << metric_a << "; ignoring the latter.";
          BlacklightWarning(message.str().c_str());
        }
        data_stream >> adiabatic_gamma;
        double temp_val;
        data_stream >> temp_val;
        data_stream >> metric_h;
        data_stream >> temp_val;
        data_stream.seekg(1, std::ios_base::cur);
        cell_data_address = data_stream.tellg();
        ConvertCoordinates();
      }
    }

    // Check coordinates
    if (first_time)
    {
      // TODO: I wantonly added fmks in here and below. Check for correctness?
      if ((simulation_coord == Coordinates::sks or simulation_coord == Coordinates::fmks)
          and x2f.n2 == 1)
      {
        bool error_low = std::abs(x2f(0,0)) > (x2f(0,1) - x2f(0,0)) * angular_domain_tolerance;
        bool error_high = std::abs(x2f(0,x2f.n1-1) - Math::pi)
            > (x2f(0,x2f.n1-1) - x2f(0,x2f.n1-2)) * angular_domain_tolerance;
        if (error_low or error_high)
        {
          std::ostringstream message;
          message.setf(std::ios_base::scientific);
          message.precision(16);
          message << "Changing theta range from [" << x2f(0,0) << ", " << x2f(0,x2f.n1-1);
          message << "] to [0, pi].";
          BlacklightWarning(message.str().c_str());
          x2f(0,0) = 0.0;
          x2f(0,x2f.n1-1) = Math::pi;
        }
      }
      if ((simulation_coord == Coordinates::sks or simulation_coord == Coordinates::fmks)
          and x3f.n2 == 1)
      {
        bool error_low = std::abs(x3f(0,0)) > (x3f(0,1) - x3f(0,0)) * angular_domain_tolerance;
        bool error_high = std::abs(x3f(0,x3f.n1-1) - 2.0 * Math::pi)
            > (x3f(0,x3f.n1-1) - x3f(0,x3f.n1-2)) * angular_domain_tolerance;
        if (error_low or error_high)
        {
          std::ostringstream message;
          message.setf(std::ios_base::scientific);
          message.precision(16);
          message << "Changing phi range from [" << x3f(0,0) << ", " << x3f(0,x3f.n1-1);
          message << "] to [0, 2*pi].";
          BlacklightWarning(message.str().c_str());
          x3f(0,0) = 0.0;
          x3f(0,x3f.n1-1) = 2.0 * Math::pi;
        }
      }
    }

    // Read cell data
    if (simulation_format == SimulationFormat::athena)
    {
      if (first_time)
      {
        VerifyVariablesAthena();
        int n5 = num_variables(ind_hydro) + num_variables(ind_bb);
        int n4 = levels.n1;
        int n3 = x3v.n1;
        int n2 = x2v.n1;
        int n1 = x1v.n1;
        for (int nn = 0; nn < num_read; nn++)
          prim[nn].Allocate(n5, n4, n3, n2, n1);
      }
      Array<float> hydro(prim[n]);
      hydro.Slice(5, 0, num_variables(ind_hydro) - 1);
      ReadHDF5FloatArray("prim", hydro);
      Array<float> bb(prim[n]);
      bb.Slice(5, num_variables(ind_hydro), num_variables(ind_hydro) + num_variables(ind_bb) - 1);
      ReadHDF5FloatArray("B", bb);
    }
    else if (simulation_format == SimulationFormat::iharm3d)
    {
      if (first_time)
      {
        // TODO support num_variables and Thetae from electron thermodynamics
        VerifyVariablesHarm();
        int n5 = num_variables(0);
        int n4 = levels.n1;
        int n3 = x3v.n1;
        int n2 = x2v.n1;
        int n1 = x1v.n1;
        for (int nn = 0; nn < num_read; nn++)
          prim[nn].Allocate(n5, n4, n3, n2, n1);
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
      ConvertPrimitives3(prim[n]);
    }
    else if (simulation_format == SimulationFormat::harm3d)
    {
      if (first_time)
      {
        int n5 = plasma_model == PlasmaModel::code_kappa ? 11 : 10;
        int n4 = levels.n1;
        int n3 = x3v.n1;
        int n2 = x2v.n1;
        int n1 = x1v.n1;
        for (int nn = 0; nn < num_read; nn++)
          prim[nn].Allocate(n5, n4, n3, n2, n1);
        prim_transpose.Allocate(n1, n2, n3, n5 + 6);
        ind_rho = 0;
        ind_pgas = 1;
        ind_kappa = 10;
        ind_u0 = 2;
        ind_uu1 = 3;
        ind_uu2 = 4;
        ind_uu3 = 5;
        ind_b0 = 6;
        ind_bb1 = 7;
        ind_bb2 = 8;
        ind_bb3 = 9;
      }
      else
        data_stream.seekg(cell_data_address);
      ReadBinary(&data_stream, prim_transpose.data, prim_transpose.n_tot);
      #pragma omp parallel
      {
        #pragma omp for schedule(static) collapse(3)
        for (int n_variable = 0; n_variable < prim[n].n5; n_variable++)
          for (int k = 0; k < x3v.n1; k++)
            for (int j = 0; j < x2v.n1; j++)
              for (int i = 0; i < x1v.n1; i++)
                prim[n](n_variable,0,k,j,i) = prim_transpose(i,j,k,n_variable+6);
        #pragma omp for schedule(static) collapse(2)
        for (int k = 0; k < x3v.n1; k++)
          for (int j = 0; j < x2v.n1; j++)
            for (int i = 0; i < x1v.n1; i++)
              prim[n](ind_pgas,0,k,j,i) *= static_cast<float>(adiabatic_gamma - 1.0);
      }
      ConvertPrimitives4(prim[n]);
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

// Function for reading initial header data in AthenaK dump
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Sets athenak_location_size, athenak_variable_size, athenak_data_offset, and athenak_time.
//   Allocates and sets variable_names and num_variable_names.
//   Changes stream pointer to beginning of input parameter section.
void SimulationReader::ReadAthenaKHeader()
{
  // Prepare buffer
  constexpr int buffer_size = 128;
  Array<char> buffer;
  buffer.Allocate(buffer_size);

  // Check version
  data_stream.seekg(0);
  buffer.Zero();
  data_stream.getline(buffer.data, buffer_size);
  data_stream.seekg(-1, std::ios_base::cur);
  if (data_stream.get() != '\n')
    throw BlacklightException("Unknown AthenaK file format.");
  if (std::strncmp(buffer.data, "Athena binary output version=1.1", buffer_size) != 0)
    throw BlacklightException("Unknown AthenaK file format.");

  // Read time
  data_stream.getline(buffer.data, buffer_size);
  buffer.Zero();
  std::ifstream::pos_type file_location = data_stream.tellg();
  data_stream.getline(buffer.data, buffer_size);
  if (std::strncmp(buffer.data, "  time=", 7) != 0)
    throw BlacklightException("Invalid AthenaK file header.");
  data_stream.seekg(file_location);
  data_stream.seekg(7, std::ios_base::cur);
  data_stream >> athenak_time;
  data_stream.seekg(1, std::ios_base::cur);

  // Read size of location
  data_stream.getline(buffer.data, buffer_size);
  buffer.Zero();
  file_location = data_stream.tellg();
  data_stream.getline(buffer.data, buffer_size);
  if (std::strncmp(buffer.data, "  size of location=", 19) != 0)
    throw BlacklightException("Invalid AthenaK file header.");
  data_stream.seekg(file_location);
  data_stream.seekg(19, std::ios_base::cur);
  data_stream >> athenak_location_size;
  if (athenak_location_size != 4 and athenak_location_size != 8)
    throw BlacklightException("Unsupported size of location.");
  data_stream.seekg(1, std::ios_base::cur);

  // Read size of variable
  buffer.Zero();
  file_location = data_stream.tellg();
  data_stream.getline(buffer.data, buffer_size);
  if (std::strncmp(buffer.data, "  size of variable=", 19) != 0)
    throw BlacklightException("Invalid AthenaK file header.");
  data_stream.seekg(file_location);
  data_stream.seekg(19, std::ios_base::cur);
  data_stream >> athenak_variable_size;
  if (athenak_location_size != 4 and athenak_location_size != 8)
    throw BlacklightException("Unsupported size of variables.");
  data_stream.seekg(1, std::ios_base::cur);

  // Read number of variables
  buffer.Zero();
  file_location = data_stream.tellg();
  data_stream.getline(buffer.data, buffer_size);
  if (std::strncmp(buffer.data, "  number of variables=", 22) != 0)
    throw BlacklightException("Invalid AthenaK file header.");
  data_stream.seekg(file_location);
  data_stream.seekg(22, std::ios_base::cur);
  data_stream >> num_variable_names;
  data_stream.seekg(1, std::ios_base::cur);

  // Read variable names
  buffer.Zero();
  data_stream.getline(buffer.data, buffer_size);
  if (std::strncmp(buffer.data, "  variables:", 12) != 0)
    throw BlacklightException("Invalid AthenaK file header.");
  variable_names = new std::string[num_variable_names];
  for (int n = 0; n < num_variable_names; n++)
  {
    if (n == 0)
      variable_names[n] = std::strtok(buffer.data + 12, " ");
    else
      variable_names[n] = std::strtok(nullptr, " ");
  }

  // Read header offset
  buffer.Zero();
  file_location = data_stream.tellg();
  data_stream.getline(buffer.data, buffer_size);
  if (std::strncmp(buffer.data, "  header offset=", 16) != 0)
    throw BlacklightException("Invalid AthenaK file header.");
  data_stream.seekg(file_location);
  data_stream.seekg(16, std::ios_base::cur);
  int header_offset;
  data_stream >> header_offset;
  data_stream.seekg(1, std::ios_base::cur);

  // Calculate data offset
  file_location = data_stream.tellg();
  data_stream.seekg(header_offset, std::ios_base::cur);
  athenak_data_offset = data_stream.tellg();
  data_stream.seekg(file_location);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for reading input parameters from Athenak dump
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Checks simulation_a, simulation_m_msun, simulation_rho_cgs, and plasma_mu for consistency with
//       user input.
//   Sets adiabatic_gamma.
//   Assumes stream pointer points to beginning of input parameter section.
//   Changes stream pointer.
void SimulationReader::ReadAthenaKInputs()
{
  // Prepare buffers and flag
  std::string buffer, section_name, variable_name;
  bool gamma_found = false;

  // Go though inputs
  while (data_stream.tellg() < athenak_data_offset)
  {
    // Read line
    std::getline(data_stream, buffer);

    // Skip comments
    if (buffer[0] == '#')
      continue;

    // Extract section name
    if (buffer.front() == '<' and buffer.back() == '>')
    {
      section_name = std::string(buffer, 1, buffer.size() - 2);
      continue;
    }

    // Extract variable name
    std::string::size_type location = buffer.find_first_of('=');
    if (location == std::string::npos)
      throw BlacklightException("Error parsing inputs in AthenaK file.");
    variable_name = buffer.substr(0, location);
    variable_name.erase(std::remove(variable_name.begin(), variable_name.end(), ' '),
        variable_name.end());

    // Extract spin
    if (section_name == "coord" and variable_name == "a")
    {
      double file_a = std::stod(buffer.substr(location + 1));
      if (file_a != simulation_a)
      {
        std::ostringstream message;
        message << "Given spin of " << simulation_a << " does not match file value of " << file_a;
        message << "; ignoring the latter.";
        BlacklightWarning(message.str().c_str());
      }
    }

    // Extract mass
    if (section_name == "units" and variable_name == "bhmass_msun")
    {
      double file_m_msun = std::stod(buffer.substr(location + 1));
      if (file_m_msun != simulation_m_msun)
      {
        std::ostringstream message;
        message << "Given mass of " << simulation_m_msun << " does not match file value of ";
        message << file_m_msun << "; ignoring the latter.";
        BlacklightWarning(message.str().c_str());
      }
    }

    // Extract density scale
    if (section_name == "units" and variable_name == "density_cgs")
    {
      double file_rho_cgs = std::stod(buffer.substr(location + 1));
      if (file_rho_cgs != simulation_rho_cgs)
      {
        std::ostringstream message;
        message << "Given density scale of " << simulation_rho_cgs;
        message << " does not match file value of " << file_rho_cgs << "; ignoring the latter.";
        BlacklightWarning(message.str().c_str());
      }
    }

    // Extract molecular weight
    if (section_name == "units" and variable_name == "mu")
    {
      double file_mu = std::stod(buffer.substr(location + 1));
      if (file_mu != plasma_mu)
      {
        std::ostringstream message;
        message << "Given density scale of " << plasma_mu << " does not match file value of ";
        message << file_mu << "; ignoring the latter.";
        BlacklightWarning(message.str().c_str());
      }
    }

    // Extract adiabatic index
    if (section_name == "mhd" and variable_name == "gamma")
    {
      adiabatic_gamma = std::stod(buffer.substr(location + 1));
      gamma_found = true;
    }
  }

  // Check for required value
  if (not gamma_found)
    throw BlacklightException("Missing adiabatic index.");
  return;
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

// Function to check that needed AthenaK variables are located as expected
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Sets indices locating specific variables in file.
//   Assumes metadata set.
void SimulationReader::VerifyVariablesAthenaK()
{
  // Check that all necessary hydrodynamical values are present
  for (athenak_ind_rho = 0; athenak_ind_rho < num_variable_names; athenak_ind_rho++)
    if (variable_names[athenak_ind_rho] == "dens")
      break;
  if (athenak_ind_rho == num_variable_names)
    throw BlacklightException("Unable to locate \"dens\" values in data file.");
  for (athenak_ind_pgas = 0; athenak_ind_pgas < num_variable_names; athenak_ind_pgas++)
    if (variable_names[athenak_ind_pgas] == "eint")
      break;
  if (athenak_ind_pgas == num_variable_names)
    throw BlacklightException("Unable to locate \"eint\" values in data file.");
  if (plasma_model == PlasmaModel::code_kappa)
  {
    for (athenak_ind_kappa = 0; athenak_ind_kappa < num_variable_names; athenak_ind_kappa++)
      if (variable_names[athenak_ind_kappa] == simulation_kappa_name)
        break;
    if (athenak_ind_kappa == num_variable_names)
      throw BlacklightException("Unable to locate electron entropy values in data file.");
  }
  for (athenak_ind_uu1 = 0; athenak_ind_uu1 < num_variable_names; athenak_ind_uu1++)
    if (variable_names[athenak_ind_uu1] == "velx")
      break;
  if (athenak_ind_uu1 == num_variable_names)
    throw BlacklightException("Unable to locate \"velx\" values in data file.");
  for (athenak_ind_uu2 = 0; athenak_ind_uu2 < num_variable_names; athenak_ind_uu2++)
    if (variable_names[athenak_ind_uu2] == "vely")
      break;
  if (athenak_ind_uu2 == num_variable_names)
    throw BlacklightException("Unable to locate \"vely\" values in data file.");
  for (athenak_ind_uu3 = 0; athenak_ind_uu3 < num_variable_names; athenak_ind_uu3++)
    if (variable_names[athenak_ind_uu3] == "velz")
      break;
  if (athenak_ind_uu3 == num_variable_names)
    throw BlacklightException("Unable to locate \"velz\" values in data file.");

  // Check that all necessary magnetic field components are present
  for (athenak_ind_bb1 = 0; athenak_ind_bb1 < num_variable_names; athenak_ind_bb1++)
    if (variable_names[athenak_ind_bb1] == "bcc1")
      break;
  if (athenak_ind_bb1 == num_variable_names)
    throw BlacklightException("Unable to locate \"bcc1\" values in data file.");
  for (athenak_ind_bb2 = 0; athenak_ind_bb2 < num_variable_names; athenak_ind_bb2++)
    if (variable_names[athenak_ind_bb2] == "bcc2")
      break;
  if (athenak_ind_bb2 == num_variable_names)
    throw BlacklightException("Unable to locate \"bcc2\" values in data file.");
  for (athenak_ind_bb3 = 0; athenak_ind_bb3 < num_variable_names; athenak_ind_bb3++)
    if (variable_names[athenak_ind_bb3] == "bcc3")
      break;
  if (athenak_ind_bb3 == num_variable_names)
    throw BlacklightException("Unable to locate \"bcc3\" values in data file.");

  // Set indices for internal arrays
  ind_rho = 0;
  ind_uu1 = 1;
  ind_uu2 = 2;
  ind_uu3 = 3;
  ind_pgas = 4;
  ind_bb1 = 5;
  ind_bb2 = 6;
  ind_bb3 = 7;
  ind_kappa = 8;
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

  // TODO make sure this is working correctly

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
