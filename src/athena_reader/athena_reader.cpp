// Blacklight Athena++ reader

// C++ headers
#include <fstream>   // ifstream
#include <ios>       // ios_base
#include <optional>  // optional
#include <string>    // string

// Blacklight headers
#include "athena_reader.hpp"
#include "../utils/array.hpp"       // Array
#include "../utils/exceptions.hpp"  // BlacklightException

//--------------------------------------------------------------------------------------------------

// Athena++ reader constructor
// Inputs:
//   p_input_reader: pointer to object containing input parameters
// Notes:
//   Opens stream for reading.
AthenaReader::AthenaReader(const InputReader *p_input_reader)
{
  // Copy general input data
  model_type = p_input_reader->model_type.value();

  // Proceed only if needed
  if (model_type == ModelType::simulation)
  {
    // Copy simulation and plasma parameters
    simulation_file = p_input_reader->simulation_file_formatted;
    plasma_model = p_input_reader->plasma_model.value();
    if (plasma_model == PlasmaModel::code_kappa)
      simulation_kappa_name = p_input_reader->simulation_kappa_name.value();

    // Open file
    data_stream = std::ifstream(simulation_file, std::ios_base::in | std::ios_base::binary);
    if (not data_stream.is_open())
      throw BlacklightException("Could not open data file.");
  }
}

//--------------------------------------------------------------------------------------------------

// Athena++ reader destructor
AthenaReader::~AthenaReader()
{
  // Free memory
  if (num_dataset_names > 0)
    delete[] dataset_names;
  if (num_variable_names > 0)
    delete[] variable_names;
  if (num_children > 0)
    delete[] children_addresses;
}

//--------------------------------------------------------------------------------------------------

// Athena++ reader read and initialize function
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Does nothing if model does not need to be read from file.
//   Initializes all member objects.
//   Implements a subset of the HDF5 standard:
//       portal.hdfgroup.org/display/HDF5/File+Format+Specification
void AthenaReader::Read()
{
  // Only proceed if needed
  if (model_type != ModelType::simulation)
    return;

  // Read basic data about file
  ReadHDF5Superblock();
  ReadHDF5RootHeap();
  ReadHDF5RootObjectHeader();
  ReadHDF5Tree();

  // Read block layout
  ReadHDF5IntArray("Levels", levels);
  ReadHDF5IntArray("LogicalLocations", locations);

  // Read coordinates
  ReadHDF5FloatArray("x1f", x1f);
  ReadHDF5FloatArray("x2f", x2f);
  ReadHDF5FloatArray("x3f", x3f);
  ReadHDF5FloatArray("x1v", x1v);
  ReadHDF5FloatArray("x2v", x2v);
  ReadHDF5FloatArray("x3v", x3v);

  // Read cell data
  ReadHDF5FloatArray("prim", prim);
  ReadHDF5FloatArray("B", bb);

  // Close data file
  data_stream.close();
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to check that needed variables located as expected
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes metadata set.
void AthenaReader::VerifyVariables()
{
  // Check that array of all primitives is present
  int ind_prim;
  int prim_offset = 0;
  for (ind_prim = 0; ind_prim < num_dataset_names; ind_prim++)
    if (dataset_names[ind_prim] == "prim")
      break;
    else
      prim_offset += num_variables(ind_prim);
  if (ind_prim == num_dataset_names)
    throw BlacklightException("Unable to locate array \"prim\" in data file.");

  // Check that all necessary primitives are present
  for (ind_rho = 0; ind_rho < num_variables(ind_prim); ind_rho++)
    if (variable_names[prim_offset+ind_rho] == "rho")
      break;
  if (ind_rho == num_variables(ind_prim))
    throw BlacklightException("Unable to locate \"rho\" slice of \"prim\" in data file.");
  if (plasma_model == PlasmaModel::ti_te_beta)
  {
    for (ind_pgas = 0; ind_pgas < num_variables(ind_prim); ind_pgas++)
      if (variable_names[prim_offset+ind_pgas] == "press")
        break;
    if (ind_pgas == num_variables(ind_prim))
      throw BlacklightException("Unable to locate \"press\" slice of \"prim\" in data file.");
  }
  if (plasma_model == PlasmaModel::code_kappa)
  {
    for (ind_kappa = 0; ind_kappa < num_variables(ind_prim); ind_kappa++)
      if (variable_names[prim_offset+ind_kappa] == simulation_kappa_name)
        break;
    if (ind_kappa == num_variables(ind_prim))
      throw
          BlacklightException("Unable to locate electron entropy slice of \"prim\" in data file.");
  }
  for (ind_uu1 = 0; ind_uu1 < num_variables(ind_prim); ind_uu1++)
    if (variable_names[prim_offset+ind_uu1] == "vel1")
      break;
  if (ind_uu1 == num_variables(ind_prim))
    throw BlacklightException("Unable to locate \"vel1\" slice of \"prim\" in data file.");
  for (ind_uu2 = 0; ind_uu2 < num_variables(ind_prim); ind_uu2++)
    if (variable_names[prim_offset+ind_uu2] == "vel2")
      break;
  if (ind_uu2 == num_variables(ind_prim))
    throw BlacklightException("Unable to locate \"vel2\" slice of \"prim\" in data file.");
  for (ind_uu3 = 0; ind_uu3 < num_variables(ind_prim); ind_uu3++)
    if (variable_names[prim_offset+ind_uu3] == "vel3")
      break;
  if (ind_uu3 == num_variables(ind_prim))
    throw BlacklightException("Unable to locate \"vel3\" slice of \"prim\" in data file.");

  // Check that array of all magnetic field components is present
  int ind_bb;
  int bb_offset = 0;
  for (ind_bb = 0; ind_bb < num_dataset_names; ind_bb++)
    if (dataset_names[ind_bb] == "B")
      break;
    else
      bb_offset += num_variables(ind_bb);
  if (ind_bb == num_dataset_names)
    throw BlacklightException("Unable to locate array \"B\" in data file.");

  // Check that all necessary magnetic field components are present
  for (ind_bb1 = 0; ind_bb1 < num_variables(ind_bb); ind_bb1++)
    if (variable_names[bb_offset+ind_bb1] == "Bcc1")
      break;
  if (ind_bb1 == num_variables(ind_bb))
    throw BlacklightException("Unable to locate \"Bcc1\" slice of \"prim\" in data file.");
  for (ind_bb2 = 0; ind_bb2 < num_variables(ind_bb); ind_bb2++)
    if (variable_names[bb_offset+ind_bb2] == "Bcc2")
      break;
  if (ind_bb2 == num_variables(ind_bb))
    throw BlacklightException("Unable to locate \"Bcc2\" slice of \"prim\" in data file.");
  for (ind_bb3 = 0; ind_bb3 < num_variables(ind_bb); ind_bb3++)
    if (variable_names[bb_offset+ind_bb3] == "Bcc3")
      break;
  if (ind_bb3 == num_variables(ind_bb))
    throw BlacklightException("Unable to locate \"Bcc3\" slice of \"prim\" in data file.");
  return;
}
