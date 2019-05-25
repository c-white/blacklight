// Ray Trace multidimensional array

// Ray Trace headers
#include "array.hpp"
#include "exceptions.hpp"  // ray_trace_exception

// Instantiations
template struct array<int>;

//--------------------------------------------------------------------------------------------------

// Multidimensional array constructor (empty)
// Inputs: (none)
template<typename type>
array<type>::array()
{
  allocated = false;
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array constructor (1D)
// Inputs: (none)
template<typename type>
array<type>::array(int n1_)
  : n1(n1_)
{
  allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array destructor
template<typename type>
array<type>::~array()
{
  if (allocated)
    delete[] data;
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (general)
// Inputs: (none)
// Outputs: (none)
template<typename type>
void array<type>::allocate()
{
  if (allocated)
    throw ray_trace_exception("Error: Attempting to reallocate array.");
  if (n1 == 0)
    throw ray_trace_exception("Error: Attempting to allocate empty array.");
  else
    data = new type[n1];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (1D)
// Inputs:
//   n1_: number of elements
// Outputs: (none)
template<typename type>
void array<type>::allocate(int n1_)
{
  if (allocated)
    throw ray_trace_exception("Error: Attempting to reallocate array.");
  n1 = n1_;
  data = new type[n1];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read accessors (1D)
// Inputs:
//   n: index
// Outputs:
//   returned value: element
template<typename type>
type array<type>::operator()(int n) const
{
  return data[n];
}
template<typename type>
type array<type>::operator()(unsigned int n) const
{
  return data[n];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read/write accessors (1D)
// Inputs:
//   n: index
// Outputs:
//   returned value: reference to element
template<typename type>
type &array<type>::operator()(int n)
{
  return data[n];
}
template<typename type>
type &array<type>::operator()(unsigned int n)
{
  return data[n];
}
