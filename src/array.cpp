// Ray Trace multidimensional array

// Ray Trace headers
#include "array.hpp"
#include "exceptions.hpp"  // ray_trace_exception

// Instantiations
template struct array<int>;
template struct array<float>;

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
// Inputs:
//   n1_: size of only dimension
template<typename type>
array<type>::array(int n1_)
  : n1(n1_),
    n2(1),
    n3(1),
    n4(1),
    n5(1)
{
  allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array constructor (2D)
// Inputs:
//   n2_: size of outermost dimension
//   n1_: size of innermost dimension
template<typename type>
array<type>::array(int n2_, int n1_)
  : n1(n1_),
    n2(n2_),
    n3(1),
    n4(1),
    n5(1)
{
  allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array constructor (3D)
// Inputs:
//   n3_: size of outermost dimension
//   n2_: size of intermediate dimension
//   n1_: size of innermost dimension
template<typename type>
array<type>::array(int n3_, int n2_, int n1_)
  : n1(n1_),
    n2(n2_),
    n3(n3_),
    n4(1),
    n5(1)
{
  allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array constructor (4D)
// Inputs:
//   n4_: size of outermost dimension
//   n3_, n2_: sizes of intermediate dimensions
//   n1_: size of innermost dimension
template<typename type>
array<type>::array(int n4_, int n3_, int n2_, int n1_)
  : n1(n1_),
    n2(n2_),
    n3(n3_),
    n4(n4_),
    n5(1)
{
  allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array constructor (5D)
// Inputs:
//   n5_: size of outermost dimension
//   n4_, n3_, n2_: sizes of intermediate dimensions
//   n1_: size of innermost dimension
template<typename type>
array<type>::array(int n5_, int n4_, int n3_, int n2_, int n1_)
  : n1(n1_),
    n2(n2_),
    n3(n3_),
    n4(n4_),
    n5(n5_)
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
  n_tot = n1 * n2 * n3 * n4 * n5;
  if (n_tot <= 0)
    throw ray_trace_exception("Error: Attempting to allocate empty array.");
  data = new type[n_tot];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (1D)
// Inputs:
//   n1_: size of only dimension
// Outputs: (none)
template<typename type>
void array<type>::allocate(int n1_)
{
  n1 = n1_;
  n2 = 1;
  n3 = 1;
  n4 = 1;
  n5 = 1;
  allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (2D)
// Inputs:
//   n2_: size of outermost dimension
//   n1_: size of innermost dimension
// Outputs: (none)
template<typename type>
void array<type>::allocate(int n2_, int n1_)
{
  n1 = n1_;
  n2 = n2_;
  n3 = 1;
  n4 = 1;
  n5 = 1;
  allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (3D)
// Inputs:
//   n3_: size of outermost dimension
//   n2_: size of intermediate dimension
//   n1_: size of innermost dimension
// Outputs: (none)
template<typename type>
void array<type>::allocate(int n3_, int n2_, int n1_)
{
  n1 = n1_;
  n2 = n2_;
  n3 = n3_;
  n4 = 1;
  n5 = 1;
  allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (4D)
// Inputs:
//   n4_: size of outermost dimension
//   n3_, n2_: sizes of intermediate dimensions
//   n1_: size of innermost dimension
// Outputs: (none)
template<typename type>
void array<type>::allocate(int n4_, int n3_, int n2_, int n1_)
{
  n1 = n1_;
  n2 = n2_;
  n3 = n3_;
  n4 = n4_;
  n5 = 1;
  allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (5D)
// Inputs:
//   n5_: size of outermost dimension
//   n4_, n3_, n2_: sizes of intermediate dimensions
//   n1_: size of innermost dimension
// Outputs: (none)
template<typename type>
void array<type>::allocate(int n5_, int n4_, int n3_, int n2_, int n1_)
{
  n1 = n1_;
  n2 = n2_;
  n3 = n3_;
  n4 = n4_;
  n5 = n5_;
  allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read accessors (1D)
// Inputs:
//   i1: only index
// Outputs:
//   returned value: element
template<typename type>
type array<type>::operator()(int i1) const
{
  return data[i1];
}
template<typename type>
type array<type>::operator()(unsigned int i1) const
{
  return data[i1];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read accessors (2D)
// Inputs:
//   i2: outermost index
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type>
type array<type>::operator()(int i2, int i1) const
{
  return data[i1 + n1 * i2];
}
template<typename type>
type array<type>::operator()(unsigned int i2, unsigned int i1) const
{
  return data[i1 + static_cast<unsigned int>(n1) * i2];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read accessors (3D)
// Inputs:
//   i3: outermost index
//   i2: intermediate index
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type>
type array<type>::operator()(int i3, int i2, int i1) const
{
  return data[i1 + n1 * (i2 + n2 * i3)];
}
template<typename type>
type array<type>::operator()(unsigned int i3, unsigned int i2, unsigned int i1) const
{
  return data[i1 + static_cast<unsigned int>(n1) * (i2 + static_cast<unsigned int>(n2) * i3)];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read accessors (4D)
// Inputs:
//   i4: outermost index
//   i3, i2: intermediate indices
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type>
type array<type>::operator()(int i4, int i3, int i2, int i1) const
{
  return data[i1 + n1 * (i2 + n2 * (i3 + n3 * i4))];
}
template<typename type>
type array<type>::operator()(unsigned int i4, unsigned int i3, unsigned int i2, unsigned int i1)
    const
{
  return data[i1 + static_cast<unsigned int>(n1) * (i2 + static_cast<unsigned int>(n2)
      * (i3 + static_cast<unsigned int>(n3) * i4))];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read accessors (5D)
// Inputs:
//   i5: outermost index
//   i4, i3, i2: intermediate indices
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type>
type array<type>::operator()(int i5, int i4, int i3, int i2, int i1) const
{
  return data[i1 + n1 * (i2 + n2 * (i3 + n3 * (i4 + n4 * i5)))];
}
template<typename type>
type array<type>::operator()(unsigned int i5, unsigned int i4, unsigned int i3, unsigned int i2,
    unsigned int i1) const
{
  return data[i1 + static_cast<unsigned int>(n1) * (i2 + static_cast<unsigned int>(n2)
      * (i3 + static_cast<unsigned int>(n3) * (i4 + static_cast<unsigned int>(n4) * i5)))];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read/write accessors (1D)
// Inputs:
//   i1: only index
// Outputs:
//   returned value: element
template<typename type>
type &array<type>::operator()(int i1)
{
  return data[i1];
}
template<typename type>
type &array<type>::operator()(unsigned int i1)
{
  return data[i1];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read/write accessors (2D)
// Inputs:
//   i2: outermost index
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type>
type &array<type>::operator()(int i2, int i1)
{
  return data[i1 + n1 * i2];
}
template<typename type>
type &array<type>::operator()(unsigned int i2, unsigned int i1)
{
  return data[i1 + static_cast<unsigned int>(n1) * i2];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read/write accessors (3D)
// Inputs:
//   i3: outermost index
//   i2: intermediate index
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type>
type &array<type>::operator()(int i3, int i2, int i1)
{
  return data[i1 + n1 * (i2 + n2 * i3)];
}
template<typename type>
type &array<type>::operator()(unsigned int i3, unsigned int i2, unsigned int i1)
{
  return data[i1 + static_cast<unsigned int>(n1) * (i2 + static_cast<unsigned int>(n2) * i3)];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read/write accessors (4D)
// Inputs:
//   i4: outermost index
//   i3, i2: intermediate indices
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type>
type &array<type>::operator()(int i4, int i3, int i2, int i1)
{
  return data[i1 + n1 * (i2 + n2 * (i3 + n3 * i4))];
}
template<typename type>
type &array<type>::operator()(unsigned int i4, unsigned int i3, unsigned int i2, unsigned int i1)
{
  return data[i1 + static_cast<unsigned int>(n1) * (i2 + static_cast<unsigned int>(n2)
      * (i3 + static_cast<unsigned int>(n3) * i4))];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read/write accessors (5D)
// Inputs:
//   i5: outermost index
//   i4, i3, i2: intermediate indices
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type>
type &array<type>::operator()(int i5, int i4, int i3, int i2, int i1)
{
  return data[i1 + n1 * (i2 + n2 * (i3 + n3 * (i4 + n4 * i5)))];
}
template<typename type>
type &array<type>::operator()(unsigned int i5, unsigned int i4, unsigned int i3, unsigned int i2,
    unsigned int i1)
{
  return data[i1 + static_cast<unsigned int>(n1) * (i2 + static_cast<unsigned int>(n2)
      * (i3 + static_cast<unsigned int>(n3) * (i4 + static_cast<unsigned int>(n4) * i5)))];
}
