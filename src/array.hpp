// Ray Trace multidimensional array header

#ifndef ARRAY_H_
#define ARRAY_H_

// C++ headers
#include <cstddef>  // size_t

//--------------------------------------------------------------------------------------------------

// Multidimensional array
template<typename type>
struct Array
{
  // Constructors and destructor
  Array();
  Array(int n1_);
  Array(int n2_, int n1_);
  Array(int n3_, int n2_, int n1_);
  Array(int n4_, int n3_, int n2_, int n1_);
  Array(int n5_, int n4_, int n3_, int n2_, int n1_);
  Array(const Array<type> &source);
  Array &operator=(const Array<type> &source);
  ~Array();

  // Data
  type *data;
  int n1, n2, n3, n4, n5;
  int n_tot;
  bool allocated = false;
  bool is_copy = false;
  static constexpr int max_dims = 5;

  // Functions - allocators and deallocator
  void Allocate();
  void Allocate(int n1_);
  void Allocate(int n2_, int n1_);
  void Allocate(int n3_, int n2_, int n1_);
  void Allocate(int n4_, int n3_, int n2_, int n1_);
  void Allocate(int n5_, int n4_, int n3_, int n2_, int n1_);
  void Deallocate();

  // Functions - read accessors
  type operator()(int i1_) const;
  type operator()(unsigned int i1_) const;
  type operator()(int i2_, int i1_) const;
  type operator()(unsigned int i2_, unsigned int i1_) const;
  type operator()(int i3_, int i2_, int i1_) const;
  type operator()(unsigned int i3_, unsigned int i2_, unsigned int i1_) const;
  type operator()(int i4_, int i3_, int i2_, int i1_) const;
  type operator()(unsigned int i4_, unsigned int i3_, unsigned int i2_, unsigned int i1_) const;
  type operator()(int i5_, int i4_, int i3_, int i2_, int i1_) const;
  type operator()(unsigned int i5_, unsigned int i4_, unsigned int i3_, unsigned int i2_,
      unsigned int i1_) const;

  // Functions - read/write accessors
  type &operator()(int i1_);
  type &operator()(unsigned int i1_);
  type &operator()(int i2_, int i1_);
  type &operator()(unsigned int i2_, unsigned int i1_);
  type &operator()(int i3_, int i2_, int i1_);
  type &operator()(unsigned int i3_, unsigned int i2_, unsigned int i1_);
  type &operator()(int i4_, int i3_, int i2_, int i1_);
  type &operator()(unsigned int i4_, unsigned int i3_, unsigned int i2_, unsigned int i1_);
  type &operator()(int i5_, int i4_, int i3_, int i2_, int i1_);
  type &operator()(unsigned int i5_, unsigned int i4_, unsigned int i3_, unsigned int i2_,
      unsigned int i1_);

  // Functions - miscellaneous
  void Slice(int dimension, int index);
  void Zero();
  std::size_t GetNumBytes() const;
};

#endif
