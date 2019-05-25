// Ray Trace multidimensional array header

#ifndef ARRAY_H_
#define ARRAY_H_

//--------------------------------------------------------------------------------------------------

// Multidimensional array
template<typename type>
struct array
{
  // Constructors and destructor
  array();
  array(const array<type> &source) = delete;
  array &operator=(const array<type> &source) = delete;
  array(int n1_);
  ~array();

  // Data
  type *data;
  int n1;
  bool allocated;

  // Functions
  void allocate();
  void allocate(int n);
  type operator()(int n) const;
  type &operator()(int n);
  type operator()(unsigned int n) const;
  type &operator()(unsigned int n);
};

#endif
