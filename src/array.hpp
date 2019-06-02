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
  array(int n1_);
  array(int n2_, int n1_);
  array(int n3_, int n2_, int n1_);
  array(int n4_, int n3_, int n2_, int n1_);
  array(int n5_, int n4_, int n3_, int n2_, int n1_);
  array(const array<type> &source);
  array &operator=(const array<type> &source);
  ~array();

  // Data
  type *data;
  int n1, n2, n3, n4, n5;
  int n_tot;
  bool allocated = false;
  bool is_copy = false;

  // Functions
  void allocate();
  void allocate(int n1_);
  void allocate(int n2_, int n1_);
  void allocate(int n3_, int n2_, int n1_);
  void allocate(int n4_, int n3_, int n2_, int n1_);
  void allocate(int n5_, int n4_, int n3_, int n2_, int n1_);
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
};

#endif
