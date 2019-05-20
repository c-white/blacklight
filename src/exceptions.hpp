// Ray Trace exceptions header

#ifndef EXCEPTION_H_
#define EXCEPTION_H_

// C++ headers
#include <stdexcept>  // runtime_error

// Ray Trace exception
struct ray_trace_exception
  : std::runtime_error
{
  explicit ray_trace_exception(const char *message)
    : std::runtime_error::runtime_error(message) {}
};

#endif
