// Ray Trace exceptions header

#ifndef EXCEPTION_H_
#define EXCEPTION_H_

// C++ headers
#include <iostream>   // cerr
#include <stdexcept>  // runtime_error
#include <string>     // string

//--------------------------------------------------------------------------------------------------

// Ray Trace exception
struct ray_trace_exception
  : std::runtime_error
{
  explicit ray_trace_exception(const char *message)
    : std::runtime_error::runtime_error(compose_message(message)) {}
  std::string compose_message(const char *message)
  {
    std::string message_str(message);
    message_str.insert(0, "Error: ");
    message_str.append("\n");
    return message_str;
  }
};

//--------------------------------------------------------------------------------------------------

// Ray Trace warning
struct ray_trace_warning
{
  explicit ray_trace_warning(const char *message)
  {
    std::string message_str(message);
    message_str.insert(0, "Warning: ");
    message_str.append("\n");
    std::cerr << message_str;
  }
};

#endif
