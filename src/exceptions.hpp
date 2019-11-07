// Ray Trace exceptions header

#ifndef EXCEPTION_H_
#define EXCEPTION_H_

// C++ headers
#include <iostream>   // cerr
#include <stdexcept>  // runtime_error
#include <string>     // string

//--------------------------------------------------------------------------------------------------

// Ray Trace exception
struct RayTraceException
  : std::runtime_error
{
  explicit RayTraceException(const char *message)
    : std::runtime_error::runtime_error(ComposeMessage(message)) {}
  std::string ComposeMessage(const char *message)
  {
    std::string message_str(message);
    message_str.insert(0, "Error: ");
    message_str.append("\n");
    return message_str;
  }
};

//--------------------------------------------------------------------------------------------------

// Ray Trace warning
struct RayTraceWarning
{
  explicit RayTraceWarning(const char *message)
  {
    std::string message_str(message);
    message_str.insert(0, "Warning: ");
    message_str.append("\n");
    std::cerr << message_str;
  }
};

#endif
