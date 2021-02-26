// Blacklight exceptions header

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

// C++ headers
#include <iostream>   // cerr
#include <stdexcept>  // runtime_error
#include <string>     // string

//--------------------------------------------------------------------------------------------------

// Blacklight exception
struct BlacklightException
  : std::runtime_error
{
  explicit BlacklightException(const char *message)
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

// Blacklight warning
struct BlacklightWarning
{
  explicit BlacklightWarning(const char *message)
  {
    std::string message_str(message);
    message_str.insert(0, "Warning: ");
    message_str.append("\n");
    std::cerr << message_str;
  }
};

#endif
