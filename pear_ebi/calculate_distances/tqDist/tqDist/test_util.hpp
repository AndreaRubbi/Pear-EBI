#ifndef TEST_UTIL_HPP
#define TEST_UTIL_HPP

#include <string>
#include <iostream>
#include <cstdlib>

template <typename T>
void assert_equal(const T &val, const T &exp, const std::string msg = "Assertion failed!")  {
  if(val != exp) {
    std::cerr << msg << " Found " << val << " but expected " << exp << std::endl;
    std::exit(-1);
  }
}

#endif
