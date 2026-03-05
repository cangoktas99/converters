#include "openfoam_utility.hpp"

#include <string>
#include <fstream>

inline void SkipFoamHeader(std::istream &in) {
  std::string line;

  /* Scan until the separator line (starts with "// *") */
  while (std::getline(in, line)) {
    if (line.rfind("// *", 0) == 0)
      return;
  }
  throw std::runtime_error("Malformed OpenFOAM file: separator not found");
}

inline std::string DetectFoamClass(const std::string& path) {
 std::ifstream 
}