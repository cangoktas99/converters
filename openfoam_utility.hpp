#ifndef OPENFOAM_UTILITY_H
#define OPENFOAM_UTILITY_H


#include <istream>
#include <fstream>
#include <vector>
#include <string>
#include <array>

#include <filesystem>
namespace fs = std::filesystem;

#include <format>
#include <iostream>

/**
 * @brief Skip FoamFile header block and the "// * * *" separator line.
 *        Positions the stream just after that line, ready to read data.
 * 
 */
void SkipFoamHeader(std::istream &in);

/**
 * @brief Read the FoamFile header and return the value of the 'class' field.
 *        Returns "" if the file cannot be opened or does not have a FoamFile
 *        header.
 * 
 */
std::string DetectFoamClass(const fs::path&);

/**
 * @brief Scan <fieldDir> for files whose FoamFile class is volScalarField or
 *  volVectorField and return their names sorted alphabetically.
 *
 */
std::vector<fs::path> DiscoverFields(const fs::path&);

/**
 * @brief Skip whitespace and comments ('//' to end-of-line)
 * 
 */
void SkipWhiteSpace(std::istream&);

/**
 * @brief Read the integer count that appears before the opening '('
 * 
 */
size_t ReadCount(std::istream&);

/**
 * @brief Expect and consume a single character
 * 
 */
void Expect(std::istream &, char);

/**
 * @brief PolyMesh utilities
 *
 */
template <typename ScalarType = double>
class PM_Util {
  using Scalar = ScalarType;

  static std::ifstream OpenFile(const fs::path &file_path, size_t& n) {
    const std::string file_path_str = file_path.string();
    std::cout << std::format("Reading the file: {}\n", file_path_str);

    std::ifstream in(file_path_str);
    if (not in)
      throw std::runtime_error(std::format("Cannot open: {}", file_path_str));
    SkipFoamHeader(in);

    n = ReadCount(in);
    if (not n)
      throw std::runtime_error(std::format("Number of entity cannot be zero: {}", file_path_str));
    Expect(in, '(');

    return in;
  }

  static void CloseFile(std::ifstream &in) {
    Expect(in, ')');
    in.close();
  }
public:
  /**
   * @brief Type for Point vector
   * 
   */
  using PntV = std::vector<std::array<Scalar, 3>>;

  /**
   * @brief Type for Face vector
   * 
   */
  using FaceV = std::vector<std::vector<size_t>>;

  using LabelV = std::vector<size_t>;

  /**
   * @brief Read points: returns vector of (x,y,z)
   *
   */
  static PntV ReadPoints(const fs::path&);

  /**
   * @brief Read faces: returns vector of face vertex lists (0-based)
   * 
   */
  static FaceV ReadFaces(const fs::path &);

  /**
   * @brief Read label list (owner or neighbour)
   * Also extracts nCells and nInternalFaces from the 'note' field if present
   *
   * @return LabelV 
   */
  static LabelV ReadLabelList(const fs::path &, size_t& n_cells_out, size_t* n_int_faces_out);
};

using DefaultScalar = double;

#endif /* OPENFOAM_UTILITY_H */
