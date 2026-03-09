#include "openfoam_utility.hpp"

#include "algorithms.hpp"

inline void SkipFoamHeader(std::istream &in) {
  std::string line;

  /* Scan until the separator line (starts with "// *") */
  while (std::getline(in, line)) {
    if (line.rfind("// *", 0) == 0)
      return;
  }
  throw std::runtime_error("Malformed OpenFOAM file: separator not found");
}

inline std::string DetectFoamClass(const fs::path& path) {
  std::ifstream in(path.string());

  if (!in)
    return "";

  std::string line;
  bool in_foam_file = false;

  while (std::getline(in, line)) {
    if (line.find("FoamFile") != std::string::npos) {
      in_foam_file = true;
      continue;
    }
    if (in_foam_file && line.find("class") != std::string::npos) {
      std::istringstream ss(line);
      std::string key, val;
      ss >> key >> val;
      if (!val.empty() && val.back() == ';')
        val.pop_back();
      return val;
    }
    if (line.rfind("// *", 0) == 0) break;
  }
  return "";
}

inline std::vector<fs::path> DiscoverFields(const fs::path& time_dir) {
  std::vector<fs::path> names;

  if (not fs::is_directory(time_dir))
    return names;

  for (auto const &dir_entry : fs::directory_iterator{time_dir}) {
    const fs::path file_path = dir_entry.path();
    if (not fs::exists(file_path) || fs::is_directory(file_path)) 
      continue;
    
    std::string foam_class = DetectFoamClass(file_path);
    if (foam_class == "volScalarField" || foam_class == "volVectorField")
      names.push_back(file_path);
  }
  std::sort(names.begin(), names.end(), CaseInsensitiveLess);
  return names;
}

inline void SkipWhiteSpace(std::istream &in) {
  while (1) {
    in >> std::ws;

    if (in.peek() != '/')
      return;

    in.get();

    if (in.peek() == '/') {
      in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    } else {
      in.unget();
      return;
    }
  }
}

inline size_t ReadCount(std::istream &in) {
  size_t n;
  SkipWhiteSpace(in);
  if (not(in >> n))
    throw std::runtime_error("Expected integer count");
  return n;
}

inline void Expect(std::istream &in, char c) {
  SkipWhiteSpace(in);
  char got;
  if (not(in >> got) || got != c) 
    throw std::runtime_error(std::format("Expected '{}'", c));
}

template <class Scalar>
inline typename PM_Util<Scalar>::PntV
PM_Util<Scalar>::ReadPoints(const fs::path &file_path) {
  size_t n{};
  auto in = OpenFile(file_path, n);

  PntV pnts(n, {0., 0., 0.});

  for (size_t i = 0; i < n; ++i) {
    Expect(in, '(');
    in >> pnts[i][0] >> pnts[i][1] >> pnts[i][2];
    Expect(in, ')');
  }
  CloseFile(in);
  return pnts;
}

template <class Scalar>
inline typename PM_Util<Scalar>::FaceV
PM_Util<Scalar>::ReadFaces(const fs::path &file_path) {
  size_t n{};
  auto in = OpenFile(file_path, n);

  FaceV faces(n);
  for (size_t i = 0; i < n; ++i) {
    SkipWhiteSpace(in);
    size_t nv;
    if (not(in >> nv))
      throw std::runtime_error(
          std::format("Expected vertex count for face {}", i));
    Expect(in, '(');
    faces[i].resize(nv);
    for (size_t j = 0; j < nv; j++) {
      in >> faces[i][j];
    }
    Expect(in, ')');
  }
  CloseFile(in);
  return faces;
}

int main() {
  fs::path dir = "LS_TestOFCase/0";
  auto vec = DiscoverFields(dir);

  if (not vec.empty()) {
    for (size_t i = 0; i < vec.size(); ++i) {
      std::cout << vec[i] << "\n";
    }
  }

  for (const char *str : {"      #1 test", "\t #2 test", "#3 test"}) {
    std::istringstream iss{std::string(str)};
    SkipWhiteSpace(iss);

    std::string rest;
    // std::getline(iss, rest);
    std::cout << iss.rdbuf() << "\n";
  }

  fs::path points_dir = "LS_TestOFCase/constant/polyMesh/points";
  auto pnts =
      PM_Util<DefaultScalar>::ReadPoints(fs::canonical(points_dir));

  std::cout << std::format("Points size: {}\n", pnts.size());

  std::cout << std::format("pnts[0][0]: {:+3.10f}, "
                           "pnts[0][1]: {:+3.10f}, "
                           "pnts[0][2]: {:+3.10f}\n",
                           pnts[0][0], pnts[0][1], pnts[0][2]);

  std::cout << std::format("pnts[-1][0]: {:+3.10f}, "
                           "pnts[-1][1]: {:+3.10f}, "
                           "pnts[-1][2]: {:+3.10f}\n",
                           pnts.back()[0], pnts.back()[1], pnts.back()[2]);

  fs::path faces_dir = "LS_TestOFCase/constant/polyMesh/faces";
  auto faces =
      PM_Util<DefaultScalar>::ReadFaces(fs::canonical(faces_dir));

  std::cout << std::format("Faces size: {}\n", faces.size());

  std::cout << std::format("faces[0][0]: {}, "
                           "faces[0][1]: {}, "
                           "faces[0][2]: {}\n",
                           faces[0][0], faces[0][1], faces[0][2]);

  std::cout << std::format("faces[-1][0]: {}, "
                           "faces[-1][1]: {}, "
                           "faces[-1][2]: {}\n",
                           faces.back()[0], faces.back()[1], faces.back()[2]);
}