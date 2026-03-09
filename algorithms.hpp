#include <algorithm>
#include <cctype>
#include <string>

inline bool CaseInsensitiveLess(const std::string &a, const std::string &b) {
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end(),
                                      [](char ac, char bc) {
                                        return std::tolower(ac) < std::tolower(bc);
                                      }
  );
}
