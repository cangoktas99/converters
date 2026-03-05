#ifndef OPENFOAM_UTILITY_H
#define OPENFOAM_UTILITY_H

#include <istream>
#include <ostream>

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
std::string DetectFoamClass(const std::string &path);



#endif /* OPENFOAM_UTILITY_H */
