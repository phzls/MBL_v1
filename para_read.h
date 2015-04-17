//
// Created by Liangsheng Zhang on 4/17/15.
//

/*
 * This file contains functions which deal with reading parameters from files and assign
 * them to parameters
 */

#ifndef MBL_V1_PARA_READ_H
#define MBL_V1_PARA_READ_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>
#include <cstdlib>
#include "parameters.h"

/*
 * Get parameters from the string vector according to keyword. The first argument is the filename,
 * second the content of the file, third the keyword, and the last is the corresponding parameter
 * member. It is assumed that for each inner vector, which corresponds to one line, the first element
 * is the name of the variable, and the second element is the value of that variable.
 */
void para_get(string, const vector<vector<string> >&, string, int&);
void para_get(string, const vector<vector<string> >&, string, bool&);
void para_get(string, const vector<vector<string> >&, string, string&);

/*
 * This functions reads content from a given filename. The string argument is the filename,
 * and the vector of strings passed in is the result. In the file anything behind "//" is
 * assumed to be a comment and ignored, and any line start with "//" is skipped.
 */
void para_file_read(string, vector<vector<string> >&);

#endif //MBL_V1_PARA_READ_H
