// Author: Mingjian Wen (wenxx151@umn.edu)

#ifndef PADDING_H_
#define PADDING_H_

#include<vector>

void set_padding(const double* cell, const int* PBC, const double cutoff,
    const int Natoms, const double* coords, const int* species,
    std::vector<double>& pad_coords, std::vector<int>& pad_species,
    std::vector<int>& pad_image);

#endif
