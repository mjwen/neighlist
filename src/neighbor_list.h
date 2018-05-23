// Author: Mingjian Wen (wenxx151@umn.edu)

#ifndef NEIGHBOR_LIST_H_
#define NEIGHBOR_LIST_H_

#include<vector>

// neighborlist structure
typedef struct
{
  int Natoms;
  double cutoff;
  int* Nneighbors;
  int* neighborList;
  int* beginIndex;
} NeighList;


void nbl_initialize(NeighList ** const nl);

// free contents of NeighList, but not NeighList itseif
void nbl_clean_content(NeighList * const nl);

// free both the contents of NeighList and itseif
void nbl_clean_all(NeighList ** const nl);

int nbl_get_neigh(NeighList const * const nl, int const request, int * const numnei,
    int ** const nei1atom);

void nbl_build(NeighList *const nl, int const Natoms, double const cutoff,
    double const * coords, int const * need_neigh);

void nbl_set_padding(int const Natoms, double const cutoff, double const * cell,
    int const * PBC, double const * coords, int const * species,
    int & Npad, std::vector<double> & pad_coords, std::vector<int> & pad_species,
    std::vector<int> & pad_image);

#endif
