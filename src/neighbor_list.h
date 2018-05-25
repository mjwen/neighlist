// Author: Mingjian Wen (wenxx151@umn.edu)

#ifndef NEIGHBOR_LIST_H_
#define NEIGHBOR_LIST_H_

#include<vector>

// neighborlist structure
typedef struct
{
  int numberOfParticles;
  double cutoff;
  int* Nneighbors;
  int* neighborList;
  int* beginIndex;
} NeighList;


void nbl_initialize(NeighList ** const nl);

void nbl_create_paddings(int const numberOfParticles, double const cutoff,
    double const * cell, int const * PBC, double const * coordinates,
    int const * speciesCode, int & numberOfPaddings,
    std::vector<double> & coordinatesOfPaddings,
    std::vector<int> & speciesCodeOfPaddings,
    std::vector<int> & masterOfPaddings);

int nbl_build(NeighList *const nl, int const numberOfParticles, double const cutoff,
    double const * coordinates, int const * needNeighbors);

int nbl_get_neigh(NeighList const * const nl, int const particleNumber,
    int * const numberOfNeighbors, int const ** const neighborsOfParticle);

int nbl_get_neigh_kim(void const * const nl, int const numberOfCutoffs,
    double const * const cutoffs, int const neighborListIndex,
    int const particleNumber, int * const numberOfNeighbors,
    int const ** const neighborsOfParticle);

void nbl_clean(NeighList ** const nl);


// helper funtion
// free contents of NeighList, but not NeighList itseif
void nbl_clean_content(NeighList * const nl);


#endif
