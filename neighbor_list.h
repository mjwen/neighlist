// Author: Mingjian Wen (wenxx151@umn.edu)

#ifndef _NEIGHBORLIST_H_
#define _NEIGHBORLIST_H_


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

void nbl_build(NeighList *const nl, int const Natoms, double const cutoff,
    double const * coords, int const * need_neigh);

// free both the contents of NeighList and itseif
void nbl_clean_all(NeighList ** const nl);

int nbl_get_neigh(NeighList const * const nl, int const request, int * const numnei,
    int ** const nei1atom);


// helper function
void coords_to_index(double const * x, int const * size, double const * max,
    double const * min, int * const index);


#endif
