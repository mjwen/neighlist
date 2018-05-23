// Author: Mingjian Wen (wenxx151@umn.edu)

#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include "neighbor_list.h"


#define DIM 3
#define TOL 1.0e-10

#define MY_ERROR(message)                                             \
  {                                                                   \
    std::cout << "* Error (Neighbor List): \"" << message << "\" : "  \
              << __LINE__ << ":" << __FILE__ << std::endl;            \
    exit(1);                                                          \
  }

#define MY_WARNING(message)                                            \
  {                                                                    \
    std::cout << "* Error (Neighbor List) : \"" << message << "\" : "  \
              << __LINE__ << ":" << __FILE__ << std::endl;             \
  }


void nbl_initialize(NeighList ** const nl)
{
  *nl = new NeighList;
  (*nl)->Natoms = -1;
  (*nl)->cutoff = -1.0;
  (*nl)->Nneighbors = 0;
  (*nl)->neighborList = 0;
  (*nl)->beginIndex = 0;
}

void nbl_clean_content(NeighList * const nl)
{
  if (nl != 0) {
    delete[] nl->Nneighbors;
    delete[] nl->neighborList;
    delete[] nl->beginIndex;
    nl->Natoms = -1;
    nl->cutoff = -1.0;
    nl->Nneighbors = 0;
    nl->neighborList = 0;
    nl->beginIndex = 0;
  }

}

void nbl_clean_all(NeighList ** const nl)
{
  // free content
  nbl_clean_content(*nl);

  // free NeighList
  delete (*nl);

  // nullify pointer
  (*nl) = 0;
}


void nbl_build(NeighList *const nl, int const Natoms, double const cutoff,
    double const * coords, int const * need_neigh)
{


  // find max and min extend of coords
  double min[DIM];
  double max[DIM];

  // init max and min of coords to that of the first atom
  for (int k=0; k<DIM; k++){
    min[k] = coords[k];
    max[k] = coords[k] + 1.0; // +1 to prevent max==min for 1D and 2D case
  }
  for (int i=0; i<Natoms; i++){
    for (int j=0; j<DIM; j++){
      if (max[j] < coords[DIM*i+j]) max[j] = coords[DIM*i+j];
      if (min[j] > coords[DIM*i+j]) min[j] = coords[DIM*i+j];
    }
  }

  // make the cell box
  int size_total = 1;
  int size[DIM];
  for (int i=0; i<DIM; i++){
    size[i] = static_cast<int> ((max[i]-min[i])/cutoff);
    size[i] = size[i] <= 0 ? 1 : size[i];
    size_total *= size[i];
  }
  if (size_total > 1000000000) {
    MY_ERROR("Cell size too large. Check if you have partilces fly away.");
  }

  // assign atoms into cells
  std::vector<std::vector<int> > cells(size_total);
  for (int i=0; i<Natoms; i++){
    int index[DIM];
    coords_to_index(&coords[DIM*i], size, max, min, index);
    int idx = index[0] + index[1]*size[0] + index[2]*size[0]*size[1];
    cells[idx].push_back(i);
  }


  // create neighbors

  // free previous neigh content first
  nbl_clean_content(nl);
  nl->Nneighbors = new int[Natoms];
  nl->beginIndex = new int[Natoms];

  // temporary neigh container
  std::vector<int> tmp_neigh;

  double cutsq = cutoff*cutoff;
  int total = 0;

  for (int i=0; i<Natoms; i++) {
    int num_neigh = 0;
    if (need_neigh[i]) {
      int index[DIM];
      coords_to_index(&coords[DIM*i], size, max, min, index);

      // loop over neighborling cells and the cell atom i resides
      for (int ii=std::max(0, index[0]-1); ii<=std::min(index[0]+1, size[0]-1); ii++)
      for (int jj=std::max(0, index[1]-1); jj<=std::min(index[1]+1, size[1]-1); jj++)
      for (int kk=std::max(0, index[2]-1); kk<=std::min(index[2]+1, size[2]-1); kk++)
      {
        int idx = ii + jj*size[0] + kk*size[0]*size[1];

        for (size_t m=0; m<cells[idx].size(); m++) {
          int n = cells[idx][m];
          if (n != i) {
            double rsq = 0.0;

            for (int k=0; k<DIM; k++) {
              double del = coords[DIM*n+k] - coords[DIM*i+k];
              rsq += del*del;
            }
            if (rsq < TOL) {
              std::ostringstream stringStream;
              stringStream <<"Collision of atoms "<<i+1<<" and "<<n+1<<". ";
              stringStream  <<"Their distance is "<<std::sqrt(rsq)<<"."<<std::endl;
              std::string my_str = stringStream.str();
              MY_ERROR(my_str);
            }
            if (rsq < cutsq) {
              tmp_neigh.push_back(n);
              num_neigh++;
            }
          }
        }
      }
    }

    nl->Nneighbors[i] = num_neigh;
    nl->beginIndex[i] = total;
    total += num_neigh;
  }

  nl->Natoms = Natoms;
  nl->cutoff = cutoff;
  nl->neighborList = new int[total];
  std::memcpy(nl->neighborList, tmp_neigh.data(), sizeof(int)*total);

}


int nbl_get_neigh(NeighList const * const nl, int const request, int * const numnei,
    int ** const nei1atom)
{

  if ((request >= nl->Natoms) || (request < 0)) {
    MY_WARNING("atom ID out of bound");
    return 1;
  }

  // number of neighbors
  *numnei = nl->Nneighbors[request];

  // neighbor list starting point
  int idx = nl->beginIndex[request];
  *nei1atom = nl->neighborList + idx;

  return 0;
}


void coords_to_index(double const * x, int const * size, double const * max,
    double const * min, int * const index)
{
  for (int i=0; i<DIM; i++) {
    index[i] = static_cast<int> (((x[i]-min[i])/(max[i]-min[i])) * size[i]);
    index[i] = std::min(index[i], size[i]-1);  // handle edge case when x[i] = max[i]
  }
}

