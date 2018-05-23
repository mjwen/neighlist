#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
//#include <stdlib.h>
//#include <stdio.h>

#include "neighbor_list.h"

#define DIM 3

void write_extendxyz(int Natoms, double const * cell, double const * coords,
    int const * code)
{

  std::fstream fs;
  fs.open("coords.xyz", std::fstream::out);

  fs << Natoms << std::endl;
  fs << "Lattice=\"";
  for (int i=0; i<9; i++) {
    fs << " " << cell[i];
  }
  fs << "\" Properties=\"species:S:1:pos:R:3 \"" <<std::endl;
  for (int i=0; i<Natoms; i++) {
    fs << code[i] <<" "<< coords[i*3] <<" "<< coords[i*3+1]
      <<" "<< coords[i*3+2] << std::endl;
  }

  fs.close();
}


int main()
{
  double alat = 2.46;
  double d = 3.35;
  double cutoff = d+0.01;
  int pbc[3] = {1,1,1};

  double cell[9] = {alat, 0, 0,   0.5*alat, sqrt(3)*0.5*alat, 0,   0, 0, 2*d};
  double* cell1 = cell;
  double* cell2 = cell+3;
  double* cell3 = cell+6;

  // initialize configurations (AB stacking graphite)
  int Natoms = 4;
  double coords[DIM*Natoms];
  int i=0;
  for (int j=0; j<3; j++) {
    coords[DIM*i+j] = 0*cell1[j] + 0*cell2[j] + 0*cell3[j];
  }
  i=1;
  for (int j=0; j<3; j++) {
    coords[DIM*i+j] = 1/3.*cell1[j] + 1/3.*cell2[j] + 0*cell3[j];
  }
  i=2;
  for (int j=0; j<3; j++) {
    coords[DIM*i+j] = 1/3.*cell1[j] + 1/3.*cell2[j] + 1/2.*cell3[j];
  }
  i=3;
  for (int j=0; j<3; j++) {
    coords[DIM*i+j] = 2/3.*cell1[j] + 2/3.*cell2[j] + 1/2.*cell3[j];
  }
  int code[Natoms];
  code[0] = 1;
  code[1] = 1;
  code[2] = 2;
  code[3] = 2;


  int Npad;
  std::vector<double> pad_coords;
  std::vector<int> pad_code;
  std::vector<int> pad_image;

  /* create padding atoms */
  nbl_set_padding(Natoms, cutoff, cell, pbc, coords, code, Npad, pad_coords,
      pad_code, pad_image);

  int total = Natoms + Npad;

  double * coords_all = new double[total*DIM];
  int * code_all = new int[total*DIM];
  std::memcpy(coords_all, coords, sizeof(double)*Natoms*DIM);
  std::memcpy(coords_all + Natoms*DIM, pad_coords.data(), sizeof(double)*Npad*DIM);
  std::memcpy(code_all, code, sizeof(int)*Natoms);
  std::memcpy(code_all + Natoms, pad_code.data(), sizeof(int)*Npad);


  /* generate neighborlist */
  int * need_neigh = new int[total];
  for (int i=0; i<Natoms; i++) {
    need_neigh[i] = 1;
  }
  for (int i=Natoms; i<total; i++) {
    need_neigh[i] = 0;
  }

  // create neighbor list
  NeighList* nl;
  nbl_initialize(&nl);
  nbl_build(nl, total, cutoff, coords_all, need_neigh);

  // use get neigh
  int request = 0;
  int numnei;
  int * nei1atom;
  nbl_get_neigh(nl, request, &numnei, &nei1atom);

  std::cout<<"get_neigh test"<<std::endl;
  std::cout<<"request = " << request << ", numnei = " << numnei << std::endl;
  for (int i=0; i<numnei; i++) {
    std::cout<<"i = " << i << ", neigh = " << nei1atom[i] << std::endl;
  }

  // print to xyz file
  write_extendxyz(total, cell,  coords_all, code_all);

  // free neighborlist
  nbl_clean_all(&nl);

  // free local data
  delete [] coords_all;
  delete [] code_all;
  delete [] need_neigh;

  return 0;
}
