#include <vector>
#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "neighbor_list.h"
#include "padding.h"

#define DIM 3

void write_extendxyz(int Natoms, double const * cell, double const * coords,
    int const * code)
{
  FILE *fp;

  fp = fopen("coords.xyz", "w");


  fprintf(fp, "%d\n",Natoms);
  fprintf(fp, "Lattice=\"");
  for (int i=0; i<9; i++) {
    fprintf(fp, " %f",cell[i]);
  }
  fprintf(fp, "\" Properties=\"species:S:1:pos:R:3 \"\n");

  for (int i=0; i<Natoms; i++) {
    fprintf(fp, "%d %f %f %f\n",code[i], coords[i*3],coords[i*3+1],coords[i*3+2]);
  }
  fclose(fp);
}


int main()
{
  double alat = 1.42;
  double cutoff = 3.01;
  int pbc[3] = {1,1,1};

  double cell[9] = {alat, 0, 0,   0.5*alat, sqrt(3)*0.5*alat, 0,   0, 0, 6};
  double* cell1 = cell;
  double* cell2 = cell+3;
  double* cell3 = cell+6;

  // initialize configurations
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


  std::vector<double> pad_coords;
  std::vector<int> pad_code;
  std::vector<int> pad_image;

  /* create padding atoms */
  set_padding(cell, pbc, cutoff, Natoms, coords, code, pad_coords, pad_code, pad_image);

  int Npad = pad_code.size();
  int total = Natoms + Npad;

  double * coords_all = new double[total*DIM];
  int * code_all = new int[total*DIM];
  std::memcpy(coords_all, coords, sizeof(double)*Natoms*DIM);
  std::memcpy(coords_all + Natoms*DIM, &pad_coords.front(), sizeof(double)*Npad*DIM);
  std::memcpy(code_all, code, sizeof(int)*Natoms);
  std::memcpy(code_all + Natoms, &pad_code.front(), sizeof(int)*Npad);


  /* generate neighborlist */
  int * need_neigh = new int[total];
  for (int i=0; i<Natoms; i++) {
    need_neigh[i] = 1;
  }
  for (int i=Natoms; i<total; i++) {
    need_neigh[i] = 0;
  }



  NeighList* nl;

  nbl_initialize(&nl);

  nbl_build(nl, total, cutoff, coords_all, need_neigh);

  int request = 0;
  int numnei;
  int * nei1atom;
  nbl_get_neigh(nl, request, &numnei, &nei1atom);

  printf("request=%d, numnei=%d\n",request,numnei);
  for (i=0; i<numnei; i++) {
    printf("i=%d, neigh=%d\n",i, nei1atom[i]);
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
