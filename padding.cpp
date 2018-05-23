// Author: Mingjian Wen (wenxx151@umn.edu)

#include <stdlib.h>
#include <cmath>
#include <vector>
#include <iostream>
#include "padding.h"


#define DIM 3
#define SMALL 1.0e-10

#define MY_ERROR(message)                                             \
  {                                                                   \
    std::cout << "* Error (Neighbor List): \"" << message << "\" : "  \
              << __LINE__ << ":" << __FILE__ << std::endl;            \
    exit(1);                                                          \
  }

#define MY_WARNING(message)                                              \
  {                                                                      \
    std::cout << "* Error (Neighbor List) : \"" << message << "\" : "  \
              << __LINE__ << ":" << __FILE__ << std::endl;               \
  }

// norm of a 3-element vector
inline double norm(const double* a) {
  return std::sqrt(a[0]*a[0]+ a[1]*a[1]+ a[2]*a[2]);
}

// dot product of two 3-element vectors
inline double dot(const double *a, const double *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// cross product of two 3-element vectors
inline void cross(const double *a, const double *b, double* axb)
{
  axb[0] = a[1]*b[2] - a[2]*b[1];
  axb[1] = a[2]*b[0] - a[0]*b[2];
  axb[2] = a[0]*b[1] - a[1]*b[0];
}

// determinant of a 3 by 3 matrix
inline double det(const double *mat)
{
  return (mat[0]*mat[4]*mat[8] - mat[0]*mat[5]*mat[7]
        - mat[1]*mat[3]*mat[8] + mat[1]*mat[5]*mat[6]
        + mat[2]*mat[3]*mat[7] - mat[2]*mat[4]*mat[6]);
}

// determinant of a 2 by 2 matrix
inline double det2(double a11, double a12, double a21, double a22)
{
  return (a11*a22) - (a12*a21);
}

// transpose of a DIM by DIM matrix
inline void transpose(const double *mat, double* trans)
{
  for (int i=0; i<DIM; i++) {
    for (int j=0; j<DIM; j++) {
      trans[DIM*i+j] = mat[DIM*j+i];
    }
  }
}

//inverse of a 3 by 3 matrix
inline void inverse(const double *mat, double *inv)
{
  inv[0] = det2(mat[4], mat[5], mat[7], mat[8]);
  inv[1] = det2(mat[2], mat[1], mat[8], mat[7]);
  inv[2] = det2(mat[1], mat[2], mat[4], mat[5]);
  inv[3] = det2(mat[5], mat[3], mat[8], mat[6]);
  inv[4] = det2(mat[0], mat[3], mat[6], mat[8]);
  inv[5] = det2(mat[2], mat[0], mat[5], mat[3]);
  inv[6] = det2(mat[3], mat[4], mat[6], mat[7]);
  inv[7] = det2(mat[1], mat[0], mat[7], mat[6]);
  inv[8] = det2(mat[0], mat[1], mat[3], mat[4]);

  double dd = det(mat);
  if(std::abs(dd) < SMALL) {
    MY_ERROR("Cannot invert cell matrix. Determinant is 0.");
  }
  for (int i=0; i<9; i++) {
    inv[i] /= dd;
  }
}


void set_padding(const double* cell, const int* PBC, const double cutoff,
    const int Natoms, const double* coords, const int* species,
    std::vector<double>& pad_coords, std::vector<int>& pad_species,
    std::vector<int>& pad_image)
{

  // transform coords into fractional coords
  double tcell[9];
  double fcell[9];
  transpose(cell, tcell);
  inverse(tcell, fcell);

  double frac_coords[DIM*Natoms];
  double min[DIM] = {1e10, 1e10, 1e10};
  double max[DIM] = {-1e10, -1e10, -1e10};
  for (int i=0; i<Natoms; i++) {
    const double* atom_coords = coords + (DIM*i);
    double x = dot(fcell, atom_coords);
    double y = dot(fcell+3, atom_coords);
    double z = dot(fcell+6, atom_coords);
    frac_coords[DIM*i+0] = x;
    frac_coords[DIM*i+1] = y;
    frac_coords[DIM*i+2] = z;
    if (x < min[0]) min[0] = x;
    if (y < min[1]) min[1] = y;
    if (z < min[2]) min[2] = z;
    if (x > max[0]) max[0] = x;
    if (y > max[1]) max[1] = y;
    if (z > max[2]) max[2] = z;
  }

  // add some extra value to deal with edge case
  for (int i=0; i<DIM; i++) {
    min[i] -= 1e-10;
    max[i] += 1e-10;
  }

    // volume of cell
  double xprod[DIM];
  cross(cell+3, cell+6, xprod);
  double volume = std::abs(dot(cell, xprod));

  // distance between parallelpiped cell faces
  double dist[DIM];
  cross(cell+3, cell+6, xprod);
  dist[0] = volume/norm(xprod);
  cross(cell+6, cell+0, xprod);
  dist[1] = volume/norm(xprod);
  cross(cell, cell+3, xprod);
  dist[2] = volume/norm(xprod);

  // number of cells in each direction
  double ratio[DIM];
  double size[DIM];
  for (int i=0; i<DIM; i++) {
    ratio[i] = cutoff/dist[i];
    size[i] = static_cast<int> (std::ceil(ratio[i]));
  }

  // creating padding atoms
  for (int i=-size[0]; i<=size[0]; i++)
  for (int j=-size[1]; j<=size[1]; j++)
  for (int k=-size[2]; k<=size[2]; k++) {

    // skip contributing atoms
    if (i==0 && j==0 && k==0) continue;

    // apply BC
    if (PBC[0]==0 && i != 0) continue;
    if (PBC[1]==0 && j != 0) continue;
    if (PBC[2]==0 && k != 0) continue;

    for (int at=0; at<Natoms; at++) {
      double x = frac_coords[DIM*at+0];
      double y = frac_coords[DIM*at+1];
      double z = frac_coords[DIM*at+2];

      // select the necessary atoms to repeate for the most outside bins
      // the follwing few lines can be easily understood when assuming size=1
      if (i == -size[0] && x - min[0] < static_cast<double>(size[0]) - ratio[0])
        continue;
      if (i ==  size[0] && max[0] - x < static_cast<double>(size[0]) - ratio[0])
        continue;
      if (j == -size[1] && y - min[1] < static_cast<double>(size[1]) - ratio[1])
        continue;
      if (j ==  size[1] && max[1] - y < static_cast<double>(size[1]) - ratio[1])
        continue;
      if (k == -size[2] && z - min[2] < static_cast<double>(size[2]) - ratio[2])
        continue;
      if (k ==  size[2] && max[2] - z < static_cast<double>(size[2]) - ratio[2])
        continue;

      // fractional coords of padding atom at
      double atom_coords[3] = {i+x, j+y, k+z};

      // absolute coords of padding atoms
      pad_coords.push_back(dot(tcell, atom_coords));
      pad_coords.push_back(dot(tcell+3, atom_coords));
      pad_coords.push_back(dot(tcell+6, atom_coords));

      // padding species code and image
      pad_species.push_back(species[at]);
      pad_image.push_back(at);

    }
  }

}
