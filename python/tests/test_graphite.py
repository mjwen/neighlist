from __future__ import absolute_import, division, print_function
import numpy as np
import neighlist as nl
from error import check_error

def create_graphite_unit_cell(alat=2.46, d=3.35):

  cell = np.zeros(9).reshape(3,3)
  cell[0][0] = alat
  cell[1][0] = 0.5*alat
  cell[1][1] = np.sqrt(3)*0.5*alat
  cell[2][2] = 2*d

  coords = []
  coords.append([0,0,0])
  coords.append(cell[0]/3. + cell[1]/3.)
  coords.append(cell[0]/3. + cell[1]/3. + cell[2]/2.)
  coords.append(cell[0]*2/3. + cell[1]*2/3. + cell[2]/2.)
  coords = np.array(coords)

  species = np.array([1,1,2,2])

  return cell, coords, species


def write_XYZ(fname, lat_vec, species, coords):

  with open (fname, 'w') as fout:

    # first line (num of atoms)
    natoms = len(species)
    fout.write('{}\n'.format(natoms))

    # second line
    # lattice
    fout.write('Lattice="')
    for line in lat_vec:
      for item in line:
        fout.write('{:.15g} '.format(item))  # g removes trailing zero
    fout.write('" ')
    # PBC
    fout.write('PBC="1 1 1" ')
    # properties
    fout.write('Properties=species:S:1:pos:R:3\n')

    # body
    for i in xrange(len(species)):
      fout.write('{:3d} '.format(species[i]))
      fout.write('{:13.5e} {:13.5e} {:13.5e}'.format(
        coords[i][0], coords[i][1], coords[i][2]))
      fout.write('\n')



def test_main():

  # create contributing atoms
  alat = 2.46
  d = 3.35
  cell, contrib_coords, contrib_species = create_graphite_unit_cell(alat, d)

  # create padding atoms
  cutoff = d+0.01
  pbc = np.array([1,1,1])
  out = nl.create_paddings(cutoff, cell, pbc, contrib_coords, contrib_species)
  pad_coords, pad_species, pad_image, error = out
  check_error(error, 'nl.create_padding')

  assert pad_coords.shape == (96, 3)
  #print('pad_coords is of shape:', pad_coords.shape)

  coords = np.concatenate((contrib_coords, pad_coords))
  species = np.concatenate((contrib_species, pad_species))
  fname = 'atoms.xyz'
  write_XYZ(fname, cell, species, coords)

  # flag to indicate wheter create neighbor list for an atom
  n_pad = pad_coords.shape[0]
  n_contrib = contrib_coords.shape[0]
  need_neigh = np.concatenate((np.ones(n_contrib), np.zeros(n_pad)))

  # create neighbor list
  neigh = nl.initialize()
  error = nl.build(neigh, cutoff, coords, need_neigh)
  check_error(error, 'nl.build')

  # build again (it will automatically empty previous neigh list)
  error = nl.build(neigh, cutoff, coords, need_neigh)
  check_error(error, 'nl.build')

  # test get neigh function
  num_neigh, neighbors, error = nl.get_neigh(neigh, 1)
  check_error(error, 'nl.get_neigh')
  assert num_neigh == 14
  #print('Atom 1 has {} neighbors:'.format(num_neigh), end=' ')
  #for i in neighbors:
  #  print(i, end=' ')

  num_neigh, neighbors, error = nl.get_neigh(neigh, 4)
  check_error(error, 'nl.get_neigh')
  assert num_neigh == 0

  # expect error message from this
  #num_neigh, neighbors, error = nl.get_neigh(neigh, n_contrib+n_pad)
  #assert error == 1

  # delete neighbor list
  nl.clean(neigh)


if __name__ == '__main__':
  test_main()


