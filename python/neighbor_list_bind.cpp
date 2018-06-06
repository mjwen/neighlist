#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdlib.h>
#include <iostream>
#include "neighbor_list.h"


#define MY_ERROR(message)                                             \
  {                                                                   \
    std::cout << "* Error (Neighbor List): \"" << message << "\" : "  \
              << __LINE__ << ":" << __FILE__ << std::endl;            \
    exit(1);                                                          \
  }

#define MY_WARNING(message)                                              \
  {                                                                      \
    std::cout << "* Error (Neighbor List) : \"" << message << "\" : "    \
              << __LINE__ << ":" << __FILE__ << std::endl;               \
  }


namespace py = pybind11;


PYBIND11_MODULE(neighlist, module) {
  module.doc() = "Python binding to neighbor list.";

  // py::nodelete, not to destroy object by python garbage collection
  py::class_<NeighList, std::unique_ptr<NeighList, py::nodelete> > (module, "NeighList")
    .def(py::init());

  module.def("initialize",
    []() {
      NeighList * nl;
      nbl_initialize(&nl);
      return nl;
    }
  );

  module.def("clean",
    [](NeighList * nl) {
      nbl_clean(&nl);
    }
  );

  module.def("build",
    [](NeighList *const nl,  double const cutoff,
      py::array_t<double> coords, py::array_t<int> need_neigh)
    {
      int Natoms_1 = coords.size()/3;
      int Natoms_2 = need_neigh.size();
      int error = Natoms_1 == Natoms_2 ? 0 : 1;
      if (error) MY_WARNING("\"coords\" size and \"need_neigh\" size does not match.");

      int Natoms = Natoms_1 <= Natoms_2 ? Natoms_1 : Natoms_2;
      auto c = coords.data();
      auto nn = need_neigh.data();
      error = error || nbl_build(nl, Natoms, cutoff, c, nn);

    return error;
    },
    py::arg("NeighList"),
    py::arg("cutoff"),
    py::arg("coords"),
    py::arg("need_neigh")
  );


  module.def("get_neigh",
    [](NeighList const * const nl, int const particleNumber) {
      int numberOfNeighbors=0;  // in case nbl_get_neigh error out
      int const * neighOfAtom;
      int error = nbl_get_neigh(nl, particleNumber, &numberOfNeighbors,
          &neighOfAtom);

      // pack as a numpy array
      auto neighborsOfParticle = py::array (py::buffer_info (
        const_cast<int*> (neighOfAtom),        // data pointer
        sizeof(int),                           // size of one element
        py::format_descriptor<int>::format(),  //Python struct-style format descriptor
        1,                                     // dimension
        {numberOfNeighbors},                            // size of each dimension
        {sizeof(int)}                          // stride of each dimension
      ));

      // short hand for the above
      //auto neighOfAtom = py::array(numberOfNeighbors, neighborsOfParticle);

      py::tuple re(3);
      re[0] = numberOfNeighbors;
      re[1] = neighborsOfParticle;
      re[2] = error;
      return re;
    },
    py::arg("NeighList"),
    py::arg("particle_number"),
    "Return(number_of_neighbors, neighbors_of_particle, error)"
  );

  module.def("get_neigh_kim",
  []() {
      // the allowed return pointer type by pybind11 is: void const *
      // so cast the function pointer to it, and we need to cast back when
      // using it
      typedef void const type;
      return (type*) &nbl_get_neigh_kim;
    }
  );


  module.def("create_paddings",
    [](double const cutoff, py::array_t<double> cell,
      py::array_t<int> PBC, py::array_t<double> coords, py::array_t<int> species)
    {

      int Natoms_1 = coords.size()/3;
      int Natoms_2 = species.size();
      int error = Natoms_1 == Natoms_2 ? 0 : 1;
      int Natoms = Natoms_1 <= Natoms_2 ? Natoms_1 : Natoms_2;
      if (error) MY_WARNING("\"coords\" size and \"species\" size does not match.");

      int Npad;
      std::vector<double> pad_coords;
      std::vector<int> pad_species;
      std::vector<int> pad_image;

      auto cell2 = cell.data();
      auto PBC2 = PBC.data();
      auto coords2 = coords.data();
      auto species2 = species.data();

      error = error || nbl_create_paddings(Natoms, cutoff, cell2, PBC2, coords2,
          species2, Npad, pad_coords, pad_species, pad_image);

    // pack as a 2D numpy array
    auto pad_coords_array = py::array (py::buffer_info (
      pad_coords.data(),
      sizeof(double),
      py::format_descriptor<double>::format(),
      2,
      {Npad, 3},
      {sizeof(double)*3, sizeof(double)}
    ));

    // pack as a numpy array
    auto pad_species_array = py::array (py::buffer_info (
      pad_species.data(),
      sizeof(int),
      py::format_descriptor<int>::format(),
      1,
      {Npad},
      {sizeof(int)}
    ));

    // pack as a numpy array
    auto pad_image_array = py::array (py::buffer_info (
      pad_image.data(),
      sizeof(int),
      py::format_descriptor<int>::format(),
      1,
      {Npad},
      {sizeof(int)}
     ));

      py::tuple re(4);
      re[0] = pad_coords_array;
      re[1] = pad_species_array;
      re[2] = pad_image_array;
      re[3] = error;
      return re;
    },
    py::arg("cutoff"),
    py::arg("cell"),
    py::arg("PBC"),
    py::arg("coordinates"),
    py::arg("species_code"),
    "Return(coordinates_of_paddings, species_code_of_paddings, \
        master_particle_of_paddings, error)"
  );

}

