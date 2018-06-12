# Neighbor List Building Library 

This C++ library creates neighbor list for a group of atoms using the binning method. That is, dividing the domain into cells and only loop over neighboring cells of the target atom. 

## Features

- Binning method to create neighbor list
- C++ API and Python bindings
- Works seamlessly with the `KIM API` 

## Installation

### C++

The C++ shared library can be created by 

```
$ git clone https://github.com/mjwen/neighlist.git
$ cd neighlist
$ make 
```

and then an example can be run by

```
$ cd examples
$ make 
$ ./graphite
```

### Python 

To use the Python API (if you want to use the ASE [kimcalculator](https://github.com/mjwen/kimcalculator), this is needed), install by

```
$ git clone https://github.com/mjwen/neighlist.git
$ cd neighlist/python
$ pip install -e .
```

Then you can run an example:

```
$ cd tests
$ python test_graphite.py
```

## Contact

 Mingjian Wen (wenxx151@umn.edu)