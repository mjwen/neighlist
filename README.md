# Neighbor List Building Library

This C++ library creates neighbor list for a group of atoms using the binning
method. That is, dividing the domain into cells and only loop over neighboring
cells of the target atom.

## Features

- C++ API and Python bindings
- The `KIM` branch works seamlessly with the `KIM API`

## Usage

### C++

The C++ shared library can be created by

```shell
$ cd neighlist
$ make
```

and then you can run an example by

```shell
$ cd examples
$ make
$ ./graphite
```

### Python

To use the Python binding, first do

```shell
$ cd neighlist/python
$ pip install -e . 
```

and then take a look at the example

```shell
$ cd tests
$ python test_graphite.py
```



Contact: Mingjian Wen (wenxx151@umn.edu)
