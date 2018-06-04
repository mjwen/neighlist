from setuptools import setup, Extension
from distutils.sysconfig import get_config_vars
import sys
import os
import subprocess


# remove `-Wstrict-prototypes' that is for C not C++
cfg_vars = get_config_vars()
for key, value in cfg_vars.items():
  if type(value) == str and '-Wstrict-prototypes' in value:
    cfg_vars[key] = value.replace('-Wstrict-prototypes', '')


def get_extra_compile_args():
  return ['-std=c++11']


class get_pybind11_includes(object):
  """Helper class to determine the pybind11 include path

  The purpose of this class is to postpone importing pybind11 until it is actually
  installed, so that the ``get_include()`` method can be invoked.

  Borrowd from: https://github.com/pybind/python_example/blob/master/setup.py
  """

  def __init__(self, user=False):
    self.user = user

  def __str__(self):
    import pybind11
    return pybind11.get_include(self.user)

def get_includes():
  neighlist_inc = ['../src/']
  pybind11_inc = [get_pybind11_includes(), get_pybind11_includes(user=True)]
  return pybind11_inc + neighlist_inc

def get_extension(module_name, sources):
  return Extension(
    module_name,
    sources = sources,
    include_dirs = get_includes(),
    extra_compile_args = get_extra_compile_args(),
    language = 'c++',
  )

neighlist = get_extension(
  'neighlist',
  ['../src/neighbor_list.cpp', 'neighbor_list_bind.cpp']
)

setup(
  name = 'neighlist',
  version = '1.0.0',
  ext_modules = [neighlist],
  install_requires = ['pybind11>=2.2', 'numpy'],

  # metadata
  author = 'Mingjian Wen',
  author_email = 'wenxx151[at]umn.edu',
  url = 'https://github.com/mjwen/neighbor_list',
  description = 'Python binding to the neighbor list building library',

  zip_safe = False,
)

