from skbuild import setup  # This line replaces 'from setuptools import setup'
import argparse

import io
import os
this_directory = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# default no pybind
with_pyroot='ON'
with_pybind='OFF'

if 'SUPERA_WITH_PYROOT' in os.environ:
    if eval(os.environ['SUPERA_WITH_PYROOT']):
        with_pyroot='ON'
    else:
        with_pyroot='OFF'

if 'SUPERA_WITH_PYBIND' in os.environ:
    if eval(os.environ['SUPERA_WITH_PYBIND']):
        with_pybind='ON'
    else:
        with_pybind='OFF'

if with_pybind=='ON' and with_pyroot=='ON':
    print('ERROR: cannot enable both PyROOT and Pybind.')
    print('       Set the shell environment variables explicitly to contro.')
    print('       $SUPERA_WITH_PYROOT=ON or OFF')
    print('       $SUPERA_WITH_PYBIND=ON or OFF')
    raise OSError

setup(
    name="supera",
    version="1.7.0",
    cmake_source_dir='src/',
    include_package_data=True,
    cmake_args=[
        #'-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON',
        #'-DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=10.9',
        '-DWITH_PYROOT={}'.format(with_pyroot),
        '-DWITH_PYBIND={}'.format(with_pybind),
    ],
    author=['Corey Adams', 'Kazuhiro Terao', 'Taritree Wongjirad', 'Marco Del Tutto', 'Jeremy Wolcott'],
    author_email='kterao@slac.stanford.edu',
    description='C++ framework to process particle physics detector simulation output for lartpc_mlreco3d machine learning data reconstruction software',
    license='MIT',
    keywords='larcv larcv3 neutrinos deep learning lartpc_mlreco3d',
    project_urls={
        'Source Code': 'https://github.com/DeepLearnPhysics/SuperaAtomic'
    },
    url='https://github.com/DeepLearnPhysics/SuperaAtomic',
    scripts=[],
    packages=['supera'],
    package_dir={'': 'python'},
    install_requires=[
        'numpy',
        'scikit-build',
        #'larcv',
    ],
    #extra_requires=[
    #    'pybind11_mkdoc',  # used to extract Python library comments from the C++ Doxygen
    #],
    long_description=long_description,
    long_description_content_type='text/markdown',
)
