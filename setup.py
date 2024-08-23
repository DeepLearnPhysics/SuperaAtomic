from skbuild import setup  # This line replaces 'from setuptools import setup'
import argparse

import io
import os
this_directory = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# default no pybind
with_python='ON'
with_pybind='OFF'

if 'SUPERA_WITHOUT_PYTHON' in os.environ and eval(os.environ['SUPERA_WITHOUT_PYTHON']):
    with_python=with_pybind='OFF'
if 'SUPERA_WITH_PYTHON' in os.environ:
    if eval(os.environ['SUPERA_WITH_PYTHON']):
        with_python='ON'
    else:
        with_python='OFF'

if 'SUPERA_WITH_PYBIND' in os.environ and eval(os.environ['SUPERA_WITH_PYBIND']):
    with_pybind='ON'
    if with_python == 'OFF':
        print('WARNING: SUPERA_WITH_PYBIND is positive and PYTHON is enabled althought it was set OFF')
        with_python = 'ON'
    


setup(
    name="supera",
    version="3.3.4",
    cmake_source_dir='src/',
    include_package_data=True,
    cmake_args=[
        #'-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON',
        '-DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=10.9',
        '-DWITH_PYTHON={}'.format(with_python),
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
