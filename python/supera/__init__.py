import os

def get_includes():

    # pip install
    data=os.path.join(os.path.dirname(__file__),"../../../../include/")
    if os.path.isdir(os.path.join(data,'supera')):
        return data
    # setup.py install
    data=os.path.join(os.path.dirname(__file__),"../include/")
    if os.path.isdir(os.path.join(data,'supera')):
        return data
    print('supera include path could not be located...')
    raise FileNotFoundError

def get_lib_dir():
    # pip install
    data=os.path.join(os.path.dirname(__file__),"../../../supera")
    if os.path.isfile(os.path.join(data,'libsupera.so')):
        return data
    if os.path.isfile(os.path.join(data,'libsupera.dylib')):
        return data
    # setup.py install
    data=os.path.join(os.path.dirname(__file__),"../lib/supera")
    if os.path.isfile(os.path.join(data,'libsupera.so')):
        return data
    if os.path.isfile(os.path.join(data,'libsupera.dylib')):
        return data
    print('supera library path could not be located...')
    raise FileNotFoundError 

__pybind,__pyroot=False,False

try:
    from . pysupera import *
    from . pysupera import test as test
    __pybind=True
except ModuleNotFoundError:
    __pybind=False

try:
    import ROOT
    lib=os.path.join(get_lib_dir(),'libsupera.so')
    ROOT.gSystem.Load(lib)
    from ROOT import supera
    # force load object definitions by instantiating one of any class
    supera.EDep
    __pyroot=True

except (ModuleNotFoundError,FileNotFoundError) as error:
    __pyroot=False

if not __pyroot and not __pybind:
    print('supera python binding could not be loaded')




