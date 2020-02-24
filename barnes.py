import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
from os import getcwd

# Functions
def array_2_pp(array):
    pp = (array.__array_interface__['data'][0]
          + np.arange(array.shape[0])
          * array.strides[0]).astype(np.uintp)
    return pp


def barnes(x_scatter,
           y_scatter,
           z_scatter,
           x_coord,
           y_coord,
           xr=1.,
           yr=1.,
           iters=1,
           gamma=0.5):

    # Load library
    # _dll = ctypes.CDLL(getcwd() + '/build/lib.linux-x86_64-3.7/barneslib.cpython-37m-x86_64-linux-gnu.so')
    _dll = ctypes.CDLL(getcwd() + '/shared/barneslib.so')

    # Info from input
    x_size, y_size = x_coord.size, y_coord.size
    x_grid, _ = np.meshgrid(x_coord, y_coord)

    # Output grid
    z_grid = np.zeros_like(x_grid)
    zpp = array_2_pp(z_grid)

    # Setup argument types
    _dll.create_grid_barnes.argtypes = [ctypes.c_double,
                                        ctypes.c_double,
                                        ctypes.c_double,
                                        ctypes.c_int,
                                        ctypes.c_int,
                                        ctypes.c_int,
                                        ctypes.c_int,
                                        ndpointer(),
                                        ndpointer(),
                                        ndpointer(dtype=np.uintp, ndim=1, flags='C'),
                                        ndpointer(),
                                        ndpointer(),
                                        ndpointer()]
    _dll.create_grid_barnes.restype = None

    # Call c function
    _dll.create_grid_barnes(ctypes.c_double(xr),
                            ctypes.c_double(yr),
                            ctypes.c_double(gamma),
                            ctypes.c_int(iters),
                            ctypes.c_int(x_size),
                            ctypes.c_int(y_size),
                            ctypes.c_int(z_scatter.size),
                            x_coord,
                            y_coord,
                            zpp,
                            x_scatter,
                            y_scatter,
                            z_scatter)

    # Output
    return z_grid
