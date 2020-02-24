"""
Python front end to the C++ libraries
"""
import ctypes
from os import getcwd
import numpy as np
from numpy.ctypeslib import ndpointer

def array_2_pp(array):
    """
    Create a C type access to the input 2D numpy array.
    """
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
    """
    Barnes interpolation of scattered data

    Wrapper to the C++ interpolation routine found in the Gri
    software.

    Mandatory inputs are:

    x_scatter:    [1D numpy array], x coordinate of scattered data
    y_scatter:    [1D numpy array], y coordinate of scattered data
    z_scatter:    [1D numpy array], scattered data
    x_grid:       [1D numpy array], x coordinate of interpolation grid
    y_grid:       [1D numpy array], y coordinate of interpolation grid

    Optional inputs are:

    xr:       [float], horizontal search radius, default 1.
    yr:       [float], vertical search radius, default 1.
    iters:    [int], number of iterations, default 1
    gamma:    [float], convergence parameter, default 0.5

    Authors:

    C++: Dan Kelley
    Python: Jean-Luc Shaw

    Reference:

    Koch, Desjardins, and Kocin (1983),
    An interactive Barnes objective map analysis scheme
    for use with satellite and conventional data,
    J. Clim. Appl. Meteorol, 22, 1487-1503
    """
    # Load library
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
