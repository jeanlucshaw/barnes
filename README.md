# Desctiption

The [Gri software](http://gri.sourceforge.net/) written by
[Dan Kelley](https://www.dal.ca/faculty/science/oceanography/people/faculty/daniel-e-kelley.html) for plotting scientific data
contains a C++ implementation of scattered data interpolation using the method described in Koch et al. (1983). This project provides
a wrapper to the C++ code, such that it can be imported as a Python function.

# Installation

Make sure numpy and matplotlib are installed, then clone this repository somewhere on your python path.

```
git clone https://github.com/jeanlucshaw/barnes.git
```

The python setup script is meant to be called from Makefile. In the same directory, type:

```
make
```

To test the installation, type:

```
python barnes_test.py
```

which should display a plot of a control function to the right, and its estimation from randomly scattered data using the Barnes
interpolation on the left.

# Usage

Once installation is complete, the import call is:

```
from barnes import barnes
```

The mandatory inputs are the scattered data and its coordinates, as well as the interpolation grid vectors. A function call could
look something like the following pseudo code:

```
x_scatter = (1D numpy array)    # horizontal coordinate of scattered data
y_scatter = (1D numpy array)    # vertical coordinate of scattered data
z_scatter = (1D numpy array)    # scattered data
x_grid = (1D numpy array)       # horizontal coordinate of interpolation grid
y_grid = (1D numpy array)       # vertical coordinate of interpolation grid

z_grid = barnes(x_scatter,
                y_scatter,
                z_scatter,
                x_grid,
                y_grid)
```

Optional parameters are the horizontal and vertical search radius values (xr, yr) defaulting to 1., the gamma parameter defaulting
to 0.5 and the number of iterations defaulting to 1.

# References

Kelley, D. E., and P. S. Galbraith (2000), Gri: A language for scientific illustration, Linux J., 75, 92–101.

Koch, S. E., Desjardins, M., and Kocin, P. J. (1983) An interactive Barnes objective map analysis scheme for use with satellite and 
conventional data, J. Clim. Appl. Meteorol., 22, 1487–1503.
