#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>


// Interpolate in box f2 f3 f0 f1 using f = f0 + f1'x + f2'y +f3'xy, where
// ()' means () - f0, and x = xx - x_at_f0, etc.
// NOTE: ii,jj is point to lower-left of desired point.  If 
// at the top or right edge, just return the edge value.
// RETURN whether point is legit.
extern "C" bool value_i_j(unsigned int ii,
                          unsigned int jj,
                          double xx,
                          double yy,
                          double *value,
                          const unsigned int xsize,
                          const unsigned int ysize,
                          const double * xmatrix,
                          const double * ymatrix,
                          double ** zmatrix)
{
	double Dx, Dy;          // width/height of domain with point
	double f0, f1, f2, f3;
	double dx;			// x - x_to_left
	double dy;			// y - y_below
	// Fiddle with dx,dy,Dx,Dy, to avoid looking past array
	dx = (ii == xsize - 1 ? 0.0 : xx - xmatrix[ii]);
	dy = (jj == ysize - 1 ? 0.0 : yy - ymatrix[jj]);
	// if (_legit_xy(ii, jj) == false
	//     || (dx != 0.0 && _legit_xy(ii + 1, jj) == false)
	//     || (dy != 0.0 && _legit_xy(ii, jj + 1) == false)
	//     || (dx != 0.0 && dy != 0.0 && _legit_xy(ii + 1, jj + 1) == false)) {
	// 	*value = gr_currentmissingvalue();
	// 	return false;
	// }
	f0 = zmatrix[ii][jj];
	f1 = dx != 0 ? zmatrix[jj][ii + 1] - f0 : 0.0;
	f2 = dy != 0 ? zmatrix[jj + 1][ii] - f0 : 0.0;
	f3 = (dx != 0 && dy != 0) ? zmatrix[jj + 1][ii + 1] - f0 - f1 - f2 : 0;
	Dx = dx != 0 ? xmatrix[ii + 1] - xmatrix[ii] : 1;
	Dy = dy != 0 ? ymatrix[jj + 1] - ymatrix[jj] : 1;
	*value = f0 + f1 * dx / Dx + f2 * dy / Dy + f3 * dx / Dx * dy / Dy;
	return true;
}

// Do interpolation search, using bisection rule on possibly irregular
// array g[].
//
// If 'x' is in the range of the grid, defined by g[0] to g[ng-1],
// then set 'b' and 'f' such that
//     x = g[b] + f * (g[b+1] - g[b])
// and return true.
//
// If 'x' is not in the range, set b to the nearest endpoint, 
// set f to the distance to the nearest endpoint and return false.
extern "C" bool nearest(double x,
                        const double g[],
                        unsigned int ng,
                        int *b,
                        double *f)
{
	int l = 0;			// left index
	int r = ng - 1;		// right index
	int m;			// middle index
	if (g[0] < g[1]) {		// ascending sequence
		if (x <= g[l])	{ *b = 0; *f = g[l] - x; return false; }
		if (g[r] <= x)	{ *b = r; *f = x - g[r]; return false; }
		m = (l + r) / 2;
		while (r - l > 1) {
			if (x < g[m])
				r = m;
			else if (g[m] < x)
				l = m;
			else { 
				*b = m;
				*f = (x - g[*b]) / (g[*b+1] - g[*b]);
				return true; 
			}
			m = (r + l) / 2;
		}
		*b = l;
		*f = (x - g[*b]) / (g[*b+1] - g[*b]);
		return true;
	} else {			// descending sequence
		if (x >= g[l])	{ *b = 0; *f = g[l] - x; return false; }
		if (g[r] >= x)	{ *b = r; *f = x - g[r]; return false; }
		m = (l + r) / 2;
		while (r - l > 1) {
			if (x > g[m])
				r = m;
			else if (g[m] > x)
				l = m;
			else {
				*b = m;
				*f = (x - g[*b]) / (g[*b+1] - g[*b]);
				return true;
			}
			m = (r + l) / 2;
		}
		*b =  l;
		*f = (x - g[*b]) / (g[*b+1] - g[*b]);
		return true;
	}
}


// Barnes-interpolate to given (xx,yy), with previously value being zz.
// 'skip' used in cross-validation studies.
extern "C" double interpolate_barnes(double xx,
                                     double yy,
                                     double zz,
                                     int n_k,
                                     const double * x,
                                     const double * y,
                                     const double * z,
                                     const std::vector<double>& z_last,
                                     double xr,
                                     double yr)
{
	// if (gr_missing(zz))
	// 	return zz;
	double sum = 0.0, sum_w = 0.0;
	for (int k = 0; k < n_k; k++) {
		double w;
                double dx = (xx - x[k]) / xr;
                dx *= dx;
                double dy = (yy - y[k]) / yr;
                dy *= dy;
                double arg = dx + dy;
                w = exp(-arg);
                sum += w * (z[k] - z_last[k]);
                sum_w += w;
	}
	if (sum_w > 0.0)
          return (zz + sum / sum_w);
	else
          return 0.0;
          // return gr_currentmissingvalue();
}

//`convert columns to grid barnes    [.xr. .yr. .gamma. .iter.]'
extern "C" bool create_grid_barnes(double xr,
                                   double yr,
                                   double gamma,
                                   unsigned int iter,
                                   const unsigned int xsize,
                                   const unsigned int ysize,
                                   int ssize,
                                   const double *xmatrix,
                                   const double *ymatrix,
                                   double **zmatrix,
                                   const double * xgood,
                                   const double * ygood,
                                   const double * zgood)
{
	// Get grid storage if it does not exist already
	// if (!_grid_exists) {
	// 	Require(allocate_grid_storage(_num_xmatrix_data, _num_ymatrix_data), 
	// 		err("Insufficient space for matrix"));
	// }
	unsigned int numgood = ssize;
	// Require(numgood > 0,
	// 	err("Cannot `convert columns to grid' since no non-missing column data"));
	// _f_xy.set_value(0.0);
	// _legit_xy.set_value(true);

	std::vector<double> z_last((size_t)numgood, 0.0);
	std::vector<double> z_last2((size_t)numgood, 0.0); 

	// bool warned = false;
	// GriTimer t;
	double xr2 = xr, yr2 = yr;
	for (unsigned int iteration = 0; iteration < iter; iteration++) {
          // Interpolate on grid
          // size_t i, j;
          for (unsigned int i = 0; i < xsize; i++) {
            for (unsigned int j = 0; j < ysize; j++) {
                          zmatrix[j][i] = interpolate_barnes(xmatrix[i],
                                                             ymatrix[j],
                                                             zmatrix[j][i],
                                                             numgood,
                                                             xgood,
                                                             ygood,
                                                             zgood,
                                                             z_last,
                                                             xr2,
                                                             yr2);
			}
			// if (!warned) {
			// 	double frac = (iteration + 1.) * (i + 1.) * (j + 1.);
			// 	frac /= iter * _num_xmatrix_data * _num_ymatrix_data;
			// 	warned = warn_if_slow(&t, frac, "convert columns to grid");
			// }
		}
		// Interpolate at data
		unsigned int k;
		for (k = 0; k < numgood; k++) {
			int ix, iy;
			double fx, fy;
			bool in_x = nearest(xgood[k], xmatrix, xsize, &ix, &fx);
			bool in_y = nearest(ygood[k], ymatrix, ysize, &iy, &fy);
			if (in_x && in_y) {
				value_i_j(ix,
                                          iy,
                                          xgood[k],
                                          ygood[k],
                                          &z_last2[k],
                                          xsize,
                                          ysize,
                                          xmatrix,
                                          ymatrix,
                                          zmatrix);
			} else {
				z_last2[k]
					= interpolate_barnes(xgood[k],
							     ygood[k],
							     z_last[k],
							     numgood,
							     xgood,
							     ygood,
							     zgood,
							     z_last,
							     xr2,
							     yr2);
			}
		}
		// Calculate change in grid
		// double rms_change = 0.0;
		// int numgood_ok = 0;
		// for (k = 0; k < numgood; k++) {
		// 	if (!gr_missing(z_last[k]) && !gr_missing(z_last2[k])) {
		// 		rms_change += (z_last[k] - z_last2[k]) * (z_last[k] - z_last2[k]);
		// 		numgood_ok++;
		// 	}
		// }
		// if (numgood_ok)
		// 	rms_change = sqrt(rms_change / numgood_ok);
		// else
		// 	rms_change = gr_currentmissingvalue();
		// if (_chatty > 0) {
		// 	sprintf(_grTempString, "  Iteration %d: lengthscales = (%g,%g); RMS(z change) = %f\n", iteration + 1, xr2, yr2, rms_change);
		// 	ShowStr(_grTempString);
		// }
		// Update z_last
		for (k = 0; k < numgood; k++)
			z_last[k] = z_last2[k];
		// Catch case of gamma=0, which means not to iterate
		if (!gamma)
			break;
		// Alter search range
		xr2 *= sqrt(gamma);
		yr2 *= sqrt(gamma);
	}				// iteration
	return true;
}
