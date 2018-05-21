#include <iostream>

#include "directfn_interface"
#include "directfn_quadratures.h"
#include "directfn_common.h"

using  std::cout;
using  std::endl;
using  std::min;

namespace Directfn {

	AbstractQuadrature::AbstractQuadrature() :
		up_points_(nullptr),
		up_weights_(nullptr), {
		}

		AbstractQuadrature::~AbstractQuadrature() {
	}

	GaussLegendreQuadrature::GaussLegendreQuadrature() :
		AbstractQuadrature() {
	}

	GaussLegendreQuadrature::~GaussLegendreQuadrature() {
	}

	bool GaussLegendreQuadrature::set_zw_N_(const size_t Nx, unique_ptr<double[]> & p_z1, unique_ptr<double[]> & p_w1) const noexcept {
		up_points_.reset(new double[Nx]);
		up_weights_.reset(new double[Nx]);

		gl_xw_1d(int(Nx), up_points_.get(), up_weights_.get());
	}

	ClenshawCurtisQuadrature::ClenshawCurtisQuadrature() :
		AbstractQuadrature() {
	}

	ClenshawCurtisQuadrature::~ClenshawCurtisQuadrature() {
	}

	bool ClenshawCurtisQuadrature::set_zw_N_(const size_t Nx, unique_ptr<double[]> & p_z1, unique_ptr<double[]> & p_w1) const noexcept {

		up_points_.reset(new double[Nx]);
		up_weights_.reset(new double[Nx]);

		ccn_compute_points_new(Nx, up_points_.get());
		nc_compute_new(Nx, -1., 1., up_points_.get(), up_weights_.get());

	}

	void ClenshawCurtisQuadrature::ccn_compute_points_new(int n, double x[]) 
	
		//****************************************************************************80
		//
		//  Purpose:
		//
		//    CCN_COMPUTE_POINTS: compute Clenshaw Curtis Nested points.
		//
		//  Discussion:
		//
		//    We want to compute the following sequence:
		//
		//    1/2,
		//    0, 1
		//    1/4, 3/4
		//    1/8, 3/8, 5/8, 7/8,
		//    1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16, and so on.
		//
		//    But we would prefer that the numbers in each row be regrouped in pairs
		//    that are symmetric about 1/2, with the number above 1/2 coming first.
		//    Thus, the last row might become:
		//    (9/16, 7/16), (11/16, 5/16), ..., (15/16, 1/16).
		//
		//    Once we have our sequence, we apply the Chebyshev transformation
		//    which maps [0,1] to [-1,+1].
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license.
		//
		//  Modified:
		//
		//    06 March 2011
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, int N, the number of elements to compute.
		//
		//    Output, double CCN_COMPUTE_POINTS_NEW[N], the elements of the sequence.
	{
		int d;
		int i;
		int k;
		int m;
		int td;
		int tu;

		//
		//  Handle first three entries specially.
		//
		if (1 <= n)
		{
			x[0] = 0.5;
		}

		if (2 <= n)
		{
			x[1] = 1.0;
		}

		if (3 <= n)
		{
			x[2] = 0.0;
		}

		m = 3;
		d = 2;

		while (m < Nx)
		{
			tu = d + 1;
			td = d - 1;

			k = min(d, n - m);

			for (i = 1; i <= k; i++)
			{
				if ((i % 2) == 1)
				{
					x[m + i - 1] = tu / 2.0 / (double)(k);
					tu = tu + 2;
				}
				else
				{
					x[m + i - 1] = td / 2.0 / (double)(k);
					td = td - 2;
				}
			}
			m = m + k;
			d = d * 2;
		}
		//
		//  Apply the Chebyshev transformation.
		//
		for (i = 0; i < n; i++)
		{
			x[i] = cos(x[i] * M_PI);
		}
		x[0] = 0.0;

		if (2 <= n)
		{
			x[1] = -1.0;
		}

		if (3 <= n)
		{
			x[2] = +1.0;
		}

	}

	void ClenshawCurtisQuadrature::nc_compute_new(int n, double x_min, double x_max, double x[], double w[])

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    NC_COMPUTE_NEW computes a Newton-Cotes quadrature rule.
		//
		//  Discussion:
		//
		//    For the interval [X_MIN,X_MAX], the Newton-Cotes quadrature rule
		//    estimates
		//
		//      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
		//
		//    using N abscissas X and weights W:
		//
		//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
		//
		//    For the CLOSED rule, the abscissas include the end points.
		//    For the OPEN rule, the abscissas do not include the end points.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license.
		//
		//  Modified:
		//
		//    17 November 2009
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, int N, the order.
		//
		//    Input, double X_MIN, X_MAX, the endpoints of the interval.
		//
		//    Input, double X[N], the abscissas.
		//
		//    Output, double NC_COMPUTE_NEW[N], the weights.
		//
	{
		unique_ptr<double[]> d;
		int i;
		int j;
		int k;
		double yvala;
		double yvalb;

		d.reset(new double[n]);

		for (i = 0; i < n; i++)
		{
			//
			//  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
			//  and zero at the other nodes.
			//
			for (j = 0; j < n; j++)
			{
				d[j] = 0.0;
			}
			d[i] = 1.0;

			for (j = 2; j <= n; j++)
			{
				for (k = j; k <= n; k++)
				{
					d[n + j - k - 1] = (d[n + j - k - 1 - 1] - d[n + j - k - 1]) / (x[n + 1 - k - 1] - x[n + j - k - 1]);
				}
			}

			for (j = 1; j <= n - 1; j++)
			{
				for (k = 1; k <= n - j; k++)
				{
					d[n - k - 1] = d[n - k - 1] - x[n - k - j] * d[n - k];
				}
			}
			//
			//  Evaluate the antiderivative of the polynomial at the left and
			//  right endpoints.
			//
			yvala = d[n - 1] / (double)(n);
			for (j = n - 2; 0 <= j; j--)
			{
				yvala = yvala * x_min + d[j] / (double)(j + 1);
			}
			yvala = yvala * x_min;

			yvalb = d[n - 1] / (double)(n);
			for (j = n - 2; 0 <= j; j--)
			{
				yvalb = yvalb * x_max + d[j] / (double)(j + 1);
			}
			yvalb = yvalb * x_max;

			w[i] = yvalb - yvala;
		}

	}

} // End of DIRECTFN namespace

		

	/****************************************************************************
	void r8mat_write(string output_filename, int m, int n, double table[])

		//****************************************************************************
		//
		//  Purpose:
		//
		//    R8MAT_WRITE writes an R8MAT file with no header.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license.
		//
		//  Modified:
		//
		//    29 June 2009
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, string OUTPUT_FILENAME, the output filename.
		//
		//    Input, int M, the spatial dimension.
		//
		//    Input, int N, the number of points.
		//
		//    Input, double TABLE[M*N], the table data.
		//
	{
		int i;
		int j;
		ofstream output;
		//
		//  Open the file.
		//
		output.open(output_filename.c_str());

		if (!output)
		{
			cerr << "\n";
			cerr << "R8MAT_WRITE - Fatal error!\n";
			cerr << "  Could not open the output file.\n";
			return;
		}
		//
		//  Write the data.
		//
		for (j = 0; j < n; j++)
		{
			for (i = 0; i < m; i++)
			{
				output << "  " << setw(24) << setprecision(16) << table[i + j*m];
			}
			output << "\n";
		}
		//
		//  Close the file.
		//
		output.close();

		return;
	}
	//****************************************************************************80

	void rescale(double a, double b, int n, double x[], double w[])

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    RESCALE rescales a Legendre quadrature rule from [-1,+1] to [A,B].
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license.
		//
		//  Modified:
		//
		//    18 October 2009
		//
		//  Author:
		//
		//    John Burkardt.
		//
		//  Reference:
		//
		//    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin,
		//    A fast algorithm for the calculation of the roots of special functions,
		//    SIAM Journal on Scientific Computing,
		//    Volume 29, Number 4, pages 1420-1438, 2007.
		//
		//  Parameters:
		//
		//    Input, double A, B, the endpoints of the new interval.
		//
		//    Input, int N, the order.
		//
		//    Input/output, double X[N], on input, the abscissas for [-1,+1].
		//    On output, the abscissas for [A,B].
		//
		//    Input/output, double W[N], on input, the weights for [-1,+1].
		//    On output, the weights for [A,B].
		//
	{
		int i;

		for (i = 0; i < n; i++)
		{
			x[i] = ((a + b) + (b - a) * x[i]) / 2.0;
		}
		for (i = 0; i < n; i++)
		{
			w[i] = (b - a) * w[i] / 2.0;
		}
		return;
	}
	//****************************************************************************80

	void rule_write(int order, string filename, double x[], double w[],
		double r[])

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    RULE_WRITE writes a quadrature rule to three files.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license.
		//
		//  Modified:
		//
		//    18 February 2010
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, int ORDER, the order of the rule.
		//
		//    Input, double A, the left endpoint.
		//
		//    Input, double B, the right endpoint.
		//
		//    Input, string FILENAME, specifies the output filenames.
		//    "filename_w.txt", "filename_x.txt", "filename_r.txt"
		//    defining weights, abscissas, and region.
		//
	{
		string filename_r;
		string filename_w;
		string filename_x;
		int i;
		int kind;

		filename_w = filename + "_w.txt";
		filename_x = filename + "_x.txt";
		filename_r = filename + "_r.txt";

		cout << "\n";
		cout << "  Creating quadrature files.\n";
		cout << "\n";
		cout << "  Root file name is     \"" << filename << "\".\n";
		cout << "\n";
		cout << "  Weight file will be   \"" << filename_w << "\".\n";
		cout << "  Abscissa file will be \"" << filename_x << "\".\n";
		cout << "  Region file will be   \"" << filename_r << "\".\n";

		r8mat_write(filename_w, 1, order, w);
		r8mat_write(filename_x, 1, order, x);
		r8mat_write(filename_r, 1, 2, r);

		return;
	}
	//****************************************************************************80

	void timestamp()

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    TIMESTAMP prints the current YMDHMS date as a time stamp.
		//
		//  Example:
		//
		//    31 May 2001 09:45:54 AM
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license.
		//
		//  Modified:
		//
		//    08 July 2009
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    None
		//
	{
# define TIME_SIZE 40

		static char time_buffer[TIME_SIZE];
		const struct std::tm *tm_ptr;
		size_t len;
		std::time_t now;

		now = std::time(NULL);
		tm_ptr = std::localtime(&now);

		len = std::strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr);

		std::cout << time_buffer << "\n";

		return;
# undef TIME_SIZE
	}
	*/

