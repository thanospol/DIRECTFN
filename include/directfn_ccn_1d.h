#ifndef _DIRECTFN_CCN_1D_H_
#define _DIRECTFN_CCN_1D_H_


using std::string;
using std::unique_ptr;
using std::shared_ptr;
using std::size_t;
using std::vector;


namespace Directfn {


	class QuadratureRule {
	public:

		QuadratureRule();
		~QuadratureRule();

		QuadratureRule(const QuadratureRule &) = delete;
		QuadratureRule(QuadratureRule &&) = delete;
		QuadratureRule & operator = (const QuadratureRule &) = delete;
		QuadratureRule & operator = (QuadratureRule &&) = delete;

		double *ccn_compute_points_new(int n) noexcept;

		int i4_min(int i1, int i2);
		double *nc_compute_new(int n, double x_min, double x_max, double x[]);

		void rescale(double a, double b, int n, double x[], double w[]);

		/*! Return precomputed value. */
		dcomplex  value() const noexcept;

		virtual void debug_print() const noexcept;

	protected:
		/*! The wave number is setup once and for all (but can be changed if needed) */
		double k0wn_;

	private:
		/*! Saves the result here */
		dcomplex precomputed_value_;

		/*! Calls the actual value */
		virtual dcomplex genuine_value_(const double R) const noexcept = 0;
	};

	double *ccn_compute_points_new(int n);
	int i4_min(int i1, int i2);
	double *nc_compute_new(int n, double x_min, double x_max, double x[]);
	//void r8mat_write(string output_filename, int m, int n, double table[]);
	void rescale(double a, double b, int n, double x[], double w[]);
	//void rule_write ( int order, string filename, double x[], double w[],
	//double r[] );
	void timestamp();
} // End of the namespace  Directfn

#endif   //  _DIRECTFN_GREENFUNC_H_

// End of the file
