#ifndef _DIRECTFN_QUADRATURES_H_
#define _DIRECTFN_QUADRATURES_H_

//#include "directfn_interface.h"
#include "directfn_common.h"

using std::string;
using std::unique_ptr;
using std::shared_ptr;
using std::size_t;
//using std::vector;


namespace Directfn {

	class AbstractQuadrature {
	public:

		AbstractQuadrature();
		virtual ~AbstractQuadrature();

		virtual bool set_zw_N(const size_t Nx, unique_ptr<double[]> & p_z1, unique_ptr<double[]> & p_w1) noexcept = 0;

	protected:
		AbstractQuadrature(const AbstractQuadrature &) = delete;
		AbstractQuadrature(AbstractQuadrature &&) = delete;
		AbstractQuadrature & operator = (const AbstractQuadrature &) = delete;
		AbstractQuadrature & operator = (AbstractQuadrature &&) = delete;
       
		
	};

	//////////////////////////////////////////////////////////////////////////////

	class GaussLegendreQuadrature: public AbstractQuadrature {
	public:

		GaussLegendreQuadrature();
		virtual ~GaussLegendreQuadrature();
		/*! Auxiliarily routine for quadrature setup */
		virtual bool set_zw_N(const size_t Nx, unique_ptr<double[]> & p_z1, unique_ptr<double[]> & p_w1) noexcept;

	protected:
		
	private:

	};

	//class ClenshawCurtisQuadrature : public AbstractQuadrature {
	//public:

	//	ClenshawCurtisQuadrature();
	//	virtual ~ClenshawCurtisQuadrature();
	//	/*! Auxiliarily routine for quadrature setup */
	//	virtual bool set_zw_N(const size_t Nx, unique_ptr<double[]> & p_z1, unique_ptr<double[]> & p_w1) noexcept;
	//protected:
	//private:
	//


	//	//void rescale(const size_t Nx, double points[], double weights[]) noexcept;
	//	bool ccn_compute_points_weights(const int Nx, double points[], double weights[]) noexcept;
	//	//bool nc_compute_new(const int Nx, const double x_min, const double x_max, double points[], double weights[]) noexcept;
	//	//virtual void rule_write ( int order, string filename, double x[], double w[], double r[] ) noexcep;
	//	//virtual void timestamp() noexcept = 0;
	//	//void r8mat_write(string output_filename, int m, int n, double table[]);

	//	//int i4_min(int i1, int i2);
	//};
	
} // End of the namespace  Directfn

#endif   //  _DIRECTFN_GREENFUNC_H_

  // End of the file
