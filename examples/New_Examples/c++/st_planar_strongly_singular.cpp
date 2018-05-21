#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
#include "directfn_contour.h"
#include "directfn_defs.h"
#include "directfn_kernel_tri.h"
#include "directfn_algorithm_st.h"
#include "directfn_algorithm_ea.h"
#include "directfn_algorithm_va.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_tri.h"

using  std::cout;
using  std::endl;
using  std::setprecision;

using Directfn::Quadrilateral_ST;
//using Directfn::Quadrilateral_EA;
//using Directfn::Quadrilateral_VA;
//using Directfn::Triangular_ST;
//using Directfn::Triangular_EA;
//using Directfn::Triangular_VA;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;
using Directfn::QuadrilateralKernel_PlanarScalar;
//using Directfn::TriangularKernel_Constant_ST;
//using Directfn::TriangularKernel_Constant_EA;


#ifndef M_PI
using Directfn::M_PI;
#endif

#ifndef DBL_EPSILON
using Directfn::DBL_EPSILON;
#endif
///////////////////////////////////////////////////////////////////////////////
 void  quad_st_strongly() noexcept {

    const double k0wn = 2 * M_PI;
    
    const double d = 0.1;
	const size_t  N_ref = 32;

	const double r1[] = { 0.0 , 0.4*d , 0.0 };
	const double r2[] = { d , 0.0 , 0.0 };
	const double r3[] = { 0.8*d,  d , 0.0 };
	const double r4[] = { 0.2*d , 0.8*d , 0.0 };

	double Error_quad;
	unsigned long T_quad;
	const dcomplex * I_quad;
	dcomplex value_quad;



    SingularContour3xn Q; // Quadrilateral
    Q.set_points(r1, r2, r3, r4);

    unique_ptr<Quadrilateral_ST<QuadrilateralKernel_PlanarScalar>> up_quad_st(new Quadrilateral_ST<QuadrilateralKernel_PlanarScalar>());

	// Setting parameters
	up_quad_st->set_wavenumber(k0wn);
	up_quad_st->set(Q);
	up_quad_st->set_Gaussian_orders_4(N_ref, N_ref, N_ref, N_ref);
	
	

	// Calculating reference value
	cout << "Computing reference values ..." << endl;
	std::chrono::steady_clock::time_point  begin_t_quad = std::chrono::steady_clock::now();
	up_quad_st->calc_Iss(); std::chrono::steady_clock::time_point  end_t_quad = std::chrono::steady_clock::now();

	T_quad = std::chrono::duration_cast<std::chrono::milliseconds>(end_t_quad - begin_t_quad).count();
	cout << "Reference value computed in " << T_quad << " milliseconds" << endl;
	

	const dcomplex * I_ref = up_quad_st->Iss();
	dcomplex ref_val = I_ref[0];
	cout << "I_ref = " << setprecision(20) << ref_val << endl;

	


	std::ofstream myfile;
	myfile.open("Results_quad_ST_SS.txt");

	const int Counter = 30;

	cout << "Convergence test starting ..." << endl;
	for (int N = 1; N <= Counter; N++)
	{
		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4( N, N, N, N_ref);

		// DIRECTFN-quad
		std::chrono::steady_clock::time_point  begin_t_quad = std::chrono::steady_clock::now();
		up_quad_st->calc_Iss();
		std::chrono::steady_clock::time_point  end_t_quad = std::chrono::steady_clock::now();

		I_quad = up_quad_st->Iss();
		value_quad = I_quad[0];
		Error_quad = fabs(abs((value_quad - ref_val)) / abs(ref_val) + DBL_EPSILON);
		T_quad = std::chrono::duration_cast<std::chrono::microseconds>(end_t_quad - begin_t_quad).count();

		cout << "N: " << N << endl;
		cout << "Runtime_quad: " << setprecision(4) << T_quad << " [mcsec]" << endl;
		cout << "I_ST_quad = " << setprecision(20) << value_quad << endl;



		// Writing results to file
		myfile << N << " ";
		myfile << setprecision(20) << Error_quad << endl;
	}

	myfile.close();
	cout << "Convergence test completed." << endl;

}

///////////////////////////////////////////////////////////////////////////////

int main(int , char *  []) {

	quad_st_strongly();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file





