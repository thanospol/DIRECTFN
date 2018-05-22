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
#include "directfn_quadratures.h"
#include "directfn_interface.h"

using  std::cout;
using  std::endl;
using  std::setprecision;

//using Directfn::Quadrilateral_ST;
//using Directfn::Quadrilateral_EA;
//using Directfn::Quadrilateral_VA;
using Directfn::Triangular_ST;
//using Directfn::Triangular_EA;
//using Directfn::Triangular_VA;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;
//using Directfn::QuadrilateralKernel_PlanarScalar;
using Directfn::TriangularKernel_Constant_ST;
//using Directfn::TriangularKernel_Constant_EA;
using Directfn::AbstractQuadrature;


#ifndef M_PI
using Directfn::M_PI;
#endif

#ifndef DBL_EPSILON
using Directfn::DBL_EPSILON;
#endif
///////////////////////////////////////////////////////////////////////////////
 void  tri_st_strongly() noexcept {

    const double k0wn = 2 * M_PI;
    
    const double d = 0.1;
	const size_t  N_ref = 30;

	const double r1[] = { 0.0 , 0.4*d , 0.0 };
	const double r2[] = { d , 0.0 , 0.0 };
	const double r3[] = { 0.8*d,  d , 0.0 };
	//const double r4[] = { 0.2*d , 0.8*d , 0.0 };

	double Error;
	unsigned long Time;
	const dcomplex * I;
	dcomplex value;



    SingularContour3xn Tr; 
    Tr.set_points(r1, r2, r3);

	unique_ptr<Triangular_ST<TriangularKernel_Constant_ST>> up_Tr(new Triangular_ST<TriangularKernel_Constant_ST>());
	// Setting parameters
	up_Tr->set_wavenumber(k0wn);
	up_Tr->set(Tr);
	up_Tr->set_Gaussian_orders_4(N_ref, N_ref, N_ref, N_ref);

	// Calculating reference value
	cout << "Computing reference values ..." << endl;
	std::chrono::steady_clock::time_point  begin_t = std::chrono::steady_clock::now();
	up_Tr->calc_Iss(); std::chrono::steady_clock::time_point  end_t = std::chrono::steady_clock::now();

	Time = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t).count();
	cout << "Reference value computed in " << Time<< " milliseconds" << endl;
	

	const dcomplex * I_ref = up_Tr->Iss();
	dcomplex ref_val = I_ref[0];
	cout << "I_ref = " << setprecision(20) << ref_val << endl;

	


	std::ofstream myfile;
	myfile.open("Results_tri_ST_weakly.txt");

	const int Counter = 30;

	cout << "Convergence test starting ..." << endl;
	for (int N = 1; N <= Counter; N++)
	{
		// Setting Gaussian orders
		up_Tr->set_Gaussian_orders_4(N, N, N, N);

		// DIRECTFN-quad
		std::chrono::steady_clock::time_point  begin_t = std::chrono::steady_clock::now();
		up_Tr->calc_Iss();
		std::chrono::steady_clock::time_point  end_t = std::chrono::steady_clock::now();

		I= up_Tr->Iss();
		value = I[0];
		Error = fabs(abs((value - ref_val)) / abs(ref_val) + DBL_EPSILON);
		Time = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();

		cout << "N: " << N << endl;
		cout << "Runtime: " << setprecision(4) << Time << " [mcsec]" << endl;
		cout << "I_ST = " << setprecision(20) << value << endl;



		// Writing results to file
		myfile << N << " ";
		myfile << setprecision(20) << Error << endl;
	}

	myfile.close();
	cout << "Convergence test completed." << endl;

}

///////////////////////////////////////////////////////////////////////////////

int main(int , char *  []) {

    tri_st_strongly();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file





