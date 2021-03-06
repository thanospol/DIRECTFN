#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
//#include "directfn_quad.h'
//#include "directfn_tri.h"
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
using Directfn::Quadrilateral_EA;
using Directfn::Quadrilateral_VA;
using Directfn::Triangular_ST;
using Directfn::Triangular_EA;
using Directfn::Triangular_VA;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;
using Directfn::QuadrilateralKernel_PlanarScalar;
using Directfn::TriangularKernel_Constant_ST;
using Directfn::TriangularKernel_Constant_EA;
using Directfn::TriangularKernel_Constant_VA;

#ifndef M_PI
using Directfn::M_PI;
#endif

#ifndef DBL_EPSILON
using Directfn::DBL_EPSILON;
#endif
///////////////////////////////////////////////////////////////////////////////
 void  tri_vs_quad_st() noexcept {

    const double k0wn = 2 * M_PI;
    
    const double d = 0.1;
	const size_t  N_ref = 25;

    double r1[] = { 0.0 , 0.0 , 0.0 };
	double r2[] = { d , 0.0 , 0.0 };
	double r3[] = { d,  d , 0.0 };
	double r4[] = { 0.0 , d , 0.0 };

	double Error_quad, Error_tri;
	unsigned long T_quad, T_tri;
	const dcomplex * I_quad;
	const dcomplex * I11, *I22, *I12, *I21;
	dcomplex value_quad, value_tri;



    SingularContour3xn Q; // Quadrilateral
	SingularContour3xn T11;  // Self-Term 11
	SingularContour3xn T22; // Self-Term 22
	SingularContour3xn T12; // Self-Term 12
	SingularContour3xn T21; // Self-Term 21
    Q.set_points(r1, r2, r3, r4);
	T11.set_points(r3, r1, r2);
	T22.set_points(r1, r3, r4);
	T12.set_points(r3, r1, r2, r4);
	T21.set_points(r1, r3, r4, r2);

    unique_ptr<Quadrilateral_ST<QuadrilateralKernel_PlanarScalar>> up_quad_st(new Quadrilateral_ST<QuadrilateralKernel_PlanarScalar>());

	// Setting parameters
	up_quad_st->set_wavenumber(k0wn);
	up_quad_st->set(Q);
	up_quad_st->set_Gaussian_orders_4(N_ref, N_ref, N_ref, N_ref);

	// Calculating reference value
	cout << "Computing reference values for ST case ..." << endl;
	std::chrono::steady_clock::time_point  begin_t_quad = std::chrono::steady_clock::now();
	up_quad_st->calc_Iss(); std::chrono::steady_clock::time_point  end_t_quad = std::chrono::steady_clock::now();

	T_quad = std::chrono::duration_cast<std::chrono::milliseconds>(end_t_quad - begin_t_quad).count();
	cout << "Reference values computed in " << T_quad << " milliseconds" << endl;
	const dcomplex * I_ref = up_quad_st->Iss();
	dcomplex ref_val = I_ref[0];


	cout << setprecision(17) << ref_val << endl;
	


	unique_ptr<Triangular_ST<TriangularKernel_Constant_ST>> up_T11(new Triangular_ST<TriangularKernel_Constant_ST>());
	// Setting parameters
	up_T11->set_wavenumber(k0wn);
	up_T11->set(T11);

	unique_ptr<Triangular_ST<TriangularKernel_Constant_ST>> up_T22(new Triangular_ST<TriangularKernel_Constant_ST>());
	// Setting parameters
	up_T22->set_wavenumber(k0wn);
	up_T22->set(T22);

	unique_ptr<Triangular_EA<TriangularKernel_Constant_EA>> up_T12(new Triangular_EA<TriangularKernel_Constant_EA>());
	// Setting parameters
	up_T12->set_wavenumber(k0wn);
	up_T12->set(T12);

	unique_ptr<Triangular_EA<TriangularKernel_Constant_EA>> up_T21(new Triangular_EA<TriangularKernel_Constant_EA>());
	// Setting parameters
	up_T21->set_wavenumber(k0wn);
	up_T21->set(T21);
	std::ofstream myfile;
	myfile.open("Results_tri_vs_quad_st.txt");

	const int Counter = 25;

	cout << "Convergence test for ST case starting ..." << endl;
	for (int N = 1; N <= Counter; N++)
	{
		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(N, N, N, N);
		up_T11->set_Gaussian_orders_4(N, N, N, N);
		up_T22->set_Gaussian_orders_4(N, N, N, N);
		up_T12->set_Gaussian_orders_4(N, N, N, N);
		up_T21->set_Gaussian_orders_4(N, N, N, N);

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


		// DIRECTFN-tri
		std::chrono::steady_clock::time_point  begin_t_tri = std::chrono::steady_clock::now();
		up_T11->calc_Iss();
		up_T22->calc_Iss();
		up_T12->calc_Iss();
		up_T21->calc_Iss();
		I11 = up_T11->Iss();
		I22 = up_T22->Iss();
		I12 = up_T12->Iss();
		I21 = up_T21->Iss();
		value_tri = I11[0] + I22[0] + I12[0] + I21[0];
		std::chrono::steady_clock::time_point  end_t_tri = std::chrono::steady_clock::now();
		Error_tri = fabs(abs((value_tri - ref_val)) / abs(ref_val) + DBL_EPSILON);
		T_tri = std::chrono::duration_cast<std::chrono::microseconds>(end_t_tri - begin_t_tri).count();

		cout << "Runtime_tri: " << setprecision(4) << T_tri << " [mcsec]" << endl;
		cout << "Error = " << setprecision(20) << Error_tri << endl;

		// Writing results to file
		myfile << N << " ";
		myfile << setprecision(20) << Error_quad << " " << setprecision(4) << T_quad << " ";
		myfile << setprecision(20) << Error_tri << " " << setprecision(4) << T_tri << " " << endl;
	}

	myfile.close();
	cout << "Convergence test for ST case completed." << endl;

}

void  tri_vs_quad_ea() noexcept {

	const double k0wn = 2 * M_PI;

	const double d = 0.1;
	const size_t  N_ref = 25;

	double r1[] = { 0.0 , 0.0 , 0.0 };
	double r2[] = { d , 0.0 , 0.0 };
	double r3[] = { d,  d , 0.0 };
	double r4[] = { 0.0 , d , 0.0 };
	double r5[] = { 0.0, 2 * d, 0.0 };
	double r6[] = { d, 2 * d, 0.0 };

	double Error_quad, Error_tri;
	unsigned long T_quad, T_tri;
	const dcomplex * I_quad;
	const dcomplex * I13, *I14, *I23, *I24;
	dcomplex value_quad, value_tri;



	SingularContour3xn Q; // Quadrilaterals EA
	SingularContour3xn T13;  // Edge Adgacent 13
	SingularContour3xn T14; // Vertex Adgacent 14
	SingularContour3xn T23; // Vertex Adgacent 23
	SingularContour3xn T24; // Vertex Adgacent 24
	Q.set_points(r1, r2, r3, r4, r5, r6);
	T13.set_points(r3, r4, r1, r5);
	T14.set_points(r3, r4, r1, r6, r5);
	T23.set_points(r3, r1, r2, r5, r4);
	T24.set_points(r3, r1, r2, r6, r5);
	

	unique_ptr<Quadrilateral_EA<QuadrilateralKernel_PlanarScalar>> up_quad_ea(new Quadrilateral_EA<QuadrilateralKernel_PlanarScalar>());
	// Setting parameters
	up_quad_ea->set_wavenumber(k0wn);
	up_quad_ea->set(Q);
	up_quad_ea->set_Gaussian_orders_4(N_ref, N_ref, N_ref, N_ref);

	// Calculating reference value
	cout << "Computing reference values for EA case..." << endl;
	std::chrono::steady_clock::time_point  begin_t_quad = std::chrono::steady_clock::now();
	up_quad_ea->calc_Iss(); std::chrono::steady_clock::time_point  end_t_quad = std::chrono::steady_clock::now();

	T_quad = std::chrono::duration_cast<std::chrono::milliseconds>(end_t_quad - begin_t_quad).count();
	cout << "Reference value computed in " << T_quad << " milliseconds" << endl;
	const dcomplex * I_ref = up_quad_ea->Iss();
	dcomplex ref_val = I_ref[0];


	cout << setprecision(17) << ref_val << endl;



	unique_ptr<Triangular_EA<TriangularKernel_Constant_EA>> up_T13(new Triangular_EA<TriangularKernel_Constant_EA>());
	// Setting parameters
	up_T13->set_wavenumber(k0wn);
	up_T13->set(T13);

	unique_ptr<Triangular_VA<TriangularKernel_Constant_VA>> up_T14(new Triangular_VA<TriangularKernel_Constant_VA>());
	// Setting parameters
	up_T14->set_wavenumber(k0wn);
	up_T14->set(T14);

	unique_ptr<Triangular_VA<TriangularKernel_Constant_VA>> up_T23(new Triangular_VA<TriangularKernel_Constant_VA>());
	// Setting parameters
	up_T23->set_wavenumber(k0wn);
	up_T23->set(T23);

	unique_ptr<Triangular_VA<TriangularKernel_Constant_VA>> up_T24(new Triangular_VA<TriangularKernel_Constant_VA>());
	// Setting parameters
	up_T24->set_wavenumber(k0wn);
	up_T24->set(T24);
	std::ofstream myfile;
	myfile.open("Results_tri_vs_quad_ea.txt");

	const int Counter = 25;
	cout << "Convergence test for EA case starting ..." << endl;
	for (int N = 1; N <= Counter; N++)
	{
		// Setting Gaussian orders
		up_quad_ea->set_Gaussian_orders_4(N, N, N, N);
		up_T13->set_Gaussian_orders_4(N, N, N, N);
		up_T14->set_Gaussian_orders_4(N, N, N, N);
		up_T23->set_Gaussian_orders_4(N, N, N, N);
		up_T24->set_Gaussian_orders_4(N, N, N, N);

		// DIRECTFN-quad
		std::chrono::steady_clock::time_point  begin_t_quad = std::chrono::steady_clock::now();
		up_quad_ea->calc_Iss();
		std::chrono::steady_clock::time_point  end_t_quad = std::chrono::steady_clock::now();

		I_quad = up_quad_ea->Iss();
		value_quad = I_quad[0];
		Error_quad = fabs(abs((value_quad - ref_val)) / abs(ref_val) + DBL_EPSILON);
		T_quad = std::chrono::duration_cast<std::chrono::microseconds>(end_t_quad - begin_t_quad).count();

		cout << "N: " << N << endl;
		cout << "Runtime_quad: " << setprecision(4) << T_quad << " [mcsec]" << endl;
		cout << "I_EA_quad = " << setprecision(20) << value_quad << endl;


		// DIRECTFN-tri
		std::chrono::steady_clock::time_point  begin_t_tri = std::chrono::steady_clock::now();
		up_T13->calc_Iss();
		up_T14->calc_Iss();
		up_T23->calc_Iss();
		up_T24->calc_Iss();
		I13 = up_T13->Iss();
		I14 = up_T14->Iss();
		I23 = up_T23->Iss();
		I24 = up_T24->Iss();
		value_tri = I13[0] + I14[0] + I23[0] + I24[0];
		std::chrono::steady_clock::time_point  end_t_tri = std::chrono::steady_clock::now();
		Error_tri = fabs(abs((value_tri - ref_val)) / abs(ref_val) + DBL_EPSILON);
		T_tri = std::chrono::duration_cast<std::chrono::microseconds>(end_t_tri - begin_t_tri).count();

		cout << "Runtime_tri: " << setprecision(4) << T_tri << " [mcsec]" << endl;
		cout << "Error = " << setprecision(20) << Error_tri << endl;

		// Writing results to file
		myfile << N << " ";
		myfile << setprecision(20) << Error_quad << " " << setprecision(4) << T_quad << " ";
		myfile << setprecision(20) << Error_tri << " " << setprecision(4) << T_tri << " " << endl;
	}

	myfile.close();
	cout << "Convergence test for EA case completed." << endl;

}

void  tri_vs_quad_va() noexcept {

	const double k0wn = 2 * M_PI;

	const double d = 0.1;
	const size_t  N_ref = 25;

	double r1[] = { 0.0 , 0.0 , 0.0 };
	double r2[] = { d , 0.0 , 0.0 };
	double r3[] = { d,  d , 0.0 };
	double r4[] = { 0.0 , d , 0.0 };
	double r5[] = { 2 * d, d, 0.0 };
	double r6[] = { 2 * d, 2 * d, 0.0 };
	double r7[] = { d, 2 * d, 0.0 };

	double Error_quad, Error_tri;
	unsigned long T_quad, T_tri;
	const dcomplex * I_quad;
	const dcomplex * I13, *I14, *I23, *I24;
	dcomplex value_quad, value_tri;



	SingularContour3xn Q; // Quadrilaterals EA
	SingularContour3xn T13;  // Vertex Adgacent 13
	SingularContour3xn T14; // Vertex Adgacent 14
	SingularContour3xn T23; // Vertex Adgacent 23
	SingularContour3xn T24; // Vertex Adgacent 24
	Q.set_points(r1, r2, r3, r4, r5, r6, r7);
	T13.set_points(r3, r4, r1, r6, r7);
	T14.set_points(r3, r4, r1, r5, r6);
	T23.set_points(r3, r1, r2, r6, r7);
	T24.set_points(r3, r1, r2, r5, r6);


	unique_ptr<Quadrilateral_VA<QuadrilateralKernel_PlanarScalar>> up_quad_va(new Quadrilateral_VA<QuadrilateralKernel_PlanarScalar>());
	// Setting parameters
	up_quad_va->set_wavenumber(k0wn);
	up_quad_va->set(Q);
	up_quad_va->set_Gaussian_orders_4(N_ref, N_ref, N_ref, N_ref);

	cout << "Computing reference values for VA case..." << endl;
	std::chrono::steady_clock::time_point  begin_t_quad = std::chrono::steady_clock::now();
	up_quad_va->calc_Iss(); std::chrono::steady_clock::time_point  end_t_quad = std::chrono::steady_clock::now();

	T_quad = std::chrono::duration_cast<std::chrono::milliseconds>(end_t_quad - begin_t_quad).count();
	cout << "Reference value computed in " << T_quad << " milliseconds" << endl;
	const dcomplex * I_ref = up_quad_va->Iss();
	dcomplex ref_val = I_ref[0];



	unique_ptr<Triangular_VA<TriangularKernel_Constant_VA>> up_T13(new Triangular_VA<TriangularKernel_Constant_VA>());
	// Setting parameters
	up_T13->set_wavenumber(k0wn);
	up_T13->set(T13);

	unique_ptr<Triangular_VA<TriangularKernel_Constant_VA>> up_T14(new Triangular_VA<TriangularKernel_Constant_VA>());
	// Setting parameters
	up_T14->set_wavenumber(k0wn);
	up_T14->set(T14);

	unique_ptr<Triangular_VA<TriangularKernel_Constant_VA>> up_T23(new Triangular_VA<TriangularKernel_Constant_VA>());
	// Setting parameters
	up_T23->set_wavenumber(k0wn);
	up_T23->set(T23);

	unique_ptr<Triangular_VA<TriangularKernel_Constant_VA>> up_T24(new Triangular_VA<TriangularKernel_Constant_VA>());
	// Setting parameters
	up_T24->set_wavenumber(k0wn);
	up_T24->set(T24);
	std::ofstream myfile;
	myfile.open("Results_tri_vs_quad_va.txt");

	const int Counter = 25;

	cout << "Convergence test for VA case starting ..." << endl;
	for (int N = 1; N <= Counter; N++)
	{
		// Setting Gaussian orders
		up_quad_va->set_Gaussian_orders_4(N, N, N, N);
		up_T13->set_Gaussian_orders_4(N, N, N, N);
		up_T14->set_Gaussian_orders_4(N, N, N, N);
		up_T23->set_Gaussian_orders_4(N, N, N, N);
		up_T24->set_Gaussian_orders_4(N, N, N, N);

		// DIRECTFN-quad
		std::chrono::steady_clock::time_point  begin_t_quad = std::chrono::steady_clock::now();
		up_quad_va->calc_Iss();
		std::chrono::steady_clock::time_point  end_t_quad = std::chrono::steady_clock::now();

		I_quad = up_quad_va->Iss();
		value_quad = I_quad[0];
		Error_quad = fabs(abs((value_quad - ref_val)) / abs(ref_val) + DBL_EPSILON);
		T_quad = std::chrono::duration_cast<std::chrono::microseconds>(end_t_quad - begin_t_quad).count();

		cout << "N: " << N << endl;
		cout << "Runtime_quad: " << setprecision(4) << T_quad << " [mcsec]" << endl;
		cout << "I_VA_quad = " << setprecision(20) << value_quad << endl;


		// DIRECTFN-tri
		std::chrono::steady_clock::time_point  begin_t_tri = std::chrono::steady_clock::now();
		up_T13->calc_Iss();
		up_T14->calc_Iss();
		up_T23->calc_Iss();
		up_T24->calc_Iss();
		I13 = up_T13->Iss();
		I14 = up_T14->Iss();
		I23 = up_T23->Iss();
		I24 = up_T24->Iss();
		value_tri = I13[0] + I14[0] + I23[0] + I24[0];
		std::chrono::steady_clock::time_point  end_t_tri = std::chrono::steady_clock::now();
		Error_tri = fabs(abs((value_tri - ref_val)) / abs(ref_val) + DBL_EPSILON);
		T_tri = std::chrono::duration_cast<std::chrono::microseconds>(end_t_tri - begin_t_tri).count();

		cout << "Runtime_tri: " << setprecision(4) << T_tri << " [mcsec]" << endl;
		cout << "Error = " << setprecision(20) << Error_tri << endl;

		// Writing results to file
		myfile << N << " ";
		myfile << setprecision(20) << Error_quad << " " << setprecision(4) << T_quad << " ";
		myfile << setprecision(20) << Error_tri << " " << setprecision(4) << T_tri << " " << endl;
	}

	myfile.close();
	cout << "Convergence test for VA case completed." << endl;

}
///////////////////////////////////////////////////////////////////////////////

int main(int , char *  []) {

	tri_vs_quad_st();
	tri_vs_quad_ea();
	tri_vs_quad_va();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file





