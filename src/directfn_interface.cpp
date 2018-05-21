#include <iostream>
#include <iomanip>
#include "directfn_interface.h"
#include "directfn_kernel_tri.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_quad_vect.h"

using std::cout;
using std::endl;
using std::setw;
using std::string;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel, typename ParticularQuadrature>
DirectfnInterface<ParticularKernel, ParticularQuadrature>::DirectfnInterface() :
up_kernel_(nullptr),
up_kerSummator_(nullptr),
up_quadrature_(nullptr),
up_w1_(nullptr),
up_z1_(nullptr),
up_w2_(nullptr),
up_z2_(nullptr),
up_w3_(nullptr),
up_z3_(nullptr),
up_w4_(nullptr),
up_z4_(nullptr),
Iss_(nullptr),
Iss_tot_(0.0, 0.0),
N1_(0),
N2_(0),
N3_(0),
N4_(0) {

    up_kernel_.reset(new ParticularKernel());
    up_kerSummator_.reset(new KernelArray<ParticularKernel>());
    up_kerSummator_->setup(up_kernel_.get());

	up_quadrature_.reset(new ParticularQuadrature());
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
DirectfnInterface<ParticularKernel, ParticularQuadrature>::~DirectfnInterface() {
}


template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnInterface<ParticularKernel, ParticularQuadrature>::set_wavenumber(const double k0_inp) noexcept {
    up_kernel_->set_wavenumber(k0_inp);
}

template <typename ParticularKernel, typename ParticularQuadrature>
bool DirectfnInterface<ParticularKernel, ParticularQuadrature>::set_Gaussian_orders_4(const size_t N1_in, const size_t N2_in,
                                                                const size_t N3_in, const size_t N4_in) noexcept {
    N1_ = N1_in;
    N2_ = N2_in;
    N3_ = N3_in;
    N4_ = N4_in;

    if (!DirectfnInterface<ParticularKernel>::set_zw_N_(N1_, up_z1_, up_w1_)) {return false;}
    if (!DirectfnInterface<ParticularKernel>::set_zw_N_(N2_, up_z2_, up_w2_)) {return false;}
    if (!DirectfnInterface<ParticularKernel>::set_zw_N_(N3_, up_z3_, up_w3_)) {return false;}
    if (!DirectfnInterface<ParticularKernel>::set_zw_N_(N4_, up_z4_, up_w4_)) {return false;}

    return true;
}

template <typename ParticularKernel, typename ParticularQuadrature>
bool DirectfnInterface<ParticularKernel, ParticularQuadrature>::set_Gaussian_orders_4(const size_t Nx[4]) noexcept {
    return DirectfnInterface::set_Gaussian_orders_4(Nx[0], Nx[1], Nx[2], Nx[3]);
}

template <typename ParticularKernel, typename ParticularQuadrature>
bool DirectfnInterface<ParticularKernel, ParticularQuadrature>::set_Gaussian_orders_4(const size_t Np) noexcept {
    return this->set_Gaussian_orders_4(Np, Np, Np, Np);
}

template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnInterface<ParticularKernel, ParticularQuadrature>::set(const SingularContour3xn &  contour_xpts) noexcept {
    up_kernel_->set(contour_xpts);
}

template <typename ParticularKernel, typename ParticularQuadrature>
size_t DirectfnInterface<ParticularKernel, ParticularQuadrature>::kernel_size() const noexcept {
    return up_kernel_->size();
}

template <typename ParticularKernel, typename ParticularQuadrature>
size_t  DirectfnInterface<ParticularKernel, ParticularQuadrature>::N1() const noexcept {
    return N1_;
}

template <typename ParticularKernel, typename ParticularQuadrature>
size_t  DirectfnInterface<ParticularKernel, ParticularQuadrature>::N2() const noexcept {
    return N2_;
}

template <typename ParticularKernel, typename ParticularQuadrature>
size_t  DirectfnInterface<ParticularKernel, ParticularQuadrature>::N3() const noexcept {
    return N3_;
}

template <typename ParticularKernel, typename ParticularQuadrature>
size_t  DirectfnInterface<ParticularKernel, ParticularQuadrature>::N4() const noexcept {
    return N4_;
}

template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnInterface<ParticularKernel, ParticularQuadrature>::calc_I_surface_surface() {
    return this->do_I_surface_surface_();
}

template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnInterface<ParticularKernel, ParticularQuadrature>::calc_Iss() {
    return this->do_I_surface_surface_();
}

template <typename ParticularKernel, typename ParticularQuadrature>
const dcomplex * DirectfnInterface<ParticularKernel, ParticularQuadrature>::Iss() const noexcept {
    return Iss_.get();
}

template <typename ParticularKernel, typename ParticularQuadrature>
const dcomplex DirectfnInterface<ParticularKernel, ParticularQuadrature>::Iss_tot() const noexcept {
    return Iss_tot_;
}

template <typename ParticularKernel, typename ParticularQuadrature>
const dcomplex DirectfnInterface<ParticularKernel, ParticularQuadrature>::Iss_arr(const size_t k) const noexcept {
    return Iss_[k];
}

template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnInterface<ParticularKernel, ParticularQuadrature>::copy_Iss_array_values_to(dcomplex * const out_array_to_be_setup) const noexcept {

    for (size_t k = 0; k < this->kernel_size(); ++k) {
        out_array_to_be_setup[k] = Iss_[k];
    }
}

template <typename ParticularKernel, typename ParticularQuadrature>
ParticularKernel * DirectfnInterface<ParticularKernel, ParticularQuadrature>::kernel_ptr() noexcept {
    return up_kernel_.get();
}

ParticularQuadrature * DirectfnInterface<ParticularKernel, ParticularQuadrature>::quadrature_ptr() noexcept {
	return up_quadrature_.get();
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnInterface<ParticularKernel, ParticularQuadrature>::debug_print() const noexcept {
    cout << "N1 = " << N1_ << "   N2 = " << N2_ << "   N3_ = " << N3_ << "   N4 = " << N4_ << endl;
}

template <typename ParticularKernel, typename ParticularQuadrature>::set_zw_N_(const size_t Nn, unique_ptr<double[]> & up_zn, unique_ptr<double[]> & up_wn) {
    //up_zn.reset(new double[Nn]);
    //up_wn.reset(new double[Nn]);

	return up_quadrature_->set_zw_N_(Nn, up_zn, up_wn);
    //gl_xw_1d(int(Nn), up_zn.get(), up_wn.get());
}


// Instantiation of Triangular elements
template class DirectfnInterface<TriangularKernel_Constant_ST, GaussLegendreQuadrature>;
template class DirectfnInterface<TriangularKernel_Constant_EA, GaussLegendreQuadrature>;
template class DirectfnInterface<TriangularKernel_Constant_VA, GaussLegendreQuadrature>;
template class DirectfnInterface<TriangularKernel_RWG_WS, GaussLegendreQuadrature>;
template class DirectfnInterface<TriangularKernel_RWG_SS, GaussLegendreQuadrature>;
template class DirectfnInterface<TriangularKernel_nxRWG_SS, GaussLegendreQuadrature>;


template class DirectfnInterface<TriangularKernel_Constant_ST, ClenshawCurtisQuadrature>;
template class DirectfnInterface<TriangularKernel_Constant_EA, ClenshawCurtisQuadrature>;
template class DirectfnInterface<TriangularKernel_Constant_VA, ClenshawCurtisQuadrature>;
template class DirectfnInterface<TriangularKernel_RWG_WS, ClenshawCurtisQuadrature>;
template class DirectfnInterface<TriangularKernel_RWG_SS, ClenshawCurtisQuadrature>;
template class DirectfnInterface<TriangularKernel_nxRWG_SS, ClenshawCurtisQuadrature>;

// Instantiation of Quadrilateral Constant kernels
template class DirectfnInterface<QuadrilateralKernel_PlanarScalar, GaussLegendreQuadrature>;
template class DirectfnInterface<QuadrilateralKernel_PlanarVectorWS, GaussLegendreQuadrature>;
template class DirectfnInterface<QuadrilateralKernel_PlanarVectorSS, GaussLegendreQuadrature>;

template class DirectfnInterface<QuadrilateralKernel_CurvilinearScalar, GaussLegendreQuadrature>;
template class DirectfnInterface<QuadrilateralKernel_CurvilinearVectorWS, GaussLegendreQuadrature>;
template class DirectfnInterface<QuadrilateralKernel_CurvilinearVectorSS, GaussLegendreQuadrature>;

// 
template class DirectfnInterface<QuadrilateralKernel_PlanarScalar, ClenshawCurtisQuadrature>;
template class DirectfnInterface<QuadrilateralKernel_PlanarVectorWS, ClenshawCurtisQuadrature>;
template class DirectfnInterface<QuadrilateralKernel_PlanarVectorSS, ClenshawCurtisQuadrature>;

template class DirectfnInterface<QuadrilateralKernel_CurvilinearScalar, ClenshawCurtisQuadrature>;
template class DirectfnInterface<QuadrilateralKernel_CurvilinearVectorWS, ClenshawCurtisQuadrature>;
template class DirectfnInterface<QuadrilateralKernel_CurvilinearVectorSS, ClenshawCurtisQuadrature>;

///////////////////////////////////////////////////////////////////////////////

}   // namespace Directfn

// End of the file


