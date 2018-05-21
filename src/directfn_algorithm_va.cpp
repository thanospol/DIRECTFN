#include <iostream>
#include <string>
#include "directfn_algorithm_va.h"
#include "directfn_kernel_tri.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_quad_vect.h"

using  std::cout;
using  std::endl;
using  std::string;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template  <typename ParticularKernel, typename ParticularQuadrature>
DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::DirectfnAlgorithm_VA():
DirectfnInterface<ParticularKernel, ParticularQuadrature>(),
pf_theta_p_lim_1_crnt_(nullptr),
pf_theta_q_lim_2_crnt_(nullptr),
pf_get_Lp_1_crnt_(nullptr),
pf_get_Lq_2_crnt_(nullptr),
pf_crnt_lam_maxlim_(nullptr),
Isub_(nullptr),
sin_Theta_p_1_(0.0),
cos_Theta_p_1_(0.0),
sin_Theta_q_2_(0.0),
cos_Theta_q_2_(0.0),
psi_a_(0.0),
psi_b_(0.0),
sin_Psi_3_(0.0),
cos_Psi_3_(0.0),
Lp_1_(0.0),
Lq_2_(0.0) {
}

template  <typename ParticularKernel, typename ParticularQuadrature>
DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::~DirectfnAlgorithm_VA() {
}

//virtual
template  <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::do_I_surface_surface_() {

    // Go throw subranges
    for (size_t m = 0; m < sub_ranges_numb_m_(); ++m) {

        // Update limits f-pointers for a given m'th subrange
        update_theta_p_lims_N1_fptr_(m);
        update_theta_q_lims_N2_fptr_(m);
        update_Lp1_fptr_(m);
        update_Lq2_fptr_(m);

        double theta_p_A = 0.0, theta_p_B = 0.0;
        pf_theta_p_lim_1_crnt_(theta_p_A, theta_p_B);
        const double half_theta_p_a_plus_b = 0.5 * (theta_p_A + theta_p_B);
        const double half_theta_p_b_mnus_a = 0.5 * (theta_p_B - theta_p_A);

        const double * const cz1 = this->up_z1_.get();
        const double * const cw1 = this->up_w1_.get();

        this->up_kerSummator_->nullify_Ipsi_1();
        for (size_t n_theta_p_1 = 0; n_theta_p_1 < this->N1(); ++n_theta_p_1) {

            const double Theta_p1 = half_theta_p_b_mnus_a * cz1[n_theta_p_1] + half_theta_p_a_plus_b;
            update_theta_p1_trigonometry_(Theta_p1);
            Lp_1_ = pf_get_Lp_1_crnt_(cos_Theta_p_1_, sin_Theta_p_1_);

            double theta_q_A = 0.0, theta_q_B = 0.0;
            pf_theta_q_lim_2_crnt_(theta_q_A, theta_q_B);
            const double half_theta_q_a_plus_b = 0.5 * (theta_q_A + theta_q_B);
            const double half_theta_q_b_mnus_a = 0.5 * (theta_q_B - theta_q_A);

            const double * const cz2 = this->up_z2_.get();
            const double * const cw2 = this->up_w2_.get();

            this->up_kerSummator_->nullify_Ieta_2();
            for (size_t n_theta_q_2 = 0; n_theta_q_2 < this->N2(); ++n_theta_q_2) {

                const double Theta_q2 = half_theta_q_b_mnus_a * cz2[n_theta_q_2] + half_theta_q_a_plus_b;
                update_theta_q2_trigonometry_(Theta_q2);
                Lq_2_ = pf_get_Lq_2_crnt_(cos_Theta_q_2_, sin_Theta_q_2_);
                // Lp1 and Lq2 are ready now:
                const double atan_Lq_to_Lp = atan(Lq_2_ / Lp_1_);
                // 1
                define_psi_a_psi_b_lamlim3_case_I_1_(atan_Lq_to_Lp);
                calc_I_psi_N3_();
                this->up_kerSummator_->accumulate_Ieta_2(cw2[n_theta_q_2]);
                // 2
                define_psi_a_psi_b_lamlim3_case_II_2_(atan_Lq_to_Lp);
                calc_I_psi_N3_();
                this->up_kerSummator_->accumulate_Ieta_2(cw2[n_theta_q_2]);
            } // N2 for
            this->up_kerSummator_->multiply_Ieta_2(half_theta_q_b_mnus_a);
            this->up_kerSummator_->accumulate_Ipsi_1(cw1[n_theta_p_1]);
        } // N1 for
        this->up_kerSummator_->multiply_Ipsi_1(half_theta_p_b_mnus_a);
        this->up_kerSummator_->assign_Ipsi_1_to( &Isub_[this->up_kernel_->size() * m] );
    } // 1 or 4  ==  m for triangles or quadrilateral
    gather_Iss_();
}

template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::gather_Iss_() noexcept {

    const size_t t_ker_sz = this->up_kernel_->size();
    // Collect sub ranges for each kernel
    this->Iss_tot_ = dcomplex(0.,0.);
    for (size_t i = 0; i < t_ker_sz; ++i) {
        dcomplex Iss_isum(0.0, 0.0);
        for (size_t m = 0; m < sub_ranges_numb_m_(); ++m) {
            Iss_isum += this->Isub_[i + t_ker_sz * m];
        }
        this->Iss_[i] = Iss_isum * this->up_kernel_->precomputed_jacobian();
        this->Iss_tot_ += this->Iss_[i];
    }
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::update_theta_p1_trigonometry_(const double Theta_p1) noexcept {

    sin_Theta_p_1_ = sin(Theta_p1);
    cos_Theta_p_1_ = cos(Theta_p1);
}

template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::update_theta_q2_trigonometry_(const double Theta_q2) noexcept {

    sin_Theta_q_2_ = sin(Theta_q2);
    cos_Theta_q_2_ = cos(Theta_q2);
}

template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::define_psi_a_psi_b_lamlim3_case_I_1_ (const double atan_Lq_to_Lp) noexcept {

    psi_a_ = 0.0;
    psi_b_ = atan_Lq_to_Lp;
    pf_crnt_lam_maxlim_ = &DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::lam_lim_b_1_;
}

template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::define_psi_a_psi_b_lamlim3_case_II_2_(const double atan_Lq_to_Lp) noexcept {

    psi_a_ = atan_Lq_to_Lp;
    psi_b_ = M_PI / 2.0;
    pf_crnt_lam_maxlim_ = &DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::lam_lim_b_2_;
}

template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::calc_I_psi_N3_() noexcept {

    const double psi_plus_ab = 0.5 * (psi_b_ + psi_a_);
    const double psi_mnus_ab = 0.5 * (psi_b_ - psi_a_);

    const double * const cz3 = this->up_z3_.get();
    const double * const cw3 = this->up_w3_.get();

    this->up_kerSummator_->nullify_Ilam_3();
    for (size_t n3_psi = 0; n3_psi < this->N3(); ++n3_psi) {

        const double t_Psi = psi_mnus_ab * cz3[n3_psi] + psi_plus_ab;
        update_psi_3_trigonometry_(t_Psi);template <typename ParticularKernel, typename ParticularQuadrature>

        calc_I_rho_4_();
        this->up_kerSummator_->accumulate_Ilam_3(cw3[n3_psi]);
    }
    this->up_kerSummator_->multiply_Ilam_3(psi_mnus_ab);
}

template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::update_psi_3_trigonometry_(const double t_Psi_3) noexcept {

    cos_Psi_3_ = cos(t_Psi_3);
    sin_Psi_3_ = sin(t_Psi_3);
}

template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::calc_I_rho_4_() noexcept {

    const double J_lam = 0.5 * (this->*pf_crnt_lam_maxlim_)();

    const double * const cz4 = this->up_z4_.get();
    const double * const cw4 = this->up_w4_.get();

    this->up_kerSummator_->nullify_Irho_4();
    for (size_t n4_lam = 0; n4_lam < this->N4(); ++n4_lam) {

        const double t_Lam = J_lam * (cz4[n4_lam] + 1.0);
        double uvxi_p[3], uvxi_q[3];
        calc_pq_simplex_(uvxi_p, uvxi_q, t_Lam); // out, out, in

        this->up_kernel_->update_rp(uvxi_p);
        this->up_kernel_->update_rq(uvxi_q);
        // Precaches the Greens function here.
        this->up_kernel_->precompute_rp_rq_data();
        this->up_kerSummator_->accumulate_Irho_4(cw4[n4_lam] * t_Lam * t_Lam * t_Lam);
    }
    this->up_kerSummator_->multiply_Irho_4(sin_Psi_3_ * cos_Psi_3_ * J_lam);
}

template <typename ParticularKernel, typename ParticularQuadrature>
double DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::lam_lim_b_1_() const noexcept {
    return Lp_1_ / cos_Psi_3_;
}

template <typename ParticularKernel, typename ParticularQuadrature>
double DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::lam_lim_b_2_() const noexcept {
    return Lq_2_ / sin_Psi_3_;
}

template <typename ParticularKernel, typename ParticularQuadrature>
void DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>::calc_pq_simplex_(double xi_p_out[3], double xi_q_out[3],
                                            const double t_Lam) noexcept {

    const double u_eta_p = t_Lam * cos_Psi_3_ * cos_Theta_p_1_ - 1.0;
    const double v_xi_p  = t_Lam * cos_Psi_3_ * sin_Theta_p_1_ - unity4_zero3_();

    va_make_simplex(xi_p_out, u_eta_p, v_xi_p);

    const double u_eta_q = t_Lam * sin_Psi_3_ * cos_Theta_q_2_ - 1.0;
    const double v_xi_q  = t_Lam * sin_Psi_3_ * sin_Theta_q_2_ - unity4_zero3_();

    va_make_simplex(xi_q_out, u_eta_q, v_xi_q);
}

// Instantiation of the Triangular Constant Kernels
template class DirectfnAlgorithm_VA<TriangularKernel_Constant_ST, GaussLegendreQuadrature>;
template class DirectfnAlgorithm_VA<TriangularKernel_Constant_EA, GaussLegendreQuadrature>;
template class DirectfnAlgorithm_VA<TriangularKernel_Constant_VA, GaussLegendreQuadrature>;

template class DirectfnAlgorithm_VA<TriangularKernel_RWG_WS, GaussLegendreQuadrature>;
template class DirectfnAlgorithm_VA<TriangularKernel_RWG_SS, GaussLegendreQuadrature>;
template class DirectfnAlgorithm_VA<TriangularKernel_nxRWG_SS, GaussLegendreQuadrature>;

// Instantiation of the Quadrilateral Planar Kernels
template class DirectfnAlgorithm_VA<QuadrilateralKernel_PlanarScalar, GaussLegendreQuadrature>;
template class DirectfnAlgorithm_VA<QuadrilateralKernel_PlanarVectorWS, GaussLegendreQuadrature>;
template class DirectfnAlgorithm_VA<QuadrilateralKernel_PlanarVectorSS, GaussLegendreQuadrature>;
// .. nxRWG?

template class DirectfnAlgorithm_VA<QuadrilateralKernel_CurvilinearScalar, GaussLegendreQuadrature>;
template class DirectfnAlgorithm_VA<QuadrilateralKernel_CurvilinearVectorWS, GaussLegendreQuadrature>;
template class DirectfnAlgorithm_VA<QuadrilateralKernel_CurvilinearVectorSS, GaussLegendreQuadrature>;

///////////////////////////////////////////////////////////////////////////////

// Instantiation of the Triangular Constant Kernels
template class DirectfnAlgorithm_VA<TriangularKernel_Constant_ST, GaussLegendreQuadrature>;
template class DirectfnAlgorithm_VA<TriangularKernel_Constant_EA, GaussLegendreQuadrature>;
template class DirectfnAlgorithm_VA<TriangularKernel_Constant_VA, GaussLegendreQuadrature>;

template class DirectfnAlgorithm_VA<TriangularKernel_RWG_WS, GaussLegendreQuadrature>;
template class DirectfnAlgorithm_VA<TriangularKernel_RWG_SS, GaussLegendreQuadrature>;
template class DirectfnAlgorithm_VA<TriangularKernel_nxRWG_SS, GaussLegendreQuadrature>;

// Instantiation of the Quadrilateral Planar Kernels
template class DirectfnAlgorithm_VA<QuadrilateralKernel_PlanarScalar, ClenshawCurtisQuadrature>;
template class DirectfnAlgorithm_VA<QuadrilateralKernel_PlanarVectorWS, ClenshawCurtisQuadrature>;
template class DirectfnAlgorithm_VA<QuadrilateralKernel_PlanarVectorSS, ClenshawCurtisQuadrature>;
// .. nxRWG?

template class DirectfnAlgorithm_VA<QuadrilateralKernel_CurvilinearScalar, ClenshawCurtisQuadrature>;
template class DirectfnAlgorithm_VA<QuadrilateralKernel_CurvilinearVectorWS, ClenshawCurtisQuadrature>;
template class DirectfnAlgorithm_VA<QuadrilateralKernel_CurvilinearVectorSS, ClenshawCurtisQuadrature>;

template <typename ParticularKernel, typename ParticularQuadrature>
Triangular_VA<ParticularKernel, ParticularQuadrature>::Triangular_VA():
DirectfnAlgorithm_VA<ParticularKernel, ParticularQuadrature>() {

    // The kernel pointer as well as the kernel-array-summator pointer
    // has been allocated in the DirectfnInterface<ParticularKernel, ParticularQuadrature>

    // Here the memory for arrays is allocated according to kernel size.
    allocate_I_vars_();
    initialize_limit_fptrs_();
}

template <typename ParticularKernel, typename ParticularQuadrature>
Triangular_VA<ParticularKernel, ParticularQuadrature>::~Triangular_VA() {
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
string Triangular_VA<ParticularKernel, ParticularQuadrature>::name() const noexcept {
    return string("Triangular_VA_Constant");
}

template <typename ParticularKernel, typename ParticularQuadrature>
void Triangular_VA<ParticularKernel, ParticularQuadrature>::allocate_I_vars_() {

    const size_t t_ker_sz = this->up_kernel_->size();
    // Allocates the memory according to the bt-type (the kernel length may vary
    // depending on type of Basis-Testing functions)
    // The kernel vector is garanteed to have finite length
    this->Iss_.reset(new dcomplex[t_ker_sz]);
    this->Isub_.reset (new dcomplex[sub_figures_1_() * t_ker_sz]);
}

template <typename ParticularKernel, typename ParticularQuadrature>
void Triangular_VA<ParticularKernel, ParticularQuadrature>::initialize_limit_fptrs_() noexcept {

    this->pf_theta_p_lim_1_crnt_ = triag_va_theta_p_limits;
    this->pf_theta_q_lim_2_crnt_ = triag_va_theta_q_limits;

    this->pf_get_Lp_1_crnt_ = triag_va_Lp1;
    this->pf_get_Lq_2_crnt_ = triag_va_Lq2;
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
size_t Triangular_VA<ParticularKernel, ParticularQuadrature>::sub_ranges_numb_m_()  const noexcept {
    return sub_figures_1_();
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
void Triangular_VA<ParticularKernel, ParticularQuadrature>::update_theta_p_lims_N1_fptr_(const size_t ) noexcept {
    // function pointer  pf_theta_p_lim_1_crnt_
    // have been already setup in the constructor. No need to reset.
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
void Triangular_VA<ParticularKernel, ParticularQuadrature>::update_theta_q_lims_N2_fptr_(const size_t ) noexcept {
    // function pointer pf_theta_q_lim_2_crnt_
    // have been already setup in the constructor. No need to reset.
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
void Triangular_VA<ParticularKernel, ParticularQuadrature>::update_Lp1_fptr_(const size_t ) noexcept {
    // function pointer pf_get_Lp_1_crnt_
    // have been already setup in the constructor. No need to reset.
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
void Triangular_VA<ParticularKernel, ParticularQuadrature>::update_Lq2_fptr_(const size_t ) noexcept {
    // function pointer pf_get_Lq_2_crnt_
    // have been already setup in the constructor. No need to reset.
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
double Triangular_VA<ParticularKernel, ParticularQuadrature>::unity4_zero3_() const noexcept {
    return 0.0;
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
void Triangular_VA<ParticularKernel, ParticularQuadrature>::va_make_simplex(double uvxi_p_out[3], const double u_eta_p,
                                                  const double v_xi_p) const noexcept {
    triag_st_make_simplex(uvxi_p_out, u_eta_p, v_xi_p);   // out[3], in, in
}


template class Triangular_VA<TriangularKernel_Constant_ST, GaussLegendreQuadrature>;
template class Triangular_VA<TriangularKernel_Constant_EA, GaussLegendreQuadrature>;
template class Triangular_VA<TriangularKernel_Constant_VA, GaussLegendreQuadrature>;

template class Triangular_VA<TriangularKernel_RWG_WS, GaussLegendreQuadrature>;
template class Triangular_VA<TriangularKernel_RWG_SS, GaussLegendreQuadrature>;
template class Triangular_VA<TriangularKernel_nxRWG_SS, GaussLegendreQuadrature>;

/////////////////////////////////////////////////////////////////////////////////

emplate class Triangular_VA<TriangularKernel_Constant_ST, ClenshawCurtisQuadrature>;
template class Triangular_VA<TriangularKernel_Constant_EA, ClenshawCurtisQuadrature>;
template class Triangular_VA<TriangularKernel_Constant_VA, ClenshawCurtisQuadrature>;

template class Triangular_VA<TriangularKernel_RWG_WS, ClenshawCurtisQuadrature>;
template class Triangular_VA<TriangularKernel_RWG_SS, ClenshawCurtisQuadrature>;
template class Triangular_VA<TriangularKernel_nxRWG_SS, ClenshawCurtisQuadrature>;

/////////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel, typename ParticularQuadrature>
Quadrilateral_VA<ParticularKernel, ParticularQuadrature>::Quadrilateral_VA():
DirectfnInterface<ParticularKernel, ParticularQuadrature>(),
pF_theta_p_lim1_arr_{nullptr, nullptr, nullptr, nullptr},
pF_theta_q_lim2_arr_{nullptr, nullptr, nullptr, nullptr},
pf_Lp_1_arr_{nullptr, nullptr, nullptr, nullptr},
pf_Lq_2_arr_{nullptr, nullptr, nullptr, nullptr} {

    // Initializtion of particular kernels is done in the
    // DirectfnInterface<ParticularKernel, ParticularQuadrature>
    allocate_I_vars_();
    initialize_limit_fptrs_();
}

template <typename ParticularKernel, typename ParticularQuadrature>
Quadrilateral_VA<ParticularKernel, ParticularQuadrature>::~Quadrilateral_VA() {
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
string Quadrilateral_VA<ParticularKernel, ParticularQuadrature>::name() const noexcept {
    return string("Quadrilateral_VA");
}

template <typename ParticularKernel, typename ParticularQuadrature>
void Quadrilateral_VA<ParticularKernel, ParticularQuadrature>::allocate_I_vars_() {

    const size_t t_ker_sz = this->up_kernel_->size();
    // Allocates the memory according to the bt-type (the kernel length may vary
    // depending on type of Basis-Testing functions)
    // The kernel vector is garanteed to have finite length
    this->Iss_.reset(new dcomplex[t_ker_sz]);
    this->Isub_.reset (new dcomplex[sub_figures_4_() * t_ker_sz]);
}

template <typename ParticularKernel, typename ParticularQuadrature>
void Quadrilateral_VA<ParticularKernel, ParticularQuadrature>::initialize_limit_fptrs_() noexcept {
    // N1
    pF_theta_p_lim1_arr_[0] = quad_va_theta_p_limits<0>;
    pF_theta_p_lim1_arr_[1] = quad_va_theta_p_limits<1>;
    pF_theta_p_lim1_arr_[2] = quad_va_theta_p_limits<2>;
    pF_theta_p_lim1_arr_[3] = quad_va_theta_p_limits<3>;
    // N2
    pF_theta_q_lim2_arr_[0] = quad_va_theta_q_limits<0>;
    pF_theta_q_lim2_arr_[1] = quad_va_theta_q_limits<1>;
    pF_theta_q_lim2_arr_[2] = quad_va_theta_q_limits<2>;
    pF_theta_q_lim2_arr_[3] = quad_va_theta_q_limits<3>;
    // N3
    pf_Lp_1_arr_[0] = quad_va_Lp1<0>;
    pf_Lp_1_arr_[1] = quad_va_Lp1<1>;
    pf_Lp_1_arr_[2] = quad_va_Lp1<2>;
    pf_Lp_1_arr_[3] = quad_va_Lp1<3>;
    // N4
    pf_Lq_2_arr_[0] = quad_va_Lq2<0>;
    pf_Lq_2_arr_[1] = quad_va_Lq2<1>;
    pf_Lq_2_arr_[2] = quad_va_Lq2<2>;
    pf_Lq_2_arr_[3] = quad_va_Lq2<3>;
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
size_t Quadrilateral_VA<ParticularKernel, ParticularQuadrature>::sub_ranges_numb_m_() const noexcept {
    return sub_figures_4_();
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
void Quadrilateral_VA<ParticularKernel, ParticularQuadrature>::update_theta_p_lims_N1_fptr_(const size_t  m) noexcept {
    this->pf_theta_p_lim_1_crnt_ = pF_theta_p_lim1_arr_[m];
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
void Quadrilateral_VA<ParticularKernel, ParticularQuadrature>::update_theta_q_lims_N2_fptr_(const size_t  m) noexcept {
    this->pf_theta_q_lim_2_crnt_ = pF_theta_q_lim2_arr_[m];
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
void Quadrilateral_VA<ParticularKernel, ParticularQuadrature>::update_Lp1_fptr_(const size_t m) noexcept {
    this->pf_get_Lp_1_crnt_ = pf_Lp_1_arr_[m];
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
void Quadrilateral_VA<ParticularKernel, ParticularQuadrature>::update_Lq2_fptr_(const size_t m) noexcept {
    this->pf_get_Lq_2_crnt_ = pf_Lq_2_arr_[m];
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
double Quadrilateral_VA<ParticularKernel, ParticularQuadrature>::unity4_zero3_() const noexcept {
    return 1.0;
}

//virtual
template <typename ParticularKernel, typename ParticularQuadrature>
void Quadrilateral_VA<ParticularKernel, ParticularQuadrature>::va_make_simplex(double uvxi_p_out[3], const double u_eta_p,
                                                               const double v_xi_p) const noexcept {
    uvxi_p_out[0] = u_eta_p;
    uvxi_p_out[1] = v_xi_p;
    uvxi_p_out[2] = 0.0;
}


// Instantiation of the Quadrilateral Constant Kernels
template class Quadrilateral_VA<QuadrilateralKernel_PlanarScalar, GaussLegendreQuadrature>;
template class Quadrilateral_VA<QuadrilateralKernel_PlanarVectorWS, GaussLegendreQuadrature>;
template class Quadrilateral_VA<QuadrilateralKernel_PlanarVectorSS, GaussLegendreQuadrature>;

template class Quadrilateral_VA<QuadrilateralKernel_CurvilinearScalar, GaussLegendreQuadrature>;
template class Quadrilateral_VA<QuadrilateralKernel_CurvilinearVectorWS, GaussLegendreQuadrature>;
template class Quadrilateral_VA<QuadrilateralKernel_CurvilinearVectorSS, GaussLegendreQuadrature>;

///////////////////////////////////////////////////////////////////////////////

// Instantiation of the Quadrilateral Constant Kernels
template class Quadrilateral_VA<QuadrilateralKernel_PlanarScalar, ClenshawCurtisQuadrature>;
template class Quadrilateral_VA<QuadrilateralKernel_PlanarVectorWS, ClenshawCurtisQuadrature>;
template class Quadrilateral_VA<QuadrilateralKernel_PlanarVectorSS, ClenshawCurtisQuadrature>;

template class Quadrilateral_VA<QuadrilateralKernel_CurvilinearScalar, ClenshawCurtisQuadrature>;
template class Quadrilateral_VA<QuadrilateralKernel_CurvilinearVectorWS, ClenshawCurtisQuadrature>;
template class Quadrilateral_VA<QuadrilateralKernel_CurvilinearVectorSS, ClenshawCurtisQuadrature>;

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

// End of the file

