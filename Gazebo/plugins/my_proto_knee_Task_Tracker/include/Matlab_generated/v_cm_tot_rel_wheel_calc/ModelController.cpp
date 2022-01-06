//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ModelController.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 11-Aug-2021 12:57:41
//

// Include Files
#include "ModelController.h"
#include <algorithm>
#include <cmath>

// Function Definitions
//
// Arguments    : void
// Return Type  : void
//
ModelController::ModelController()
{
  static const double dv[15]{0.06,
                             0.128,
                             0.22749999999999998,
                             0.118,
                             0.0363,
                             0.3,
                             0.29484,
                             0.165888,
                             0.152928,
                             10.0,
                             0.001,
                             0.0050866042499999988,
                             0.000905969664,
                             0.00070978982399999991,
                             0.1};
  static const double dv1[7]{
      0.1, 0.145, 0.236, 0.256, 0.355, 0.78539816339744828, 2.5132741228718345};
  std::copy(&dv[0], &dv[15], &this->robot_inertial_par[0]);
  for (int i{0}; i < 7; i++) {
    this->robot_dim[i] = dv1[i];
  }
}

//
// Arguments    : void
// Return Type  : void
//
ModelController::~ModelController()
{
  // (no terminate code required)
}

//
// V_CM_TOT_REL_WHEEL_CALC
//     V_CM_TOT_REL_WHEEL =
//     V_CM_TOT_REL_WHEEL_CALC(ALPHA,BETA,BETA_DOT,GAMMA,H1,H2,L_A,L_D,L_CM_C,L_CM_D,L_CM_AH,L_CM_BODY,L_H,M_AH,M_C,M_BODY,M_D,M_W,PHI_R,PHI_R_DOT,THETA_ST,THETA_ST_BODY)
//
// Arguments    : double alpha
//                double b_gamma
//                double h1
//                double h2
//                double phi_r
//                double beta
//                double phi_r_dot
//                double beta_dot
//                double v_cm_tot_rel_wheel[2]
// Return Type  : void
//
void ModelController::v_cm_tot_rel_wheel_calc(double alpha, double b_gamma,
                                              double h1, double h2,
                                              double phi_r, double beta,
                                              double phi_r_dot, double beta_dot,
                                              double v_cm_tot_rel_wheel[2])
{
  double b_v_cm_tot_rel_wheel_tmp;
  double c_v_cm_tot_rel_wheel_tmp;
  double d_v_cm_tot_rel_wheel_tmp;
  double e_v_cm_tot_rel_wheel_tmp;
  double f_v_cm_tot_rel_wheel_tmp;
  double g_v_cm_tot_rel_wheel_tmp;
  double h_v_cm_tot_rel_wheel_tmp;
  double i_v_cm_tot_rel_wheel_tmp;
  double j_v_cm_tot_rel_wheel_tmp;
  double k_v_cm_tot_rel_wheel_tmp;
  double l_v_cm_tot_rel_wheel_tmp;
  double m_v_cm_tot_rel_wheel_tmp;
  double n_v_cm_tot_rel_wheel_tmp;
  double o_v_cm_tot_rel_wheel_tmp;
  double t10;
  double t11;
  double t12;
  double t13;
  double t14;
  double t16;
  double t19;
  double t2;
  double t3;
  double t4;
  double t6;
  double t7;
  double v_cm_tot_rel_wheel_tmp;
  //     This function was generated by the Symbolic Math Toolbox version 8.7.
  //     11-Aug-2021 12:49:03
  t2 = this->robot_inertial_par[6] * 2.0;
  t3 = this->robot_inertial_par[8] * 2.0;
  t4 = this->robot_inertial_par[7] * 2.0;
  t6 = (beta + phi_r) + this->robot_dim[5];
  t7 = (b_gamma + phi_r) + this->robot_dim[5];
  t10 = std::cos(t6);
  t11 = std::cos(t7);
  t12 = std::sin(t6);
  t13 = std::sin(t7);
  t14 = phi_r + -this->robot_dim[6];
  t6 = (phi_r + -alpha) + this->robot_dim[5];
  t16 = std::cos(t6);
  t7 = std::sin(t6);
  t19 = 1.0 / ((((this->robot_inertial_par[9] + t2) + t3) + t4) +
               this->robot_inertial_par[5] * 2.0);
  v_cm_tot_rel_wheel_tmp = beta_dot * h1;
  b_v_cm_tot_rel_wheel_tmp = v_cm_tot_rel_wheel_tmp * this->robot_dim[4];
  c_v_cm_tot_rel_wheel_tmp =
      beta_dot * this->robot_dim[3] * this->robot_inertial_par[9];
  d_v_cm_tot_rel_wheel_tmp = beta_dot * this->robot_inertial_par[1] * t4;
  e_v_cm_tot_rel_wheel_tmp =
      this->robot_dim[0] * this->robot_inertial_par[8] * phi_r_dot;
  f_v_cm_tot_rel_wheel_tmp =
      this->robot_inertial_par[2] * this->robot_inertial_par[6] * phi_r_dot;
  g_v_cm_tot_rel_wheel_tmp =
      this->robot_dim[3] * this->robot_inertial_par[9] * phi_r_dot;
  h_v_cm_tot_rel_wheel_tmp =
      this->robot_dim[4] * this->robot_inertial_par[8] * phi_r_dot;
  i_v_cm_tot_rel_wheel_tmp =
      this->robot_dim[4] * this->robot_inertial_par[9] * phi_r_dot;
  j_v_cm_tot_rel_wheel_tmp =
      this->robot_dim[4] * this->robot_inertial_par[7] * phi_r_dot;
  k_v_cm_tot_rel_wheel_tmp = this->robot_inertial_par[1] * phi_r_dot * t4;
  l_v_cm_tot_rel_wheel_tmp =
      this->robot_inertial_par[4] * this->robot_inertial_par[9] * phi_r_dot;
  m_v_cm_tot_rel_wheel_tmp =
      b_v_cm_tot_rel_wheel_tmp * this->robot_inertial_par[9];
  n_v_cm_tot_rel_wheel_tmp = v_cm_tot_rel_wheel_tmp * this->robot_dim[0] * t3;
  o_v_cm_tot_rel_wheel_tmp = beta_dot * h2 * this->robot_inertial_par[3];
  v_cm_tot_rel_wheel_tmp =
      v_cm_tot_rel_wheel_tmp * this->robot_inertial_par[2] * t2;
  t6 = b_v_cm_tot_rel_wheel_tmp * t3;
  b_v_cm_tot_rel_wheel_tmp *= t4;
  v_cm_tot_rel_wheel[0] =
      t19 *
      ((((((((((((((((c_v_cm_tot_rel_wheel_tmp * t10 +
                      d_v_cm_tot_rel_wheel_tmp * t10) -
                     e_v_cm_tot_rel_wheel_tmp * t16 * 2.0) -
                    f_v_cm_tot_rel_wheel_tmp * t16 * 2.0) -
                   this->robot_inertial_par[3] * this->robot_inertial_par[8] *
                       phi_r_dot * t13 * 2.0) +
                  g_v_cm_tot_rel_wheel_tmp * t10) -
                 h_v_cm_tot_rel_wheel_tmp * t16 * 2.0) -
                i_v_cm_tot_rel_wheel_tmp * t16) -
               j_v_cm_tot_rel_wheel_tmp * t16 * 2.0) +
              k_v_cm_tot_rel_wheel_tmp * t10) +
             l_v_cm_tot_rel_wheel_tmp * std::sin(t14)) -
            o_v_cm_tot_rel_wheel_tmp * this->robot_inertial_par[8] * t13 *
                2.0) +
           m_v_cm_tot_rel_wheel_tmp * t16) +
          n_v_cm_tot_rel_wheel_tmp * t16) +
         v_cm_tot_rel_wheel_tmp * t16) +
        t6 * t16) +
       b_v_cm_tot_rel_wheel_tmp * t16);
  v_cm_tot_rel_wheel[1] =
      t19 * ((((((((((((((((c_v_cm_tot_rel_wheel_tmp * t12 +
                            d_v_cm_tot_rel_wheel_tmp * t12) -
                           e_v_cm_tot_rel_wheel_tmp * t7 * 2.0) -
                          f_v_cm_tot_rel_wheel_tmp * t7 * 2.0) +
                         g_v_cm_tot_rel_wheel_tmp * t12) -
                        h_v_cm_tot_rel_wheel_tmp * t7 * 2.0) -
                       i_v_cm_tot_rel_wheel_tmp * t7) -
                      j_v_cm_tot_rel_wheel_tmp * t7 * 2.0) +
                     this->robot_inertial_par[3] * phi_r_dot * t3 * t11) +
                    k_v_cm_tot_rel_wheel_tmp * t12) -
                   l_v_cm_tot_rel_wheel_tmp * std::cos(t14)) +
                  m_v_cm_tot_rel_wheel_tmp * t7) +
                 n_v_cm_tot_rel_wheel_tmp * t7) +
                o_v_cm_tot_rel_wheel_tmp * t3 * t11) +
               v_cm_tot_rel_wheel_tmp * t7) +
              t6 * t7) +
             b_v_cm_tot_rel_wheel_tmp * t7);
}

//
// File trailer for ModelController.cpp
//
// [EOF]
//
