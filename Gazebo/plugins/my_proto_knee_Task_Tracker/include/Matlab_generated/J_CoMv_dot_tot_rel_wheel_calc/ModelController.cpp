//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ModelController.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 06-Aug-2021 09:38:00
//

// Include Files
#include "ModelController.h"
#include <algorithm>
#include <cmath>

// Function Definitions
//
// Arguments    : double alpha
//                double b_gamma
//                double h1
//                double h2
//                double phi_r
//                double beta
//                double phi_r_dot
//                double beta_dot
//                double J_CoMv_dot_tot_rel_wheel[10]
// Return Type  : void
//
void ModelController::J_CoMv_dot_tot_rel_wheel_calc(
    double alpha, double b_gamma, double h1, double h2, double phi_r,
    double beta, double phi_r_dot, double beta_dot,
    double J_CoMv_dot_tot_rel_wheel[10])
{
  double J_CoMv_dot_tot_rel_wheel_tmp;
  double b_J_CoMv_dot_tot_rel_wheel_tmp;
  double t11;
  double t12;
  double t14;
  double t15;
  double t17;
  double t18;
  double t19;
  double t2;
  double t20;
  double t23;
  double t24;
  double t28;
  double t29;
  double t3;
  double t30;
  double t35;
  double t37;
  double t4;
  double t5;
  double t7;
  double t8;
  t2 = this->robot_dim[0] + this->robot_dim[4];
  t3 = this->robot_inertial_par[6] * 2.0;
  t4 = this->robot_inertial_par[8] * 2.0;
  t5 = this->robot_inertial_par[7] * 2.0;
  t7 = (beta + phi_r) + this->robot_dim[5];
  t8 = (b_gamma + phi_r) + this->robot_dim[5];
  t11 = std::cos(t7);
  t12 = std::cos(t8);
  t7 = std::sin(t7);
  t14 = std::sin(t8);
  t15 = phi_r + -this->robot_dim[6];
  t8 = (phi_r + -alpha) + this->robot_dim[5];
  t17 = this->robot_dim[3] * t11;
  t18 = this->robot_inertial_par[1] * t11;
  t19 = this->robot_dim[3] * t7;
  t20 = this->robot_inertial_par[1] * t7;
  t23 = std::cos(t8);
  t24 = std::sin(t8);
  t8 = 1.0 / ((((this->robot_inertial_par[9] + t3) + t4) + t5) +
              this->robot_inertial_par[5] * 2.0);
  t11 = this->robot_dim[4] * t23;
  t28 = this->robot_dim[4] * t24;
  t29 = h1 * t11;
  t30 = h1 * t28;
  t7 = beta_dot * t8;
  t35 = t7 * (this->robot_inertial_par[9] * t17 + t5 * t18);
  t37 = -(t7 * (this->robot_inertial_par[9] * t19 + t5 * t20));
  J_CoMv_dot_tot_rel_wheel[0] = 0.0;
  J_CoMv_dot_tot_rel_wheel[1] = 0.0;
  J_CoMv_dot_tot_rel_wheel[2] = 0.0;
  J_CoMv_dot_tot_rel_wheel[3] = 0.0;
  J_CoMv_dot_tot_rel_wheel[4] = 0.0;
  J_CoMv_dot_tot_rel_wheel[5] = 0.0;
  J_CoMv_dot_tot_rel_wheel_tmp = this->robot_inertial_par[2] * t3;
  b_J_CoMv_dot_tot_rel_wheel_tmp = phi_r_dot * t8;
  J_CoMv_dot_tot_rel_wheel[6] =
      t37 +
      b_J_CoMv_dot_tot_rel_wheel_tmp *
          (((this->robot_inertial_par[9] *
                 ((-t19 + t28) + this->robot_inertial_par[4] * std::cos(t15)) -
             t5 * (t20 - t28)) -
            t4 * (this->robot_inertial_par[3] * t12 - t2 * t24)) +
           J_CoMv_dot_tot_rel_wheel_tmp * t24);
  J_CoMv_dot_tot_rel_wheel[7] =
      t35 +
      b_J_CoMv_dot_tot_rel_wheel_tmp *
          (((this->robot_inertial_par[9] *
                 ((t17 - t11) + this->robot_inertial_par[4] * std::sin(t15)) +
             t5 * (t18 - t11)) -
            t4 * (this->robot_inertial_par[3] * t14 + t2 * t23)) -
           J_CoMv_dot_tot_rel_wheel_tmp * t23);
  J_CoMv_dot_tot_rel_wheel_tmp = h2 * this->robot_inertial_par[3];
  t8 = h1 * t2;
  t7 = h1 * this->robot_inertial_par[2] * t3;
  J_CoMv_dot_tot_rel_wheel[8] =
      t37 -
      b_J_CoMv_dot_tot_rel_wheel_tmp *
          (((this->robot_inertial_par[9] * (t19 + t30) + t5 * (t20 + t30)) +
            t4 * (J_CoMv_dot_tot_rel_wheel_tmp * t12 + t8 * t24)) +
           t7 * t24);
  J_CoMv_dot_tot_rel_wheel[9] =
      t35 +
      b_J_CoMv_dot_tot_rel_wheel_tmp *
          (((this->robot_inertial_par[9] * (t17 + t29) + t5 * (t18 + t29)) -
            this->robot_inertial_par[8] *
                (J_CoMv_dot_tot_rel_wheel_tmp * t14 - t8 * t23) * 2.0) +
           t7 * t23);
}

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
// File trailer for ModelController.cpp
//
// [EOF]
//
