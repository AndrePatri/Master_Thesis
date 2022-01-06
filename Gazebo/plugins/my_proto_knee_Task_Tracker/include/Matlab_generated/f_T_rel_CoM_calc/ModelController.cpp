//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ModelController.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 06-Aug-2021 19:05:40
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
// F_T_CALC
//     F_T =
//     F_T_CALC(X_REF1,X_REF2,ALPHA,BETA,BETA_DOT,GAMMA,H1,H2,L_A,L_D,L_CM_C,L_CM_D,L_CM_AH,L_CM_BODY,L_H,M_AH,M_C,M_BODY,M_D,M_W,PHI_R,PHI_R_DOT,THETA_ST,THETA_ST_BODY)
//
// Arguments    : double alpha
//                double b_gamma
//                double h1
//                double h2
//                const double Chi_command_ddot[2]
//                double phi_r
//                double beta
//                double phi_r_dot
//                double beta_dot
//                double f_T_rel_CoM[5]
// Return Type  : void
//
void ModelController::f_T_rel_CoM_calc(double alpha, double b_gamma, double h1,
                                       double h2,
                                       const double Chi_command_ddot[2],
                                       double phi_r, double beta,
                                       double phi_r_dot, double beta_dot,
                                       double f_T_rel_CoM[5])
{
  double t11;
  double t15;
  double t16;
  double t19;
  double t2;
  double t20;
  double t21;
  double t22;
  double t23;
  double t24;
  double t3;
  double t30;
  double t39;
  double t4;
  double t40;
  double t43;
  double t46;
  double t5;
  double t54;
  double t7;
  double t8;
  double t83;
  double t84;
  //     This function was generated by the Symbolic Math Toolbox version 8.7.
  //     29-Jul-2021 16:00:59
  t2 = this->robot_dim[0] + this->robot_dim[4];
  t3 = this->robot_inertial_par[6] * 2.0;
  t4 = this->robot_inertial_par[8] * 2.0;
  t5 = this->robot_inertial_par[7] * 2.0;
  t7 = (beta + phi_r) + this->robot_dim[5];
  t8 = (b_gamma + phi_r) + this->robot_dim[5];
  t11 = std::cos(t7);
  t7 = std::sin(t7);
  t15 = phi_r + -this->robot_dim[6];
  t16 = (phi_r + -alpha) + this->robot_dim[5];
  t19 = this->robot_dim[3] * t11;
  t20 = this->robot_inertial_par[1] * t11;
  t21 = this->robot_inertial_par[3] * std::cos(t8);
  t22 = this->robot_dim[3] * t7;
  t23 = this->robot_inertial_par[1] * t7;
  t24 = this->robot_inertial_par[3] * std::sin(t8);
  t11 = std::cos(t16);
  t30 = std::sin(t16);
  t54 = 1.0 / ((((this->robot_inertial_par[9] + t3) + t4) + t5) +
               this->robot_inertial_par[5] * 2.0);
  t39 = this->robot_dim[4] * t11;
  t40 = this->robot_dim[4] * t30;
  t43 = t2 * t11;
  t2 *= t30;
  t7 = this->robot_inertial_par[2] * t3;
  t46 = t7 * t11;
  t11 = t7 * t30;
  t16 = h1 * t39;
  t8 = h1 * t40;
  t7 = beta_dot * t54;
  t83 = t7 * (this->robot_inertial_par[9] * t19 + t5 * t20);
  t84 = t7 * (this->robot_inertial_par[9] * t22 + t5 * t23);
  t40 = ((t11 + -t5 * (t23 - t40)) + -t4 * (t21 - t2)) +
        this->robot_inertial_par[9] *
            ((this->robot_inertial_par[4] * std::cos(t15) + -t22) + t40);
  t3 = ((h1 * t11 + this->robot_inertial_par[9] * (t22 + t8)) +
        t5 * (t23 + t8)) +
       t4 * (h2 * t21 + h1 * t2);
  t2 = ((h1 * t46 + this->robot_inertial_par[9] * (t19 + t16)) +
        t5 * (t20 + t16)) +
       this->robot_inertial_par[8] * (h2 * t24 + -(h1 * t43)) * -2.0;
  t30 = ((t46 + t4 * (t24 + t43)) + -t5 * (t20 + -t39)) +
        -(this->robot_inertial_par[9] *
          ((t19 + this->robot_inertial_par[4] * std::sin(t15)) + -t39));
  t7 = phi_r_dot * t54;
  t16 = beta_dot * (t84 + t7 * t3) + phi_r_dot * (t84 + -(t7 * t40));
  t7 = beta_dot * (t83 + t7 * t2) + phi_r_dot * (t83 + -(t7 * t30));
  f_T_rel_CoM[0] = 0.0;
  f_T_rel_CoM[1] = 0.0;
  f_T_rel_CoM[2] = 0.0;
  t11 = Chi_command_ddot[0] * t54;
  t8 = Chi_command_ddot[1] * t54;
  f_T_rel_CoM[3] =
      ((t11 * t30 * 2.0 + t8 * t40 * 2.0) + t54 * t30 * t16 * 2.0) -
      t54 * t40 * t7 * 2.0;
  f_T_rel_CoM[4] = ((t11 * t2 * -2.0 - t8 * t3 * 2.0) + t54 * t3 * t7 * 2.0) -
                   t54 * t2 * t16 * 2.0;
}

//
// File trailer for ModelController.cpp
//
// [EOF]
//
