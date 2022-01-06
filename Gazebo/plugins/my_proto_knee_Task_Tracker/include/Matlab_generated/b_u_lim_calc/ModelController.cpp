//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ModelController.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 06-Aug-2021 10:29:52
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
  static const double dv2[7]{20.0,
                             0.0,
                             0.0,
                             0.0,
                             1.6406094968746698,
                             0.29575137404340573,
                             0.28873817557208159};
  std::copy(&dv[0], &dv[15], &this->robot_inertial_par[0]);
  for (int i{0}; i < 7; i++) {
    this->robot_dim[i] = dv1[i];
    this->robot_knee_el_par[i] = dv2[i];
  }
  this->g = 9.81;
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
// B_U_LIM_CALC
//     B_U_LIM =
//     B_U_LIM_CALC(C_LB_KNEE,C_LB_WHEEL,C_UB_KNEE,C_UB_WHEEL,ALPHA,ALPHA0,BETA0,BETA,BETA_DOT,G,GAMMA0,GAMMA,H1,H2,K1T,K2T,K3T,K4T,L_A,L_D,L_CM_C,L_CM_D,L_CM_AH,L_CM_BODY,L_H,M_AH,M_C,M_BODY,M_D,PHI_R,PHI_R_DOT,THETA_ST,THETA_ST_BODY)
//
// Arguments    : double alpha
//                double b_gamma
//                double h1
//                double h2
//                const double u_lim[4]
//                double phi_r
//                double beta
//                double phi_r_dot
//                double beta_dot
//                double b_u_lim[8]
// Return Type  : void
//
void ModelController::b_u_lim_calc(double alpha, double b_gamma, double h1,
                                   double h2, const double u_lim[4],
                                   double phi_r, double beta, double phi_r_dot,
                                   double beta_dot, double b_u_lim[8])
{
  double t10;
  double t101;
  double t101_tmp;
  double t102;
  double t102_tmp;
  double t105;
  double t106;
  double t108;
  double t11;
  double t111;
  double t12;
  double t123;
  double t125;
  double t13;
  double t133;
  double t14;
  double t15;
  double t151;
  double t152;
  double t153;
  double t154;
  double t155;
  double t156;
  double t157;
  double t159;
  double t16;
  double t160;
  double t161;
  double t162;
  double t163;
  double t164;
  double t165;
  double t166;
  double t167;
  double t17;
  double t172;
  double t18;
  double t182;
  double t183;
  double t19;
  double t194;
  double t20;
  double t202;
  double t21;
  double t22;
  double t23;
  double t24;
  double t25;
  double t26;
  double t27;
  double t28;
  double t29;
  double t37;
  double t39;
  double t4;
  double t5;
  double t6;
  double t64;
  double t65;
  double t66;
  double t67;
  double t68;
  double t69;
  double t7;
  double t70;
  double t73;
  double t74;
  double t75;
  double t76;
  double t77;
  double t78;
  double t79;
  double t8;
  double t82;
  double t83;
  double t87;
  double t88;
  double t89;
  double t9;
  double t90;
  double t91;
  double t92;
  double t93;
  double t94;
  double t95;
  double t96;
  double t97;
  double t97_tmp;
  double t98;
  double t98_tmp;
  //     This function was generated by the Symbolic Math Toolbox version 8.7.
  //     29-Jul-2021 16:01:07
  t4 = alpha * this->robot_knee_el_par[1];
  t5 = this->robot_knee_el_par[1] * this->robot_knee_el_par[5];
  t6 = this->robot_knee_el_par[1] * this->robot_knee_el_par[4];
  t7 = this->robot_knee_el_par[3] * this->robot_knee_el_par[4];
  t8 = beta * this->robot_knee_el_par[1];
  t9 = beta * this->robot_knee_el_par[3];
  t10 = beta_dot * beta_dot;
  t11 = h1 * h1;
  t12 = h2 * h2;
  t13 = phi_r_dot * phi_r_dot;
  t15 = alpha * h1 * this->robot_knee_el_par[0];
  t16 = this->robot_knee_el_par[5] * h1 * this->robot_knee_el_par[0];
  t18 = alpha * h2 * this->robot_knee_el_par[0];
  t20 = this->robot_knee_el_par[5] * h2 * this->robot_knee_el_par[0];
  t24 = this->robot_knee_el_par[6] * h1 * this->robot_knee_el_par[0];
  t39 = this->robot_knee_el_par[6] * h2;
  t25 = t39 * this->robot_knee_el_par[0];
  t26 = t39 * this->robot_knee_el_par[2];
  t27 = b_gamma * h1 * this->robot_knee_el_par[0];
  t39 = b_gamma * h2;
  t28 = t39 * this->robot_knee_el_par[0];
  t29 = t39 * this->robot_knee_el_par[2];
  t14 = std::cos(alpha + b_gamma);
  t17 = h1 * t4;
  t19 = h1 * t5;
  t21 = h1 * t6;
  t22 = h1 * t8;
  t23 = std::sin(alpha + beta);
  t37 = std::cos((beta + this->robot_dim[5]) + this->robot_dim[6]);
  t39 = std::sin((beta + phi_r) + this->robot_dim[5]);
  t64 = std::cos((-alpha + this->robot_dim[5]) + this->robot_dim[6]);
  t65 = std::sin((phi_r + -alpha) + this->robot_dim[5]);
  t66 = this->g * this->robot_inertial_par[3] * this->robot_inertial_par[8] *
        std::cos((b_gamma + phi_r) + this->robot_dim[5]);
  t67 = this->g * this->robot_dim[3] * this->robot_inertial_par[9] * t39;
  t68 =
      this->g * this->robot_inertial_par[1] * this->robot_inertial_par[7] * t39;
  t194 = beta_dot * this->robot_dim[3];
  t69 =
      t194 * this->robot_dim[4] * this->robot_inertial_par[9] * phi_r_dot * t23;
  t202 = this->robot_inertial_par[1] * this->robot_dim[4] *
         this->robot_inertial_par[7];
  t73 = t202 * t13 * t23;
  t39 = beta_dot * this->robot_inertial_par[1] * this->robot_dim[4] *
        this->robot_inertial_par[7] * phi_r_dot * t23;
  t75 = t39 * 2.0;
  t76 = t39 * 4.0;
  t78 = t194 * this->robot_inertial_par[4] * this->robot_inertial_par[9] *
        phi_r_dot * t37;
  t194 = this->robot_dim[3] * this->robot_dim[4] * this->robot_inertial_par[9];
  t82 = t194 * t10 * t23;
  t83 = t202 * t10 * t23;
  t93 = h1 * this->robot_dim[0] * this->robot_inertial_par[3] *
        this->robot_inertial_par[8] * t13 * t14;
  t94 = h2 * this->robot_dim[0] * this->robot_inertial_par[3] *
        this->robot_inertial_par[8] * t13 * t14;
  t95 = h1 * this->robot_inertial_par[3] * this->robot_dim[4] *
        this->robot_inertial_par[8] * t13 * t14;
  t96 = h2 * this->robot_inertial_par[3] * this->robot_dim[4] *
        this->robot_inertial_par[8] * t13 * t14;
  t125 = beta_dot * h1;
  t97_tmp = t125 * this->robot_dim[0] * this->robot_inertial_par[3] *
            this->robot_inertial_par[8] * phi_r_dot * t14;
  t97 = t97_tmp * 2.0;
  t39 = beta_dot * h2;
  t98_tmp = t39 * this->robot_dim[0] * this->robot_inertial_par[3] *
            this->robot_inertial_par[8] * phi_r_dot * t14;
  t98 = t98_tmp * 2.0;
  t101_tmp = t125 * this->robot_inertial_par[3] * this->robot_dim[4] *
             this->robot_inertial_par[8] * phi_r_dot * t14;
  t101 = t101_tmp * 2.0;
  t102_tmp = t39 * this->robot_inertial_par[3] * this->robot_dim[4] *
             this->robot_inertial_par[8] * phi_r_dot * t14;
  t102 = t102_tmp * 2.0;
  t202 = this->robot_dim[3] * this->robot_inertial_par[4] *
         this->robot_inertial_par[9];
  t106 = t202 * t10 * t37;
  t123 = t125 * this->robot_inertial_par[1] * this->robot_dim[4] *
         this->robot_inertial_par[7] * phi_r_dot * t23 * -2.0;
  t133 = t194 * t13 * t23 / 2.0;
  t39 = this->robot_dim[0] * this->robot_inertial_par[3] *
        this->robot_inertial_par[8] * t10;
  t151 = t39 * t11 * t14;
  t152 = t39 * t12 * t14;
  t39 = this->robot_inertial_par[3] * this->robot_dim[4] *
        this->robot_inertial_par[8] * t10;
  t153 = t39 * t11 * t14;
  t154 = t39 * t12 * t14;
  t172 = t202 * t13 * t37 / 2.0;
  t183 = h1 * this->robot_dim[3] * this->robot_dim[4] *
         this->robot_inertial_par[9] * t13 * t23 * -0.5;
  t70 = h2 * t66;
  t74 = t69 * 2.0;
  t77 = h1 * t69;
  t79 = this->g * this->robot_inertial_par[4] * this->robot_inertial_par[9] *
        std::cos(phi_r + -this->robot_dim[6]);
  t87 = this->g * this->robot_dim[0] * this->robot_inertial_par[8] * t65;
  t88 =
      this->g * this->robot_inertial_par[2] * this->robot_inertial_par[6] * t65;
  t39 = this->g * this->robot_dim[4];
  t89 = t39 * this->robot_inertial_par[8] * t65;
  t90 = t39 * this->robot_inertial_par[9] * t65;
  t91 = t39 * this->robot_inertial_par[7] * t65;
  t92 = h1 * t83;
  t105 = h1 * t73;
  t108 = h1 * t75;
  t111 = t67 / 2.0;
  t155 = t11 * t82;
  t156 = t11 * t83;
  t157 = t82 / 2.0;
  t160 = h1 * t152;
  t161 = h2 * t151;
  t162 = h1 * t154;
  t163 = h2 * t153;
  t164 = t125 * this->robot_inertial_par[4] * this->robot_dim[4] *
         this->robot_inertial_par[9] * phi_r_dot * t64;
  t166 = h1 * t133;
  t167 = t106 / 2.0;
  t182 = h1 * t82 * -0.5;
  t39 = this->robot_inertial_par[4] * this->robot_dim[4] *
        this->robot_inertial_par[9] * t10 * t11 * t64;
  t23 = h1 * this->robot_inertial_par[4] * this->robot_dim[4] *
        this->robot_inertial_par[9] * t13 * t64 / 2.0;
  t65 = h1 * t87;
  t125 = h1 * t88;
  t64 = h1 * t89;
  t10 = h1 * t91;
  t13 = t79 / 2.0;
  t159 = t90 / 2.0;
  t165 = h1 * t157;
  t194 = t11 * -t83;
  t202 = h1 * t90 * -0.5;
  t14 = t155 / 2.0;
  t37 = t39 / 2.0;
  t12 = h1 * t159;
  t39 =
      h1 *
      ((((((((((((((((((((((((((((t67 + t66 * 2.0) + t68 * 2.0) + t74) + t76) +
                              t82) +
                             t97_tmp * 4.0) +
                            t98_tmp * 4.0) +
                           t101_tmp * 4.0) +
                          t102_tmp * 4.0) +
                         h1 * t74) +
                        h1 * t76) +
                       t83 * 2.0) +
                      -(t78 * 2.0)) +
                     -t79) +
                    -(t87 * 2.0)) +
                   -(t88 * 2.0)) +
                  -(t89 * 2.0)) +
                 -t90) +
                -(t91 * 2.0)) +
               -t106) +
              t152 * 2.0) +
             t154 * 2.0) +
            -(t151 * 2.0)) +
           -(t153 * 2.0)) +
          -t155) +
         t156 * -2.0) +
        -(t164 * 2.0)) +
       t39) /
      2.0;
  b_u_lim[0] =
      ((((((((((((((((((((((((((((u_lim[1] + t78) + -t66) + -t68) + -t69) +
                              -t75) +
                             t87) +
                            t88) +
                           t89) +
                          t91) +
                         -t97) +
                        -t98) +
                       -t101) +
                      -t102) +
                     -t77) +
                    t123) +
                   -t111) +
                  -t83) +
                 t13) +
                t151) +
               t153) +
              t156) +
             t159) +
            t164) +
           t167) +
          -t157) +
         -t152) +
        -t154) +
       t14) +
      -t37;
  b_u_lim[1] =
      ((((((((((((((((((((((((((((((((((((((((((((((u_lim[3] + t5) + t6) + t7) +
                                                 t16) +
                                                t19) +
                                               t20) +
                                              t21) +
                                             t24) +
                                            t25) +
                                           t26) +
                                          -t4) +
                                         -t8) +
                                        -t9) +
                                       -t15) +
                                      -t17) +
                                     -t18) +
                                    -t22) +
                                   -t27) +
                                  -t28) +
                                 -t29) +
                                t73) +
                               -t68) +
                              -t70) +
                             t92) +
                            t93) +
                           t94) +
                          t95) +
                         t96) +
                        t105) +
                       -t111) +
                      t133) +
                     t156) +
                    t160) +
                   t161) +
                  t162) +
                 t163) +
                t165) +
               t166) +
              -t65) +
             -t125) +
            -t64) +
           -t10) +
          -t172) +
         t202) +
        t14) +
       -t23) +
      -t39;
  b_u_lim[2] =
      ((((((((((((((((((((((((((((u_lim[1] + t66) + t68) + t69) + t75) + t77) +
                             t83) +
                            t97) +
                           t98) +
                          t101) +
                         t102) +
                        t108) +
                       t111) +
                      -t78) +
                     -t87) +
                    -t88) +
                   -t89) +
                  -t91) +
                 t152) +
                t154) +
               t157) +
              -t13) +
             -t159) +
            -t167) +
           -t151) +
          -t153) +
         t194) +
        -t164) +
       -t14) +
      t37;
  b_u_lim[3] =
      ((((((((((((((((((((((((((((((((((((((((((((((u_lim[3] + t4) + t8) + t9) +
                                                 t15) +
                                                t17) +
                                               t18) +
                                              t22) +
                                             t27) +
                                            t28) +
                                           t29) +
                                          -t5) +
                                         -t6) +
                                        -t7) +
                                       -t16) +
                                      -t19) +
                                     -t20) +
                                    -t21) +
                                   -t24) +
                                  -t25) +
                                 -t26) +
                                t68) +
                               t70) +
                              t111) +
                             -t73) +
                            t65) +
                           t125) +
                          t64) +
                         t10) +
                        -t92) +
                       -t93) +
                      -t94) +
                     -t95) +
                    -t96) +
                   -t105) +
                  -t133) +
                 t172) +
                t182) +
               t183) +
              t12) +
             t194) +
            -t160) +
           -t161) +
          -t162) +
         -t163) +
        -t14) +
       t23) +
      t39;
  b_u_lim[4] =
      ((((((((((((((((((((((((((((-u_lim[0] + t66) + t68) + t69) + t75) + t77) +
                             t83) +
                            t97) +
                           t98) +
                          t101) +
                         t102) +
                        t108) +
                       t111) +
                      -t78) +
                     -t87) +
                    -t88) +
                   -t89) +
                  -t91) +
                 t152) +
                t154) +
               t157) +
              -t13) +
             -t159) +
            -t167) +
           -t151) +
          -t153) +
         t194) +
        -t164) +
       -t14) +
      t37;
  b_u_lim[5] =
      ((((((((((((((((((((((((((((((((((((((((((((((t4 + t8) + t9) + t15) +
                                                 t17) +
                                                t18) +
                                               t22) +
                                              t27) +
                                             t28) +
                                            t29) +
                                           -u_lim[2]) +
                                          -t5) +
                                         -t6) +
                                        -t7) +
                                       -t16) +
                                      -t19) +
                                     -t20) +
                                    -t21) +
                                   -t24) +
                                  -t25) +
                                 -t26) +
                                t68) +
                               t70) +
                              t111) +
                             -t73) +
                            t65) +
                           t125) +
                          t64) +
                         t10) +
                        -t92) +
                       -t93) +
                      -t94) +
                     -t95) +
                    -t96) +
                   -t105) +
                  -t133) +
                 t172) +
                t182) +
               t183) +
              t12) +
             t194) +
            -t160) +
           -t161) +
          -t162) +
         -t163) +
        -t14) +
       t23) +
      t39;
  b_u_lim[6] =
      ((((((((((((((((((((((((((((-u_lim[0] + t78) + -t66) + -t68) + -t69) +
                              -t75) +
                             t87) +
                            t88) +
                           t89) +
                          t91) +
                         -t97) +
                        -t98) +
                       -t101) +
                      -t102) +
                     -t77) +
                    t123) +
                   -t111) +
                  -t83) +
                 t13) +
                t151) +
               t153) +
              t156) +
             t159) +
            t164) +
           t167) +
          -t157) +
         -t152) +
        -t154) +
       t14) +
      -t37;
  b_u_lim[7] =
      ((((((((((((((((((((((((((((((((((((((((((((((t5 + t6) + t7) + t16) +
                                                 t19) +
                                                t20) +
                                               t21) +
                                              t24) +
                                             t25) +
                                            t26) +
                                           -u_lim[2]) +
                                          -t4) +
                                         -t8) +
                                        -t9) +
                                       -t15) +
                                      -t17) +
                                     -t18) +
                                    -t22) +
                                   -t27) +
                                  -t28) +
                                 -t29) +
                                t73) +
                               -t68) +
                              -t70) +
                             t92) +
                            t93) +
                           t94) +
                          t95) +
                         t96) +
                        t105) +
                       -t111) +
                      t133) +
                     t156) +
                    t160) +
                   t161) +
                  t162) +
                 t163) +
                t165) +
               t166) +
              -t65) +
             -t125) +
            -t64) +
           -t10) +
          -t172) +
         t202) +
        t14) +
       -t23) +
      -t39;
}

//
// File trailer for ModelController.cpp
//
// [EOF]
//
