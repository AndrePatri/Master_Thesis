//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ModelController.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 11-Aug-2021 12:59:32
//

#ifndef MODELCONTROLLER_H
#define MODELCONTROLLER_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Definitions
class ModelController {
public:
  ModelController();
  ~ModelController();
  void p_cm_tot_rel_wheel_calc(double alpha, double b_gamma, double phi_r,
                               double beta, double p_cm_tot_rel_wheel[2]);
  double robot_dim[7];
  double robot_inertial_par[15];
};

#endif
//
// File trailer for ModelController.h
//
// [EOF]
//
