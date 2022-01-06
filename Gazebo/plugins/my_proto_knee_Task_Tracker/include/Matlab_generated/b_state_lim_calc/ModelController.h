//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ModelController.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 06-Aug-2021 10:18:26
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
  void b_state_lim_calc(double delta, const double phi_r_lim[2],
                        const double beta_lim[2], double phi_r, double beta,
                        double phi_r_dot, double beta_dot,
                        double b_state_lim[4]);
};

#endif
//
// File trailer for ModelController.h
//
// [EOF]
//
