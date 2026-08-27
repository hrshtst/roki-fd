/* RoKi-FD - Robot Kinetics library: forward dynamics extention
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rkfd_property - properties for forward dynamics simulation
 * contributer: 2014- Naoki Wakisaka
 */

#include <roki_fd/rkfd_property.h>

bool rkFDPrpInit(rkFDPrp *prp)
{
  prp->dt              = RK_FD_DT_DEFAULT;
  prp->pyramid         = RK_FD_FRICTION_PYRAMID_ORDER_DEFAULT;
  prp->friction_weight = RK_FD_KINETIC_FRICTION_WEIGHT_DEFAULT;
  prp->max_iter        = RK_FD_MAX_ITER_DEFAULT;
  prp->vel_eps         = RK_FD_VEL_EPSILON_DEFAULT;
  prp->joint_limit_k     = RK_FD_JOINT_LIMIT_K_DEFAULT;
  prp->joint_limit_c     = RK_FD_JOINT_LIMIT_C_DEFAULT;
  prp->joint_limit_width = RK_FD_JOINT_LIMIT_WIDTH_DEFAULT;
  return true;
}

void rkFDPrpDestroy(rkFDPrp *prp){}
