/* RoKi-FD - Robot Kinetics library: forward dynamics extention
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rkfd_defs - definitions of default parameters
 * contributer: 2014- Naoki Wakisaka
 */

#ifndef __RKFD_DEFS_H__
#define __RKFD_DEFS_H__

#include <zeda/zeda.h>

#include <roki_fd/roki_fd_export.h>

#define RK_FD_DT_DEFAULT                        0.001
#define RK_FD_CDT_DEFAULT                       0.002
#define RK_FD_JOINT_COMP_K_DEFAULT            100
#define RK_FD_JOINT_COMP_L_DEFAULT              0.01
#define RK_FD_KINETIC_FRICTION_WEIGHT_DEFAULT 100
#define RK_FD_FRICTION_PYRAMID_ORDER_DEFAULT    8

#define RK_FD_MAX_ITER_DEFAULT                 10
#define RK_FD_VEL_EPSILON_DEFAULT              (1.0e-8)

/* joint range limits as compliant stops: stiffness [Nm/rad] or [N/m],
 * damping [Nms/rad] or [Ns/m], and the largest penetration of the
 * integrated state beyond a limit [rad] or [m], beyond which the state
 * is projected back onto the limit */
#define RK_FD_JOINT_LIMIT_K_DEFAULT           100
#define RK_FD_JOINT_LIMIT_C_DEFAULT             1
#define RK_FD_JOINT_LIMIT_WIDTH_DEFAULT         0.0872664626 /* 5 degrees */

#endif /* __RKFD_DEFS_H__ */
