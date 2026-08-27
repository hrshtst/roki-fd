/* RoKi-FD - Robot Kinetics library: forward dynamics extention
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rkfd_property - properties for forward dynamics simulation
 * contributer: 2014- Naoki Wakisaka
 */

#ifndef __RKFD_PROPERTY_H__
#define __RKFD_PROPERTY_H__

#include <roki_fd/rkfd_defs.h>

__BEGIN_DECLS

typedef struct{
  double dt;
  int pyramid;
  double friction_weight;

  int max_iter;
  double vel_eps;
  double joint_limit_k;     /* stiffness of the joint range limits */
  double joint_limit_c;     /* damping of the joint range limits */
  double joint_limit_width; /* largest penetration of a joint range limit */
} rkFDPrp;

#define rkFDPrpDT(p)             (p)->dt
#define rkFDPrpPyramid(p)        (p)->pyramid
#define rkFDPrpFrictionWeight(p) (p)->friction_weight
#define rkFDPrpMaxIter(p)        (p)->max_iter
#define rkFDPrpVelEps(p)         (p)->vel_eps
#define rkFDPrpJointLimitK(p)     (p)->joint_limit_k
#define rkFDPrpJointLimitC(p)     (p)->joint_limit_c
#define rkFDPrpJointLimitWidth(p) (p)->joint_limit_width

#define rkFDPrpSetDT(f,t)             ( (f)->prp.dt = (t) )
#define rkFDPrpSetPyramid(f,n)        ( (f)->prp.pyramid = (n) )
#define rkFDPrpSetFrictionWeight(f,w) ( (f)->prp.friction_weight = (w) )
#define rkFDPrpSetMaxIter(f,i)        ( (f)->prp.max_iter = (i) )
#define rkFDPrpSetVelEps(f,e)         ( (f)->prp.vel_eps = (e) )
#define rkFDPrpSetJointLimitK(f,k)     ( (f)->prp.joint_limit_k = (k) )
#define rkFDPrpSetJointLimitC(f,c)     ( (f)->prp.joint_limit_c = (c) )
#define rkFDPrpSetJointLimitWidth(f,w) ( (f)->prp.joint_limit_width = (w) )

__ROKI_FD_EXPORT bool rkFDPrpInit(rkFDPrp *prp);
__ROKI_FD_EXPORT void rkFDPrpDestroy(rkFDPrp *prp);

__END_DECLS

#endif /* __RKFD_PROPERTY_H__ */
