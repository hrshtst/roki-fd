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

/* ZTK */

/*! \brief joint range limit stops in a ZTK file.
 *
 * The compliant stops of the joint range limits are configured by a
 * section
 *   [roki_fd::jointlimit]
 *   stiffness : 100
 *   damping : 1
 *   width : 5
 * where the width is in degrees (or meters for prismatic joints), as the
 * joint range limits in chain files; a zero stiffness disables the
 * compliant stop, leaving the perfectly inelastic one at the limit.
 * Keys not given keep their values.  rkFDPrpFromZTK() evaluates such a
 * section of a ZTK format processor \a ztk into \a prp, and
 * rkFDPrpReadZTK() reads it from a file \a filename; a file without the
 * section leaves \a prp untouched.  rkFDContactInfoScanFile() evaluates
 * the section of the contact information file as well.
 * \return
 * rkFDPrpFromZTK() and rkFDPrpReadZTK() return \a prp.
 */
#define ZTK_TAG_ROKI_FD_JOINTLIMIT           "roki_fd::jointlimit"
#define ZTK_KEY_ROKI_FD_JOINTLIMIT_STIFFNESS "stiffness"
#define ZTK_KEY_ROKI_FD_JOINTLIMIT_DAMPING   "damping"
#define ZTK_KEY_ROKI_FD_JOINTLIMIT_WIDTH     "width"

__ROKI_FD_EXPORT rkFDPrp *rkFDPrpFromZTK(rkFDPrp *prp, ZTK *ztk);
__ROKI_FD_EXPORT rkFDPrp *rkFDPrpReadZTK(rkFDPrp *prp, const char filename[]);

__END_DECLS

#endif /* __RKFD_PROPERTY_H__ */
