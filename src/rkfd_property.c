/* RoKi-FD - Robot Kinetics library: forward dynamics extention
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rkfd_property - properties for forward dynamics simulation
 * contributer: 2014- Naoki Wakisaka
 */

#include <roki_fd/rkfd_property.h>
#include <zm/zm_misc.h>

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

/* ZTK */

static void *_rkFDPrpJointLimitStiffnessFromZTK(void *obj, int i, void *arg, ZTK *ztk){
  ((rkFDPrp *)obj)->joint_limit_k = ZTKDouble(ztk);
  return obj;
}
static void *_rkFDPrpJointLimitDampingFromZTK(void *obj, int i, void *arg, ZTK *ztk){
  ((rkFDPrp *)obj)->joint_limit_c = ZTKDouble(ztk);
  return obj;
}
static void *_rkFDPrpJointLimitWidthFromZTK(void *obj, int i, void *arg, ZTK *ztk){
  ((rkFDPrp *)obj)->joint_limit_width = zDeg2Rad( ZTKDouble(ztk) );
  return obj;
}

static const ZTKPrp __ztk_prp_rkfd_jointlimit[] = {
  { ZTK_KEY_ROKI_FD_JOINTLIMIT_STIFFNESS, 1, _rkFDPrpJointLimitStiffnessFromZTK, NULL },
  { ZTK_KEY_ROKI_FD_JOINTLIMIT_DAMPING,   1, _rkFDPrpJointLimitDampingFromZTK, NULL },
  { ZTK_KEY_ROKI_FD_JOINTLIMIT_WIDTH,     1, _rkFDPrpJointLimitWidthFromZTK, NULL },
};

static void *_rkFDPrpJointLimitFromZTK(void *obj, int i, void *arg, ZTK *ztk)
{
  if( !_ZTKEvalKey( obj, arg, ztk, __ztk_prp_rkfd_jointlimit ) ) return NULL;
  return obj;
}

static const ZTKPrp __ztk_prp_tag_rkfd_jointlimit[] = {
  { ZTK_TAG_ROKI_FD_JOINTLIMIT, 1, _rkFDPrpJointLimitFromZTK, NULL },
};

rkFDPrp *rkFDPrpFromZTK(rkFDPrp *prp, ZTK *ztk)
{
  if( ZTKCountTag( ztk, ZTK_TAG_ROKI_FD_JOINTLIMIT ) == 0 ) return prp;
  _ZTKEvalTag( prp, NULL, ztk, __ztk_prp_tag_rkfd_jointlimit );
  return prp;
}

rkFDPrp *rkFDPrpReadZTK(rkFDPrp *prp, const char filename[])
{
  ZTK ztk;

  ZTKInit( &ztk );
  if( ZTKParse( &ztk, filename ) )
    prp = rkFDPrpFromZTK( prp, &ztk );
  else
    prp = NULL;
  ZTKDestroy( &ztk );
  return prp;
}
