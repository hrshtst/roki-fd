/* RoKi-FD - Robot Kinetics library: forward dynamics extention
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rkfd_sim - forward dynamics simulation
 * contributer: 2014- Naoki Wakisaka
 */

#ifndef __RKFD_SIM_H__
#define __RKFD_SIM_H__

#include <roki_fd/rkfd_property.h>
#include <roki_fd/rkfd_chain.h>
#include <roki_fd/rkfd_cd.h>
#include <roki_fd/rkfd_solver.h>

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: rkFD
   forward dynamics class
 * ********************************************************** */
typedef struct{
  rkFDChain fc;    /* derived by rkChain */
  int _offset;     /* offset */
  zVecStruct _dis; /* displacement */
  zVecStruct _vel; /* velocity */
  zVecStruct _acc; /* acceleration */
} rkFDCellDat;

ZEDA_DEF_LIST_CLASS( rkFDCellList, rkFDCell, rkFDCellDat );

#define rkFDCellDatChain(d) rkFDChainBase(&(d)->fc)
#define rkFDCellChain(c)    rkFDCellDatChain(&(c)->data)

/* forward dynamics simulator */
typedef struct _rkFD{
  double t;
  rkFDPrp prp;
  rkFDSolver solver;     /* solver */
  rkFDCellList list;     /* chain list */
  rkContactInfoArray ci; /* contact information list */
  rkContactInfo cidef;   /* default contact infomation */
  rkFDCD cd;             /* contact maneger */

  zODE2 ode;
  int ode_step;
  zVec dis, vel;        /* total joint state */
  zVec acc;
  int size;             /* total joint size */
} rkFD;

#define rkFDTime(f)   (f)->t
#define rkFDDT(f)     (f)->prp.dt
#define rkFDGetPrp(f) ( &(f)->prp )

__ROKI_FD_EXPORT rkFD *rkFDCreate(rkFD *fd);
__ROKI_FD_EXPORT void rkFDDestroy(rkFD *fd);

__ROKI_FD_EXPORT rkFDCell *rkFDChainReg(rkFD *fd, rkChain *chain);
__ROKI_FD_EXPORT rkFDCell *rkFDChainRegFile(rkFD *fd, char filename[]);
__ROKI_FD_EXPORT bool rkFDChainUnreg(rkFD *fd, rkFDCell *cell);

__ROKI_FD_EXPORT void rkFDFK(rkFD *fd, zVec dis);
__ROKI_FD_EXPORT void rkFDUpdateRate(rkFD *fd, zVec vel, zVec acc);
__ROKI_FD_EXPORT void rkFDUpdateFKRate(rkFD *fd);

/*! \brief joint range limits as compliant stops.
 *
 * rkJointSetDis() clamps the displacement of a joint at its range limit
 * when the state is written to the chain, while the state the forward
 * dynamics integrates is constrained by nothing.  The forward dynamics
 * therefore treats a limit as a compliant stop, in the manner of the
 * penalty contact: the integrated state may penetrate the limit, which
 * the chain does not show, and the penetration produces a reaction
 * torque that enters the driving torque of the joint.
 *
 * rkFDJointLimitTrq() is that torque for a displacement \a dis and a
 * velocity \a vel of a joint ranging from \a min to \a max: a spring of
 * stiffness \a k and a damper of coefficient \a c on the penetration,
 * which only pushes the joint back into the range and never pulls it
 * toward the limit.  rkFDChainJointLimitTrq() stores it for every
 * one-degree-of-freedom joint of \a chain in \a trq, indexed by joint
 * identifier; \a k = 0 disables the stops.
 *
 * rkFDChainStopJointsAtLimits() bounds the penetration: it projects the
 * joint state \a dis and \a vel of \a chain onto the range limits widened
 * by \a width, setting a displacement beyond a bound to the bound and
 * zeroing the velocity that points beyond it, i.e. the stop becomes a
 * perfectly inelastic collision there.  With \a width = 0 the stop is that
 * collision at the limit itself.
 *
 * The forward dynamics applies both to every state it connects to the
 * chains; the stiffness, the damping and the width are properties of
 * rkFDPrp (rkFDPrpSetJointLimitK(), rkFDPrpSetJointLimitC() and
 * rkFDPrpSetJointLimitWidth()).  So that the chain shows the
 * penetration, and its positions and velocities agree, the joint range
 * limits are widened by the width while the state is written to the
 * chain and restored right after: rkFDChainJointLimitStore() stores the
 * limits of the one-degree-of-freedom joints of \a chain in \a min and
 * \a max at registration (a limit changed afterwards is not seen), and
 * rkFDChainJointLimitWiden() sets them widened by \a width (restored
 * with \a width = 0).
 * \return
 * rkFDJointLimitTrq() returns the torque.
 * rkFDChainJointLimitTrq() returns \a trq.
 * rkFDChainStopJointsAtLimits() returns the number of joints whose
 * velocity is stopped.
 */
__ROKI_FD_EXPORT int rkFDChainStopJointsAtLimits(rkChain *chain, zVec dis, zVec vel, double width);
__ROKI_FD_EXPORT double rkFDJointLimitTrq(double dis, double vel, double min, double max, double k, double c);
__ROKI_FD_EXPORT zVec rkFDChainJointLimitTrq(rkChain *chain, zVec dis, zVec vel, double k, double c, zVec trq);
__ROKI_FD_EXPORT void rkFDChainJointLimitStore(rkChain *chain, zVec min, zVec max);
__ROKI_FD_EXPORT void rkFDChainJointLimitWiden(rkChain *chain, zVec min, zVec max, double width);

__ROKI_FD_EXPORT void rkFDChainSetDis(rkFDCell *lc, zVec dis);
__ROKI_FD_EXPORT void rkFDChainSetVel(rkFDCell *lc, zVec vel);
__ROKI_FD_EXPORT bool rkFDContactInfoScanFile(rkFD *fd, char filename[]);

/* for a fake-crawler */
__ROKI_FD_EXPORT void rkFDCDCellSetSlideMode(rkCDCell *cell, bool mode);
__ROKI_FD_EXPORT void rkFDCDCellSetSlideVel(rkCDCell *cell, double vel);
__ROKI_FD_EXPORT void rkFDCDCellSetSlideAxis(rkCDCell *cell, zVec3D *axis);
__ROKI_FD_EXPORT rkCDCell *rkFDShape3DGetCDCell(rkFD *fd, zShape3D *shape);
__ROKI_FD_EXPORT rkCDCell *rkFDShape3DSetSlideMode(rkFD *fd, zShape3D *shape, bool mode);
__ROKI_FD_EXPORT rkCDCell *rkFDShape3DSetSlideVel(rkFD *fd, zShape3D *shape, double vel);
__ROKI_FD_EXPORT rkCDCell *rkFDShape3DSetSlideAxis(rkFD *fd, zShape3D *shape, zVec3D *axis);

/* ODE updater */
__ROKI_FD_EXPORT zVec rkFDODECatDefault(zVec x, double k, zVec v, zVec xnew, void *util);
__ROKI_FD_EXPORT zVec rkFDODESubDefault(zVec x1, zVec x2, zVec dx, void *util);
#define rkFDODE2Assign(f,t)        zODE2Assign( &(f)->ode, t, rkFDODECatDefault, NULL, rkFDODESubDefault, NULL )
#define rkFDODE2AssignRegular(f,t) zODE2AssignRegular( &(f)->ode, t )

/* solver */
#define rkFDSetSolver(f,type) do{                                 \
    rkFDSolverReset( &(f)->solver );                              \
    rkFDSolverCreate_##type( &(f)->solver );                      \
    rkFDSolverGetDefaultContactInfo( &(f)->solver, &(f)->cidef ); \
  } while(0)

__ROKI_FD_EXPORT rkFD *rkFDSolve(rkFD *fd);
__ROKI_FD_EXPORT void rkFDUpdateInit(rkFD *fd);
__ROKI_FD_EXPORT rkFD *rkFDUpdate(rkFD *fd);
__ROKI_FD_EXPORT void rkFDUpdateDestroy(rkFD *fd);
#define rkFDSolveContact(f,b) (f)->_solve_contact(f,b)

/* for debug */
__ROKI_FD_EXPORT void rkFDPrint(rkFD *fd);

__END_DECLS

#endif /* __RKFD_SIM_H__ */
