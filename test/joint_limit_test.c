#include <roki_fd/roki_fd.h>

/* Joint range limits as compliant stops.
 *
 * A pendulum (a rod hinged at its end, pointing horizontally at rest)
 * falls under gravity against the lower range limit of its joint.
 * roki clamps the joint displacement at the limit when the state is
 * written to the chain; the forward dynamics must follow the
 * unconstrained motion until the limit, must stop the joint there
 * without letting velocity or displacement run away in the integrated
 * state, must react with the torque of the stop, and must release the
 * joint when it moves away.  A double pendulum whose first joint is
 * locked by a two-sided stop shows that the reaction torque of a stop
 * is transmitted through the chain. */

#define DT       0.001
#define LIMIT    zDeg2Rad(-30.0)
#define TOL      1.0e-9
#define G        9.8

typedef struct{
  rkFD fd;
  rkFDCell *cell;
  zVec dis, vel;
} pendulum_t;

static bool pendulum_create(pendulum_t *p, const char *filename, double k)
{
  rkFDCreate( &p->fd );
  if( !( p->cell = rkFDChainRegFile( &p->fd, (char *)filename ) ) ) return false;
  p->dis = zVecAlloc( rkChainJointSize( rkFDCellChain(p->cell) ) );
  p->vel = zVecAlloc( rkChainJointSize( rkFDCellChain(p->cell) ) );
  rkFDChainSetDis( p->cell, p->dis );
  rkFDChainSetVel( p->cell, p->vel );
  rkFDODE2Assign( &p->fd, Regular );
  rkFDODE2AssignRegular( &p->fd, RK4 );
  rkFDPrpSetDT( &p->fd, DT );
  rkFDPrpSetJointLimitK( &p->fd, k );
  rkFDSetSolver( &p->fd, Vert );
  rkFDUpdateInit( &p->fd );
  return true;
}

static void pendulum_destroy(pendulum_t *p)
{
  rkFDUpdateDestroy( &p->fd );
  zVecFree( p->dis );
  zVecFree( p->vel );
  rkFDDestroy( &p->fd );
}

/* the chain's joint state (the displacement clamped at the limits) */
static void pendulum_step(pendulum_t *p)
{
  rkFDUpdate( &p->fd );
  rkChainGetJointDisAll( rkFDCellChain(p->cell), p->dis );
  rkChainGetJointVelAll( rkFDCellChain(p->cell), p->vel );
}

/* the integrated displacement of a joint, which may penetrate a limit */
static double pendulum_state_dis(pendulum_t *p, int i)
{
  return zVecElem( p->fd.dis, i );
}

/* the projection onto the bounds */
static void assert_stop_function(void)
{
  pendulum_t p;
  zVec dis, vel;
  int n;

  pendulum_create( &p, "pendulum.ztk", 0 );
  dis = zVecAlloc( 1 );
  vel = zVecAlloc( 1 );
  /* inside the range: untouched */
  zVecSetElem( dis, 0, 0.0 );
  zVecSetElem( vel, 0, -1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel, 0 );
  zAssert( rkFDChainStopJointsAtLimits (inside the range), n == 0 && zVecElem(vel,0) == -1.0 );
  /* at the lower limit, moving beyond it: stopped */
  zVecSetElem( dis, 0, LIMIT );
  zVecSetElem( vel, 0, -1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel, 0 );
  zAssert( rkFDChainStopJointsAtLimits (at the lower limit, moving beyond), n == 1 && zVecElem(vel,0) == 0.0 );
  /* at the lower limit, moving away from it: untouched */
  zVecSetElem( vel, 0, 1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel, 0 );
  zAssert( rkFDChainStopJointsAtLimits (at the lower limit, moving away), n == 0 && zVecElem(vel,0) == 1.0 );
  /* at the upper limit, moving beyond it: stopped */
  zVecSetElem( dis, 0, zDeg2Rad(90.0) );
  zVecSetElem( vel, 0, 1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel, 0 );
  zAssert( rkFDChainStopJointsAtLimits (at the upper limit, moving beyond), n == 1 && zVecElem(vel,0) == 0.0 );
  /* beyond the lower limit, moving beyond it: projected onto the limit and stopped */
  zVecSetElem( dis, 0, LIMIT - 0.1 );
  zVecSetElem( vel, 0, -1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel, 0 );
  zAssert( rkFDChainStopJointsAtLimits (beyond the lower limit, moving beyond), n == 1 && zVecElem(dis,0) == LIMIT && zVecElem(vel,0) == 0.0 );
  /* beyond the lower limit, moving away from it: projected onto the limit, velocity kept */
  zVecSetElem( dis, 0, LIMIT - 0.1 );
  zVecSetElem( vel, 0, 1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel, 0 );
  zAssert( rkFDChainStopJointsAtLimits (beyond the lower limit, moving away), n == 0 && zVecElem(dis,0) == LIMIT && zVecElem(vel,0) == 1.0 );
  /* within the width of the stop: untouched */
  zVecSetElem( dis, 0, LIMIT - 0.05 );
  zVecSetElem( vel, 0, -1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel, 0.1 );
  zAssert( rkFDChainStopJointsAtLimits (within the width of the stop), n == 0 && zVecElem(dis,0) == LIMIT - 0.05 && zVecElem(vel,0) == -1.0 );
  /* beyond the width of the stop: projected onto the bound and stopped */
  zVecSetElem( dis, 0, LIMIT - 0.2 );
  zVecSetElem( vel, 0, -1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel, 0.1 );
  zAssert( rkFDChainStopJointsAtLimits (beyond the width of the stop), n == 1 && zVecElem(dis,0) == LIMIT - 0.1 && zVecElem(vel,0) == 0.0 );
  zVecFree( dis );
  zVecFree( vel );
  pendulum_destroy( &p );
}

/* the reaction torque of a stop */
static void assert_limit_trq(void)
{
  double k = 100, c = 1, min = -0.5, max = 0.5;

  zAssert( rkFDJointLimitTrq (inside the range), rkFDJointLimitTrq( 0.0, -1.0, min, max, k, c ) == 0.0 );
  zAssert( rkFDJointLimitTrq (beyond the lower limit, approaching), zIsTiny( rkFDJointLimitTrq( min - 0.01, -1.0, min, max, k, c ) - ( k * 0.01 + c * 1.0 ) ) );
  zAssert( rkFDJointLimitTrq (beyond the lower limit, receding slowly), zIsTiny( rkFDJointLimitTrq( min - 0.01, 0.5, min, max, k, c ) - ( k * 0.01 - c * 0.5 ) ) );
  zAssert( rkFDJointLimitTrq (beyond the lower limit, receding fast: never pulls), rkFDJointLimitTrq( min - 0.01, 10.0, min, max, k, c ) == 0.0 );
  zAssert( rkFDJointLimitTrq (beyond the upper limit, approaching), zIsTiny( rkFDJointLimitTrq( max + 0.01, 1.0, min, max, k, c ) + ( k * 0.01 + c * 1.0 ) ) );
  zAssert( rkFDJointLimitTrq (beyond the upper limit, receding fast: never pulls), rkFDJointLimitTrq( max + 0.01, -10.0, min, max, k, c ) == 0.0 );
}

/* the pendulum on the compliant stop (the default) */
static void assert_pendulum(void)
{
  pendulum_t p, q;
  double t, t_hit = -1, max_dev = 0, v_hit = 0, max_vel_after = 0, max_pen = 0, d, k, trq;
  int i;

  pendulum_create( &p, "pendulum.ztk", RK_FD_JOINT_LIMIT_K_DEFAULT );
  pendulum_create( &q, "pendulum_free.ztk", RK_FD_JOINT_LIMIT_K_DEFAULT );
  k = rkFDPrpJointLimitK( &p.fd.prp );
  for( t=0, i=0; t<2.0; t+=DT, i++ ){
    pendulum_step( &p );
    pendulum_step( &q );
    if( t_hit < 0 ){
      /* until the limit is reached, the two motions coincide */
      d = fabs( zVecElem(p.dis,0) - zVecElem(q.dis,0) );
      if( zVecElem(p.dis,0) <= LIMIT + TOL ){
        t_hit = t;
        v_hit = fabs( zVecElem(p.vel,0) );
      } else if( d > max_dev ){
        max_dev = d;
      }
    } else {
      /* the stop is compliant: the rod still accelerates for a moment
         until the spring balances gravity, then it is caught */
      if( fabs( zVecElem(p.vel,0) ) > max_vel_after ) max_vel_after = fabs( zVecElem(p.vel,0) );
      if( LIMIT - pendulum_state_dis( &p, 0 ) > max_pen ) max_pen = LIMIT - pendulum_state_dis( &p, 0 );
    }
  }
  zAssert( compliant stop (the limit is reached), t_hit > 0 && t_hit < 0.5 );
  zAssert( compliant stop (the free motion is unaltered before the limit), max_dev < TOL );
  zAssert( compliant stop (no velocity accumulates at the stop), max_vel_after < 1.5 * v_hit );
  zAssert( compliant stop (the penetration stays within the width), max_pen < rkFDPrpJointLimitWidth( &p.fd.prp ) );
  zAssert( compliant stop (the joint comes to rest on the stop), zVecElem(p.dis,0) < LIMIT && zVecElem(p.dis,0) > LIMIT - rkFDPrpJointLimitWidth( &p.fd.prp ) && fabs( zVecElem(p.vel,0) ) < 1.0e-6 );
  /* the chain shows the penetration */
  zAssert( compliant stop (the chain shows the penetration), fabs( zVecElem(p.dis,0) - pendulum_state_dis( &p, 0 ) ) < TOL );
  /* the penetration carries the gravity torque of the rod at its rest angle */
  trq = k * ( LIMIT - pendulum_state_dis( &p, 0 ) );
  zAssert( compliant stop (the stop carries the load), fabs( trq - G * 0.1 * cos( zVecElem(p.dis,0) ) ) < 1.0e-3 * G * 0.1 );
  zAssert( compliant stop (the penetration is within the width), LIMIT - pendulum_state_dis( &p, 0 ) < rkFDPrpJointLimitWidth( &p.fd.prp ) );
  /* the unconstrained pendulum keeps swinging far beyond the limit */
  zAssert( compliant stop (the free pendulum swings beyond the limit), zVecElem(q.dis,0) < LIMIT - 0.1 || fabs(zVecElem(q.vel,0)) > 0.1 );
  /* release: a velocity away from the stop moves the joint off the limit */
  zVecSetElem( p.vel, 0, 2.0 );
  rkFDChainSetVel( p.cell, p.vel );
  for( i=0; i<20; i++ ) pendulum_step( &p );
  zAssert( compliant stop (the joint leaves the limit when it moves away), zVecElem(p.dis,0) > LIMIT + 0.01 );

  pendulum_destroy( &p );
  pendulum_destroy( &q );
}

/* the pendulum on the inelastic stop (zero stiffness) */
static void assert_hard_stop(void)
{
  pendulum_t p;
  double t, t_hit = -1, max_vel_after = 0;
  int i, n_off = 0;

  pendulum_create( &p, "pendulum.ztk", 0 );
  for( t=0, i=0; t<1.0; t+=DT, i++ ){
    pendulum_step( &p );
    if( t_hit < 0 ){
      if( zVecElem(p.dis,0) <= LIMIT + TOL ) t_hit = t;
    } else {
      if( fabs( zVecElem(p.dis,0) - LIMIT ) > TOL ) n_off++;
      if( fabs( zVecElem(p.vel,0) ) > max_vel_after ) max_vel_after = fabs( zVecElem(p.vel,0) );
    }
  }
  zAssert( inelastic stop (the joint stays at the limit), t_hit > 0 && n_off == 0 );
  zAssert( inelastic stop (no velocity accumulates at the stop), max_vel_after < TOL );
  zAssert( inelastic stop (the integrated state stays at the limit), fabs( pendulum_state_dis( &p, 0 ) - LIMIT ) < TOL );
  pendulum_destroy( &p );
}

/* the reaction of a stop is transmitted through the chain: a double
 * pendulum whose first joint is locked horizontally by a two-sided
 * stop and whose second link falls onto its lower stop.  At rest the
 * penetration of each stop carries the gravity torque about its
 * joint: that of the second link about the second joint, and that of
 * both links about the first, which includes the moment the second
 * stop transmits. */
static void assert_reaction(void)
{
  pendulum_t p;
  double t, k, trq1, trq2, trq1_expected, trq2_expected;

  pendulum_create( &p, "double_pendulum.ztk", RK_FD_JOINT_LIMIT_K_DEFAULT );
  k = rkFDPrpJointLimitK( &p.fd.prp );
  for( t=0; t<3.0; t+=DT ) pendulum_step( &p );
  zAssert( transmitted reaction (the chain comes to rest on the stops), zVecElem(p.dis,0) < 0 && zVecElem(p.dis,1) < LIMIT && fabs( zVecElem(p.vel,0) ) < 1.0e-6 && fabs( zVecElem(p.vel,1) ) < 1.0e-6 );
  zAssert( transmitted reaction (the chain shows the penetrations), fabs( zVecElem(p.dis,0) - pendulum_state_dis( &p, 0 ) ) < TOL && fabs( zVecElem(p.dis,1) - pendulum_state_dis( &p, 1 ) ) < TOL );
  /* the gravity torques at the rest angles: about the second joint, the
     second link; about the first, both links, the second one through the
     moment the second stop transmits */
  trq2 = k * ( LIMIT - pendulum_state_dis( &p, 1 ) );
  trq2_expected = G * 0.1 * cos( zVecElem(p.dis,0) + zVecElem(p.dis,1) );
  trq1 = k * ( 0.0 - pendulum_state_dis( &p, 0 ) );
  trq1_expected = G * 0.1 * cos( zVecElem(p.dis,0) ) + G * ( 0.2 * cos( zVecElem(p.dis,0) ) + 0.1 * cos( zVecElem(p.dis,0) + zVecElem(p.dis,1) ) );
  zAssert( transmitted reaction (the second stop carries the second link), fabs( trq2 - trq2_expected ) < 1.0e-3 * trq2_expected );
  zAssert( transmitted reaction (the first stop carries both links), fabs( trq1 - trq1_expected ) < 1.0e-3 * trq1_expected );
  pendulum_destroy( &p );
}

int main(int argc, char *argv[])
{
  assert_stop_function();
  assert_limit_trq();
  assert_pendulum();
  assert_hard_stop();
  assert_reaction();
  return 0;
}
