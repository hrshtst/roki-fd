#include <roki_fd/roki_fd.h>

/* A pendulum (a rod hinged at its end, pointing horizontally at rest)
 * falls under gravity against the lower range limit of its joint.
 * roki clamps the joint displacement at the limit when the state is
 * composed; the forward dynamics must stop the velocity there as well,
 * must not alter the motion before the limit is reached, and must not
 * hold the joint at the limit when it moves away from it. */

#define DT       0.001
#define T_END    1.0
#define LIMIT    zDeg2Rad(-30.0)
#define TOL      1.0e-9

typedef struct{
  rkFD fd;
  rkFDCell *cell;
  zVec dis, vel;
} pendulum_t;

static bool pendulum_create(pendulum_t *p, const char *filename)
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

static void pendulum_step(pendulum_t *p)
{
  rkFDUpdate( &p->fd );
  rkChainGetJointDisAll( rkFDCellChain(p->cell), p->dis );
  rkChainGetJointVelAll( rkFDCellChain(p->cell), p->vel );
}

/* the function itself */
static void assert_stop_function(void)
{
  pendulum_t p;
  zVec dis, vel;
  int n;

  pendulum_create( &p, "pendulum.ztk" );
  dis = zVecAlloc( 1 );
  vel = zVecAlloc( 1 );
  /* inside the range: untouched */
  zVecSetElem( dis, 0, 0.0 );
  zVecSetElem( vel, 0, -1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel );
  zAssert( rkFDChainStopJointsAtLimits (inside the range), n == 0 && zVecElem(vel,0) == -1.0 );
  /* at the lower limit, moving beyond it: stopped */
  zVecSetElem( dis, 0, LIMIT );
  zVecSetElem( vel, 0, -1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel );
  zAssert( rkFDChainStopJointsAtLimits (at the lower limit, moving beyond), n == 1 && zVecElem(vel,0) == 0.0 );
  /* at the lower limit, moving away from it: untouched */
  zVecSetElem( vel, 0, 1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel );
  zAssert( rkFDChainStopJointsAtLimits (at the lower limit, moving away), n == 0 && zVecElem(vel,0) == 1.0 );
  /* at the upper limit, moving beyond it: stopped */
  zVecSetElem( dis, 0, zDeg2Rad(90.0) );
  zVecSetElem( vel, 0, 1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel );
  zAssert( rkFDChainStopJointsAtLimits (at the upper limit, moving beyond), n == 1 && zVecElem(vel,0) == 0.0 );
  /* beyond the lower limit, moving beyond it: projected onto the limit and stopped */
  zVecSetElem( dis, 0, LIMIT - 0.1 );
  zVecSetElem( vel, 0, -1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel );
  zAssert( rkFDChainStopJointsAtLimits (beyond the lower limit, moving beyond), n == 1 && zVecElem(dis,0) == LIMIT && zVecElem(vel,0) == 0.0 );
  /* beyond the lower limit, moving away from it: projected onto the limit, velocity kept */
  zVecSetElem( dis, 0, LIMIT - 0.1 );
  zVecSetElem( vel, 0, 1.0 );
  n = rkFDChainStopJointsAtLimits( rkFDCellChain(p.cell), dis, vel );
  zAssert( rkFDChainStopJointsAtLimits (beyond the lower limit, moving away), n == 0 && zVecElem(dis,0) == LIMIT && zVecElem(vel,0) == 1.0 );
  zVecFree( dis );
  zVecFree( vel );
  pendulum_destroy( &p );
}

/* the pendulum */
static void assert_pendulum(void)
{
  pendulum_t p, q;
  double t, t_hit = -1, max_dev = 0, max_vel_after = 0, d;
  int i, n_hit = 0;

  pendulum_create( &p, "pendulum.ztk" );
  pendulum_create( &q, "pendulum_free.ztk" );
  for( t=0, i=0; t<T_END; t+=DT, i++ ){
    pendulum_step( &p );
    pendulum_step( &q );
    if( t_hit < 0 ){
      /* until the limit is reached, the two motions coincide */
      d = fabs( zVecElem(p.dis,0) - zVecElem(q.dis,0) );
      if( zVecElem(p.dis,0) <= LIMIT + TOL ){
        t_hit = t;
      } else if( d > max_dev ){
        max_dev = d;
      }
    } else {
      /* on the stop: pinned at the limit, and no velocity accumulates */
      if( zVecElem(p.dis,0) > LIMIT + TOL || zVecElem(p.dis,0) < LIMIT - TOL ) n_hit++;
      if( fabs( zVecElem(p.vel,0) ) > max_vel_after ) max_vel_after = fabs( zVecElem(p.vel,0) );
    }
  }
  zAssert( joint limit (the limit is reached), t_hit > 0 && t_hit < 0.5 );
  zAssert( joint limit (the free motion is unaltered before the limit), max_dev < TOL );
  zAssert( joint limit (the joint stays at the limit), n_hit == 0 );
  zAssert( joint limit (no velocity accumulates at the stop), max_vel_after < TOL );
  /* the unconstrained pendulum keeps swinging far beyond the limit */
  zAssert( joint limit (the free pendulum swings beyond the limit), zVecElem(q.dis,0) < LIMIT - 0.1 || fabs(zVecElem(q.vel,0)) > 0.1 );

  /* release: a velocity away from the stop moves the joint off the limit */
  zVecSetElem( p.vel, 0, 2.0 );
  rkFDChainSetVel( p.cell, p.vel );
  for( i=0; i<20; i++ ) pendulum_step( &p );
  zAssert( joint limit (the joint leaves the limit when it moves away), zVecElem(p.dis,0) > LIMIT + 0.01 );

  pendulum_destroy( &p );
  pendulum_destroy( &q );
}

int main(int argc, char *argv[])
{
  assert_stop_function();
  assert_pendulum();
  return 0;
}
