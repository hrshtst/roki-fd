#include <roki_fd/roki_fd.h>

/* joint range limits as compliant stops: a pendulum falling onto the
 * lower limit of its joint, and a double pendulum whose first joint is
 * locked by a two-sided stop while its second link falls onto its
 * lower limit (the models of test/joint_limit_test.c, given shapes for
 * the display).  The joint displacements are recorded for rk_anim; the
 * pendulums swing in the x-z plane about the y axis:
 *   rk_anim ../model/pendulum.ztk pendulum.zvs -coord 0.05 -y 1.5 -pan 90
 *   rk_anim ../model/double_pendulum.ztk double_pendulum.zvs -coord 0.05 -y 1.5 -pan 90
 * The penetration of a stop is visible as the joint resting a little
 * beyond its limit; with a zero stiffness (the second argument, default
 * RK_FD_JOINT_LIMIT_K_DEFAULT) the stop is a perfectly inelastic
 * collision at the limit. */

#define T   3.0
#define DT  0.001
#define DTC 0.01

static void record(const char *model, const char *zvsfile, double k)
{
  rkFD fd;
  rkFDCell *cell;
  zVec dis;
  FILE *fp;
  double t, t_cnt;

  if( !( fp = fopen( zvsfile, "w" ) ) ) return;
  rkFDCreate( &fd );
  cell = rkFDChainRegFile( &fd, (char *)model );
  /* the shapes are for the display: no contact between the links */
  rkCDPairChainUnreg( rkFDCDBase(&fd.cd), rkFDCellChain(cell) );
  dis = zVecAlloc( rkChainJointSize( rkFDCellChain(cell) ) );
  rkFDChainSetDis( cell, dis );
  rkFDChainSetVel( cell, dis );
  rkFDODE2Assign( &fd, Regular );
  rkFDODE2AssignRegular( &fd, RK4 );
  rkFDPrpSetDT( &fd, DT );
  rkFDPrpSetJointLimitK( &fd, k );
  rkFDSetSolver( &fd, Vert );
  rkFDUpdateInit( &fd );
  for( t=0, t_cnt=0; t<T; t+=DT ){
    if( t >= t_cnt ){
      rkChainGetJointDisAll( rkFDCellChain(cell), dis );
      fprintf( fp, "%g ", DTC );
      zVecFPrint( fp, dis );
      t_cnt += DTC;
    }
    rkFDUpdate( &fd );
  }
  rkChainGetJointDisAll( rkFDCellChain(cell), dis );
  eprintf( "%s: final joint displacements [deg]:", model );
  for( t=0; t<zVecSizeNC(dis); t++ ) eprintf( " %.3f", zRad2Deg( zVecElem( dis, (int)t ) ) );
  eprintf( "\n" );
  rkFDUpdateDestroy( &fd );
  zVecFree( dis );
  rkFDDestroy( &fd );
  fclose( fp );
}

int main(int argc, char *argv[])
{
  double k = argc > 1 ? atof( argv[1] ) : RK_FD_JOINT_LIMIT_K_DEFAULT;

  record( "../model/pendulum.ztk", "pendulum.zvs", k );
  record( "../model/double_pendulum.ztk", "double_pendulum.zvs", k );
  return 0;
}
