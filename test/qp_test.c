#include <roki_fd/roki_fd.h>
#include <roki_fd/rkfd_opt_qp.h>

/* The contact QP solvers.
 *
 * rkFDQPSolveGI must return the exact optimum of
 *   min 1/2 f^T Q f + c^T f  subject to  A f >= b
 * for a positive-definite Q, which the Karush=Kuhn=Tucker conditions
 * certify: the answer is feasible, the gradient Q f + c is a
 * nonnegative combination of the rows of the constraints marked
 * active, and those constraints hold with equality.  Besides random
 * problems, it is tested on the contact QP of the Vert solver for a
 * rigid body resting on four coplanar vertices, whose Hessian is
 * singular but for the relaxation and whose friction pyramids are
 * redundant: the problem the projected Gauss=Seidel method converges
 * on only slowly. */

#define N_TRIAL   1000
#define TOL       1.0e-8
#define TOL_RIGID 1.0e-7 /* the Hessian's condition number reaches 1e9 */

#define NV        4    /* contact vertices */
#define NC        8    /* facets of a friction pyramid */
#define MU        0.5  /* friction coefficient */
#define RELAX     1.0e-4
#define MASS      1.0
#define GRAVITY   9.8
#define DT        0.001

/* the Karush=Kuhn=Tucker conditions of the answer f with the active
 * set idx; returns the worst violation relative to the scale of the
 * problem */
static double kkt_violation(zMat q, zVec c, zMat a, zVec b, zVec f, zIndex idx)
{
  zMat aa, gram;
  zVec g, ag, lam, res;
  int n, mb, na, i, j, k;
  double scale, viol = 0, d;

  n = zVecSizeNC(f);
  mb = zVecSizeNC(b);
  for( na=0, i=0; i<mb; i++ )
    if( zIndexElemNC(idx,i) ) na++;
  g = zVecAlloc( n );
  res = zVecAlloc( n );
  zMulMatVec( q, f, g );
  zVecAddDRC( g, c );
  scale = 1 + zVecNorm( c ) + zVecNorm( g );
  /* primal feasibility, and equality of the active constraints */
  for( i=0; i<mb; i++ ){
    d = ( zRawVecInnerProd( zMatRowBufNC(a,i), zVecBufNC(f), n ) - zVecElemNC(b,i) ) / scale;
    viol = zMax( viol, -d );
    if( zIndexElemNC(idx,i) ) viol = zMax( viol, fabs( d ) );
  }
  zVecCopy( g, res );
  if( na > 0 ){
    /* multipliers: the least-squares fit of the gradient by the
       active rows, which are linearly independent */
    aa = zMatAlloc( na, n );
    gram = zMatAllocSqr( na );
    ag = zVecAlloc( na );
    lam = zVecAlloc( na );
    for( k=0, i=0; i<mb; i++ )
      if( zIndexElemNC(idx,i) ){
        for( j=0; j<n; j++ )
          zMatElemNC(aa,k,j) = zMatElemNC(a,i,j);
        k++;
      }
    zMulMatMatT( aa, aa, gram );
    zMulMatVec( aa, g, ag );
    zLESolveGauss( gram, ag, lam );
    for( k=0; k<na; k++ )
      viol = zMax( viol, -zVecElemNC(lam,k) / scale ); /* dual feasibility */
    for( k=0; k<na; k++ )
      for( j=0; j<n; j++ )
        zVecElemNC(res,j) -= zVecElemNC(lam,k) * zMatElemNC(aa,k,j);
    zMatFreeAtOnce( 2, aa, gram );
    zVecFreeAtOnce( 2, ag, lam );
  }
  viol = zMax( viol, zVecNorm( res ) / scale ); /* stationarity */
  zVecFreeAtOnce( 2, g, res );
  return viol;
}

/* a random problem with a known feasible point */
static void assert_random(void)
{
  zMat q, gen, a;
  zVec c, b, x0, f;
  zIndex idx;
  int trial, n, mb, i, j, nfail = 0;
  double viol, viol_max = 0;

  for( trial=0; trial<N_TRIAL; trial++ ){
    n = zRandI( 2, 12 );
    mb = zRandI( 1, 3*n );
    q = zMatAllocSqr( n );
    gen = zMatAllocSqr( n );
    a = zMatAlloc( mb, n );
    c = zVecAlloc( n );
    b = zVecAlloc( mb );
    x0 = zVecAlloc( n );
    f = zVecAlloc( n );
    idx = zIndexCreate( mb );
    zMatRandUniform( gen, -1, 1 );
    zMulMatTMat( gen, gen, q );
    for( i=0; i<n; i++ ) zMatElemNC(q,i,i) += 0.1;
    zVecRandUniform( c, -1, 1 );
    zMatRandUniform( a, -1, 1 );
    zVecRandUniform( x0, -1, 1 );
    for( i=0; i<mb; i++ ){
      for( zVecElemNC(b,i)=-zRandF(0,1), j=0; j<n; j++ )
        zVecElemNC(b,i) += zMatElemNC(a,i,j) * zVecElemNC(x0,j);
    }
    if( !rkFDQPSolveGI( q, c, a, b, f, idx ) ) nfail++;
    else if( ( viol = kkt_violation( q, c, a, b, f, idx ) ) > viol_max ) viol_max = viol;
    zMatFreeAtOnce( 3, q, gen, a );
    zVecFreeAtOnce( 4, c, b, x0, f );
    zIndexFree( idx );
  }
  zAssert( rkFDQPSolveGI (random problems solved), nfail == 0 );
  zAssert( rkFDQPSolveGI (random problems: KKT conditions), viol_max < TOL );
}

/* the Vert solver's QP of a rigid body on NV coplanar vertices, the
 * first nv of which touch: the relative velocity error of the vertices
 * after a step, w = W f + w0 with the Delassus operator W, is minimized
 * in the least-squares sense with the relaxation, in the friction
 * pyramids; each vertex force is (normal, tangential, tangential) */
static void rigid_body_qp(int nv, zVec3D *vel, zMat q, zVec c, zMat a, zVec b)
{
  zMat jac, w;
  zVec w0;
  zVec3D p[NV];
  double inertia[3] = { 0.002, 0.003, 0.004 };
  int i, j, k, n;

  n = 3 * nv;
  zVec3DCreate( &p[0],  0.05,  0.03, -0.02 );
  zVec3DCreate( &p[1], -0.05,  0.03, -0.02 );
  zVec3DCreate( &p[2], -0.05, -0.03, -0.02 );
  zVec3DCreate( &p[3],  0.05, -0.03, -0.02 );
  /* the velocity of each vertex along (z, x, y) from the body twist
     (v, omega): v + omega x p */
  jac = zMatAlloc( n, 6 );
  for( i=0; i<nv; i++ ){
    zMatElemNC(jac,3*i  ,2) = 1;
    zMatElemNC(jac,3*i  ,3) = p[i].c.y;
    zMatElemNC(jac,3*i  ,4) = -p[i].c.x;
    zMatElemNC(jac,3*i+1,0) = 1;
    zMatElemNC(jac,3*i+1,4) = p[i].c.z;
    zMatElemNC(jac,3*i+1,5) = -p[i].c.y;
    zMatElemNC(jac,3*i+2,1) = 1;
    zMatElemNC(jac,3*i+2,3) = -p[i].c.z;
    zMatElemNC(jac,3*i+2,5) = p[i].c.x;
  }
  /* W = J M^-1 J^T */
  w = zMatAllocSqr( n );
  for( i=0; i<n; i++ )
    for( j=0; j<n; j++ ){
      for( k=0; k<3; k++ )
        zMatElemNC(w,i,j) += zMatElemNC(jac,i,k) * zMatElemNC(jac,j,k) / MASS
                           + zMatElemNC(jac,i,k+3) * zMatElemNC(jac,j,k+3) / inertia[k];
    }
  /* the free relative velocity: gravity over a step plus a slip */
  w0 = zVecAlloc( n );
  for( i=0; i<nv; i++ ){
    zVecElemNC(w0,3*i  ) = -GRAVITY * DT + vel[i].c.z;
    zVecElemNC(w0,3*i+1) = vel[i].c.x;
    zVecElemNC(w0,3*i+2) = vel[i].c.y;
  }
  zMulMatTMat( w, w, q );
  for( i=0; i<n; i++ ) zMatElemNC(q,i,i) += RELAX;
  zMulMatTVec( w, w0, c );
  /* the friction pyramids: MU f_n + cos(th) f_t1 + sin(th) f_t2 >= 0 */
  zMatZero( a );
  zVecZero( b );
  for( i=0; i<nv; i++ )
    for( k=0; k<NC; k++ ){
      zMatElemNC(a,NC*i+k,3*i  ) = MU;
      zMatElemNC(a,NC*i+k,3*i+1) = cos( zPIx2 * k / NC );
      zMatElemNC(a,NC*i+k,3*i+2) = sin( zPIx2 * k / NC );
    }
  zMatFreeAtOnce( 2, jac, w );
  zVecFree( w0 );
}

static void assert_rigid_body(void)
{
  zMat q, a;
  zVec c, b, f;
  zIndex idx;
  zVec3D vel[NV];
  int trial, nv, i, nfail = 0;
  double viol, viol_max = 0, fn_min = 0, fn_sum, rest_err = 0;

  for( trial=0; trial<N_TRIAL; trial++ ){
    nv = trial % 2 == 0 ? NV : zRandI( 1, NV ); /* full sole, or edges and corners */
    q = zMatAllocSqr( 3*nv );
    a = zMatAlloc( NC*nv, 3*nv );
    c = zVecAlloc( 3*nv );
    b = zVecAlloc( NC*nv );
    f = zVecAlloc( 3*nv );
    idx = zIndexCreate( NC*nv );
    for( i=0; i<nv; i++ ) zVec3DZero( &vel[i] );
    if( trial % 4 >= 2 ){ /* a body landing and slipping */
      for( i=0; i<nv; i++ )
        zVec3DCreate( &vel[i], zRandF(-0.5,0.5), zRandF(-0.5,0.5), zRandF(-1,0) );
    }
    rigid_body_qp( nv, vel, q, c, a, b );
    if( !rkFDQPSolveGI( q, c, a, b, f, idx ) ){
      nfail++;
    } else{
      if( ( viol = kkt_violation( q, c, a, b, f, idx ) ) > viol_max ) viol_max = viol;
      for( fn_sum=0, i=0; i<nv; i++ ){
        fn_min = zMin( fn_min, zVecElemNC(f,3*i) );
        fn_sum += zVecElemNC(f,3*i);
      }
      /* resting on the full sole, the normal impulses carry the weight */
      if( nv == NV && trial % 4 == 0 )
        rest_err = zMax( rest_err, fabs( fn_sum / ( MASS * GRAVITY * DT ) - 1 ) );
    }
    zMatFreeAtOnce( 2, q, a );
    zVecFreeAtOnce( 3, c, b, f );
    zIndexFree( idx );
  }
  zAssert( rkFDQPSolveGI (rigid body contact problems solved), nfail == 0 );
  zAssert( rkFDQPSolveGI (rigid body contact: KKT conditions), viol_max < TOL_RIGID );
  zAssert( rkFDQPSolveGI (rigid body contact: no pulling normal force), fn_min > -TOL_RIGID );
  zAssert( rkFDQPSolveGI (rigid body at rest: the weight is carried), rest_err < 1.0e-3 );
}

/* an infeasible problem: f >= 1 and -f >= 0 */
static void assert_infeasible(void)
{
  zMat q, a;
  zVec c, b, f;
  zIndex idx;

  q = zMatAllocSqr( 1 );
  a = zMatAlloc( 2, 1 );
  c = zVecAlloc( 1 );
  b = zVecAlloc( 2 );
  f = zVecAlloc( 1 );
  idx = zIndexCreate( 2 );
  zMatElemNC(q,0,0) = 1;
  zMatElemNC(a,0,0) = 1;
  zMatElemNC(a,1,0) = -1;
  zVecElemNC(b,0) = 1;
  zAssert( rkFDQPSolveGI (an infeasible problem is reported), !rkFDQPSolveGI( q, c, a, b, f, idx ) );
  zMatFreeAtOnce( 2, q, a );
  zVecFreeAtOnce( 3, c, b, f );
  zIndexFree( idx );
}

int main(int argc, char *argv[])
{
  zRandInit();
  assert_random();
  assert_rigid_body();
  assert_infeasible();
  return 0;
}
