/* RoKi-FD - Robot Kinetics library: forward dynamics extention
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * rkfd_opt_qp - optimization tools: quadratic programming.
 * contributer: 2014- Naoki Wakisaka
 */

#include <roki_fd/rkfd_opt_qp.h>
#include <zm/zm_opt.h>

typedef struct{ /* list of active set indices */
  zIndex idx;
  double min;
} _rkFDQPASMIndexData;
ZEDA_DEF_LIST_CLASS( _rkFDQPASMIndexList, _rkFDQPASMIndex, _rkFDQPASMIndexData );

static zVec _rkFDQPSolveASMInitDefault(zMat a, zVec b, zVec ans, void *util)
{
  return zVecZero( ans );
}

static double _rkFDQPSolveASMConditionDefault(zMat a, zVec ans, int i, void *util)
{
  return zRawVecInnerProd(zMatRowBuf(a,i),zVecBuf(ans),zVecSizeNC(ans));
}

static uint _rkFDQPSolveASMInitIndex(zIndex idx, zMat a, zVec b, zVec ans, void *util, double cond(zMat,zVec,int,void*))
{
  int m, i;

  /* initialize the active set of constraints */
  for( m=0, i=0; i<zVecSizeNC(b); i++ ){
    if( zEqual( cond( a, ans, i, util ), zVecElemNC(b,i), zTOL ) ){
      zIndexElemNC(idx,i) = 1;
      m++;
    } else
      zIndexElemNC(idx,i) = 0;
  }
  return m;
}

#define RKFD_OPT_QP_ASM_TOL ( 1.0e-8 )
bool rkFDQPSolveASM(zMat q, zVec c, zMat a, zVec b, zVec ans, zIndex idx, zVec init(zMat,zVec,zVec,void*), double cond(zMat,zVec,int,void*), void *util)
{
  zMat qinv, aact, tm, s;
  zVec xy, bact, rhs, lam, tmp, xv, d;
  double tempd, tempd2, objv;
  _rkFDQPASMIndexList ilist;
  _rkFDQPASMIndex *idata;
  bool endflag, ret = true;
  int tempi;
  int n, mb, m;
  int i, j, k;

  /* function settings */
  if( !init ) init = _rkFDQPSolveASMInitDefault;
  if( !cond ) cond = _rkFDQPSolveASMConditionDefault;

  /* initalize */
  init( a, b, ans, util );

  /* allocate workspace */
  n = zVecSizeNC(ans);
  mb = zVecSizeNC(b);
  zListInit( &ilist );
  d = zVecAlloc(n);
  qinv = zMatAllocSqr( n );
  aact = zMatAlloc( mb, n );
  tm = zMatAlloc( mb, n );
  s = zMatAllocSqr( mb );
  xy = zVecAlloc( n + mb );
  bact = zVecAlloc( mb );
  rhs = zVecAlloc( mb );
  lam = zVecAlloc( mb );
  tmp = zVecAlloc( n );
  xv = zVecAlloc( n );
  if( !d || !qinv || !aact || !tm || !s || !xy || !bact || !rhs || !lam || !tmp || !xv ){
    ret = false;
    goto RET;
  }
  /* the Hessian is positive definite by construction (a Gram matrix plus
     the relaxation regularizer), so each active-set subproblem is solved
     through its inverse and the small Schur complement of the active
     rows, instead of one big indefinite KKT matrix whose rank-revealing
     solution is numerically fragile under degenerate active sets */
  if( !zMatInv( q, qinv ) ){
    ret = false;
    goto RET;
  }

  /* solve */
  m = _rkFDQPSolveASMInitIndex( idx, a, b, ans, util, cond );
  do{
    zVecSetSize( xy, n + m );
    if( m > 0 ){
      zMatSetSize( aact, m, n );
      zVecSetSize( bact, m );
      for( k=0, i=0; i<m; i++ ){
        while( !zIndexElemNC(idx,k) && k < mb ) k++;
        for( j=0; j<n; j++ )
          zMatElemNC(aact,i,j) = zMatElemNC(a,k,j);
        zVecElemNC(bact,i) = zVecElemNC(b,k);
        k++;
      }
      zMatSetSize( tm, m, n );
      zMatSetSize( s, m, m );
      zVecSetSize( rhs, m );
      zVecSetSize( lam, m );
      /* S lam = b_act + A_act Q^-1 c,  x = Q^-1 ( A_act^T lam - c ) */
      zMulMatMat( aact, qinv, tm );
      zMulMatMatT( tm, aact, s );
      zMulMatVec( tm, c, rhs );
      zVecAddDRC( rhs, bact );
      zVecZero( lam );
      zLESolveMP( s, rhs, NULL, NULL, lam );
      zMulMatTVec( aact, lam, tmp );
      zVecSubDRC( tmp, c );
      zMulMatVec( qinv, tmp, xv );
      for( i=0; i<m; i++ )
        zVecElemNC(xy,n+i) = zVecElemNC(lam,i);
    } else{
      /* unconstrained minimum */
      zMulMatVec( qinv, c, xv );
      zVecRevDRC( xv );
    }
    for( i=0; i<n; i++ )
      zVecElemNC(xy,i) = zVecElemNC(xv,i);

    /* progress check: the KKT solve of a degenerate active set only
       reproduces the pinned solution to its own accuracy, so compare
       against the active-set tolerance, not zTOL, lest numerical noise
       be taken for a feasible direction */
    for( i=0; i<n; i++ )
      if( !zEqual( zVecElemNC(xy,i), zVecElemNC(ans,i), RKFD_OPT_QP_ASM_TOL ) ) goto STEP2;
    for( i=0; i<n; i++ )
      zVecElemNC(ans,i) = zVecElemNC(xy,i);
    for( i=0; i<m; i++ )
      if( zVecElemNC(xy,n+i) < 0 ) goto NEXT;
    goto RET; /* found the optimal solution */

   NEXT:
    tempd = zVecElemNC(xy,n);
    for( i=1; i<m; i++ )
      if( zVecElemNC(xy,n+i) < tempd )
        tempd = zVecElemNC(xy,n+i);
    tempi = 0;
    for( i=0; i<zVecSizeNC(b); i++ )
      if( zIndexElemNC(idx,i) != 0 ){
        if( fabs( zVecElemNC(xy,tempi+n) - tempd ) < RKFD_OPT_QP_ASM_TOL ){
          zIndexElemNC(idx,i) = 0;
          m--;
        }
        tempi++;
      }
    continue;

   STEP2:
    /* find a new feasible direction */
    for( i=0; i<n; i++ )
      zVecElemNC(d,i) = zVecElemNC(xy,i) - zVecElemNC(ans,i);
    tempd = 1.0;
    for( i=0; i<zVecSizeNC(b); i++ ){
      tempd2 = zRawVecInnerProd( zMatRowBuf(a,i), zVecBuf(d), n );
      if( zIndexElemNC(idx,i) == 0 && tempd2 < 0 ){
        tempd2 = ( zVecElemNC(b,i) - cond( a, ans, i, util ) ) / tempd2;
        if( tempd2 < tempd ) tempd = tempd2;
      }
    }
    if( tempd < 0 ) tempd = 0; /* a blocking constraint is already tight */
    for( i=0; i<n; i++ )
      zVecElemNC(ans,i) += tempd * zVecElemNC(d,i);
    for( i=0; i<zVecSizeNC(b); i++ )
      if( zIndexElemNC(idx,i) == 0 && zEqual( cond( a, ans, i, util ), zVecElemNC(b,i), zTOL ) ){
        zIndexElemNC(idx,i) = 1;
        m++;
      }
    /* check if circulation happens due to degeneracy */
    endflag = false;
    objv = zQuadraticValue( q, c, ans );
    zListForEach( &ilist, idata ){
      for( i=0; i<zVecSizeNC(b);i++ )
        if( zIndexElemNC(idx,i) != zIndexElemNC(idata->data.idx,i) ) goto CONTINUE;
      if( fabs( idata->data.min/objv - 1.0 ) > RKFD_OPT_QP_ASM_TOL ) goto CONTINUE;
      endflag = true;
      break;
     CONTINUE: ;
    }
    if( endflag == true ) goto RET;

    /* register to the list of bases */
    idata = zAlloc( _rkFDQPASMIndex, 1 );
    idata->data.idx = zIndexCreate( zVecSizeNC(b) );
    for( i=0; i<zVecSizeNC(b); i++ )
      zIndexElemNC(idata->data.idx,i) = zIndexElemNC(idx,i);
    idata->data.min = objv;
    zListInsertTail( &ilist ,idata );
  } while( 1 );

 RET:
  zListForEach( &ilist, idata )
    zIndexFree( idata->data.idx );
  zListDestroy( _rkFDQPASMIndex, &ilist );
  zMatFreeAtOnce( 4, qinv, aact, tm, s );
  zVecFreeAtOnce( 7, xy, bact, rhs, lam, tmp, xv, d );
  return ret;
}

/* solve a quadratic program with a positive-definite Hessian and
 * inequality constraints through its dual linear complementarity
 * problem by the projected Gauss=Seidel method.  Unlike the active-set
 * method above, no rank-revealing linear solve is involved, so
 * degenerate (redundant) constraint sets such as friction pyramids of
 * coplanar contact vertices are handled without numerical hazard.
 * The primal answer satisfies the constraints only at convergence, and
 * impacts and edge contacts of a sole may take far more than 10^4
 * sweeps; stopping earlier leaves pulling normal forces. */
#define RKFD_OPT_QP_PGS_ITER_MAX 200000
#define RKFD_OPT_QP_PGS_TOL ( 1.0e-12 )
bool rkFDQPSolvePGS(zMat q, zVec c, zMat a, zVec b, zVec ans, zIndex idx)
{
  zMat qinv, t, m;
  zVec pfree, q0, lam;
  double w, lnew, dmax, lmax;
  bool ret = true;
  int mb, n, i, j, iter;

  n = zVecSizeNC(ans);
  mb = zVecSizeNC(b);
  qinv = zMatAllocSqr( n );
  t = zMatAlloc( mb, n );
  m = zMatAllocSqr( mb );
  pfree = zVecAlloc( n );
  q0 = zVecAlloc( mb );
  lam = zVecAlloc( mb );
  if( !qinv || !t || !m || !pfree || !q0 || !lam ){
    ret = false;
    goto RET;
  }
  /* the Hessian is positive definite by construction: a Gram matrix
     plus the relaxation regularizer on the diagonal */
  if( !zMatInv( q, qinv ) ){
    zVecZero( ans );
    ret = false;
    goto RET;
  }
  /* dual LCP: 0 <= lam  complementary to  M lam + q0 >= 0, where
     M = A Q^-1 A^T and q0 = A pfree - b with pfree = -Q^-1 c;
     the primal optimum is ans = pfree + Q^-1 A^T lam */
  zMulMatVec( qinv, c, pfree );
  zVecRevDRC( pfree );
  zMulMatMat( a, qinv, t );
  zMulMatMatT( t, a, m );
  zMulMatVec( a, pfree, q0 );
  zVecSubDRC( q0, b );
  for( iter=0; iter<RKFD_OPT_QP_PGS_ITER_MAX; iter++ ){
    dmax = lmax = 0;
    for( i=0; i<mb; i++ ){
      if( zMatElemNC(m,i,i) < zTOL ) continue; /* degenerate row */
      w = zVecElemNC(q0,i);
      for( j=0; j<mb; j++ )
        w += zMatElemNC(m,i,j) * zVecElemNC(lam,j);
      lnew = zVecElemNC(lam,i) - w / zMatElemNC(m,i,i);
      if( lnew < 0 ) lnew = 0;
      if( fabs( lnew - zVecElemNC(lam,i) ) > dmax )
        dmax = fabs( lnew - zVecElemNC(lam,i) );
      zVecElemNC(lam,i) = lnew;
      if( lnew > lmax ) lmax = lnew;
    }
    if( dmax <= RKFD_OPT_QP_PGS_TOL * zMax( 1, lmax ) ) break;
  }
  /* primal solution: t^T = Q^-1 A^T since Q is symmetric */
  zMulMatTVec( t, lam, ans );
  zVecAddDRC( ans, pfree );
  for( i=0; i<mb; i++ )
    zIndexElemNC(idx,i) = zVecElemNC(lam,i) > 0 ? 1 : 0;

 RET:
  zMatFreeAtOnce( 3, qinv, t, m );
  zVecFreeAtOnce( 3, pfree, q0, lam );
  return ret;
}

/* solve s x = r for x in place of the first n components of r by the
 * Cholesky decomposition of the leading n x n block of a symmetric
 * positive-definite matrix s, which is overwritten by its factor. */
static bool _rkFDQPCholeskySolveDST(zMat s, zVec r, int n)
{
  int i, j, k;
  double d;

  for( j=0; j<n; j++ ){
    d = zMatElemNC(s,j,j);
    for( k=0; k<j; k++ )
      d -= zMatElemNC(s,j,k) * zMatElemNC(s,j,k);
    if( d <= 0 ) return false;
    zMatElemNC(s,j,j) = sqrt( d );
    for( i=j+1; i<n; i++ ){
      d = zMatElemNC(s,i,j);
      for( k=0; k<j; k++ )
        d -= zMatElemNC(s,i,k) * zMatElemNC(s,j,k);
      zMatElemNC(s,i,j) = d / zMatElemNC(s,j,j);
    }
  }
  for( i=0; i<n; i++ ){ /* forward substitution */
    d = zVecElemNC(r,i);
    for( k=0; k<i; k++ )
      d -= zMatElemNC(s,i,k) * zVecElemNC(r,k);
    zVecElemNC(r,i) = d / zMatElemNC(s,i,i);
  }
  for( i=n-1; i>=0; i-- ){ /* backward substitution */
    d = zVecElemNC(r,i);
    for( k=i+1; k<n; k++ )
      d -= zMatElemNC(s,k,i) * zVecElemNC(r,k);
    zVecElemNC(r,i) = d / zMatElemNC(s,i,i);
  }
  return true;
}

/* solve a quadratic program with a positive-definite Hessian and
 * inequality constraints by the dual active-set method of Goldfarb and
 * Idnani (Math. Program. 27:1-33, 1983).  Starting from the
 * unconstrained minimum, the method adds the most violated constraint
 * one at a time and solves the problem exactly on the active set,
 * dropping on the way any active constraint whose multiplier would turn
 * negative.  The cost increases at every step, so that the iteration
 * ends in a finite number of steps at the exact optimum; unlike the
 * projected Gauss=Seidel method above, an ill-conditioned problem costs
 * accuracy, not iterations.  A constraint linearly dependent on the
 * active ones, as a fourth facet of a friction pyramid is, never enters
 * the active set: a pure dual step drops an active one instead, so that
 * the active constraints stay linearly independent and their Gram
 * matrix positive definite.  It returns false, leaving ans and idx
 * untouched, if the constraints are infeasible, the Gram matrix breaks
 * down numerically, or the iteration does not end. */
#define RKFD_OPT_QP_GI_ITER_MAX 10000
#define RKFD_OPT_QP_GI_TOL      ( 1.0e-12 )
#define RKFD_OPT_QP_GI_DEP_TOL  ( 1.0e-10 )
bool rkFDQPSolveGI(zMat q, zVec c, zMat a, zVec b, zVec ans, zIndex idx)
{
  zMat qinv, e, s;
  zVec f, u, z, r, qn, rn;
  zIndex act;
  double val, fmax, worst, an, zn, rmax, t, t1, t2, up;
  bool ret = false;
  int n, mb, na = 0, i, j, k, l, p, k1, iter = 0;

  n = zVecSizeNC(ans);
  mb = zVecSizeNC(b);
  qinv = zMatAllocSqr( n );
  e = zMatAllocSqr( n ); /* Q^-1 a_k^T of the active constraints, row-wise */
  s = zMatAllocSqr( n );
  f = zVecAlloc( n );
  u = zVecAlloc( n ); /* multipliers of the active constraints */
  z = zVecAlloc( n );
  r = zVecAlloc( n );
  qn = zVecAlloc( n );
  rn = zVecAlloc( mb );
  act = zIndexCreate( n );
  if( !qinv || !e || !s || !f || !u || !z || !r || !qn || !rn || !act )
    goto RET;
  if( !zMatInv( q, qinv ) ) goto RET;
  /* the unconstrained minimum */
  zMulMatVec( qinv, c, f );
  zVecRevDRC( f );
  for( i=0; i<mb; i++ )
    zVecElemNC(rn,i) = sqrt( zRawVecInnerProd( zMatRowBufNC(a,i), zMatRowBufNC(a,i), n ) );
  while( 1 ){
    /* the most violated inactive constraint, normalized by its row */
    for( fmax=0, j=0; j<n; j++ )
      fmax = zMax( fmax, fabs( zVecElemNC(f,j) ) );
    for( p=-1, worst=0, i=0; i<mb; i++ ){
      if( zVecElemNC(rn,i) < zTOL ) continue; /* degenerate row */
      for( k=0; k<na; k++ )
        if( zIndexElemNC(act,k) == i ) break;
      if( k < na ) continue;
      val = ( zRawVecInnerProd( zMatRowBufNC(a,i), zVecBufNC(f), n ) - zVecElemNC(b,i) ) / zVecElemNC(rn,i);
      if( val < -RKFD_OPT_QP_GI_TOL * ( fmax + RKFD_OPT_QP_GI_TOL ) && val < worst ){
        worst = val;
        p = i;
      }
    }
    if( p < 0 ) break; /* optimal */
    for( j=0; j<n; j++ ) /* Q^-1 a_p^T, Q being symmetric */
      zVecElemNC(qn,j) = zRawVecInnerProd( zMatRowBufNC(qinv,j), zMatRowBufNC(a,p), n );
    an = zRawVecInnerProd( zMatRowBufNC(a,p), zVecBufNC(qn), n );
    for( up=0; ; ){
      if( ++iter > RKFD_OPT_QP_GI_ITER_MAX ) goto RET;
      /* r = ( A Q^-1 A^T )^-1 A Q^-1 a_p^T over the active rows A,
         z = Q^-1 a_p^T - Q^-1 A^T r: the primal step direction */
      for( k=0; k<na; k++ ){
        for( l=0; l<na; l++ )
          zMatElemNC(s,k,l) = zRawVecInnerProd( zMatRowBufNC(a,zIndexElemNC(act,k)), zMatRowBufNC(e,l), n );
        zVecElemNC(r,k) = zRawVecInnerProd( zMatRowBufNC(a,zIndexElemNC(act,k)), zVecBufNC(qn), n );
      }
      if( na > 0 && !_rkFDQPCholeskySolveDST( s, r, na ) ) goto RET;
      for( j=0; j<n; j++ ){
        for( val=zVecElemNC(qn,j), k=0; k<na; k++ )
          val -= zVecElemNC(r,k) * zMatElemNC(e,k,j);
        zVecElemNC(z,j) = val;
      }
      /* dual step length: the first active multiplier to vanish */
      for( rmax=0, k=0; k<na; k++ )
        rmax = zMax( rmax, fabs( zVecElemNC(r,k) ) );
      for( t1=HUGE_VAL, k1=-1, k=0; k<na; k++ )
        if( zVecElemNC(r,k) > RKFD_OPT_QP_GI_TOL * rmax &&
            zVecElemNC(u,k) / zVecElemNC(r,k) < t1 ){
          t1 = zVecElemNC(u,k) / zVecElemNC(r,k);
          k1 = k;
        }
      /* primal step length: the added constraint satisfied; none if
         it depends linearly on the active ones, as it always does on
         n active ones, which round-off may hide in an ill-conditioned
         problem */
      zn = zRawVecInnerProd( zVecBufNC(z), zMatRowBufNC(a,p), n );
      t2 = na < n && zn > RKFD_OPT_QP_GI_DEP_TOL * an ?
        ( zVecElemNC(b,p) - zRawVecInnerProd( zMatRowBufNC(a,p), zVecBufNC(f), n ) ) / zn : HUGE_VAL;
      if( t1 == HUGE_VAL && t2 == HUGE_VAL ) goto RET; /* infeasible */
      t = zMin( t1, t2 );
      if( t2 < HUGE_VAL )
        for( j=0; j<n; j++ )
          zVecElemNC(f,j) += t * zVecElemNC(z,j);
      for( k=0; k<na; k++ )
        zVecElemNC(u,k) -= t * zVecElemNC(r,k);
      up += t;
      if( t2 <= t1 ){ /* full step: add the constraint */
        zIndexElemNC(act,na) = p;
        zVecElemNC(u,na) = up;
        zRawVecCopy( zVecBufNC(qn), zMatRowBufNC(e,na), n );
        na++;
        break;
      }
      /* partial step: drop the active constraint whose multiplier vanished */
      for( k=k1; k<na-1; k++ ){
        zIndexElemNC(act,k) = zIndexElemNC(act,k+1);
        zVecElemNC(u,k) = zVecElemNC(u,k+1);
        zRawVecCopy( zMatRowBufNC(e,k+1), zMatRowBufNC(e,k), n );
      }
      na--;
    }
  }
  zVecCopy( f, ans );
  for( i=0; i<mb; i++ )
    zIndexElemNC(idx,i) = 0;
  for( k=0; k<na; k++ )
    zIndexElemNC(idx,zIndexElemNC(act,k)) = zVecElemNC(u,k) > 0 ? 1 : 0;
  ret = true;

 RET:
  zMatFreeAtOnce( 3, qinv, e, s );
  zVecFreeAtOnce( 6, f, u, z, r, qn, rn );
  zIndexFree( act );
  return ret;
}
