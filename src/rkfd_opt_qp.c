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
