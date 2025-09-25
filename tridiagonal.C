// Routines for solving tridiagonal linear systems
#include "A++.h"

typedef double Real;
typedef doubleSerialArray RealArray;

// define macros for array references
#define AX(i,j) Ax_p[(i)+3*(j)]
#define RHS(i) rhs_p[(i)]

// =====================================================================================
// Factor a tridiagonal matrix -- the factorization is stored in Ax
// =====================================================================================
int factorTridiagonalMatrix( RealArray & Ax )
{
    const int iax=Ax.getBase(1), ibx=Ax.getBound(1);
    Real *Ax_p = Ax.getDataPointer();

    // Factor: (no pivoting)
    for( int i1=iax+1; i1<=ibx; i1++ )
    {
        Real d = -AX(0,i1)/AX(1,i1-1);   // -a[i1]/b[i1-1]
        AX(1,i1) += d*AX(2,i1-1);
        AX(0,i1) = d; // save d here
    }

    return 0;
}

// =====================================================================================
// Solve the tridiagonal matrix problem given the factored matrix Ax
// =====================================================================================
int solveTridiagonal( RealArray & Ax, RealArray & rhs )
{
    const int iax=Ax.getBase(1), ibx=Ax.getBound(1);

    const Real *Ax_p = Ax.getDataPointer();
    Real *rhs_p = rhs.getDataPointer();

    // --- forward elimination ---
    for( int i1=iax+1; i1<=ibx; i1++ )
        RHS(i1) += AX(0,i1)*RHS(i1-1);

    // --- back-substitution ---
    RHS(ibx) = RHS(ibx)/AX(1,ibx);
    for( int i1=ibx-1; i1>=iax; i1-- )
        RHS(i1) = (RHS(i1) - AX(2,i1)*RHS(i1+1) )/AX(1,i1);

    return 0;
}
