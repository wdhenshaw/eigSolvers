// This file automatically generated from residual.bC with bpp.
#include "mpi.h"
#include "Overture.h"
#include "ParallelUtility.h"
#include "display.h"
#include "CompositeGridOperators.h"
#include "gridFunctionNorms.h" 

// Put this last
#include "Ogev.h"

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )


// ==========================================================================================
///  \brief Compute the relative residual in the eigenvalue equation:
///      || L_h v + lambda^2 v ||_max / | lambda |
/// 
// ==========================================================================================
Real Ogev::getEigenPairResidual( Real lambda, realCompositeGridFunction & v,
                                                                realCompositeGridFunction & res,  
                                                                CompositeGridOperators & operators, 
                                                                int component /* =0 */ )
{

  // const Real & c  = dbase.get<real>("c");
  // CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");
    

  // realCompositeGridFunction res(cg);   // ***** do this for now ... is there a work space we can use instead?

    CompositeGrid & cg = *v.getCompositeGrid();

    Index I1,I2,I3;
    Index D1,D2,D3;
    int i1,i2,i3;

  // Compute   lap = L v 
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
        OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
        OV_GET_SERIAL_ARRAY(Real,res[grid],resLocal);


    // Compute at active interior points **CHECK ME ***  avoid interp points
        int extra=-1; // ** check me ** fix for periodic or interp boundaries
        getIndex(cg[grid].gridIndexRange(),I1,I2,I3,extra);

        getIndex(cg[grid].dimension(),D1,D2,D3);  
        bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,D1,D2,D3);
        RealArray lapLocal(D1,D2,D3);

        resLocal(D1,D2,D3,component)=0.;

        ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
        if( ok )
        {
            operators[grid].derivative(MappedGridOperators::laplacianOperator,vLocal,lapLocal,I1,I2,I3,component);
      // res = c^2 Lap (v) + lambda^2 * v
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( maskLocal(i1,i2,i3)>0 )
                {
                    resLocal(i1,i2,i3,component) = lapLocal(i1,i2,i3) + lambda * vLocal(i1,i2,i3,component);
                }
                else
                {
                    resLocal(i1,i2,i3,component)=0.;
                }
            }
        }
    }
    Real maxRes;

    maxRes = maxNorm( res,component )/fabs(lambda);


    return maxRes;


}

