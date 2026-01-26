// This file automatically generated from coarseToFine.bC with bpp.
#include "mpi.h"
#include "Overture.h"
#include "ParallelUtility.h"
#include "display.h"
#include "CompositeGridOperators.h"
#include "gridFunctionNorms.h"
#include "GridStatistics.h"

#include "Integrate.h"

#include "Ogev.h"

// Boundary conditions:
const int periodic=-1, interpolation=0, displacement=1, traction=2, dirichlet=1, neumann=2;

#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(int i3=I3Base; i3<=I3Bound; i3++) for(int i2=I2Base; i2<=I2Bound; i2++) for(int i1=I1Base; i1<=I1Bound; i1++)

// --------------------------------------------------------------------------------------
//   Macro: return the index's for possible active points
//            
//  NOTE: This macro appears in solveSLEPc.bC and eigenModes.bC 
// --------------------------------------------------------------------------------------
    

// ==========================================================================================
///  \brief Evaluate the Rayleigh Quotient for a component of the grid function v
// ==========================================================================================
Real Ogev::getRayleighQuotient( realCompositeGridFunction & v, int component, CompositeGrid & cg, CompositeGridOperators & cgop )
{
    
  // printF("Ogev::getRayleighQuotient:: estimate the eigenvalue using the Rayleigh quotient...\n");

  // CompositeGrid & cg = *v.getCompositeGrid();

    Real lambda=0.;

  // CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");
    CompositeGridOperators & operators = cgop;
    
  // useDiscreteInnerProduct=true  : use the Integrate class to compute integrals
  //                        =false : sum over active interior points 

  // bool useDiscreteInnerProduct = false; // true;  // These options need to be compared
    const int & useAccurateInnerProduct = dbase.get<int>("useAccurateInnerProduct");
    if( !useAccurateInnerProduct )
    {
    // ----- Compute the RQ by summing over active points ------
    // This should be ALMOST AS GOOD as using an inner product when V_j is close to an eigenvector
    //     L_h V_j = \lambda^2 V_j 


        realCompositeGridFunction lap(cg);   // ***** do this for now ... is there a work space we can use instead?

        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
        int i1,i2,i3;
        int iab[2];    

    // Compute   lap = L v 
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            const IntegerArray & gid = mg.gridIndexRange();

            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);

                Iv[2]=Range(0,0);
                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                {
                    for( int side=0; side<=1; side++ )
                    {
                        int is = 1-2*side;
                        iab[side]=gid(side,axis);
                        const int bc = mg.boundaryCondition(side,axis);
                        if( bc==dirichlet )
                        {
                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                        }
                        else if( bc==neumann )
                        {
              // include boundary
                        }
                        else if( bc>0 )
                        {
                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                            OV_ABORT("error");
                        }
                        else if( bc<0 )
                        {
              // periodic -- include left end
                            if( side==1 )
                                iab[side] += is; 
                        }
                        else
                        {
              // interpolation boundary : include end 
                        }
                    }
                    Iv[axis] = Range(iab[0],iab[1]);
                }
          
            operators[grid].derivative(MappedGridOperators::laplacianOperator,vLocal,lapLocal,I1,I2,I3,component);

        }

    // Compute   ( v, vLap ) and (v, v )
        Real vDotLap=0., vDotv=0.;
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            const IntegerArray & gid = mg.gridIndexRange();

            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);  

                Iv[2]=Range(0,0);
                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                {
                    for( int side=0; side<=1; side++ )
                    {
                        int is = 1-2*side;
                        iab[side]=gid(side,axis);
                        const int bc = mg.boundaryCondition(side,axis);
                        if( bc==dirichlet )
                        {
                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                        }
                        else if( bc==neumann )
                        {
              // include boundary
                        }
                        else if( bc>0 )
                        {
                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                            OV_ABORT("error");
                        }
                        else if( bc<0 )
                        {
              // periodic -- include left end
                            if( side==1 )
                                iab[side] += is; 
                        }
                        else
                        {
              // interpolation boundary : include end 
                        }
                    }
                    Iv[axis] = Range(iab[0],iab[1]);
                }
            bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
            
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( maskLocal(i1,i2,i3) > 0 )
                {    
                    vDotLap += vLocal(i1,i2,i3,component)*lapLocal(i1,i2,i3);
                    vDotv   += vLocal(i1,i2,i3,component)*vLocal(i1,i2,i3,component);
                }
            }

        }
        ParallelUtility::getSum(vDotLap);
        ParallelUtility::getSum(vDotv);

        lambda = vDotLap/vDotv;  // this is really -lambda^2 .. sqrt taken below 

    }
    else
    {
    // Compute the RQ by using discrete approximations to the L2 inner product

        if( !dbase.has_key("integrate") )
        {
            Integrate & integrate = dbase.put<Integrate>("integrate");
            integrate.updateToMatchGrid(cg);
        }
        Integrate & integrate = dbase.get<Integrate>("integrate");

        realCompositeGridFunction lap(cg);   // ***** do this for now ... is there a work space we can use instead?

        Index I1,I2,I3;

    // Compute   lap = L v 
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);

            getIndex(cg[grid].dimension(),I1,I2,I3);
            bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);

            operators[grid].derivative(MappedGridOperators::laplacianOperator,vLocal,lapLocal,I1,I2,I3,component);

        }

    // Compute ( v, vLap )
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            getIndex(mg.dimension(),I1,I2,I3);
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);      
            bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
            if( ok )
                lapLocal(I1,I2,I3) = vLocal(I1,I2,I3,component)*lapLocal(I1,I2,I3);
        }

    // *** vLocal should have zero BC's -- so no need to change lap I think ****
        lap.periodicUpdate(); // is this needed ? interpolate may do this

        lap.interpolate();   // Probably do NOT want to do this if v is not smooth


        Real vDotLap = integrate.volumeIntegral(lap);    

    // Compute ( v,v )
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            getIndex(mg.dimension(),I1,I2,I3);
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);  

            bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
            if( ok )
                lapLocal(I1,I2,I3) = vLocal(I1,I2,I3,component)*vLocal(I1,I2,I3,component); // store in lapLocal
        }

        Real vDotv = integrate.volumeIntegral(lap);  

        lambda = vDotLap/vDotv;  // this is really -lambda^2 .. sqrt taken below 

    }

    lambda = -lambda; // to match Ogev
  // if( lambda<0 )
  // {
  //   lambda=sqrt(-lambda);
  //   if( false )
  //     printF("Ogev::getRayleighQuotient: component=%d: lambda(RQ)=%16.8e (Lap(phi) = -lambda^2 phi)\n",
  //          component,lambda);
  // }
  // else
  // {
  //   printF("Ogev::getRayleighQuotient::ERROR: rayleigh quotient for Laplacian is postive!? lambda=%12.4e\n",lambda);
  //   OV_ABORT("error");
  // }

    return lambda;
}


// ================================================================================================
/// \brief Convert coarse grid eigenpairs to fine grid eigenpairs
///
// ================================================================================================
int Ogev::coarseToFine( const aString & eigCase, const int numberOfComponents,
                                                                  int orderOfAccuracy, int & numEigenValues, int & numEigenVectors, 
                                                                  RealArray & eig, realCompositeGridFunction & uev, 
                                                                  IntegerArray & eigMultiplicity, IntegerArray & eigStartIndex, 
                                                                  CompositeGrid & cg, CompositeGridOperators & cgop  )
{

    const int & interpWidth = dbase.get<int>("interpWidth");
    int interpolationWidth = interpWidth>0 ? interpWidth : orderOfAccuracy+1;
    const int & useAccurateInnerProduct = dbase.get<int>("useAccurateInnerProduct");

    printF("\n @@@@@@@ ENTERING coarseToFine: interpolationWidth=%d useAccurateInnerProduct=%d @@@@@@@@@@ \n",interpolationWidth,useAccurateInnerProduct);


    if( !dbase.has_key("integrate") )
    {
        dbase.put<Integrate>("integrate");
    }
    Integrate & integrate = dbase.get<Integrate>("integrate");

    printF("==== Ogev::coarseToFine: compute integration weights...\n");
    integrate.updateToMatchGrid(cg);

    if( true )
    {
        Real volume;
        volume = integrate.volume();
        printf("... computed volume of the domain = %9.2e\n",volume);
    }

    Real lx=1., ly=1., lz=1.;
    int numEigsTrue = numEigenValues + 20; // compute a few more in case we have multiplicityTrue > 1
    RealArray eigsTrue(numEigsTrue); 

    int discreteEigenvalues=1; // get true discrete eigs if known
    bool eigenValuesAreKnown, eigenVectorsAreKnown; 

    getTrueEigenValues( eigCase, cg, numEigsTrue, eigsTrue, eigenValuesAreKnown, eigenVectorsAreKnown, NULL, discreteEigenvalues );

    RealArray lambdaRQ(numEigenVectors); // holds Rayleigh Quoitient


    for( int ie=0; ie<numEigenVectors; ie++ )
    {  
        lambdaRQ(ie) = getRayleighQuotient( uev, ie, cg, cgop );
        int je=ie; 
        for( int j=1; j<numEigsTrue; j++ )
        {
            if( fabs( lambdaRQ(ie) - eigsTrue(j) ) < fabs( lambdaRQ(ie) - eigsTrue(je) ) )
                je=j;
        }
        Real lambdaTrue = eigsTrue(je);
        Real relErr    =fabs(   eig(0,ie)-lambdaTrue)/lambdaTrue;
        Real relErrNew =fabs(lambdaRQ(ie)-lambdaTrue)/lambdaTrue;
        printF("ie=%3d: eig=%12.5e (rel-err=%9.2e), lambdaRQ=%12.5e (rel-err=%9.2e), err-ratio=%8.2e, true=%12.5e\n",
                  ie,eig(0,ie),relErr,lambdaRQ(ie),relErrNew, relErr/relErrNew,lambdaTrue );
    }

  // ----- Use the RQ value ----
    for( int ie=0; ie<numEigenVectors; ie++ )
    {
        eig(0,ie) = lambdaRQ(ie);
    }


    return 0;
}
