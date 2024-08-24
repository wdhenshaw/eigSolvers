// This file automatically generated from fillInterpolationCoefficients.bC with bpp.
#include "mpi.h"
#include "Overture.h"
#include "ParallelUtility.h"

// Put this last
#include "Ogev.h"



#define initExplicitInterp EXTERN_C_NAME(initexplicitinterp)
extern "C"
{
    void initExplicitInterp(const int&ndc1, const int&ndc2, const int&ndc3, const int&ndci,
                    const int&ipar, Real&coeff, const Real&ci, Real&pr, Real&ps, Real&pt,
                    const Real&gridSpacing, const int&indexStart,
                    const int&variableInterpolationWidth, const int&interpoleeLocation, const int&interpoleeGrid);
}


// // do this for now:
// int numberOfComponents=1, numberOfProcessors=1;
// int *pnab, *pnoffset;

#define nab(side,axis,p,grid) pnab[(side)+2*( (axis) + 3*( (p) + numberOfProcessors*( (grid) ) ) )]
#define ndab(axis,p,grid) (nab(1,axis,p,grid)-nab(0,axis,p,grid)+1)
#define noffset(p,grid) pnoffset[(p)+numberOfProcessors*(grid)]


int Ogev::
getGlobalIndex( int n, int *iv, int grid, int p ) 
// ===============================================================================
/// \brief: Return the global index (equation number) given the point, grid and processor
/// \note: These equation numbers are base=0.
// ===============================================================================
{
  // printF("getGlobalIndex: n=%d, (i1,i2,i3)=(%d,%d,%d) grid=%d, p=%d\n",n,iv[0],iv[1],iv[2],grid,p);

    return  n + numberOfComponents*(
                  (iv[0]-nab(0,axis1,p,grid))+ndab(0,p,grid)*(
                    iv[1]-nab(0,axis2,p,grid) +ndab(1,p,grid)*(
                    iv[2]-nab(0,axis3,p,grid))) + noffset(p,grid) );
}


int Ogev::
getGlobalIndex( int n, int i1, int i2, int i3, int grid, int p ) 
// ===============================================================================
/// \brief: Return the global index (equation number) given the point, grid and processor
/// \note: These equation numbers are base=0.
// ===============================================================================
{
  // printF("getGlobalIndex: n=%d, (i1,i2,i3)=(%d,%d,%d) grid=%d, p=%d\n",n,i1,i2,i3,grid,p);

    return  n + numberOfComponents*(
                  (i1-nab(0,axis1,p,grid))+ndab(0,p,grid)*(
                    i2-nab(0,axis2,p,grid) +ndab(1,p,grid)*(
                    i3-nab(0,axis3,p,grid))) + noffset(p,grid) );
}


int Ogev::
buildGlobalIndexing(CompositeGrid & cg )
// ============================================================================================
/// \brief  Build the arrays needed for defining the global indexing (equation numbers) 
///
///   nab(side,axis,p,grid) : bounds of grid on processor p.
// ===========================================================================================
{
    debug = 7; // *****


    bool & globalIndexingIsComputed = dbase.get<bool>("globalIndexingIsComputed");

    printF("\n ++++ Ogev::buildGlobalIndexing : globalIndexingIsComputed=%d\n",(int)globalIndexingIsComputed);
  
    if( globalIndexingIsComputed )
        return 0;
    
    globalIndexingIsComputed=true;  // *wdh* What happens if the number of grids or grid points changes??

    const int & processorIsActive = dbase.get<bool>("processorIsActive");

    if( !processorIsActive )
        return 0;

    const int myid=max(0,Communication_Manager::My_Process_Number);

    int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2]; 
    int jv[3], &j1=jv[0], &j2=jv[1], &j3=jv[2]; 

    const int nabDimension=2*3*numberOfProcessors*cg.numberOfComponentGrids();
    delete [] pnab;  // *wdh* 091128 -- for refactor
    pnab = new int[nabDimension];  
    delete [] pnoffset; // *wdh* 091128 -- for refactor
    pnoffset = new int [numberOfProcessors*cg.numberOfComponentGrids()]; 


    numberOfGridPoints=0;
    numberOfGridPointsThisProcessor=0;

    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        intArray & ug= cg[grid].mask();
    // realArray & ug= uu[grid];
            
        int nd1a = ug.getBase(0), nd1b=ug.getBound(0);
        int nd2a = ug.getBase(1), nd2b=ug.getBound(1);
        int nd3a = ug.getBase(2), nd3b=ug.getBound(2);

        IndexBox uBox;
        for( int p=0; p<numberOfProcessors; p++ )
        {
            CopyArray::getLocalArrayBox( p, ug, uBox ); // find the array bounds on proc. p (no ghost)
            for( int axis=0; axis<3; axis++ )
            {
                nab(0,axis,p,grid)=uBox.base(axis);
                nab(1,axis,p,grid)=uBox.bound(axis);
            }
            if( debug & 1 )
                printF("Ogev::buildGlobalIndexing: grid=%i p=%i local bounds=[%i,%i][%i,%i][%i,%i] \n",grid,p,
                              nab(0,axis1,p,grid),nab(1,axis1,p,grid),
                              nab(0,axis2,p,grid),nab(1,axis2,p,grid),
                              nab(0,axis3,p,grid),nab(1,axis3,p,grid));
        }
        
        numberOfGridPointsThisProcessor+=ndab(0,myid,grid)*ndab(1,myid,grid)*ndab(2,myid,grid);
        numberOfGridPoints+= (nd1b-nd1a+1)*(nd2b-nd2a+1)*(nd3b-nd3a+1);   
    } // end for grid


    int offset=0;
    for( int p=0; p<numberOfProcessors; p++ )
    {
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            noffset(p,grid)=offset;
            offset+=ndab(0,p,grid)*ndab(1,p,grid)*ndab(2,p,grid);  // add number of pts on this local array
        }
    }
    
    if( debug & 4 )
    {
        printf(" Ogev::buildGlobalIndexing >>>>>>>>>>> myid=%i numberOfGridPointsThisProcessor=%i, numberOfGridPoints=%i<<<<<<<<<<<<<<<<<<<\n",
                      myid,numberOfGridPointsThisProcessor,numberOfGridPoints);
        fflush(0);
    }
    
    return 0;
}


int Ogev::
fillInterpolationCoefficients( Mat & A, realCompositeGridFunction & uu )
//===================================================================
// /Description:
//   Fill the matrix with the interpolation coefficients
// (this routine comes from PETScSolver::fillfillInterpolationCoefficients)
//===================================================================
{
    
    CompositeGrid & cg = *uu.getCompositeGrid();
    const int myid=Communication_Manager::My_Process_Number;

    if( debug & 4 )
    {
        printF("fillInterpolationCoefficients: max(cg.numberOfInterpolationPoints)=%i\n",max(cg.numberOfInterpolationPoints));
    }
    

    if( cg.numberOfBaseGrids() ==0 || max(cg.numberOfInterpolationPoints) <= 0 )
        return 0;

    Real time0=getCPU();

    const int numberOfDimensions=cg.numberOfDimensions();

  // for now we use only one width per grid
    int axis,grid;
    IntegerArray width(3,cg.numberOfComponentGrids()); width=1;
    Range NG(0,cg.numberOfComponentGrids()-1);
    for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    for( axis=axis1; axis<numberOfDimensions; axis++ ) 
        width(axis,grid)=max(width(axis,grid),max(cg.interpolationWidth(axis,grid,NG)));

    const int maxWidth=max(width);

    IntegerArray indexStart(3,cg.numberOfComponentGrids());
    RealArray gridSpacing(3,cg.numberOfComponentGrids());
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        for( axis=0; axis<3; axis++ )
        {
            indexStart(axis,grid)=mg.gridIndexRange(0,axis);
            gridSpacing(axis,grid)=mg.gridSpacing(axis);
        }
            
    }
        
    int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2]; 
    int ivd[3], &i1d=ivd[0], &i2d=ivd[1], &i3d=ivd[2];  // for donor position

    realSerialArray coeff;
    
    int ierr;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        const realArray & ug = uu[grid];

        int ni=cg.numberOfInterpolationPoints(grid);
        if( ni==0 ) continue;
                
        intArray & interpoleeLocation = cg.interpoleeLocation[grid];
        intSerialArray il; 
        
    // *wdh* 091129 -- fix to use local interpolation arrays if they are there.
        if( ( grid<cg.numberOfBaseGrids() && 
                    cg->localInterpolationDataState==CompositeGridData::localInterpolationDataForAMR ) || 
                    cg->localInterpolationDataState==CompositeGridData::noLocalInterpolationData )
        {

      // use the interpolation data in the parallel arrays
            getLocalArrayWithGhostBoundaries(cg.interpoleeLocation[grid],il);
        }
        else
        {
      // use the interpolation data in the serial arrays (for now these are refinement grids)
            il.reference(cg->interpoleeLocationLocal[grid]);
        }

        int n1a = il.getBase(0) +interpoleeLocation.getGhostBoundaryWidth(0), 
                n1b = il.getBound(0)-interpoleeLocation.getGhostBoundaryWidth(0);

        ni = n1b-n1a+1;  // number of interpolation points on this processor
        
        if( ni==0 ) continue;
        

        Range R(n1a,n1b);

        coeff.redim(R,width(axis1,grid),width(axis2,grid),width(axis3,grid));

        intSerialArray ip; 
        intSerialArray interpoleeGrid; 
        intSerialArray viw; 
        realSerialArray ci; 

        if( ( grid<cg.numberOfBaseGrids() && 
                    cg->localInterpolationDataState==CompositeGridData::localInterpolationDataForAMR ) || 
                    cg->localInterpolationDataState==CompositeGridData::noLocalInterpolationData )
        {

      // use the interpolation data in the parallel arrays
            getLocalArrayWithGhostBoundaries(cg.interpolationPoint[grid],ip);
            getLocalArrayWithGhostBoundaries(cg.interpoleeGrid[grid],interpoleeGrid);
            getLocalArrayWithGhostBoundaries(cg.variableInterpolationWidth[grid],viw);
            getLocalArrayWithGhostBoundaries(cg.interpolationCoordinates[grid],ci);

        }
        else
        {
      // use the interpolation data in the serial arrays (for now these are refinement grids)
      // printf("PETScSolver::USE LOCAL INTERP ARRAYS\n");
            
            ip.reference( cg->interpolationPointLocal[grid]);
            il.reference( cg->interpoleeLocationLocal[grid]);
            interpoleeGrid.reference( cg->interpoleeGridLocal[grid]);
            viw.reference( cg->variableInterpolationWidthLocal[grid]);
            ci.reference(cg->interpolationCoordinatesLocal[grid]);
        }
            
    // int debug = 0; // Oges::debug;
        if( debug & 4 )
        {
            printf("Ogev::fillInterp: myid=%i: grid=%i [n1a,n1b]=[%i,%i] numberOfComponents=%i interp-width=%i\n",myid,grid,n1a,n1b,
                          numberOfComponents,width(axis1,grid));
      // for( int i=n1a; i<=n1b; i++ )
      // {
      //        printf(" grid=%i i=%i ip=(%i,%i) il=(%i,%i) viw=%i\n",grid,i,ip(i,0),ip(i,1),il(i,0),il(i,1),viw(i));
      // }
        }

        int ipar[7]={numberOfDimensions,
                                  grid,
                                  ni,
                                  mg.isCellCentered(0),
                                  (int)Interpolant::precomputeAllCoefficients,
                                  maxWidth,
                                  0}; //  This last position is saved for a return value of useVariableWidthInterpolation

        realSerialArray & cc = coeff;

    // ************ warning -- watch out for indexing of local arrays -- base and bound : n1a,..,n1b    
    //   use interpoleeGrid(n1a,0)
        RealArray pr(R),ps(R),pt(R);
        initExplicitInterp(cc.getLength(0),cc.getLength(1),cc.getLength(2),il.getLength(0),
                                              ipar[0], 
                                              *cc.getDataPointer(),
                                              ci(n1a,0),pr(n1a),ps(n1a),pt(n1a),gridSpacing(0,0),indexStart(0,0),
                                              viw(n1a),il(n1a,0),interpoleeGrid(n1a,0));
                
        int useVariableWidthInterpolation=ipar[6]; 

    // printf("Ogev::useVariableWidthInterpolation=%i\n",useVariableWidthInterpolation);

    // cc.display("coeff after initExplicitInterp");

        const Real epsForInterpCoeff=REAL_EPSILON*100.; // neglect interpolation coeff smaller than this

        assert( width(0,grid)==width(1,grid) );

        i3=0; i3d=0;
        int n=0;
        for( int i=n1a; i<=n1b; i++ )
        {
            i1=ip(i,0);
            i2=ip(i,1);
            if( numberOfDimensions==3 ) i3=ip(i,2);
            
            const int gridi = interpoleeGrid(i);
            const int iw = useVariableWidthInterpolation ? viw(i) : width(0,grid);  // *wdh* 100113 -- added support for VIW
            const int iw3 = numberOfDimensions==2 ? 1 : iw;

            #ifdef USE_PPP
                int p= ug.Array_Descriptor.findProcNum( iv );  // processor number
            #else
                int p=0; 
            #endif
            for( int n=0; n<numberOfComponents; n++ )
            {

        // fill in value -1 for grid,(i1,i2,i3)
                const int ig=getGlobalIndex( n, iv, grid, p );  // get the global index
                Real v=-1.; 
                ierr = MatSetValues(A,1,&ig,1,&ig,&v,INSERT_VALUES);CHKERRQ(ierr);

        // if( debug & 1 ) printf("interp: grid=%i p=%i i=%i ip=(%i,%i,%i) n=%i coeff[m1,m2]:\n",grid,p,i,i1,i2,i3,n);


                for( int m3=0; m3<iw3; m3++ )
                {
                    if( numberOfDimensions==3 ) i3d=il(i,2)+m3;
                    for( int m2=0; m2<iw; m2++ )
                    {
                        i2d=il(i,1)+m2;
                        for( int m1=0; m1<iw; m1++ )
                        {
                            i1d=il(i,0)+m1;
                
                            #ifdef USE_PPP        
                                int p= uu[gridi].Array_Descriptor.findProcNum( ivd );  // processor number
                            #else
                                int p=0;
                            #endif
                            int jg=getGlobalIndex( n, ivd, gridi, p );  // get the global index

              // fill in value  grid, (i1,i2,i3) coeff(i,m1,m2,m3)
                            v=coeff(i,m1,m2,m3);

                            if( fabs(v)>epsForInterpCoeff )
                            {
                // if( debug & 1 ) printf(" (ig,jg)=(%i,%i) [m1=%i,m2=%i,m3=%i](gridi=%i,p=%i)=%4.2f \n",ig,jg,m1,m2,m3,gridi,p,v);
                                ierr = MatSetValues(A,1,&ig,1,&jg,&v,INSERT_VALUES);CHKERRQ(ierr);

                            }
                        
                        }
                    }
                }
            }
        }

    }
    Real time=getCPU()-time0;
    time0=getCPU();
    if( debug & 1 ) printF("*** fillInterpolationCoefficients: cpu time = %8.2e\n",time);

        
    return 0;

}