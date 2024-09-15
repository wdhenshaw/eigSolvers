// This file automatically generated from Ogev.bC with bpp.
//===========================================================================
//
//  Overlapping Grid EigenValue/EigenVector solver
//
//===========================================================================

#include "Ogev.h"
#include "display.h"
#include <numeric>      // std::iota

#define rjbesl EXTERN_C_NAME(rjbesl)
extern "C"
{
    void rjbesl( const Real & ka, const Real & alpha, const int & nb, Real & jnka, int & ncalc);
}

#define FOR_3D(i1,i2,i3,I1,I2,I3)                                       int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(); int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(int i3=I3Base; i3<=I3Bound; i3++)                                       for(int i2=I2Base; i2<=I2Bound; i2++)                                     for(int i1=I1Base; i1<=I1Bound; i1++)

// --------------------------------------------------------------------
// Macro to evaluate the spherical Bessel function near the origin
// Input: nn,x
// Output : R 
// -------------------------------------------------------------------

// --------------------------------------------------------------------
// Macro to evaluate the associated Legendre function
//  n = mPhi   : degree
//  k = mTheta : order 
// *wdh* Formulae from Wolfram 
// -------------------------------------------------------------------


    

//===========================================================================
//
/// \brief Constructor for the Overlapping Grid EigenValue/EigenVector solver
//
//===========================================================================
Ogev::Ogev()
{
    dbase.put<bool>("globalIndexingIsComputed") = false;
    dbase.put<bool>("processorIsActive")        = true;

    debug = 0;
    numberOfComponents=1;
    numberOfProcessors = max(1,Communication_Manager::Number_Of_Processors);

    numberOfGridPoints = 0 ;
    numberOfGridPointsThisProcessor = 0;

  // --- arrays for global indexing -----
    pnab = NULL;  
    pnoffset = NULL;

    dbase.put<int>("orderOfAccuracy")=2; 

    dbase.put<Real>("eigenValueTolForMultiplicity")=1e-5; // tolerance for checking whether eigenvalues are the same 
    dbase.put<int>("numberOfEigenvectors")=0;

}

//===========================================================================
//
/// \brief Destructor for the Overlapping Grid EigenValue/EigenVector solver
//
//===========================================================================
Ogev::~Ogev()
{
    delete [] pnab;
    delete [] pnoffset;
}

// ==============================================================================
/// \brief Return the name of the boundary condition. **put this here for now**
// ==============================================================================
aString Ogev::bcName( int bc )
{
    if( bc==-1 )
        return "p";
    else if( bc==0 )
        return "i"; // interp
    else if( bc==1 )
        return "d"; // displacement
    else if( bc==2 )
        return "t"; // traction
    else
        return "u"; // unknown

}


#include <vector>
#include <algorithm>       // STL algorithms class library


// ========================================================
// \brief Return the symbol an p-th order accurate approximation
// to minus the second derivative
// ========================================================
Real Ogev::getDiscreteSymbol( const Real modeNumber, const Real dx ) const
{
    const int & orderOfAccuracy= dbase.get<int>("orderOfAccuracy");

    Real symbol;
    if( orderOfAccuracy==2 )
    {

        symbol = SQR( sin(modeNumber*Pi*dx/2.)/(dx/2.) );
    }
    else if( orderOfAccuracy==4 )
    {
        Real sine = sin(modeNumber*Pi*dx/2.);
        symbol = SQR(sine/(dx/2.)) * ( 1. + sine*sine/3. );
    }
    else
    {
        OV_ABORT("finish me");
    }
    return symbol;
}

//===========================================================================
//
/// \brief Return the eigenvalues (and optionally the eigenvectors) of the minus Laplacian
///       - Delta( u ) = lambda^2 u 
///  for a line, square or box, with given boundatu conditions.
//
/// \param numEigs (input) : number of eigenvalues requested
/// \param eigs(i)   : eigenvalues from smallest to largest
/// \param numberOfDimensions (input) : 1,2 or 3
/// \param mg (input) : mg.boundayCondition(side,axis) defines boundary conditions, 1=dirichlet, 2=neumann
/// \param lx,ly,lz : box dimensions
/// \param discreteEigenvalues Input) : if true return the eigenvalues for the discretized problem
// 
//===========================================================================
int Ogev::
getEigenvaluesBox( int numEigs, RealArray & eigs, CompositeGrid & cg , 
                                      Real lx  /* =1.0 */ , Real ly  /* =1.0 */ , Real lz /* =1.0 */,
                                      RealCompositeGridFunction *eigenvector /* = NULL */,
                                      const bool discreteEigenvalues /* =false */  )
{
    const int & orderOfAccuracy= dbase.get<int>("orderOfAccuracy");

    MappedGrid & mg = cg[0]; // for BC's

    Real dx[3]={1.,1.,1.}, xab[2][3]={0.,0.,0.,0.,0.,0.};
    mg.getRectangularGridParameters( dx, xab );

    const IntegerArray & gid = mg.gridIndexRange();

    Real lv[3] = { lx,ly,lz };
    for( int axis=0; axis<3; axis++ )
    {
        lv[axis] = xab[1][axis] - xab[0][axis];
    }

  // mg.getDeltaX(dx); 

    if( discreteEigenvalues )
    {
        assert( mg.isRectangular() );
    } 

    const bool computeEigenvectors = eigenvector!=NULL; 

    const int dirichlet=1, neumann=2; // **FIX ME**

    int nv[3], &nx=nv[0], &ny=nv[1], &nz=nv[2];
    nv[0]=nv[1]=nv[2]=1;


    const int numberOfDimensions = mg.numberOfDimensions();
    const IntegerArray & bc = mg.boundaryCondition();

    if( numberOfDimensions==1 )
    {
        nx=numEigs; ny=1; nz=1; 
    }
    else if( numberOfDimensions==2 )
    {
    // we need to create enough potential eignvalues so that we have the correct numEigs after sorting
        nx=numEigs; ny=numEigs; // this is too many 
    }
    else
    {
        nx=numEigs; 
        ny=numEigs; 
        nz=numEigs; // this is too many 
    }
  // No sense looking for more eigs than are on the grid
    int maxEigs = 1;  // max number of eigs in out list 
    for( int axis=0; axis<numberOfDimensions; axis++ )
    {
        int maxModes = max( 1, gid(1,axis)-gid(0,axis)+1); // max modes this direction
        maxEigs *= maxModes;
        if( mg.isPeriodic(axis) )
            maxModes = max(1, maxModes/2);  // in periodic directions we will double

        nv[axis] = min( nv[axis],maxModes );
    }
    maxEigs += 10; // is this needed?

  // if( bcOpt !=0 ) // finish me
  // {
  //   printF("Ogev::getEigenvaluesBox: bcOpt=%d is not supported yet\n",bcOpt);
  //   OV_ABORT("error");
  // }

  // Make a vector of values and then sort to find smallest
    std::vector<Real> ev;

    RealArray modeNumber(3,maxEigs); // holds mx,my,mz for each eigenMode
    IntegerArray phaseNumber(3,maxEigs); // for periodic directions we have a sine and cosine
    phaseNumber=0;
    Range all;

    int mv[3], &m1=mv[0], &m2=mv[1], &m3=mv[2];
    bool periodicEig[3]={false,false,false};  
    int eigCount=0;  // counts eigenvectors
    for( m3=1; m3<=nz; m3++ )
    for( m2=1; m2<=ny; m2++ )
    for( m1=1; m1<=nx; m1++ )
    {
        Real eigValue = 0.;

        for( int axis=0; axis<numberOfDimensions; axis++ )
        {
            const int m = mv[axis]; 
            assert( eigCount < maxEigs );
            periodicEig[axis]=false; // set to true if this is a periodic double eigenvalue 

            modeNumber(axis,eigCount)=m;
            phaseNumber(axis,eigCount)=0; 

            const Real length = axis==0 ? lx : axis==1 ? ly : lz;
            if( bc(0,axis)==dirichlet && bc(1,axis)==dirichlet )
            {
                if( discreteEigenvalues )
          // eigValue += SQR( sin(m*Pi*dx[axis]/2.)/(dx[axis]/2.) ); 
                    eigValue += getDiscreteSymbol( 1.*m/lv[axis], dx[axis] );
                else
                    eigValue +=       SQR(m*Pi/length);
            }
            else if( (bc(0,axis)==neumann   && bc(1,axis)==dirichlet) ||
                              (bc(0,axis)==dirichlet && bc(1,axis)==neumann  )  ) 
            {
                if( discreteEigenvalues )
          // eigValue += SQR( sin( .5*(2*m-1)*Pi*dx[axis]/2.)/(dx[axis]/2.) );  // check me 
                    eigValue += getDiscreteSymbol( .5*(2*m-1)/lv[axis], dx[axis] );
                else
                    eigValue +=       SQR(.5*(2*m-1)*Pi/length);
            }
            else if( bc(0,axis)==neumann && bc(1,axis)==neumann )
            {
                if( discreteEigenvalues )
          // eigValue += SQR( sin( (m-1)*Pi*dx[axis]/2.)/(dx[axis]/2.) );  // check me 
                    eigValue += getDiscreteSymbol( (m-1.)/lv[axis], dx[axis] );
                else
                    eigValue +=      SQR( (m-1)*Pi/length);
            }
            else if( bc(0,axis)<0 && bc(1,axis)<0 )
            {
        // periodic 
                if( discreteEigenvalues )
          // eigValue += SQR( sin( 2.*(m-1)*Pi*dx[axis]/2.)/(dx[axis]/2.) );  // check me 
                    eigValue += getDiscreteSymbol( 2.*(m-1.)/lv[axis], dx[axis] );
                else
                    eigValue +=      SQR( 2.*(m-1)*Pi/length);

                if( m>1 )
                    periodicEig[axis]=true; // sin and cos 
            }      
            else if( bc(0,axis)>0 || bc(1,axis)>0 )
            {
                printF("Ogev::getEigenvaluesBox:error: unknown bc");
                ::display(bc,"bc");
                OV_ABORT("Ogev::getEigenvaluesBox:error: unknown bc");
            }

        } // end for axis 

        if( false )
        {
            const int mxe = modeNumber(0,eigCount);
            const int mye = modeNumber(1,eigCount);
            const int mze = modeNumber(2,eigCount);
              printF("boxEigs: Add true eig: eigCount=%d, [mx,my,mz]=[%d,%d,%d], eigValue=%14.6e \n",eigCount,mxe,mye,mze,eigValue);     
        }
        ev.push_back(eigValue); 
        eigCount++;

    // -- set the phase for periodic directions ---
        for( int p3=0; p3<=periodicEig[2]; p3++ )
        for( int p2=0; p2<=periodicEig[1]; p2++ )
        for( int p1=0; p1<=periodicEig[0]; p1++ )
        {
            if( p1>0 || p2>0 || p3>0 )
            {
                ev.push_back(eigValue);

                assert( eigCount < maxEigs );
                modeNumber(all,eigCount)=modeNumber(all,eigCount-1);

                phaseNumber(0,eigCount)=p1;  // set phase = 0 or 1 (sine or cosine)
                phaseNumber(1,eigCount)=p2;  
                phaseNumber(2,eigCount)=p3;  
                eigCount++; 
            }     
        }
    // for( int axis=0; axis<numberOfDimensions; axis++ )
    // { 
    //   if( periodicEig[axis] )
    //   {
    //     ev.push_back(eigValue);

    //     assert( eigCount < maxEigs );
    //     modeNumber(all,eigCount)=modeNumber(all,eigCount-1);

    //     phaseNumber(all,eigCount)=0;  
    //     phaseNumber(axis,eigCount)=1;   // this mode is a cosine
    //     eigCount++;
    //   }
    // }
    }
  // OV_ABORT("stop here for now");

    if( !computeEigenvectors )
    {
    // sort eigenvalues from smallest to largest
        sort(ev.begin(),ev.end());

        eigs.redim(numEigs);
        for( int i=0; i<numEigs; i++ )
        {
            eigs(i) = ev[i];
        }
    }
    else
    {
    // ---Return eignvalues AND eigenvectors ---
        printF("\n >>> Ogev::getEigenvaluesBox: return eigenvalues and eigenvectors...\n");


    // --- sort an index array myIndex according to values in ev 
        assert( eigCount>0 );

        vector<int> myIndex(eigCount);
        std::iota(myIndex.begin(),myIndex.end(),0); //Initializing
        sort( myIndex.begin(),myIndex.end(), [&](int i,int j){return ev[i]<ev[j];} );    

    // --- first fill in eigenvalues ---
        IntegerArray phaseNumberNew(3,maxEigs);
        for( int i=0; i<numEigs; i++ )
        {
            const int i0 = myIndex[i]; // sorted indicies
            eigs(i) = ev[i0];

            const int mx = modeNumber(0,i0);
            const int my = modeNumber(1,i0);
            const int mz = modeNumber(2,i0);
            if( debug & 8 )
                printF("boxEigs: i=%3d, i0=%5d: eigs=%12.4e, [lx,ly,lz]=[%g,%g,%g], [mx,my,mz]=[%d,%d,%d] phase=[%d,%d,%d]\n",
                    i,i0,eigs(i),lx,ly,lz,mx,my,mz,phaseNumber(0,i0),phaseNumber(1,i0),phaseNumber(2,i0));      
        }

    // --- fill in eigenvectors ---
        RealCompositeGridFunction & eva = *eigenvector;
        eva.updateToMatchGrid(cg,all,all,all,numEigs);
        eva=1.;    

        Index I1,I2,I3;
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];  
                      
            mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
            OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);
            getIndex(mg.dimension(),I1,I2,I3);

            for( int i=0; i<numEigs; i++ )
            {
                const int i0 = myIndex[i]; // sorted indicies

                for( int axis=0; axis<numberOfDimensions; axis++ )
                {
          // const Real length = axis==0 ? lx : axis==1 ? ly : lz; 
                    const Real length = lv[axis];
                    if( bc(0,axis)==dirichlet && bc(1,axis)==dirichlet )
                    {
                        const Real freq = modeNumber(axis,i0)*Pi/length;        
                        eva[grid](I1,I2,I3,i) *= sin(freq*(xLocal(I1,I2,I3,axis)-xab[0][axis]));
                    }
                    else if( bc(0,axis)==neumann && bc(1,axis)==neumann )
                    {
                        const Real freq = (modeNumber(axis,i0)-1)*Pi/length;   
                        eva[grid](I1,I2,I3,i) *= cos(freq*(xLocal(I1,I2,I3,axis)-xab[0][axis]));            
                    }
                    else if( bc(0,axis)==dirichlet && bc(1,axis)==neumann )
                    {
                        const Real freq = 0.5*(2*modeNumber(axis,i0)-1)*Pi/length;   
                        eva[grid](I1,I2,I3,i) *= sin(freq*(xLocal(I1,I2,I3,axis)-xab[0][axis]));            
                    }
                    else if( bc(0,axis)==neumann && bc(1,axis)==dirichlet )
                    {
                        const Real freq = 0.5*(2*modeNumber(axis,i0)-1)*Pi/length;   
                        eva[grid](I1,I2,I3,i) *= cos(freq*(xLocal(I1,I2,I3,axis)-xab[0][axis]));            
                    }                    
                    else if( bc(0,axis)<0 && bc(1,axis)<0 )
                    {
            // periodic -- could be sin or cos 
                        const Real freq = (modeNumber(axis,i0)-1)*2.*Pi/length;  
                        if( phaseNumber(axis,i0)==0 )
                            eva[grid](I1,I2,I3,i) *= sin(freq*(xLocal(I1,I2,I3,axis)-xab[0][axis]));    
                        else
                            eva[grid](I1,I2,I3,i) *= cos(freq*(xLocal(I1,I2,I3,axis)-xab[0][axis]));         
                    }          
                    else
                    {
                        OV_ABORT("Ogev::getEigenvaluesBox: fill eigenvectors: FINISH ME FOR other BCs");
                    }
                }
            }
        }

    }

  // Real eigTrue = i==0 ? SQR(   Pi/lx) + SQR(   Pi/ly) : 
  //                i==1 ? SQR(   Pi/lx) + SQR(2.*Pi/ly) :
  //                i==2 ? SQR(   Pi/lx) + SQR(2.*Pi/ly) : 
  //                i==3 ? SQR(2.*Pi/lx) + SQR(2.*Pi/ly) :
  //                i==4 ? SQR(3.*Pi/lx) + SQR(1.*Pi/ly) :
  //                i==5 ? SQR(3.*Pi/lx) + SQR(1.*Pi/ly) :
  //                       SQR(2.*Pi/lx) + SQR(3.*Pi/ly) ;


    return 0;
}


//===========================================================================
//
/// \brief Return the eigenvalues, and optionally eigenvectors, of minus Laplacian
///       - Delta( u ) = lambda^2 u 
///  for an annulus, disk, hollow cylinder or cylinder
//
/// \param numEigs (input) : number of eigenvalues requested
/// \param eigs(i)   : eigenvalues from smallest to largest
/// \param numberOfDimensions (input) : 1,2 or 3
/// \param bcOpt (input) : 0=Dirichlet BC's 
/// \param ra,rb,za,zb : cylinder bounds
//===========================================================================
int Ogev::
getEigenvaluesCylinder( int numEigs, RealArray & eigs, CompositeGrid & cg, 
                                                Real ra  /* = 0.5 */ , Real rb  /* =1.0 */ , Real za /* = 0.0 */, Real zb /* = 1.0 */,
                                                RealCompositeGridFunction *eigenvector /* = NULL */  )
{

    const bool computeEigenvectors = eigenvector!=NULL; 
  // if( computeEigenvectors )
  // {
  //   OV_ABORT("getEigenvaluesCylinder: FINISH ME");
  // }

    MappedGrid & mg = cg[0]; // for BC's

    const int dirichlet=1, neumann=2; // **FIX ME**

    const int numberOfDimensions = mg.numberOfDimensions();
    const IntegerArray & bc = mg.boundaryCondition();


  // Make a vector of values and then sort to find smallest
    std::vector<Real> ev;

    if( ra>0.0 && numberOfDimensions==2 )
    {
    // --- annulus  ----

        if( bc(0,1)==dirichlet && bc(1,1)==dirichlet )
        {
            #include "annulusEigenvaluesDirichletDirichlet.h"

            for( int i=0; i<numRoot; i++ )
            {
                Real eigValue = SQR( lambda[i] ); // NOTE SQUARE
                ev.push_back(eigValue); 
            }
        }
        else if( bc(0,1)==neumann && bc(1,1)==neumann )
        {
            #include "annulusEigenvaluesNeumannNeumann.h"

      // Note all but the first root lambda=0 is a double root I think 
            for( int i=0; i<numRoot; i++ )
            {
                Real eigValue = SQR( lambda[i] ); // NOTE SQUARE
                ev.push_back(eigValue); 
                if( i>0 )
                    ev.push_back(eigValue); // double root
            }
        }
        else
        {
            printF("Ogev::getEigenvaluesCylinder: These boundary conditions not implemented yet\n");
            ::display(bc,"bc");
            OV_ABORT("error");
        }


    }
  // else if( numberOfDimensions==2 && ra==0. )
    else if( ra==0. )
    {
    // ------ DISK OR CYLINDER -------

        printf("Get true eigs for a disk or cyl, numberOfDimensions=%d, ra=%g, rb=%g, za=%g, zb=%g\n",numberOfDimensions,ra,rb,za,zb); 

    // DIRICHLET BC
            
    // File generated by mx/codes/bessel.maple
    // besselZeros[ndbz][mdbz]
        #include "besselZeros.h" 

        const Real Lr = rb-ra; 
        const Real Lz = zb-za; 
        
    // modeNumber(0,eigCount) = n 
    // modeNumber(1,eigCount) = m
    // modeNumber(2,eigCount) = l  ( for z direction in 3D)
    // modeNumber(3,eigCount) = 1,2  for multiple eigenvalues

        const int ndlz   = numberOfDimensions==2 ? 1 : ndbz; // max modes in z-direction
        const int lStart = numberOfDimensions==2 ? 0 : 1;   
        const int lEnd   = numberOfDimensions==2 ? 0 : ndlz;

        RealArray modeNumber(4,ndbz*mdbz*ndlz*2); // holds n,m,for each eigenvalue so we can form the eigenvector
        int eigCount=0; 
        for( int n=0; n<ndbz; n++ )
        {
            for( int m=0; m<mdbz; m++ )
            {

                for( int l=lStart; l<=lEnd; l++ ) // z-direction 
                {        

                    Real eigValue = SQR( besselZeros[n][m]/rb )  + SQR( l*Pi/Lz); // NOTE SQUARE

                    modeNumber(0,eigCount)=n; 
                    modeNumber(1,eigCount)=m; 
                    modeNumber(2,eigCount)=l; 
                    modeNumber(3,eigCount)=1; eigCount++;

                    ev.push_back(eigValue); 
          // Double root occurs if n >0 : cos(n*theta) and sin(n*theta)
                    if( n>0 )
                    {
                        ev.push_back(eigValue); // double root
                        modeNumber(0,eigCount)=n; 
                        modeNumber(1,eigCount)=m;  
                        modeNumber(2,eigCount)=l;  
                        modeNumber(3,eigCount)=2;  eigCount++;
                    }
                }
            }
        }


    // const Real jzmn = besselZeros[n][m];  // m'th zero of Jn
    // lambda=jzmn/a; 
    //  ua(i1,i2,i3,nn) = jn(n,lambda*r)*cos(n*theta)   

        if( !computeEigenvectors )
        {
            sort(ev.begin(),ev.end());

            eigs.redim(numEigs);
            eigs=0;
            int maxEigs = min(numEigs,ev.size()); 
            for( int i=0; i<maxEigs; i++ )
            {
                eigs(i) = ev[i];
            }
        }
        else
        {
      // ---Return eignvalues AND eigenvectors ---

            printF("\n >>> Ogev::getEigenvaluesCylinder: return eigenvalues and eigenvectors...\n");

      // --- sort an index array myIndex according to values in ev 
            vector<int> myIndex(eigCount);
            std::iota(myIndex.begin(),myIndex.end(),0); //Initializing
            sort( myIndex.begin(),myIndex.end(), [&](int i,int j){return ev[i]<ev[j];} );    

      // --- first fill in eigenvalues ---
            eigs.redim(numEigs);
            eigs=0;
            int maxEigs = min(numEigs,ev.size()); 
            for( int i=0; i<maxEigs; i++ )
            {
                const int i0 = myIndex[i]; // sorted indicies      
                eigs(i) = ev[i0];
            }

      // --- fill in eigenvectors ---
            RealCompositeGridFunction & eva = *eigenvector;
            Range all;
            eva.updateToMatchGrid(cg,all,all,all,numEigs);
            eva=0.;    

            Index I1,I2,I3;
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];  
                          
                mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
                OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);
                getIndex(mg.dimension(),I1,I2,I3);

                for( int i=0; i<numEigs; i++ )
                {
                    const int i0 = myIndex[i]; // sorted indicies
                    const int n       = modeNumber(0,i0);
                    const int m       = modeNumber(1,i0);
                    const int l       = modeNumber(2,i0);
                    const int eigMult = modeNumber(3,i0); // multiplicity number: 1 or 2
                    const Real lamz = l*Pi/Lz; 
                    const Real lambda = sqrt(eigs(i) - lamz*lamz ); 

                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        const Real xd = xLocal(i1,i2,i3,0), yd = xLocal(i1,i2,i3,1);
                        const Real theta = atan2(yd,xd);
                        const Real r = sqrt( xd*xd + yd*yd );
                        const Real zFact = numberOfDimensions==2 ? 1.0 : sin(lamz*(xLocal(i1,i2,i3,2)-za));             
                        if( eigMult==1 )
                            eva[grid](i1,i2,i3,i) = jn(n,lambda*r)*cos(n*theta)*zFact; 
                        else        
                            eva[grid](i1,i2,i3,i) = jn(n,lambda*r)*sin(n*theta)*zFact; 

                    }        
                }
            }

        }

    }  
    else
    {
        printF("Ogev::getEigenvaluesCylinder: Error: finish me for numberOfDimensions=%d, ra=%g, rb=%g\n",
                      numberOfDimensions,ra,rb);
        OV_ABORT("error");
    }

    return 0;
}


//===========================================================================
//
/// \brief Return the eigenvalues, and optionally eigenvectors, of minus Laplacian
///       - Delta( u ) = lambda^2 u 
///  for a SPHERE or SPHERICAL SHELL
//
/// \param numEigs (input) : number of eigenvalues requested
/// \param eigs(i)   : eigenvalues from smallest to largest
/// \param numberOfDimensions (input) : 1,2 or 3
/// \param bcOpt (input) : 0=Dirichlet BC's 
/// \param ra,rb : sphere inner and outer radii
//===========================================================================


int Ogev::
getEigenvaluesSphere( int numEigs, RealArray & eigs, CompositeGrid & cg, 
                                                Real ra /* =0.0 */, Real rb /* =1.0 */, 
                                                RealCompositeGridFunction *eigenvector /* = NULL  */ )
{

    const bool computeEigenvectors = eigenvector!=NULL; 
  // if( computeEigenvectors )
  // {
  //   OV_ABORT("getEigenvaluesCylinder: FINISH ME");
  // }

    MappedGrid & mg = cg[0]; // for BC's

    const int dirichlet=1, neumann=2; // **FIX ME**

    const int numberOfDimensions = mg.numberOfDimensions();
    const IntegerArray & bc = mg.boundaryCondition();


  // Make a vector of values and then sort to find smallest
    std::vector<Real> ev;


    if( ra==0. )
    {
    // ------ SOLID SPHERE -------

        printf("Get true eigs for a SPHERE, numberOfDimensions=%d, ra=%g, rb=%g\n",numberOfDimensions,ra,rb); 


      
    // DIRICHLET BC
    // #include "src/sphereEigenvalueDirichletRoots.h"    

    // ***** Values from cgwave/maple/sphereEigenvalus.mpl *****

        const int mPhiMax =5, mrMax=4, mThetaMax=2*mPhiMax+1;
    // n = mPhi, m=mr 
        Real omegaArray[]={
        // mPhi=0, mr=1,2,..,4
                3.1415926535897932e+00, 6.2831853071795865e+00, 9.4247779607693797e+00, 1.2566370614359173e+01,
        // mPhi=1, mr=1,2,..,4
                4.4934094579090642e+00, 7.7252518369377072e+00, 1.0904121659428900e+01, 1.4066193912831473e+01,
        // mPhi=2, mr=1,2,..,4
                5.7634591968945498e+00, 9.0950113304763552e+00, 1.2322940970566582e+01, 1.5514603010886748e+01,
        // mPhi=3, mr=1,2,..,4
                6.9879320005005200e+00, 1.0417118547379365e+01, 1.3698023153249249e+01, 1.6923621285213840e+01,
        // mPhi=4, mr=1,2,..,4
                8.1825614525712427e+00, 1.1704907154570391e+01, 1.5039664707616521e+01, 1.8301255959541990e+01,        
        };
        #define omegaRoot(mPhi,mr) omegaArray[(mr-1)+mrMax*(mPhi-0)]

        RealArray modeNumber(4,mPhiMax*mThetaMax*mrMax); // holds (mPhi,mr,mTheta) each eigenvalue so we can form the eigenvector
        int eigCount=0; 
        for( int mPhi=0; mPhi<mPhiMax; mPhi++ )
        {
      // for( int mr=0; mr<mrMax; mr++ ) // r-direction : different zeros of J(mPhi+1/2, xi ) 
            for( int mr=1; mr<=mrMax; mr++ )   // r-direction : different zeros of J(mPhi+1/2, xi )  *NOTE* mr starts at "1"
            {        
                Real eigValue= SQR( omegaRoot(mPhi,mr)/rb );

                for( int mTheta=-mPhi; mTheta<=mPhi; mTheta++ ) // multiple eigs have different values for mTheta
                {
                    modeNumber(0,eigCount) = mPhi;
                    modeNumber(1,eigCount) = mr; 
                    modeNumber(2,eigCount) = mTheta; 
                    modeNumber(3,eigCount) = 2*mPhi+1; // multiplicity of this eig. Is this used ??
                    eigCount++;
                    ev.push_back(eigValue); 
                }
            }
        }

        if( !computeEigenvectors )
        {
            sort(ev.begin(),ev.end());

            eigs.redim(numEigs);
            eigs=0;
            int maxEigs = min(numEigs,ev.size()); 
            for( int i=0; i<maxEigs; i++ )
            {
                eigs(i) = ev[i];
            }
        }
        else
        {
      // ---Return eignvalues AND eigenvectors ---

            printF("\n >>> Ogev::getEigenvaluesSphere: return eigenvalues and eigenvectors...\n");

      // --- sort an index array myIndex according to values in ev 
            vector<int> myIndex(eigCount);
            std::iota(myIndex.begin(),myIndex.end(),0); //Initializing
            sort( myIndex.begin(),myIndex.end(), [&](int i,int j){return ev[i]<ev[j];} );    

      // --- first fill in eigenvalues ---
            eigs.redim(numEigs);
            eigs=0;
            int maxEigs = min(numEigs,ev.size()); 
            for( int i=0; i<maxEigs; i++ )
            {
                const int i0 = myIndex[i]; // sorted indicies      
                eigs(i) = ev[i0];
            }

      // --- fill in eigenvectors ---
            RealCompositeGridFunction & eva = *eigenvector;
            Range all;
            eva.updateToMatchGrid(cg,all,all,all,numEigs);
            eva=0.;    

            const Real rTol = 1e-3;                 // change formula for r smaller than this  
            const Real phiTol= REAL_EPSILON*1000.;  // what should this be ?
            const Real epsx= REAL_EPSILON*100.;
            
            const Real alpha=0.5;    // fractional part of Bessel order
            int ncalc;         // output flag
            Real sphericalBessel;   


            const int maxBesselOrder=mPhiMax+1;
            Real jnka[maxBesselOrder];           // holds evaluation of Bessel function   

            Index I1,I2,I3;
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];  
                          
                mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
                OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);
                OV_GET_SERIAL_ARRAY(Real,eva[grid],evLocal);

                getIndex(mg.dimension(),I1,I2,I3);

                for( int i=0; i<numEigs; i++ )
                {
                    const int i0 = myIndex[i]; // sorted indicies

                    const int mPhi    = modeNumber(0,i0);  
                    const int mr      = modeNumber(1,i0);
                    const int mTheta  = modeNumber(2,i0);    // *** CAN WE MOVE THE loop over mTheta inside the loop below ?
                    const int mThetaAbs = abs(mTheta); 

                    const Real lambda = sqrt(eigs(i));  // note: eigs are sorted, use "i" as the index 

                    const int nb = mPhi+1; // nterm+1  ! eval J0, J1, ... J(nb)  --           

                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    { 
              
                        Real xd = xLocal(i1,i2,i3,0), yd = xLocal(i1,i2,i3,1), zd = xLocal(i1,i2,i3,2);
                        if( abs(xd)<epsx && abs(yd)<epsx )
                        {
                            xd=epsx;  //  avoid atan2(0,0)
                            yd=epsx;
                        }
                        const Real theta = atan2(yd,xd);
                        Real r = sqrt( xd*xd + yd*yd + zd*zd );

            // const Real cosphi = (abs(r)<tol) ? (zd/tol) : (zd/r);
                        const Real cosphi = fabs(r)<phiTol ? 1.0 : zd/r;

                        Real Pnk; 
                            if( mPhi == 0 )
                            {
                                Pnk = 1.;
                            }  
                            else if( mPhi == 1 )
                            {
                                if( mThetaAbs==0 )
                                    Pnk = cosphi;
                                else
                                    Pnk = -sqrt(1.-cosphi*cosphi);
                            }
                            else if( mPhi == 2 )
                            {
                                if( mThetaAbs==0 )
                                    Pnk = .5*( 3.0*cosphi*cosphi - 1.0 );
                                else if( mThetaAbs==1 )
                                    Pnk = - 3.0*cosphi*sqrt(1. - cosphi*cosphi);
                                else
                                    Pnk = 3.0*( 1.-cosphi*cosphi ); 
                            }
                            else if( mPhi == 3 )
                            {
                                if( mThetaAbs==0 )
                                    Pnk = .5*cosphi*( 5.0*cosphi*cosphi - 3.0 );
                                else if( mThetaAbs==1 )
                                    Pnk = 1.5*(1. - 5.*cosphi*cosphi)*sqrt(1 - cosphi*cosphi);
                                else if( mThetaAbs==2 )
                                    Pnk = 15.0*cosphi*(1.0-cosphi*cosphi);
                                else
                                    Pnk = -15.*pow( 1.0 -cosphi*cosphi, 1.5 );
                            }
                            else if( mPhi == 4 )
                            {
                                if( mThetaAbs==0 )
                                    Pnk = (1./8.)*( 3. + (cosphi*cosphi)*( -30. + (cosphi*cosphi)*35. ) );
                                else if( mThetaAbs==1 )
                                    Pnk = 2.5*cosphi*( 3.- 7.*cosphi*cosphi)*sqrt(1.0 -cosphi*cosphi);
                                else if( mThetaAbs==2 )
                                    Pnk = (15./2.)*( 7.*cosphi*cosphi-1. )*(1.-cosphi*cosphi);
                                else if( mThetaAbs==3 )
                                    Pnk = -105.*cosphi * pow( 1.0 -cosphi*cosphi, 1.5 );
                                else
                                    Pnk = 105.*SQR( 1. - cosphi*cosphi );
                            }
                            else
                            {
                                printF("associatedLegendreFunctions ERROR: mPhi=%d, mThetaAbs=%d not supported.\n",mPhi,mThetaAbs);
                                OV_ABORT("ERROR");
                            }

                        Real kr=lambda*r;  // argument of the Bessel function

                        if( fabs(r)<rTol )
                        {
              // Use Small r approximations for the spherical Bessel: 
                                if( mPhi == 0 )
                                {
                                    sphericalBessel = 1.0 - (1.0/6.0)*(kr*kr) + (1.0/120.0)*(kr*kr*kr*kr);
                                }
                                else if( mPhi == 1 )
                                {
                                    sphericalBessel = (1.0/3.0)*kr - (1.0/30.0)*(kr*kr*kr) + (1.0/840.0)*(kr*kr*kr*kr*kr);
                                }
                                else if( mPhi == 2 )
                                {
                                    sphericalBessel = (1.0/15.0)*(kr*kr) - (1.0/210.0)*(kr*kr*kr*kr);
                                }
                                else if( mPhi == 3 )
                                {
                                    sphericalBessel = (1.0/105.0)*(kr*kr*kr) - (1.0/1890.0)*(kr*kr*kr*kr*kr);
                                }
                                else
                                {
                                    printF("sphericalBesselFunctions ERROR: mPhi > 4 is not defined\n");
                                    OV_ABORT("ERROR");
                                }
                        }
                        else
                        {
              // Note: To get J_(mPhi+1/2, kr)  we want we have to evaluate previous J_( m, kr ) m=0,1,...,mPhi-1
                            rjbesl(kr, alpha, nb, jnka[0], ncalc); // this computes jn_0, jn_1, ..., jn_{nb-1}
                            sphericalBessel = sqrt(Pi/(2.0*r))*jnka[mPhi]; 
                        }
                        Real cosSineTheta = mTheta<=0 ? cos(mTheta*theta) : sin(mTheta*theta);

                        evLocal(i1,i2,i3,i) = sphericalBessel*cosSineTheta*Pnk;

                    } // end for 3d 


                } // end for i 
            } // end for grid 

        } // end save eigs and eigenvectors 

    }  // end if ra==0.
    else
    {
        printF("Ogev::getEigenvaluesSphere: Error: finish me for a spherical shell: ra=%g, rb=%g\n",ra,rb);
        OV_ABORT("error");
    }


    return 0;
}







