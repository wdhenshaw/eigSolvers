//===========================================================================
//
//  Overlapping Grid EigenValue/EigenVector solver
//
//===========================================================================

#include "Ogev.h"
#include "display.h"
#include <numeric>      // std::iota

#define FOR_3D(i1,i2,i3,I1,I2,I3)                                       \
int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(); \
int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
for(int i3=I3Base; i3<=I3Bound; i3++)                                       \
  for(int i2=I2Base; i2<=I2Bound; i2++)                                     \
    for(int i1=I1Base; i1<=I1Bound; i1++)

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
  else if( numberOfDimensions==2 && ra==0. )
  {
    // ------ DISK -------
    // if( bcOpt !=0 ) // finish me 
    // {
    //   printF("Ogev::getEigenvaluesCylinder: bcOpt=%d is not supported yet (Neumann BC for disk) \n",bcOpt);
    //   OV_ABORT("error");
    // }

    // DIRICHLET BC
      
    // File generated by mx/codes/bessel.maple
    // besselZeros[ndbz][mdbz]
    #include "besselZeros.h" 
    
    // modeNumber(0,eigCount) = n 
    // modeNumber(1,eigCount) = m
    // modeNumber(2,eigCount) = 1,2  for multiple eigenvalues
    RealArray modeNumber(3,ndbz*mdbz*2); // holds n,m,for each eigenvalue so we can form the eigenvector
    int eigCount=0; 
    for( int n=0; n<ndbz; n++ )
    {
      for( int m=0; m<mdbz; m++ )
      {
        Real eigValue = SQR( besselZeros[n][m] ); // NOTE SQUARE
        modeNumber(0,eigCount)=n; modeNumber(1,eigCount)=m; modeNumber(2,eigCount)=1; eigCount++;
        ev.push_back(eigValue); 
        // Double root occurs if n >0 : cos(n*theta) and sin(n*theta)
        if( n>0 )
        {
          ev.push_back(eigValue); // double root
          modeNumber(0,eigCount)=n; modeNumber(1,eigCount)=m;  modeNumber(2,eigCount)=2;  eigCount++;
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
          const int n = modeNumber(0,i0);
          const int m = modeNumber(1,i0);
          const int eigMult = modeNumber(2,i0); // multiplicity number: 1 or 2
          const Real lambda = sqrt(eigs(i))/rb; 

          FOR_3D(i1,i2,i3,I1,I2,I3)
          {
            const Real xd = xLocal(i1,i2,i3,0), yd = xLocal(i1,i2,i3,1);
            const Real theta = atan2(yd,xd);
            const Real r = sqrt( xd*xd + yd*yd );            
            if( eigMult==1 )
              eva[grid](i1,i2,i3,i) = jn(n,lambda*r)*cos(n*theta); 
            else        
              eva[grid](i1,i2,i3,i) = jn(n,lambda*r)*sin(n*theta); 
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



  // else if( numberOfDimensions==2 && ra==0. )
  // {
  //   // ***** OLD VERSION ***
  //   // ---- DISK -------

  //   // % ---- Eigenvalues for eigenfunctions of the heat equation on an disk ----
  //   // % File written by cg/ad/codes/diskEigenvalues.maple 
  //   // % numBesselOrderDBC  : Bessel orders are m=0,1,2,...,numBesselOrderDBC-1 
  //   // % numRootDBC         : number of roots
  //   // % ----------- DIRCHLET BCs at r=rb ----- 
  //   // % Solution is u = cJ*Jn(lambda*r)
  //   // rb=1.00000000000000e+00; 
  //   // numBesselOrder=30; numRoot=15;
  //      const int numRoot = 30*15; 
  //      Real lambda [] ={
  //      2.40482555769577e+00,
  //      5.52007811028631e+00,
  //      8.65372791291101e+00,
  //      1.17915344390143e+01,
  //      1.49309177084878e+01,
  //      1.80710639679109e+01,
  //      2.12116366298793e+01,
  //      2.43524715307493e+01,
  //      2.74934791320403e+01,
  //      3.06346064684320e+01,
  //      3.37758202135736e+01,
  //      3.69170983536640e+01,
  //      4.00584257646282e+01,
  //      4.31997917131767e+01,
  //      4.63411883716618e+01,
  //      3.83170597020751e+00,
  //      3.83170597020751e+00,  // double 
  //      7.01558666981562e+00,
  //      1.01734681350627e+01,
  //      1.33236919363142e+01,
  //      1.64706300508776e+01,
  //      1.96158585104682e+01,
  //      2.27600843805928e+01,
  //      2.59036720876184e+01,
  //      2.90468285349169e+01,
  //      3.21896799109744e+01,
  //      3.53323075500839e+01,
  //      3.84747662347716e+01,
  //      4.16170942128145e+01,
  //      4.47593189976528e+01,
  //      4.79014608871854e+01,
  //      5.13562230184068e+00,
  //      5.13562230184068e+00,  // double 
  //      8.41724414039986e+00,
  //      1.16198411721491e+01,
  //      1.47959517823513e+01,
  //      1.79598194949878e+01,
  //      2.11169970530218e+01,
  //      2.42701123135731e+01,
  //      2.74205735499846e+01,
  //      3.05692044955164e+01,
  //      3.37165195092227e+01,
  //      3.68628565112838e+01,
  //      4.00084467334782e+01,
  //      4.31534537783715e+01,
  //      4.62979966772369e+01,
  //      4.94421641104169e+01,
  //      6.38016189592398e+00,
  //      6.38016189592398e+00,  // double 
  //      9.76102312998167e+00,
  //      1.30152007216984e+01,
  //      1.62234661603188e+01,
  //      1.94094152264350e+01,
  //      2.25827295931044e+01,
  //      2.57481666992950e+01,
  //      2.89083507809218e+01,
  //      3.20648524070977e+01,
  //      3.52186707386101e+01,
  //      3.83704724347569e+01,
  //      4.15207196704068e+01,
  //      4.46697431166173e+01,
  //      4.78177856915333e+01,
  //      5.09650299062052e+01,
  //      7.58834243450380e+00,
  //      1.10647094885012e+01,
  //      1.43725366716176e+01,
  //      1.76159660498048e+01,
  //      2.08269329569624e+01,
  //      2.40190195247711e+01,
  //      2.71990877659813e+01,
  //      3.03710076671172e+01,
  //      3.35371377118192e+01,
  //      3.66990011287446e+01,
  //      3.98576273021809e+01,
  //      4.30137377233544e+01,
  //      4.61678535129244e+01,
  //      4.93203606863903e+01,
  //      5.24715513984580e+01,
  //      8.77148381595995e+00,
  //      1.23386041974669e+01,
  //      1.57001740797117e+01,
  //      1.89801338751799e+01,
  //      2.22177998965613e+01,
  //      2.54303411542227e+01,
  //      2.86266183072911e+01,
  //      3.18117167240478e+01,
  //      3.49887812945593e+01,
  //      3.81598685619671e+01,
  //      4.13263832540474e+01,
  //      4.44893191232197e+01,
  //      4.76493998066971e+01,
  //      5.08071652030063e+01,
  //      5.39630265583781e+01,
  //      9.93610952421768e+00,
  //      1.35892901705412e+01,
  //      1.70038196678160e+01,
  //      2.03207892135665e+01,
  //      2.35860844355814e+01,
  //      2.68201519834114e+01,
  //      3.00337223865705e+01,
  //      3.32330417628471e+01,
  //      3.64220196682585e+01,
  //      3.96032394160754e+01,
  //      4.27784816131995e+01,
  //      4.59490159980426e+01,
  //      4.91157737247643e+01,
  //      5.22794539036011e+01,
  //      5.54405920688531e+01,
  //      1.10863700192451e+01,
  //      1.48212687270132e+01,
  //      1.82875828324817e+01,
  //      2.16415410198484e+01,
  //      2.49349278876730e+01,
  //      2.81911884594832e+01,
  //      3.14227941922656e+01,
  //      3.46370893520693e+01,
  //      3.78387173828536e+01,
  //      4.10307736915855e+01,
  //      4.42154085052613e+01,
  //      4.73941657555705e+01,
  //      5.05681846797956e+01,
  //      5.37383253719633e+01,
  //      5.69052499919788e+01,
  //      1.22250922640047e+01,
  //      1.60377741908877e+01,
  //      1.95545364309971e+01,
  //      2.29451731318746e+01,
  //      2.62668146411766e+01,
  //      2.95456596709985e+01,
  //      3.27958000373415e+01,
  //      3.60256150638696e+01,
  //      3.92404479951781e+01,
  //      4.24438877432736e+01,
  //      4.56384441821991e+01,
  //      4.88259303815539e+01,
  //      5.20076914566869e+01,
  //      5.51847479392890e+01,
  //      5.83578890252697e+01,
  //      1.33543004774353e+01,
  //      1.72412203824891e+01,
  //      2.08070477892641e+01,
  //      2.42338852577506e+01,
  //      2.75837489635730e+01,
  //      3.08853789676967e+01,
  //      3.41543779238551e+01,
  //      3.74000999771566e+01,
  //      4.06285537189645e+01,
  //      4.38438014203373e+01,
  //      4.70487007376540e+01,
  //      5.02453269553054e+01,
  //      5.34352271570421e+01,
  //      5.66195802665084e+01,
  //      5.97993016309602e+01,
  //      1.44755006865545e+01,
  //      1.84334636669666e+01,
  //      2.20469853646978e+01,
  //      2.55094505541828e+01,
  //      2.88873750635305e+01,
  //      3.22118561997127e+01,
  //      3.54999092053739e+01,
  //      3.87618070178817e+01,
  //      4.20041902366718e+01,
  //      4.52315741035350e+01,
  //      4.84471513872694e+01,
  //      5.16532516681659e+01,
  //      5.48516190759633e+01,
  //      5.80435879282325e+01,
  //      6.12301979772927e+01,
  //      1.55898478844555e+01,
  //      1.96159669039669e+01,
  //      2.32758537262634e+01,
  //      2.67733225455095e+01,
  //      3.01790611787849e+01,
  //      3.35263640755886e+01,
  //      3.68335713418949e+01,
  //      4.01118232709542e+01,
  //      4.33683609475217e+01,
  //      4.66081326762749e+01,
  //      4.98346535103967e+01,
  //      5.30504989591351e+01,
  //      5.62576047151145e+01,
  //      5.94574569083880e+01,
  //      6.26512173882029e+01,
  //      1.66982499338482e+01,
  //      2.07899063600784e+01,
  //      2.44948850438814e+01,
  //      2.80267099499731e+01,
  //      3.14599600353180e+01,
  //      3.48299869902902e+01,
  //      3.81563775046814e+01,
  //      4.14510923079397e+01,
  //      4.47219435431911e+01,
  //      4.79742935312690e+01,
  //      5.12119670041011e+01,
  //      5.44377769283251e+01,
  //      5.76538448119069e+01,
  //      6.08618046824805e+01,
  //      6.40629378248501e+01,
  //      1.78014351532824e+01,
  //      2.19562440678363e+01,
  //      2.57051030539247e+01,
  //      2.92706304418748e+01,
  //      3.27310533109784e+01,
  //      3.61236576664488e+01,
  //      3.94692068252439e+01,
  //      4.27804392654472e+01,
  //      4.60657109115756e+01,
  //      4.93307800964435e+01,
  //      5.25797690643834e+01,
  //      5.58157198763058e+01,
  //      5.90409340372493e+01,
  //      6.22571893937317e+01,
  //      6.54658837972321e+01,
  //      1.88999979531740e+01,
  //      2.31157783472528e+01,
  //      2.69073689761821e+01,
  //      3.05059501638960e+01,
  //      3.39931849847815e+01,
  //      3.74081851286397e+01,
  //      4.07728278535019e+01,
  //      4.41005905657983e+01,
  //      4.74003477805432e+01,
  //      5.06782369464799e+01,
  //      5.39386662091269e+01,
  //      5.71848985981193e+01,
  //      6.04194098521303e+01,
  //      6.36441175089623e+01,
  //      6.68605330122601e+01,
  //      1.99944306298164e+01,
  //      2.42691800262089e+01,
  //      2.81024152316678e+01,
  //      3.17334133443745e+01,
  //      3.52470867857933e+01,
  //      3.86842763892896e+01,
  //      4.20679169986568e+01,
  //      4.54121896147331e+01,
  //      4.87264641162407e+01,
  //      5.20172412788816e+01,
  //      5.52892041465600e+01,
  //      5.85458289043851e+01,
  //      6.17897598959450e+01,
  //      6.50230502510422e+01,
  //      6.82473219964208e+01,
  //      2.10851461130647e+01,
  //      2.54170190063428e+01,
  //      2.92908706963134e+01,
  //      3.29536648850688e+01,
  //      3.64933979124465e+01,
  //      3.99525534902212e+01,
  //      4.33550732038518e+01,
  //      4.67158094350258e+01,
  //      5.00446060168425e+01,
  //      5.33483123318529e+01,
  //      5.66318759428130e+01,
  //      5.98989787287768e+01,
  //      6.31524281933803e+01,
  //      6.63944090385882e+01,
  //      6.96266508798697e+01,
  //      2.21724946188263e+01,
  //      2.65597841380254e+01,
  //      3.04732799463335e+01,
  //      3.41672678538405e+01,
  //      3.77326805220540e+01,
  //      4.12135670591350e+01,
  //      4.46348297531127e+01,
  //      4.80119629360655e+01,
  //      5.13552646507517e+01,
  //      5.46719191775795e+01,
  //      5.79671288334661e+01,
  //      6.12447740981704e+01,
  //      6.45078204027352e+01,
  //      6.77585801138869e+01,
  //      7.09988874898542e+01,
  //      2.32567760851100e+01,
  //      2.76978983508553e+01,
  //      3.16501181518570e+01,
  //      3.53747172179548e+01,
  //      3.89654320476540e+01,
  //      4.24678072133308e+01,
  //      4.59076638663650e+01,
  //      4.93011113380306e+01,
  //      5.26588836515202e+01,
  //      5.59884872205543e+01,
  //      5.92953699442867e+01,
  //      6.25836041801356e+01,
  //      6.58563082805149e+01,
  //      6.91159185022860e+01,
  //      7.23643708714945e+01,
  //      2.43382496234072e+01,
  //      2.88317303513009e+01,
  //      3.28218027618736e+01,
  //      3.65764507589996e+01,
  //      4.01920951008371e+01,
  //      4.37157124179554e+01,
  //      4.71740045682365e+01,
  //      5.05836711400240e+01,
  //      5.39558652828627e+01,
  //      5.72984036543732e+01,
  //      6.06169711271735e+01,
  //      6.39158255761563e+01,
  //      6.71982335006616e+01,
  //      7.04667514173550e+01,
  //      7.37234143308358e+01,
  //      2.54171408140725e+01,
  //      2.99616037916252e+01,
  //      3.39887027852352e+01,
  //      3.77728578443991e+01,
  //      4.14130655138926e+01,
  //      4.49576767484217e+01,
  //      4.84342391952057e+01,
  //      5.18600199280746e+01,
  //      5.52465756146393e+01,
  //      5.86020220738467e+01,
  //      6.19322730728826e+01,
  //      6.52417659938166e+01,
  //      6.85339109388210e+01,
  //      7.18113812037196e+01,
  //      7.50763080770358e+01,
  //      2.64936474160190e+01,
  //      3.10878045460303e+01,
  //      3.51511462445783e+01,
  //      3.89642865475951e+01,
  //      4.26286989308185e+01,
  //      4.61940558943890e+01,
  //      4.96887188180811e+01,
  //      5.31305012504038e+01,
  //      5.65313488968685e+01,
  //      5.98996663967787e+01,
  //      6.32415888283659e+01,
  //      6.65617274042516e+01,
  //      6.98636315104012e+01,
  //      7.31500878919949e+01,
  //      7.64233215263568e+01,
  //      2.75679438912622e+01,
  //      3.22105865494821e+01,
  //      3.63094262234586e+01,
  //      4.01510494809325e+01,
  //      4.38393162543919e+01,
  //      4.74251721615751e+01,
  //      5.09377627926220e+01,
  //      5.43954287365206e+01,
  //      5.78104912784283e+01,
  //      6.11916342175167e+01,
  //      6.45452068206938e+01,
  //      6.78759887703120e+01,
  //      7.11876646342920e+01,
  //      7.44831314264882e+01,
  //      7.77647053193715e+01,
  //      2.86401850307640e+01,
  //      3.33301765307659e+01,
  //      3.74638058174514e+01,
  //      4.13334286141939e+01,
  //      4.50452081825569e+01,
  //      4.86513186682398e+01,
  //      5.21816626034961e+01,
  //      5.56550895985435e+01,
  //      5.90842839864422e+01,
  //      6.24781996896588e+01,
  //      6.58433934695246e+01,
  //      6.91848084146234e+01,
  //      7.25062603808247e+01,
  //      7.58107536154390e+01,
  //      7.91006930938100e+01,
  //      2.97105088898112e+01,
  //      3.44467778847176e+01,
  //      3.86145222220815e+01,
  //      4.25116792859419e+01,
  //      4.62466390132317e+01,
  //      4.98727628890996e+01,
  //      5.34206851310378e+01,
  //      5.69097476241952e+01,
  //      6.03529860588510e+01,
  //      6.37596160161843e+01,
  //      6.71363954544691e+01,
  //      7.04884260839277e+01,
  //      7.38196513520235e+01,
  //      7.71331798448655e+01,
  //      8.04315030482500e+01,
  //      3.07790391865673e+01,
  //      3.55605738670345e+01,
  //      3.97617901342240e+01,
  //      4.36860335666763e+01,
  //      4.74438498564916e+01,
  //      5.10897496663600e+01,
  //      5.46550754431627e+01,
  //      5.81596457487303e+01,
  //      6.16168367044990e+01,
  //      6.50361176103562e+01,
  //      6.84244416914216e+01,
  //      7.17870647559507e+01,
  //      7.51280543337913e+01,
  //      7.84506205900872e+01,
  //      8.17573393260169e+01,
  //      3.18458872786873e+01,
  //      3.66717302506885e+01,
  //      4.09058046024923e+01,
  //      4.48567030971160e+01,
  //      4.86370613629595e+01,
  //      5.23025037824672e+01,
  //      5.58850591940292e+01,
  //      5.94050082759006e+01,
  //      6.28760573508886e+01,
  //      6.63079219790168e+01,
  //      6.97077450622004e+01,
  //      7.30809322269462e+01,
  //      7.64316717529225e+01,
  //      7.97632727533480e+01,
  //      8.30783932439071e+01,
  //      3.29111538049842e+01,
  //      3.77803975505171e+01,
  //      4.20467434316519e+01,
  //      4.60238814981534e+01,
  //      4.98264760544082e+01,
  //      5.35112321699736e+01,
  //      5.71108446949870e+01,
  //      6.06460428046863e+01,
  //      6.41308534286104e+01,
  //      6.75752313706000e+01,
  //      7.09865039340716e+01,
  //      7.43702225115652e+01,
  //      7.77306929663708e+01,
  //      8.10713208513562e+01,
  //      8.43948443864676e+01,
  //      3.39749300587487e+01,
  //      3.88867128985445e+01,
  //      4.31847692232648e+01,
  //      4.71877464287473e+01,
  //      5.10122803250912e+01,
  //      5.47161258190625e+01,
  //      5.83326247108322e+01,
  //      6.18829419079998e+01,
  //      6.53814159326284e+01,
  //      6.88382342227700e+01,
  //      7.22609034994972e+01,
  //      7.56551170814429e+01,
  //      7.90252954057647e+01,
  //      8.23749380728051e+01,
  //      8.57068615833657e+01,
  //      3.50372991442602e+01,
  //      3.99908016345971e+01,
  //      4.43200311174680e+01,
  //      4.83484613524947e+01,
  //      5.21946461688142e+01,
  //      5.59173614317029e+01,
  //      5.95505780245508e+01,
  //      6.31158846012404e+01,
  //      6.66279227935537e+01,
  //      7.00971064386397e+01,
  //      7.35311169614047e+01,
  //      7.69357859646205e+01,
  //      8.03156455965276e+01,
  //      8.36742872228401e+01,
  //      8.70146037847473e+01
  //     };
  //     for( int i=0; i<numRoot; i++ )
  //     {
  //       Real eigValue = SQR( lambda[i] ); // NOTE SQUARE
  //       ev.push_back(eigValue); 
  //       // ** WE ARE MISSING SOME DOUBLE ROOTS ***  **FIX ME***
  //       if( eigValue>49. )  // for some reason we are missing double roots above about 49, but not all! fix me
  //         ev.push_back(eigValue); // double root
  //     }
  // }
