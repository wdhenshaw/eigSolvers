// This file automatically generated from computeEigenvalues.bC with bpp.
// ==================================================================================
// Compute the eigenvalues and eigenvectors using SLEPc
// ==================================================================================

#include "mpi.h"
#include "Overture.h"
#include "ParallelUtility.h"
#include "display.h"
#include "CompositeGridOperators.h"
#include "SparseRep.h" 

#include "Ogev.h"
// put last to avoid conflicts with "real"
#include <slepceps.h>

#define FOR3N(i1,i2,i3,n,n1a,n1b,n2a,n2b,n3a,n3b)       for( i3=n3a; i3<=n3b; i3++ )                        for( i2=n2a; i2<=n2b; i2++ )                      for( i1=n1a; i1<=n1b; i1++ )                    for( n=0; n<numberOfComponents; n++ )

#define ForBoundary(side,axis)   for( int axis=0; axis<numberOfDimensions; axis++ ) for( int side=0; side<=1; side++ )  

// ============================================================================
// Compute the global matrix index ig from the grid-function index (i1,i2,i3,n)
// ===========================================================================

// =================================================================================
// Compute the local grid-function index (i1,i2,i3,n) from the global matrix index ig
// ================================================================================
// #beginMacro getLocalIndex( ig,i1,i2,i3,n )
//    i3 = 0;
//    n = ig/(nd1a*nd2a); 
//    i2 = n2a + (ig- nd1a*nd2a*(n) )/nd1a;
//    i1 = ig + n1a - nd1a*( (i2)-n2a + nd2a*(n) );
// #endMacro



// ==================================================================================
/// \brief Compute the eigenvalues using SLEPc
///
/// \param numberOfEigenvalues (input) :   
/// \param numberOfEigenvectors(input) :   
/// \param eig (output) : holds eigenvalues
/// \param ucg (output) : holds eigevectors.
///
// ==================================================================================
int Ogev::
computeEigenvalues( const aString & problem, const int numberOfComponents,
                                        int orderOfAccuracy, int & numEigenValues, int & numEigenVectors, 
                                        RealArray & eig, CompositeGrid & cg, realCompositeGridFunction & ucg, 
                                        CompositeGridOperators & cgop, 
                                        Real tol, int eigOption, int maxIterations,
                                        const int setTargetEigenvalue, const Real targetEigenvalue,
                                        IntegerArray & bc, int numGhost, int saveMatlab, int useWideStencils,
                                        int maximumProjectedDimension /* = -1 */  )
{
    Real cpu0=getCPU();

    Real & shift = dbase.get<Real>("shift");

    const Real lambdaShift = setTargetEigenvalue ? targetEigenvalue : 0.;
    const int numberOfDimensions = cg.numberOfDimensions();

  // --- Determine if we have radiation boundary conditions ---
    int & complexProblem = dbase.get<int>("complexProblem");

    bool isSingular=true; // true if the matrix A is singular (e.g. all Neumann/Radiation BCs)
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        ForBoundary(side,axis)
        {
            if( mg.boundaryCondition(side,axis)==characteristic )
            {
                complexProblem=true;
                break;
            }
            else if(  mg.boundaryCondition(side,axis)==dirichlet )
            {
                isSingular=false;
            }
        }
    }

    if( isSingular )
    {
        if( shift==0.0 )
        {
            shift =1.; 
            printF("++++ computeEigenvalues: INFO the boundary conditions lead to a singular matrix A,\n"
                          "++++ a shift will be added to avoid this, A -> A + shift*B. Setting shift=%g\n",shift);
        }
    }

    printF(" === computeEigenvalues: problem=%s, orderOfAccuracy=%d, numberOfComponents=%d isSingular=%d shift=%g complexProblem=%d, eigOption=%d (0: Ax=lam*B*x, 1:Bx=(1/lam)A*x)=====\n",
              (const char*)problem,orderOfAccuracy, numberOfComponents,(int)isSingular,shift,(int)complexProblem, eigOption);


    const int myid=max(0,Communication_Manager::My_Process_Number);

    dbase.get<int>("orderOfAccuracy")=orderOfAccuracy; 

  // realMappedGridFunction & u = ucg[0]; // do this for now 

    Mat            A,B;             /* matrices */
    EPS            eps;             /* eigenproblem solver context */
    EPSType        type;
    Vec            xr,xi,*Iv,*Cv;
    PetscInt       nev,maxit,i,its,lits,nconv,nini=0,ncon=0;
    char           filename[PETSC_MAX_PATH_LEN];
    PetscViewer    viewer;
    PetscBool      flg,evecs,ishermitian;
    PetscErrorCode ierr;

  // SlepcInitialize(&argc,&argv,(char*)0,help);


    PetscBool      flag;
    int j,N,m,n,Istart,Iend,II;

  // ierr = PetscOptionsGetInt(NULL,"-n",&n,NULL);CHKERRQ(ierr);
  // ierr = PetscOptionsGetInt(NULL,"-m",&m,&flag);CHKERRQ(ierr);
  // if (!flag) m=n;


  // CompositeGrid & cg = *ucg.getCompositeGrid();

  // --- build info needed for global indexing ----
    buildGlobalIndexing( cg );
  

  // // ---- Count the total number of grid points -----
  // int totalNumberOfGridPoints=0; 
  // for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  // {
  //   MappedGrid & mg = cg[grid];
  //   const IntegerArray & gid = mg.gridIndexRange();
  //   const int nd1a = gid(1,0)-gid(0,0)+1 + 2*numGhost;
  //   const int nd2a = gid(1,1)-gid(0,1)+1 + 2*numGhost;
  //   const int n1a = gid(0,0)-numGhost;
  //   const int n2a = gid(0,1)-numGhost;
  //   const int n1b = gid(1,0)+numGhost;
  //   const int n2b = gid(1,1)+numGhost;  

  //   totalNumberOfGridPoints += nd1a*nd2a*numberOfComponents; // total number of grid points     


  // }
    printF("computeEigenvalues: numberOfComponentGrids=%d, numberOfGridPoints=%d\n",
                  cg.numberOfComponentGrids(),numberOfGridPoints);

  // MappedGrid & mg = *ucg.getMappedGrid();
  // const bool isRectangular = mg.isRectangular();
  // real dx[3]={1.,1.,1.};
  // if( isRectangular )
  //   mg.getDeltaX(dx);

  // assert( isRectangular );

  // const IntegerArray & gid = mg.gridIndexRange();
  // n = gid(1,0)-gid(0,0)+1 - 2;
  // m = gid(1,1)-gid(0,1)+1 - 2;

  // ----- Number of ghost points ----
  // const int numGhost= orderOfAccuracy/2;


  // const int nd1a = gid(1,0)-gid(0,0)+1 + 2*numGhost;
  // const int nd2a = gid(1,1)-gid(0,1)+1 + 2*numGhost;
  // const int n1a = gid(0,0)-numGhost;
  // const int n2a = gid(0,1)-numGhost;
  // const int n1b = gid(1,0)+numGhost;
  // const int n2b = gid(1,1)+numGhost;  

    int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2]; 
    int ige;


    bool useNew=true;


  // N = n*m;
  // N = nd1a*nd2a*numberOfComponents; // total number of grid points

    N = numberOfGridPoints*numberOfComponents; // total number of grid points

    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n Generalized Eigenvalue Problem Ax = k Bx, N=%D\n\n",N);CHKERRQ(ierr);



   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   //      Load the matrices that define the eigensystem, Ax=kBx
   //   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//   ierr = PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized eigenproblem stored in file.\n\n");CHKERRQ(ierr);
//   ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-f1",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
//   if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate a file name for matrix A with the -f1 option");

// #if defined(PETSC_USE_COMPLEX)
//   ierr = PetscPrintf(PETSC_COMM_WORLD," Reading COMPLEX matrices from binary files...\n");CHKERRQ(ierr);
// #else
//   ierr = PetscPrintf(PETSC_COMM_WORLD," Reading REAL matrices from binary files...\n");CHKERRQ(ierr);
// #endif
//   ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
//   ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
//   ierr = MatSetFromOptions(A);CHKERRQ(ierr);
//   ierr = MatLoad(A,viewer);CHKERRQ(ierr);
//   ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

//   ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-f2",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
//   if (flg) {
//     ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
//     ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
//     ierr = MatSetFromOptions(B);CHKERRQ(ierr);
//     ierr = MatLoad(B,viewer);CHKERRQ(ierr);
//     ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
//   } else {
//     ierr = PetscPrintf(PETSC_COMM_WORLD," Matrix B was not provided, setting B=I\n\n");CHKERRQ(ierr);
//     B = NULL;
//   }


  // ---- Fill Matrix A and B for the minus Laplacian ----
    if( problem=="laplace" )
    {
        if( complexProblem )
        {
      // radiation BCs --> leads to a quadratic eigenvalue problem
            fillMatrixLaplacianComplex( orderOfAccuracy, ucg, cgop, A, B, numGhost, useNew, tol, eigOption,bc, saveMatlab, lambdaShift );
        }
        else
        {
            fillMatrixLaplacian( orderOfAccuracy, ucg, cgop, A, B, numGhost, useNew, tol, eigOption,bc, saveMatlab, lambdaShift );
        }
    }
    else if( problem=="ile" )
    {
        fillMatrixIncompressibleElasticity( numberOfComponents, orderOfAccuracy, cg[0], A, B, numGhost, useNew, tol, 
                                                                                eigOption,bc,saveMatlab, useWideStencils );
    }
    else
    {
        OV_ABORT("ERROR: unknown problem");
    }

    
    if( shift!=0.0 )
    {
    // Add a shift to avoid a singular "A" that could occur with all Neumann/Radiation BCs

    // A x = lam * B * x 
    // Shift:
    //  (A + shift B) x = (lam + shift) B x 

    // str - either SAME_NONZERO_PATTERN, DIFFERENT_NONZERO_PATTERN, UNKNOWN_NONZERO_PATTERN, or SUBSET_NONZERO_PATTERN (nonzeros of X is a subset of Y’s)
        printF("SHIFT A -> A + shift*B ...\n");    
        ierr= MatAXPY(A, shift, B, DIFFERENT_NONZERO_PATTERN ); 
        printF("... done shifting A\n");  
    }



    ierr = MatGetVecs(A,NULL,&xr);CHKERRQ(ierr);
    ierr = MatGetVecs(A,NULL,&xi);CHKERRQ(ierr);

    /*
          Read user constraints if available
    */
  // ierr = PetscOptionsGetInt(NULL,"-nconstr",&ncon,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-nconstr",&ncon,&flg);CHKERRQ(ierr);
    if (flg) {
        if (ncon<=0) SETERRQ(PETSC_COMM_WORLD,1,"The number of constraints must be >0");
        ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-fconstr",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
        if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must specify the name of the file storing the constraints");
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
        ierr = VecDuplicateVecs(xr,ncon,&Cv);CHKERRQ(ierr);
        for (i=0;i<ncon;i++) {
            ierr = VecLoad(Cv[i],viewer);CHKERRQ(ierr);
        }
        ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }

    /*
          Read initial guesses if available
    */
  // ierr = PetscOptionsGetInt(NULL,"-ninitial",&nini,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-ninitial",&nini,&flg);CHKERRQ(ierr);
    if (flg) {
        if (nini<=0) SETERRQ(PETSC_COMM_WORLD,1,"The number of initial vectors must be >0");
        ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-finitial",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
        if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must specify the name of the file containing the initial vectors");
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
        ierr = VecDuplicateVecs(xr,nini,&Iv);CHKERRQ(ierr);
        for (i=0;i<nini;i++) {
            ierr = VecLoad(Iv[i],viewer);CHKERRQ(ierr);
        }
        ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                Create the eigensolver and set various options
          - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /*
          Create eigensolver context
    */
    ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);




    /*
          Set operators. In this case, it is a generalized eigenvalue problem
    */
    if( eigOption==0 )
    {
     // Solve A x = k B x
        ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
    }
    else
    {
     // Solve B x = k A x
        ierr = EPSSetOperators(eps,B,A);CHKERRQ(ierr);
    }

    /*
          If the user provided initial guesses or constraints, pass them here
    */
    ierr = EPSSetInitialSpace(eps,nini,Iv);CHKERRQ(ierr);
    ierr = EPSSetDeflationSpace(eps,ncon,Cv);CHKERRQ(ierr);



  // ---- choose eigenvalues to find ----


    if( 1==0 &&                       // DID NOT SEEM TO WORK ???  TURN OFF ...
            setTargetEigenvalue == 1 )
    {
    // Set a Target eigenvalue 
        PetscScalar target = eigOption==0 ? targetEigenvalue : 1/targetEigenvalue;
        printF("\n @@@@@@@@@ computeEigenvalues: targetEigenvalue=%e, SLEPC target=%e (for SLEPSc, target=1/targetEig for eigOption=1) @@@@@@@@@@\n\n",
                  targetEigenvalue,target );

    // EPSGetST(eps,&st);  
    // STSetType(st,STSINVERT); 
    // STSetShift( st,target );

        ierr = EPSSetTarget(eps,target ); CHKERRQ(ierr);
        ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE); CHKERRQ(ierr);

    
  
    }  
    else if( eigOption==0 )
    {
        ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE); CHKERRQ(ierr);
    }
    else
    {
        ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE); CHKERRQ(ierr);
    }

    
  // set eigenSolver type
    const aString & eigenSolver = dbase.get<aString>("eigenSolver");
    if( eigenSolver=="KrylovSchur" )
    {
        printF("computeEigenvalues: set SLEPSc eigenSolver: ESPSetType=[%s]\n",(const char*)eigenSolver);
        EPSSetType( eps, EPSKRYLOVSCHUR );
    }
    else if( eigenSolver=="ARPACK" )
    {
        printF("computeEigenvalues: set SLEPSc eigenSolver: ESPSetType=[%s]\n",(const char*)eigenSolver);
        EPSSetType( eps, EPSARPACK );
    }
    else if( eigenSolver=="LAPACK" )
    {
        printF("computeEigenvalues: set SLEPSc eigenSolver: ESPSetType=[%s]\n",(const char*)eigenSolver);
        EPSSetType( eps, EPSLAPACK );
    }
    else
    {
        printF("computeEigenValues:ERROR: unknown eigenSolver=[%s]\n",(const char*)eigenSolver );
        OV_ABORT("ERROR");
    }

  // PetscInt maxIt = 20000; 
    ierr = EPSSetTolerances(eps,tol,maxIterations); CHKERRQ(ierr);

    /*
          Set solver parameters at runtime
    */
    ierr = EPSSetFromOptions(eps);CHKERRQ(ierr); // *wdh* added Aug 5, 2024

    PetscInt mpd = PETSC_DEFAULT; // numEigenVectors; // maximum projected dimension, decrease to save space
    if( maximumProjectedDimension>0 )
        mpd = maximumProjectedDimension;

    printF("Setting numEigenValues=%d, numEigenVectors=%d, maximumProjectedDimension=%d, maxIterations=%d\n",
                numEigenValues,numEigenVectors,maximumProjectedDimension,maxIterations );
    PetscInt ncv = PETSC_DEFAULT; // numEigenValues*2+1; // size of column space 
    ierr = EPSSetDimensions(eps,numEigenValues,ncv,mpd); CHKERRQ(ierr);

  // ----------- get info about ksp object ----------
    ST st; // spectral transform
    EPSGetST(eps,&st); 
    KSP ksp;
    STGetKSP(st,&ksp); 
    aString name="";
    const int maxLen=100;
    if( ksp !=NULL ) // *new* way June 20, 2017
    {
        KSPType type;
        ierr=KSPGetType(ksp,&type);CHKERRQ(ierr); 
        name = name + "ksp[" + type + "]";
    }
    
    PC pc;
    ierr=KSPGetPC(ksp,&pc); CHKERRQ(ierr);  
    if( pc!=NULL )
    {
        PCType type;
        PCGetType(pc,&type);
        if( type == PCILU )
        {
      // trouble here if -st_pc_type bjacobi
            PetscInt levels=0;
      // printF("computeEigenvalues: call PCFactorGetLevels...\n"); 
            ierr = PCFactorGetLevels(pc, &levels); CHKERRQ(ierr);
      // printF("computeEigenvalues: ... done call PCFactorGetLevels\n"); 

            name = name + " pc[" + type + sPrintF("(%d)",levels) + "]";
        }
        else
        {
              name = name + " pc[" + type + "]";
        }
    } 
    printF("****** computeEigenvalues: KSP name=%s\n",(const char *)name);


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                            Solve the eigensystem
          - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    ierr = EPSSolve(eps);     CHKERRQ(ierr);



    Real cpuTotal = getCPU()-cpu0;
    printF("--------------------------------------------------------------\n");
    printF("--------------- Ogev : computeEigenvalues --------------------\n");
    printF("   numberOfGridPoints=%d, numEigenvalues=%d, numEigenVectors=%d\n",
                    numberOfGridPoints,numEigenValues,numEigenVectors);
    printF("      total cpu = %9.2e (s) (includes build matrix)\n",cpuTotal);
    printF(" cpu/eigenvalue = %9.2e (s)\n",cpuTotal/numEigenValues);


    /*
          Optional: Get some information from the solver and display it
    */
    ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);CHKERRQ(ierr);
  // old: ierr = EPSGetOperationCounters(eps,NULL,NULL,&lits);CHKERRQ(ierr);
  // new for v18.2: 
    ierr = EPSGetIterationNumber(eps,&lits);CHKERRQ(ierr);  
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %d\n",lits);CHKERRQ(ierr);
    ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);

    ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %d\n",nev);CHKERRQ(ierr);

    PetscCall(EPSGetConverged(eps,&nconv));
    printF(">>Number of converged eigenvalues= %d\n",nconv);

  // nev = numEigenValues;
    nev = max(numEigenValues,nconv);    // keep all converged
  // nev = min(numEigenValues,nconv);
    numEigenValues = nev;
    numEigenVectors= nev;
    printF(">> Setting numEigenValues=numEigenVectors=%d\n",numEigenValues);

    ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",tol,maxit);CHKERRQ(ierr);

        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        Display solution and clean up
          - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    ierr = EPSPrintSolution(eps,NULL);CHKERRQ(ierr);

    printF("--------------------------------------------------------------\n");
    
    /*
          Save eigenvectors, if requested
    */
  // ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-evecs",filename,PETSC_MAX_PATH_LEN,&evecs);CHKERRQ(ierr);

  // ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
    
  // printF("Number of converged=%d\n",nconv);

    const int numComp = complexProblem ? 2*numberOfComponents : numberOfComponents;

    if( nconv>0 && nev>0 ) 
    {
        Range all;
        ucg.updateToMatchGrid(cg,all,all,all,numComp*numEigenVectors);
        ucg.setName("phi");                 // give names to grid function ...
        for( int i=0, k=0; i<numEigenVectors; i++ )
        {
            if( numberOfComponents==1 )
            {
                if( complexProblem )
                {
                    ucg.setName(sPrintF("psir%d",i), k); k++;  // real part of psi
                    ucg.setName(sPrintF("psii%d",i), k); k++;  // imag part of psi
                }
                else
                    ucg.setName(sPrintF("psi%d",i), i);
            }
            else
            {
                assert( numberOfComponents==2 );
                ucg.setName(sPrintF("psi%d",i), k); k++;  // first component is psi 
                ucg.setName(sPrintF("p%d",i), k); k++;  // 2nd component is p 
            }
        }

    // ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    // ierr = EPSIsHermitian(eps,&ishermitian);CHKERRQ(ierr);
    // for (i=0;i<nconv;i++) 
        eig.redim(2,numEigenVectors);

    // complexEigFound = 0  : this is not a complex eig
    //                 = 1  : first vector in a complex eig
    //                 = 2  : second vector in a complex eig
        int complexEigFound=0; 
        for( int i=0; i<numEigenVectors; i++ ) 
        {

            PetscScalar kr, ki;
      // printF("Get eigenpair i=%d\n",i);

      // --- Get (complex) eigenvalue and complex eigenvector
            ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi); CHKERRQ(ierr);


            if( ierr!=0 )
                OV_ABORT("Error from EPSGetEigenpair");


            if( eigOption==0 )
            {
                eig(0,i) = kr + lambdaShift - shift;
                eig(1,i) = ki;
            }
            else
            {
        // --- Here we have computed the inverse of the eigenvalues we want ----
        // 1/z = zBar/ |z|^2 
                eig(0,i) =  kr/( kr*kr + ki*ki ) + lambdaShift - shift;
                eig(1,i) = -ki/( kr*kr + ki*ki );
            }

      // printF("Eigenvalue %d : k=%18.14e + %18.14e I \n",i,kr,ki);

      // ierr = EPSGetEigenvector(eps,i,xr,xi);CHKERRQ(ierr);

      // ierr = VecView(xr,viewer);CHKERRQ(ierr);
      // if (!ishermitian) { ierr = VecView(xi,viewer);CHKERRQ(ierr); }

            if( i<numEigenVectors )
            {
        // ---- Save the eigenvector ----
        // printF("Save eigenvector %d to the grid function.\n",i);

                PetscScalar *xrv, *xiv;
                VecGetArray(xr,&xrv);  // get the local array from Petsc
                ierr = VecGetOwnershipRange(xr,&Istart,&Iend);CHKERRQ(ierr);
                if( ierr!=0 )
                    OV_ABORT("PETSC ERROR");

        // If the imaginary part of k is non-zero then
        // the eigenvectors appear as complex conjugates
        //     x = xr + I * xi 
        //     x = xr - I * xi 
                if( ki!=0 || complexProblem )
                {
                    if( complexEigFound==2 )
                      complexEigFound=0;  // reset -- this must be a new complex eig 
                    
                    complexEigFound++;
                    assert( complexEigFound<=2 );

                    VecGetArray(xi,&xiv); // get 
                }
                else
                {
                    complexEigFound=0; 
                }

                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {  
                    realArray & ug= ucg[grid];
                    realSerialArray uLocal; getLocalArrayWithGhostBoundaries(ug,uLocal); 

                    int n1a = uLocal.getBase(0) +ug.getGhostBoundaryWidth(0), 
                            n1b = uLocal.getBound(0)-ug.getGhostBoundaryWidth(0);
                    int n2a = uLocal.getBase(1) +ug.getGhostBoundaryWidth(1), 
                            n2b = uLocal.getBound(1)-ug.getGhostBoundaryWidth(1);
                    int n3a = uLocal.getBase(2) +ug.getGhostBoundaryWidth(2), 
                            n3b = uLocal.getBound(2)-ug.getGhostBoundaryWidth(2);
                        
                    if( false )
                        printf("Ogev: myid=%i local array bounds = [%i,%i][%i,%i][%i,%i]\n",myid,n1a,n1b,n2a,n2b,n3a,n3b);

                    i1=n1a, i2=n2a, i3=n3a;
                    int n=0;
                    int ig=getGlobalIndex( n, iv, grid, myid );  // get the global index for the first point

                    if( complexProblem==0 )
                    {

                        const int nb=uLocal.getBase(3);
                        const int nComp = nb + n + numberOfComponents*(i); // fill in this component
                        FOR3N(i1,i2,i3,n,n1a,n1b,n2a,n2b,n3a,n3b)
                        {

              // ******** NOTE we can probably just increment ig by 1 if we start correctly
              // int ig=getGlobalIndex( iv, grid, myid );  // get the global index

                            if( ig>=Istart && ig<=Iend )
                            {
                // if( false ) printf(" myid=%i: i1,i2=%i,%i, ig=%i xrv[ig]=%6.4f\n",myid,i1,i2,ig,xrv[ig-Istart]);
                                if( complexEigFound<=1 )
                                    uLocal(i1,i2,i3,nComp)=xrv[ig-Istart];  // use real part of eigenvector 
                                else
                                    uLocal(i1,i2,i3,nComp)=xiv[ig-Istart];  // use imaginary part of eigenvector 
                            }
                            else
                            {
                                int p=myid;
                                printf("Ogev::ERROR: myid=%i, i1,i2=%i,%i, ig=%i Istart,Iend=[%i,%i]\n", myid,i1,i2,ig,Istart,Iend);
                            }
                            ig++;
                        }
                    }
                    else
                    {
            //  --- complex case: fill in Real and Imag parts ---

                        const int nb=uLocal.getBase(3);
                        const int nComp = nb + n + numComp*(i); 
                        FOR3N(i1,i2,i3,n,n1a,n1b,n2a,n2b,n3a,n3b)
                        {
                            if( ig>=Istart && ig<=Iend )
                            {
                                uLocal(i1,i2,i3,nComp  )=xrv[ig-Istart];  // real part of eigenvector       Re(phi)
                                uLocal(i1,i2,i3,nComp+1)=xiv[ig-Istart];  // imaginary part of eigenvector  Im(phi) 
                            }
                            else
                            {
                                int p=myid;
                                printf("Ogev::ERROR: myid=%i, i1,i2=%i,%i, ig=%i Istart,Iend=[%i,%i]\n", myid,i1,i2,ig,Istart,Iend);
                            }
                            ig++;
                        }


                    }
                    
                }
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {  
                    ucg[grid].periodicUpdate();  // added for periodic box April 14, 2023
                    ucg[grid].updateGhostBoundaries();
                    if( debug & 8 )
                        display(ucg[grid],sPrintF("Eigenvectors: ucg[%i]",grid),"%6.3f ");
                }          

            }
        }
    // ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }

    /*
          Free work space
    */
    ierr = EPSDestroy(&eps);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = MatDestroy(&B);CHKERRQ(ierr);
    ierr = VecDestroy(&xr);CHKERRQ(ierr);
    ierr = VecDestroy(&xi);CHKERRQ(ierr);
    if (nini > 0) {
        ierr = VecDestroyVecs(nini,&Iv);CHKERRQ(ierr);
    }
    if (ncon > 0) {
        ierr = VecDestroyVecs(ncon,&Cv);CHKERRQ(ierr);
    }
  // ierr = SlepcFinalize();





    return 0;

}// end compute Eigenvalues 

