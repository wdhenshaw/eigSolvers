// This file automatically generated from genEigs.bC with bpp.
// ==================================================================================================
// Solve a generalized eigenvalue problem with SLEPc
//            A x = k B x
// 
// Started from SLEPc examples:
//      src/eps/examples/tutorials/ex2.c  : Laplace
//      src/eps/examples/tutorials/ex7.c  : generalized eig problem, read matrices from files
// ==================================================================================================



static char help[] = "Compute some  eigenvalues\n";


// **TESTING : March 19, 2023 

#include "mpi.h"
#include "Overture.h"

#include "display.h"
#include "PlotStuff.h"  
#include "SquareMapping.h" 
#include "Ogshow.h"

#include "CompositeGridOperators.h"
#include "SparseRep.h" 
#include "Oges.h"
#include "ParallelUtility.h"
#include "gridFunctionNorms.h"

#include "Integrate.h"


#include "HDF_DataBase.h"

// Put this last
#include "Ogev.h"



// lapack routines
#ifdef OV_USE_DOUBLE
    #define GETRF EXTERN_C_NAME(dgetrf)
    #define GETRI EXTERN_C_NAME(dgetri)
    #define GETRS EXTERN_C_NAME(dgetrs)
    #define GECON EXTERN_C_NAME(dgecon)
    #define LANGE EXTERN_C_NAME(dlange)
  // #define GEEV  EXTERN_C_NAME(dgeev)
#else
    #define GETRF EXTERN_C_NAME(sgetrf)
    #define GETRI EXTERN_C_NAME(sgetri)
    #define GETRS EXTERN_C_NAME(sgetrs)
    #define GECON EXTERN_C_NAME(sgecon)
    #define LANGE EXTERN_C_NAME(slange)
  // #define GEEV  EXTERN_C_NAME(sgeev)
#endif

extern "C"
{
  // PA = LU factor
    void GETRF( int & m, int & n, Real & a, const int & lda, int & ipvt, int & info );
  // Solve given LU
    void GETRS( char *trans, int & n, int & nhrs, Real & a, const int & lda, int & ipvt, Real & b, const int & ldb, int & info );

  // compute inverse:
    void GETRI( int & n, Real & a, const int & lda, const int & ipvt, Real & work, const int & iwork, int & info );

    void GECON( char *norm, int & n, Real & a, const int & lda, Real & anorm, Real & rcond, Real & work, int & iwork, int & info );
    Real LANGE( char *norm, int & m, int & n, Real & a, const int & lda, Real & work );

  // void GEEV( char *jobvl, char* jobvr, int & n, Real & a, const int & lda,
  //             Real & wr, Real & wi, Real &vl, int & ldvl, Real & vr, int & ldvr, Real & work, int & lwork, int & info );

}

// #include "CgSolverUtil.h"

// Boundary conditions:
const int periodic=-1, interpolation=0, displacement=1, traction=2, dirichlet=1, neumann=2;

#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(int i3=I3Base; i3<=I3Bound; i3++) for(int i2=I2Base; i2<=I2Bound; i2++) for(int i1=I1Base; i1<=I1Bound; i1++)  

#define FOR_3(i1,i2,i3,I1,I2,I3) I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(int i3=I3Base; i3<=I3Bound; i3++) for(int i2=I2Base; i2<=I2Bound; i2++) for(int i1=I1Base; i1<=I1Bound; i1++)  

#define ForBoundary(side,axis)   for( int axis=0; axis<numberOfDimensions; axis++ ) for( int side=0; side<=1; side++ )  

// // ============================================================================
// // Compute the global matrix index ig from the grid-function index (i1,i2,i3,n)
// // ===========================================================================
// #beginMacro getGlobalIndex(i1,i2,i3,n,ig)
//   ig = (i1)-n1a + nd1a*( (i2)-n2a + nd2a*(n) );
//   assert( ig>=0 && ig<N );
// #endMacro

// // =================================================================================
// // Compute the local grid-function index (i1,i2,i3,n) from the global matrix index ig
// // ================================================================================
// #beginMacro getLocalIndex( ig,i1,i2,i3,n )
//    i3 = 0;
//    n = ig/(nd1a*nd2a); 
//    i2 = n2a + (ig- nd1a*nd2a*(n) )/nd1a;
//    i1 = ig + n1a - nd1a*( (i2)-n2a + nd2a*(n) );
// #endMacro

// #beginMacro printMatrixEntry(ig0,ig,val,label)
//   if( printEntries )
//   {
//     printF(" ig0=%4d ig=%4d value=%10.3f (label)\n",ig0,ig,val);
//   }
// #endMacro







// =================================================================
/// \brief Return the name of the boundary condition.
// =================================================================
aString bcName( int bc )
{
    if( bc==-1 )
        return "p";
    else if( bc==0 )
        return "i"; // interp
    else if( bc==1 )
        return "d"; // displacement or Dirichlet 
    else if( bc==2 )
        return "n"; // traction or Neumann
    else
        return "u"; // unknown

}



// =========================================================
// Macro: set a BC from a command arg argument
// =========================================================

// =========================================================
// Macro: Assign a bcNumber, later used to set BCs
// =========================================================


// // ======================================================================
// // Macro: Compute the inner product of u(.,.,.,ii) and v(.,.,.,jj) 
// //   w = temp grid function to hold product 
// //   
// //   **DUPLICATED and changed FROM eveSolver.bC **** FIX ME 
// // ======================================================================
// #beginMacro innerProductor(ii,jj,u,v,w,dotProduct)   
//   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
//   {
//     MappedGrid & mg = cg[grid];
//     getIndex(mg.dimension(),I1,I2,I3);
//     OV_GET_SERIAL_ARRAY(Real,u[grid],uuLocal);
//     OV_GET_SERIAL_ARRAY(Real,v[grid],vvLocal);
//     OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);

//     wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ii)*vvLocal(I1,I2,I3,jj);
//   }

//   dotProduct = integrate.volumeIntegral(w);
// #endMacro


// --------------------------------------------------------------------------------------
//   Macro: return the index's for possible active points
//            
//  NOTE: This macro appears in solveSLEPc.bC and eigenModes.bC 
// --------------------------------------------------------------------------------------
    

// ===============================================================
// Macro: Compute the inner product of two grid functions          
// ===============================================================

// =================================================================================
// Macro: Compute the errors between the true (continuous) and discrete eigenvectors
// =================================================================================


// =================================================================================
// Macro: Compute multiplicities for the true eigenvalues
// =================================================================================


// ==================================================================================
// ================================ MAIN ============================================
// ==================================================================================
int main(int argc,char **argv)
{

    Overture::start(argc,argv);  // initialize Overture

  // This macro will initialize the PETSc solver if OVERTURE_USE_PETSC is defined.
    INIT_PETSC_SOLVER();

    SlepcInitialize(&argc,&argv,(char*)0,help);  

    aString commandFileName="";

    aString problem = "laplace"; 
    aString eigCase = "square"; // for getting true eigenvaliues
    bool plotOption=true;
    aString nameOfShowFile=""; 
    int flushFrequency = 100; // number of solutions per show sub-file
    aString nameOfGridFile="square32.order2.hdf";

    aString tableFileName = "genEigsTable"; // name of LaTeX file with table of results 


    bool loadBalance=true; // false;
    int numberOfParallelGhost=2;


    int debug=0;
    int numPoints=101;
    int nx=52, ny=52; 
    int numEigenValues  = 1;            // number of eigenvalues to compute 
    int numEigenVectors = -1;           // number of eigenvectors to save
    int maximumProjectedDimension = -1; // reduce storage when solving for many eigs (default=ncv)
    int discreteEigenvalues=false;      // if true, use exact discrete eigs' instead of exact continuous eigs (if available)

    int orthogonalize=true;            // if true, count multiplicities, orthogonalize eigenvectors, choose a basis of EVs for multiple eigs *new* Sept 3, 2024

    bool useAccurateInnerProduct=false; // if true, use an accurate inner product
  // int orderOfAccuracy = 2 ;
    Real tol=1.e-8; 

    Real rho=1., mu=1.;

    int includePressure=0;  // solve for pressure as well as psi

    int saveMatlab=0;       // 1 = save A and B matrices to matlab format

    Real lx=1., ly=1., lz=1.0;      // Bounds on the square domain

    IntegerArray bc(2,3);
    bc = displacement; 

  // eigOption = 0 : =solve Ax =kBx for smallest k, 
  //           = 1 : solve Bx = k A x for largest k 

    int eigOption = 1;

    Real eigSign = 1.;  // use to change sign of the eigenfunction
    int maxIterations = 200000; // for PETSc solves to invert the shifted matrix 

    int setTargetEigenvalue=0;  // 1 = look for eigenvalues near the target: targetEigenvalue
    Real targetEigenvalue=0; 
    int useWideStencils = 0; // use wide stencils when they fit.
    
    int eigc=0; // which eigen-vector to choose **make this an option***

    Real xa=0., xb=1., ya=0., yb=1.;

    const int numberOfBoundaryConditions=200;
    IntegerArray bcNumber(numberOfBoundaryConditions+1); 
    bcNumber = dirichlet;   // default boundary condition

    aString matlabFileName = "genEigs"; 

    if( argc<=1 )
    {
          printF(" ------------------ genEigs: Compute eigenvalues and eigenvectors --------------------------\n"
          "Usage:\n"
          "  genEigs [-noplot] eigs.cmd -problem=<s> -eigCase=<s> -g=<s> -numEigenValues=<i> -tol=<f> \n"
          "          -bc[123456]=[d|n] -show=<s> -matlab=<s> -table=<s> -orthogonalize=[0|1] \n"
          "          -discreteEigenValues=[0|1] -go=<s>\n"
          "  \n"
          "   -noplot : run without graphics\n"
          "   -problem=[laplace|ile] : laplace = (negative) laplacian\n"
          "                           :ile = incompressible elasticity\n"
          "   -eigCase=[square|disk|sphere] : compare answers to the known eigenpairs for a square (or box),\n"
          "         disk (or annulus or cylinder) or sphere.\n"
          "   -g=gridName : name of overset grid generated by ogen, example -g=square32.order2.hdf\n"
          "   -numEigenValues=i : number of eigenvalues to compute (and eigenvectors)\n"
          "   -tol=f : tolerance for the eigenvalues, e.g. -tol=1e-12.\n"
          "   -bc[1234567]=[d|n] : set boundary condition to Dirichlet or Neumann on a boundary with bc flag 1,2,3,4,6.\n"
          "               e.g. -bc1=d sets boundaries with bc flag=1 to Dirichlet.\n"
          "   -orthogonalize=[0|1] : 1=orthogonalize and eigenvectors corresponding to multiple eigenvalues.\n"
          "   -discreteEigenValues=[0|1] : 1=compare to true discrete eigenvalues (square or box only)\n"
          "   -matlab=matlabFileName : name of Matlab output file holding eigenvalues and errors,\n"
          "                            e.g. -matlab=diskG4Order2\n"
          "   -table=nameOfTableFile : name of file holding a LaTeX table of results, \n"
          "                            e.g. -table=square32Order4Table.\n"
          "   -show=showFileName : save results to a show file with this name. e.g. -show=myShowFile.show \n"
          "                        (use plotStuff to display results from the show file).\n" 
          "   -go=[go|og] : -go=go : run and exit. -go=og (open graphics) : when running with -noplot,\n" 
          "         open graphics windows after commands have been read.   \n" 
          "---------------------------------------------------------------------------------------------------------\n");  
    }
    
    int len=0;
    if( argc > 1 )
    { 
        for( int i=1; i<argc; i++ )
        {
            aString arg = argv[i];
            if( arg=="-noplot" || arg=="noplot" )
            {
                plotOption=false;
            }
            else if( (len=arg.matches("-problem="))  )
            {
                problem = arg(len,arg.length()-1);
                printF("Setting problem=[%s]\n",(const char*)problem);
            }
            else if( (len=arg.matches("-eigCase="))  )
            {
                eigCase = arg(len,arg.length()-1);
                printF("Setting eigCase=[%s]\n",(const char*)eigCase);
            }      
            else if( (len=arg.matches("-show="))  )
            {
                nameOfShowFile = arg(len,arg.length()-1);
                printF("Setting nameOfShowFile=[%s]\n",(const char*)nameOfShowFile);
            }      
            else if( arg(0,6)=="-debug=" )
            {
                sScanF(arg(7,arg.length()-1),"%i",&debug);
                printF("Setting debug=%i\n",debug);
            }
            else if( (len=arg.matches("-matlab="))  )
            {
                matlabFileName = arg(len,arg.length()-1);
                printF("Setting matlabFileName=[%s.m]\n",(const char*)matlabFileName);
            } 
            else if( (len=arg.matches("-table="))  )
            {
                tableFileName = arg(len,arg.length()-1);
                printF("Setting tableFileName=[%s.tex] (name of LateX table output file)\n",(const char*)tableFileName);
            }            
            else if( (len=arg.matches("-includePressure="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&includePressure);
                printF("Setting includePressure=%d\n",includePressure);
            }  
            else if( (len=arg.matches("-discreteEigenvalues="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&discreteEigenvalues);
                printF("Setting discreteEigenvalues=%d\n",discreteEigenvalues);
            }  
            else if( (len=arg.matches("-useAccurateInnerProduct="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&useAccurateInnerProduct);
                printF("Setting useAccurateInnerProduct=%d\n",useAccurateInnerProduct);
            } 
            else if( (len=arg.matches("-orthogonalize="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&orthogonalize);
                printF("Setting orthogonalize=%d\n",orthogonalize);
            }

            else if( (len=arg.matches("-eigOption="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&eigOption);
                printF("Setting eigOption=%d\n",eigOption);
            } 

            else if( (len=arg.matches("-setTargetEigenvalue="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&setTargetEigenvalue);
                printF("Setting setTargetEigenvalue=%d\n",setTargetEigenvalue);
            }
            
            else if( (len=arg.matches("-targetEigenvalue="))  )
            {
                sScanF(arg(len,arg.length()-1),"%e",&targetEigenvalue);
                printF("Setting targetEigenvalue=%e\n",targetEigenvalue);
            }               

            else if( (len=arg.matches("-maxIterations="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&maxIterations);
                printF("Setting maxIterations=%d\n",maxIterations);
            }  

            else if( (len=arg.matches("-eigc="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&eigc);
                printF("Setting eigc=%d\n",eigc);
            }     
            else if( (len=arg.matches("-saveMatlab="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&saveMatlab);
                printF("Setting saveMatlab=%d\n",saveMatlab);
            }                        
            else if( (len=arg.matches("-useWideStencils="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&useWideStencils);
                printF("Setting useWideStencils=%d\n",useWideStencils);
            }

            else if( (len=arg.matches("-bc1="))  )
            {
                    aString bcChar = arg(len,arg.length()-1);
                    if( bcChar=="d" )
                    {
                        bc(0,0)=displacement;
                        bcNumber(1)=dirichlet; // new way 
                    }
                    else if( bcChar=="t" || bcChar=="n" )
                    {
                        bc(0,0)=traction;
                        bcNumber(1)=neumann; // new way 
                    }
                    else
                    {
                        printF("setBC: Unknown bc=[%s]\n",(const char*)bcChar);
                        OV_ABORT("error");
                    }
                    printF("Setting bc(0,0)=%s\n",(const char*)bcName(bc(0,0)));
            }
            else if( (len=arg.matches("-bc2="))  )
            {
                    aString bcChar = arg(len,arg.length()-1);
                    if( bcChar=="d" )
                    {
                        bc(1,0)=displacement;
                        bcNumber(2)=dirichlet; // new way 
                    }
                    else if( bcChar=="t" || bcChar=="n" )
                    {
                        bc(1,0)=traction;
                        bcNumber(2)=neumann; // new way 
                    }
                    else
                    {
                        printF("setBC: Unknown bc=[%s]\n",(const char*)bcChar);
                        OV_ABORT("error");
                    }
                    printF("Setting bc(1,0)=%s\n",(const char*)bcName(bc(1,0)));
            }
            else if( (len=arg.matches("-bc3="))  )
            {
                    aString bcChar = arg(len,arg.length()-1);
                    if( bcChar=="d" )
                    {
                        bc(0,1)=displacement;
                        bcNumber(3)=dirichlet; // new way 
                    }
                    else if( bcChar=="t" || bcChar=="n" )
                    {
                        bc(0,1)=traction;
                        bcNumber(3)=neumann; // new way 
                    }
                    else
                    {
                        printF("setBC: Unknown bc=[%s]\n",(const char*)bcChar);
                        OV_ABORT("error");
                    }
                    printF("Setting bc(0,1)=%s\n",(const char*)bcName(bc(0,1)));
            }
            else if( (len=arg.matches("-bc4="))  )
            {
                    aString bcChar = arg(len,arg.length()-1);
                    if( bcChar=="d" )
                    {
                        bc(1,1)=displacement;
                        bcNumber(4)=dirichlet; // new way 
                    }
                    else if( bcChar=="t" || bcChar=="n" )
                    {
                        bc(1,1)=traction;
                        bcNumber(4)=neumann; // new way 
                    }
                    else
                    {
                        printF("setBC: Unknown bc=[%s]\n",(const char*)bcChar);
                        OV_ABORT("error");
                    }
                    printF("Setting bc(1,1)=%s\n",(const char*)bcName(bc(1,1)));
            }
            else if( (len=arg.matches("-bc5="))  )
            {
        // **FIX ME** Could do the general case
                aString bcChar = arg(len,arg.length()-1);
                    if( bcChar=="d" )
                    {
                        bcNumber(5)=dirichlet; 
                    }
                    else if( bcChar=="t" || bcChar=="n" )
                    {
                        bcNumber(5)=neumann; 
                    }
                    else
                    {
                        printF("setBC: Unknown bc=[%s]\n",(const char*)bcChar);
                        OV_ABORT("error");
                    }
                    printF("Setting bcNumber%d = %s\n",5,(const char*)bcChar);
            }
          
            else if( (len=arg.matches("-bc6="))  )
            {
        // **FIX ME** Could do the general case
                aString bcChar = arg(len,arg.length()-1);
                    if( bcChar=="d" )
                    {
                        bcNumber(6)=dirichlet; 
                    }
                    else if( bcChar=="t" || bcChar=="n" )
                    {
                        bcNumber(6)=neumann; 
                    }
                    else
                    {
                        printF("setBC: Unknown bc=[%s]\n",(const char*)bcChar);
                        OV_ABORT("error");
                    }
                    printF("Setting bcNumber%d = %s\n",6,(const char*)bcChar);
            } 
            else if( (len=arg.matches("-bc7="))  )
            {
        // **FIX ME** Could do the general case
                aString bcChar = arg(len,arg.length()-1);
                    if( bcChar=="d" )
                    {
                        bcNumber(7)=dirichlet; 
                    }
                    else if( bcChar=="t" || bcChar=="n" )
                    {
                        bcNumber(7)=neumann; 
                    }
                    else
                    {
                        printF("setBC: Unknown bc=[%s]\n",(const char*)bcChar);
                        OV_ABORT("error");
                    }
                    printF("Setting bcNumber%d = %s\n",7,(const char*)bcChar);
            }             
            else if( (len=arg.matches("-bc100="))  )
            {
        // **FIX ME** Could do the general case
                aString bcChar = arg(len,arg.length()-1);
                    if( bcChar=="d" )
                    {
                        bcNumber(100)=dirichlet; 
                    }
                    else if( bcChar=="t" || bcChar=="n" )
                    {
                        bcNumber(100)=neumann; 
                    }
                    else
                    {
                        printF("setBC: Unknown bc=[%s]\n",(const char*)bcChar);
                        OV_ABORT("error");
                    }
                    printF("Setting bcNumber%d = %s\n",100,(const char*)bcChar);
            }            

            else if( len=arg.matches("-g=") )
            {
                nameOfGridFile=arg(len,arg.length()-1);
                printF("setting nameOfGridFile=[%s]\n",(const char*)nameOfGridFile);
            }

            else if( arg=="loadBalance" || arg=="-loadBalance" )
            {
                loadBalance=true;
            }
            else if( len=arg.matches("-numberOfParallelGhost=") )
            {
                sScanF(arg(len,arg.length()-1),"%i",&numberOfParallelGhost);
                if( numberOfParallelGhost<0 || numberOfParallelGhost>10 )
                {
                    printF("ERROR: numberOfParallelGhost=%i is no valid!\n",numberOfParallelGhost);
                    OV_ABORT("error");
                }
                printF("Setting numberOfParallelGhost=%i\n",numberOfParallelGhost);
            }       

            else if( (len=arg.matches("-nx="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&nx);
                printF("Setting nx=%d\n",nx);
            }

            else if( (len=arg.matches("-tol="))  )
            {
                sScanF(arg(len,arg.length()-1),"%e",&tol);
                printF("Setting tol=%e\n",tol);
            }

            else if( (len=arg.matches("-lx="))  )
            {
                sScanF(arg(len,arg.length()-1),"%e",&lx);
                printF("Setting lx=%e\n",lx);
            } 
            else if( (len=arg.matches("-ly="))  )
            {
                sScanF(arg(len,arg.length()-1),"%e",&ly);
                printF("Setting ly=%e\n",ly);
            }
            else if( (len=arg.matches("-lz="))  )
            {
                sScanF(arg(len,arg.length()-1),"%e",&lz);
                printF("Setting lz=%e\n",lz);
            }      

            else if( (len=arg.matches("-xb="))  )
            {
                sScanF(arg(len,arg.length()-1),"%e",&xb);
                printF("Setting xb=%e\n",xb);
            }      
            else if( (len=arg.matches("-yb="))  )
            {
                sScanF(arg(len,arg.length()-1),"%e",&yb);
                printF("Setting yb=%e\n",yb);
            }  

            else if( (len=arg.matches("-eigSign="))  )
            {
                sScanF(arg(len,arg.length()-1),"%e",&eigSign);
                printF("Setting eigSign=%e\n",eigSign);
            }

            else if( (len=arg.matches("-ny="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&ny);
                printF("Setting ny=%d\n",ny);
            }
            else if( (len=arg.matches("-numEigenValues="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&numEigenValues);
                printF("Setting numEigenValues=%d\n",numEigenValues);
            }          
            else if( (len=arg.matches("-numEigenVectors="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&numEigenVectors);
                printF("Setting numEigenVectors=%d\n",numEigenVectors);
            } 
            else if( (len=arg.matches("-mpd="))  )
            {
                sScanF(arg(len,arg.length()-1),"%i",&maximumProjectedDimension);
                printF("Setting maximumProjectedDimension=%d\n",maximumProjectedDimension);
            }               
            else if( commandFileName=="" )
            {
                commandFileName=arg;    
                printf("genEigs: reading commands from file [%s]\n",(const char*)commandFileName);
            }
        }
    }
    else
    {
        printF("Usage: `genEigs [-noplot] [-g=<gridName>] [file.cmd] [-debug=<value>] ' \n");
    }
    
    if( numEigenVectors<0 )
        numEigenVectors=numEigenValues; // May 4, 2023

    GenericGraphicsInterface & gi = *Overture::getGraphicsInterface("genEigs",plotOption,argc,argv);
    PlotStuffParameters psp;

  // By default start saving the command file called:
    aString logFile="genEigs.cmd";
    gi.saveCommandFile(logFile);
    printF("User commands are being saved in the file `%s'\n",(const char *)logFile);

    aString outputFileName="genEigs.log";
    FILE *outFile = NULL;
    

  // read from a command file if given
    if( commandFileName!="" )
    {
        printF("read command file =[%s].\n",(const char*)commandFileName);
        gi.readCommandFile(commandFileName);
    }


    CompositeGrid cg;
  // bool loadBalance=true; // turn on or off the load balancer
    getFromADataBase(cg,nameOfGridFile,loadBalance);

    int firstChar=0; 
    for( int i=nameOfGridFile.length()-1; i>=0; i-- )
    {
        if( nameOfGridFile[i]=='/' ){ firstChar=i+1; break; } // start from end, work backwards and look for a directory symbol
    }
    aString gridNameNoPrefix = nameOfGridFile(firstChar,nameOfGridFile.length()-1);

    const int numberOfDimensions= cg.numberOfDimensions();

  // --- Set boundary conditions ----


    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        const IntegerArray & bc = mg.boundaryCondition();
        ForBoundary(side,axis)
        {
            int ibc = bc(side,axis);
            if( ibc>0 && ibc<=numberOfBoundaryConditions )
            {
                bc(side,axis) = bcNumber(ibc);
                printF("Setting bc(%d,%d) = %d on grid=%d\n",side,axis,bc(side,axis),grid);
            }
        }


    }

    int minDiscretizationWidth=INT_MAX;
    Range R=cg.numberOfDimensions();
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        const IntegerArray & dw = mg.discretizationWidth();
        minDiscretizationWidth=min(minDiscretizationWidth,min(dw(R)));
    }

    int orderOfAccuracy = minDiscretizationWidth-1;

    if( orderOfAccuracy % 2 ==1 )
        orderOfAccuracy--;   // must be even


    int numGhost = orderOfAccuracy/2;

    if( (problem=="ile" &&  min(abs(bc-traction))==0) || useWideStencils )
    {
    // includePressure=1;  // traction BCs need the pressure --> no more
        numGhost++;         // and one more ghost arg
    }
    else
    {
    // numGhost++; // TEST
    }

    printF("genEigs: minDiscretizationWidth=%i, Setting orderOfAccuracy=%d, numGhost=%d.\n",
        minDiscretizationWidth,orderOfAccuracy,numGhost);

    int numberOfComponents=1;
    if( problem=="laplace" )
    {
        numberOfComponents=1;
    }
    else if( problem=="ile" )
    {
    // Incompressible elasticity
        if( includePressure )
        {
            numberOfComponents=2; 
        }
    }
    else
    {
      OV_ABORT("ERROR: unknown problem");
    }


  // Overlapping grid eigenvalue solver:
    Ogev ogev;


    cg.update(MappedGrid::THEmask | MappedGrid::THEvertex | MappedGrid::THEcenter);  


  // MappedGrid & mg = cg[0]; // *** DO THIS FOR NOW ***

    
    Range all;       
    int totalNumberOfComponents = numEigenVectors*numberOfComponents;

  // realCompositeGridFunction u(cg,all,all,all,totalNumberOfComponents);  // create a grid function

    realCompositeGridFunction u(cg);  // holds eigenvectors -- now allocated in computeEigenvalues

  
    Index I1,I2,I3; 
  // getIndex(mg.dimension(),I1,I2,I3);                        // assign I1,I2,I3 from dimension

    u=0.;
  // u(I1,I2,I3,0)=sin(Pi*mg.vertex()(I1,I2,I3,axis1))         // component 0 : sin(pi*x)*cos(pi*y)
  //              *cos(Pi*mg.vertex()(I1,I2,I3,axis2));        
  // u(I1,I2,I3,1)=cos(Pi*mg.vertex()(I1,I2,I3,axis1))         // component 1 : cos(pi*x)*sin(pi*y)
  //              *sin(Pi*mg.vertex()(I1,I2,I3,axis2));       

  

    CompositeGridOperators cgop(cg);
    cgop.setOrderOfAccuracy(orderOfAccuracy);
    MappedGridOperators & op = cgop[0];

  // u.setOperators(op);
    u.setOperators(cgop);

    const int u1c=0, u2c=1, pc=2, psic=3;
    realCompositeGridFunction uv(cg,all,all,all,numberOfDimensions+2); // ** FIX ME ***
    uv.setName("u1",u1c);  
    uv.setName("u2",u2c);  
    uv.setName("p",pc); 
    uv.setName("psi",psic); 
    uv.setOperators(cgop);
    realCompositeGridFunction p(cg,all,all,all);
    p.setName("p",0);
    p.setOperators(cgop);

    RealCompositeGridFunction evTrue; // holds true eigenvectors
    RealCompositeGridFunction w(cg,all,all,all,1); // hold temporary variables

  // Integrate integrate;  // for computing inner products for errors in eigenvectors

    
    Real time=0.; // plot curve at this time
    
  // ========== create the GUI and dialog ================
    GUIState dialog;
    dialog.setWindowTitle("Eigenvalues Code");
    dialog.setExitCommand("exit", "exit");

    aString cmds[] = {"compute",
                                        "plot eigenvectors",
                                        "plot true eigenvectors",
                    // "contour",
                                        "check",
                                        "erase",
                                        "save show file",
                                        "save matlab file",
                                        "" };
    int numberOfPushButtons=0;  // number of entries in cmds
    while( cmds[numberOfPushButtons]!="" ){numberOfPushButtons++;}; // 
    int numRows=(numberOfPushButtons+1)/2;
    dialog.setPushButtons( cmds, cmds, numRows ); 

    const int numberOfTextStrings=15;  // max number allowed
    aString textLabels[numberOfTextStrings];
    aString textStrings[numberOfTextStrings];

    int nt=0;
    textLabels[nt] = "time:";             sPrintF(textStrings[nt],"%g",time);  nt++; 
    textLabels[nt] = "numPoints:";        sPrintF(textStrings[nt],"%i",numPoints);  nt++; 
    textLabels[nt] = "show file:";        sPrintF(textStrings[nt],"%s",(const char*)nameOfShowFile);  nt++; 
    textLabels[nt] = "flush frequency:";  sPrintF(textStrings[nt],"%i",flushFrequency);  nt++; 
  // null strings terminal list
    textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
    dialog.setTextBoxes(textLabels, textLabels, textStrings);

    gi.pushGUI(dialog);
    gi.appendToTheDefaultPrompt("genEigs>");
    psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

    RealArray eig; 

    RealArray eigErr, evectErr;
    IntegerArray multiplicityTrue, multIndexTrue;
    bool eigenVectorsAreKnown=false;
    bool eigenValuesAreKnown=false;

    IntegerArray eigMultiplicity, eigStartIndex;  
  
    aString answer,buff;  
    for( ;; )
    {
        gi.getAnswer(answer,"");  
  
        if( answer=="continue" )
        {
            break;
        }
        else if( answer=="exit" || answer=="done" )
        {
            break;
        }
        else if( answer=="compute" )
        {
      // --- compute eigenvalues and eigenvectors :
      //  eig : holds computed eigenvalues
      //  u : holds computed eigenvectors
            ogev.computeEigenvalues( problem,numberOfComponents, orderOfAccuracy, numEigenValues,numEigenVectors, eig, 
                                                              cg, u, cgop, tol, eigOption, maxIterations, setTargetEigenvalue, targetEigenvalue,
                                                              bc, numGhost, saveMatlab, useWideStencils, maximumProjectedDimension );

          
            if( orthogonalize && numEigenVectors>0 )  
            {
        //  Count multiplicities and orthogonalize eigenvectors for multiple eigenvalues

                ogev.orthogonalizeEigenvectors( problem, numberOfComponents, orderOfAccuracy, numEigenValues, numEigenVectors, eig,  u, eigMultiplicity, eigStartIndex );

            }
            else if( numEigenVectors>0 )
            {
                eigMultiplicity.redim(numEigenVectors);
                eigMultiplicity=1; 
            }


      // ----- compute residuals ----
            printF("\n >>>> genEigs:: compute residuals ||  A - lambda u ||/lambda ....\n");
            realCompositeGridFunction res(cg,all,all,all,numEigenVectors);
            RealArray resMax(numEigenVectors);
            for( int ie=0; ie<numEigenVectors; ie++ )
            {
                Real lambda = eig(0,ie);
                resMax(ie) = ogev.getEigenPairResidual( lambda, u, res, cgop, ie );
        // printF(" ie=%4d  max norm resdiual ||  A - lambda u ||/lambda =%9.2e\n",resMax(ie));
            }
            printF("... done compute residuals\n");


            tableFileName = tableFileName + ".tex";

      // FILE *outFile = fopen("genEigsTable.tex","w" );     // Save some tex output here 
            FILE *outFile = fopen((const char*)tableFileName,"w" );     // Save some tex output here 

            printF("\n======================== GenEigs problem=%s =======================\n",(const char*)problem);
            printF(" grid=%s\n"
                          " orderOfAccuracy=%d, numberOfComponents=%d numGhost=%d \n",
                        (const char*)nameOfGridFile,orderOfAccuracy,numberOfComponents,numGhost);
            printF(" useAccurateInnerProduct=%d, discreteEigenvalues=%d\n",useAccurateInnerProduct,discreteEigenvalues);
            printF(" includePressure=%d, useWideStencils=%d, eigCase=%s\n",includePressure,useWideStencils,(const char*)eigCase);
            fPrintF(outFile,"%% grid=%s\n",(const char*)nameOfGridFile);
            fPrintF(outFile,"%% orderOfAccuracy=%d, \n",orderOfAccuracy);
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                const IntegerArray & bc = mg.boundaryCondition();

                if( numberOfDimensions==2 )
                {
                    printF(" bc = [%s,%s,%s,%s] for grid=%i (%s)\n",
                                  (const char*)bcName(bc(0,0)),(const char*)bcName(bc(1,0)),
                                  (const char*)bcName(bc(0,1)),(const char*)bcName(bc(1,1)),
                                  grid,(const char*)mg.getName());
                    fPrintF(outFile,
                                  "%% bc = [%s,%s,%s,%s] for grid=%i (%s)\n",
                                  (const char*)bcName(bc(0,0)),(const char*)bcName(bc(1,0)),
                                  (const char*)bcName(bc(0,1)),(const char*)bcName(bc(1,1)),
                                  grid,(const char*)mg.getName());
                }
                else
                {
                    printF(" bc = [%s,%s,%s,%s,%s,%s] for grid=%i (%s)\n",
                                  (const char*)bcName(bc(0,0)),(const char*)bcName(bc(1,0)),
                                  (const char*)bcName(bc(0,1)),(const char*)bcName(bc(1,1)),
                                  (const char*)bcName(bc(0,2)),(const char*)bcName(bc(1,2)),
                                  grid,(const char*)mg.getName());
                    fPrintF(outFile,
                                  "%% bc = [%s,%s,%s,%s,%s,%s] for grid=%i (%s)\n",
                                  (const char*)bcName(bc(0,0)),(const char*)bcName(bc(1,0)),
                                  (const char*)bcName(bc(0,1)),(const char*)bcName(bc(1,1)), 
                                  (const char*)bcName(bc(0,2)),(const char*)bcName(bc(1,2)),
                                grid,(const char*)mg.getName());          
                }        
            }
            printF("===================================================================\n");

      // --- In some cases we know the true eigenvalues ----
            int numEigs = numEigenVectors; 
            int numEigsTrue = numEigs + 20; // compute a few more in case we have multiplicityTrue > 1
            if( setTargetEigenvalue )
            {
                numEigsTrue = numEigsTrue + 1000;   // add extra to account for setting a target value
            }

            RealArray eigsTrue(numEigsTrue); 

            eigErr.redim(numEigs);   eigErr=0.;
            evectErr.redim(numEigs); evectErr=0.;

            eigsTrue=1.;

            int bcOpt=0; // Dirichlet BC's 
            if( eigCase == "square" || eigCase == "box" )
            {
                eigenVectorsAreKnown=true;
                eigenValuesAreKnown=true;
                ogev.getEigenvaluesBox( numEigsTrue, eigsTrue, cg, lx,ly,lz, &evTrue, discreteEigenvalues );
            }
            else if( eigCase == "annulus" || eigCase == "disk" )
            {
                eigenValuesAreKnown=true;
                if( eigCase == "disk" )
                    eigenVectorsAreKnown=true;

                Real ra=.5, rb=1., za=0., zb=1.;
                if( eigCase == "disk" ) 
                {
                    ra=0.;
                    if( numberOfDimensions==3 )
                        rb=.5; 
                }
                if( eigenVectorsAreKnown )
                    ogev.getEigenvaluesCylinder( numEigsTrue, eigsTrue, cg, ra,rb, za, zb, &evTrue );
                else
                    ogev.getEigenvaluesCylinder( numEigsTrue, eigsTrue, cg, ra,rb, za, zb );
            }
            else if( eigCase == "sphere"  )
            {
                eigenValuesAreKnown=true;
                eigenVectorsAreKnown=true;

                Real ra=0., rb=1.;
                if( eigenVectorsAreKnown )
                    ogev.getEigenvaluesSphere( numEigsTrue, eigsTrue, cg, ra,rb,  &evTrue );
                else
                    ogev.getEigenvaluesSphere( numEigsTrue, eigsTrue, cg, ra,rb );
            }      
            else
            {
                printF("\n **** INFO: genEigs: unknown eigCase=[%s]\n\n",(const char*)eigCase);
            }

      // ---- COMPUTE MULTIPLICITIES ----
      // Save multiplicities to the show file
      //   -> true eigs are sorted at this point 

      // multIndexTrue(ie) = first true eig with given multiplicityTrue
            multiplicityTrue.redim(numEigsTrue);
            multIndexTrue.redim(numEigsTrue);
            multiplicityTrue=1;
            if( eigenVectorsAreKnown )
            {
                    int i=0; 
                    while( i<numEigsTrue )
                    {
                          int mult=1; 
                          for( int ie=i; ie<numEigsTrue-1; ie++ )
                          {
                              if( fabs(eigsTrue(ie+1)-eigsTrue(ie)) < tol*(1.+fabs(eigsTrue(ie))) ) 
                              {
                                  mult++;
                              }
                              else
                                break;
                          }
                          for( int j=0; j<mult; j++ )
                          {
                              multiplicityTrue(i+j)=mult;
                              multIndexTrue(i+j)=i;
               // printF("eig=%3d,  multiplicityTrue=%d\n",i+j,multiplicityTrue(i+j));
                          }
                          i += mult;
                    }
            }

            Integrate myIntegrate;
            Integrate & integrate = orthogonalize ? ogev.dbase.get<Integrate>("integrate") : myIntegrate;

            if( eigenVectorsAreKnown && useAccurateInnerProduct && !orthogonalize )
            {
                integrate.updateToMatchGrid(cg);
                w=1;
                Real volume;
                volume = integrate.volumeIntegral(w);
                printF("Computed volume of the domain is %12.4e\n",volume);
            }

            fPrintF(outFile,"\\begin{table}[H]\\tableFont % you should set \\tableFont to \\footnotesize or other size\n");
            fPrintF(outFile,"\\begin{center}\n");
            if( eigenValuesAreKnown )
            {
                fPrintF(outFile,"\\begin{tabular}{|c|c|c|c|c|c|c|}  \\hline\n");
                fPrintF(outFile,"\\multicolumn{7}{|c|}{%s, order=%d} \\\\ \\hline\n",(const char*)gridNameNoPrefix,orderOfAccuracy);
                fPrintF(outFile,"   j    &         $\\lambda_j$        & $\\lambda_j$-err  & $\\phi_j$-err  & multe & multc & $\\| A\\phi - \\lambda\\phi\\|/\\lambda$     \\\\ \\hline\n");
            }
            else
            {
                fPrintF(outFile,"\\begin{tabular}{|c|c|c|c|}  \\hline\n");
                fPrintF(outFile,"\\multicolumn{4}{|c|}{%s, order=%d} \\\\ \\hline\n",(const char*)gridNameNoPrefix,orderOfAccuracy);
                fPrintF(outFile,"   j    &         $\\lambda_j$      & multc & $\\| A\\phi - \\lambda\\phi\\|/\\lambda$     \\\\ \\hline\n");
            }

      // --------------- START LOOP -----------
      // Real maxRelErr=0., maxEvectErr=0.;

            printF("  resid = || A phi - lambda phi ||_max / |lambda| \n");
            printF("  VECT max-err = max-norm( closest vector in 2 norm to eigenpsace )/( max-norm of phi)  \n");
            printF("  VECT l2-err  =      L2h( closest vector in 2 norm to eigenpsace )/( L2h norm of phi)  \n");
            for( int i=0; i<numEigenVectors; i++ )
            {
        // printF("Eigenvalue %d : k=%18.14e + %18.14e I",i,eig(0,i),eig(1,i));
                printF("Eigenvalue %3d : k=%14.8f + (%+6.1e) i",i,eig(0,i),eig(1,i));
                fPrintF(outFile,"  %4d  &  %9.3f + (%+6.1e) i  ",i,eig(0,i),eig(1,i));
                if( eigenValuesAreKnown )
                {
  
          // -- first find closest true eigenvalue ---
                    Real diffMin = REAL_MAX*.1;
                    int eigIndex=0;  
                    for( int ii=0; ii<numEigsTrue; ii++ )
                    {
                      Real diff = fabs( eig(0,i) - eigsTrue(ii) ) ;
                      if( diff < diffMin )
                      {
                            eigIndex=ii; diffMin=diff;
                      }
                    }
                    const int mult = multiplicityTrue(eigIndex);
                    const int ie = multIndexTrue(eigIndex);    // index for exact solution, and FIRST member of a multiple EIG

                    const Real eigTrue = eigsTrue(ie);
                    const Real absErr = fabs(eig(0,i)-eigTrue); 
                    const Real relErr = eigTrue !=0. ? absErr/eigTrue : absErr;

                    eigErr(i)=relErr; // save 

          // printF(", true=%14.8f, err=%8.2e, rel-err=%8.2e",eigTrue,absErr,relErr);
                    printF(", eig=%4d, true=%14.8f, rel-err=%8.2e",eigIndex,eigTrue,relErr);
                    fPrintF(outFile,"&    %8.2e      ",relErr);

                    if( eigenVectorsAreKnown )
                    {

              // -- check errors in eigenvectors  ---
                            int md=mult; 
                            RealArray a(md,md), ai(md,md), b(md), alpha(md);
              // *new* way 
              // find 
              //   u[i] approx  SUM_j alpha_j evTrue[j]
              // ( evTrue[ii], u[i] ) = SUM_j alpha_j ( evTrue[ii], evTrue[j] )
              //   A0 alpha = b0  
              // Solve
              //   A0^T A0 alpha = A^T b0
              //    A alpha = b 
                            for( int i1=0; i1<md; i1++ )
                            {
                // int ii = i1==0 ? i : j;
                                int ii = ie + i1;
                                    if( useAccurateInnerProduct )
                                    {
                                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                        {
                                            getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                            OV_GET_SERIAL_ARRAY(Real,evTrue[grid],uuLocal);
                                            OV_GET_SERIAL_ARRAY(Real,u[grid],vvLocal);
                                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                            wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ii)*vvLocal(I1,I2,I3,i);
                                        }
                                        b(i1) = integrate.volumeIntegral(w);
                                    }
                                    else
                                    {
                                        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                                        int iab[2]; 
                                        int ii1,ii2,ii3;  
                                        b(i1)=0.; 
                                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                        {
                                            MappedGrid & mg = cg[grid];
                                            const IntegerArray & gid = mg.gridIndexRange();
                                            OV_GET_SERIAL_ARRAY(Real,evTrue[grid],uuLocal);
                                            OV_GET_SERIAL_ARRAY(Real,u[grid],vvLocal);
                                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
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
                                            FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                            {
                        // wLocal(ii1,ii2,ii3) = uuLocal(ii1,ii2,ii3,ii)*vvLocal(ii1,ii2,ii3,i); // ** TEMP***
                                                if( maskLocal(ii1,ii2,ii3)>0 )
                                                    b(i1) += uuLocal(ii1,ii2,ii3,ii)*vvLocal(ii1,ii2,ii3,i);
                                            }
                                        }
                                    }  
                // innerProductor( ii, i, evTrue,u, w, b(i1) );  // ue[ii] . u[i] 
                                for( int i2=0; i2<md; i2++ )
                                {
                  // int jj = i2==0 ? i : j;
                                    int jj = ie + i2;
                                        if( useAccurateInnerProduct )
                                        {
                                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                            {
                                                getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                                OV_GET_SERIAL_ARRAY(Real,evTrue[grid],uuLocal);
                                                OV_GET_SERIAL_ARRAY(Real,evTrue[grid],vvLocal);
                                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                                wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ii)*vvLocal(I1,I2,I3,jj);
                                            }
                                            a(i1,i2) = integrate.volumeIntegral(w);
                                        }
                                        else
                                        {
                                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                                            int iab[2]; 
                                            int ii1,ii2,ii3;  
                                            a(i1,i2)=0.; 
                                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                            {
                                                MappedGrid & mg = cg[grid];
                                                const IntegerArray & gid = mg.gridIndexRange();
                                                OV_GET_SERIAL_ARRAY(Real,evTrue[grid],uuLocal);
                                                OV_GET_SERIAL_ARRAY(Real,evTrue[grid],vvLocal);
                                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
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
                                                FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                                {
                          // wLocal(ii1,ii2,ii3) = uuLocal(ii1,ii2,ii3,ii)*vvLocal(ii1,ii2,ii3,jj); // ** TEMP***
                                                    if( maskLocal(ii1,ii2,ii3)>0 )
                                                        a(i1,i2) += uuLocal(ii1,ii2,ii3,ii)*vvLocal(ii1,ii2,ii3,jj);
                                                }
                                            }
                                        }  
                  // innerProductor( ii,jj, evTrue,evTrue, w,a(i1,i2));    // ue[ii] . ue[jj]
                                }
                            }  
              // Compute the 1-norm of "a" as needed below for the condition number 
                            Real *rwork = new Real [4*md];
                            int *iwork = new int [md];
                            Real anorm = LANGE( "1", md, md, a(0,0), md, rwork[0] ); // "1" = 1-norm 
              // PA = LU factor
                            IntegerArray ipvt(md);
                            int info;
                            GETRF( md,md,a(0,0), md, ipvt(0), info );                
                            if( info!=0 )
                            {
                                printF("ERROR return from GETRF, info=%d\n",info);
                                OV_ABORT("error");
                            }
              // ----- check the condition number ---
                            Real rcond;
                            GECON( "1", md, a(0,0), md, anorm, rcond, rwork[0], iwork[0], info ); // "1" = 1-norm 
                            if( info!=0 )
                                OV_ABORT("GECON: info !=0 ");
                            if( rcond < 1e-6 )
                                printF("\n **WARNING** : 1-norm inverse condition number of A = %9.2e, anorm=%9.2e\n\n",rcond,anorm);
              // Solve given LU
                            int nrhs=1;
                            GETRS( "N", md, nrhs, a(0,0), md, ipvt(0), b(0), md, info );
                            alpha=b;
              // ::display(alpha,"alpha");
                            if( info!=0 )
                            {
                                printF("ERROR return from GETRS, info=%d\n",info);
                                OV_ABORT("error");
                            }  
                            delete [] rwork;
                            delete [] iwork;   
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                // general case 
                                w[grid] = u[grid](all,all,all,i); 
                                for( int ii=0; ii<mult; ii++ )
                                    w[grid] -= alpha(ii)*evTrue[grid](all,all,all,ie+ii);
                            }     
                            const Real phiMax = maxNorm( u, i );
                            const Real phiL2 =   l2Norm( u, i );
                            const Real errMax = maxNorm( w )/phiMax;
                            const Real errL2  =  l2Norm( w )/phiL2;

            // printF("i=%3d: max-err=%8.2e, l2-err=%8.2e\n",i,errMax,errL2);
                        printF(", multe=%d mult=%d VECT: max-err=%8.2e, l2-err=%8.2e, resid=%8.2e",multiplicityTrue(ie),eigMultiplicity(i),errMax,errL2,resMax(i));
                        
                        if( errMax > .05 ) printF(" @EV@");
                        if( resMax(i) > 1.0e-10 ) printF(" *RES*");
                        if( orthogonalize &&  multiplicityTrue(ie) != eigMultiplicity(i) ){ printF(" *mult differ*"); }

                        fPrintF(outFile,"&   %8.2e    &   %d    &  %d   &    %8.2e ",errMax,multiplicityTrue(ie),eigMultiplicity(i),resMax(i));

                        evectErr(i)=errMax;

                        if( eig(1,i)!=0. )
                            printF(" **COMPLEX** ");
                        
                    }
                    else
                    {
                        fPrintF(outFile,"&    "); // leave blank
                    }
                    printF("\n");

                }
                else // eigenvalues not known
                {
                    printF(", mult=%d, resid=%8.2e\n",eigMultiplicity(i),resMax(i));
                    fPrintF(outFile," &   %d     &    %8.2e   ",eigMultiplicity(i),resMax(i));
                }
                fPrintF(outFile,"\\\\\n");

            }

            Real maxRelErr=0.,maxEvectErr=0.; 
            if( eigenValuesAreKnown )
            {
                  maxRelErr=max(eigErr);
                  maxEvectErr=max(evectErr);
                  printF(" max-rel-err=%8.2e, max-evect-err=%8.2e",maxRelErr,maxEvectErr);
            }
            Real maxResMax = max(resMax);      

            fPrintF(outFile,"\\hline\n");
            fPrintF(outFile,"\\end{tabular}\n");
            fPrintF(outFile,"\\caption{Computed eigenvalues, relative-error in the eigenvalues, and relative error in the eigenvectors. orthogonalize=%d",orthogonalize);
            if( eigenValuesAreKnown )
                fPrintF(outFile,", max-rel-err=%8.2e, max-evect-err=%8.2e, max-residual=%8.2e\n",maxRelErr,maxEvectErr,maxResMax);
            else
                fPrintF(outFile,", max-residual=%8.2e\n",maxResMax);

            fPrintF(outFile,"}\\label{table:genEigs%s}\n",(const char*)gridNameNoPrefix);
            fPrintF(outFile,"\\end{center}\n");
            fPrintF(outFile,"\\end{table}\n");    
            fclose(outFile);


            printF("\n ==== SUMMARY:");
            printF(" max-residual = %8.2e\n",maxResMax);

            printF("\n Wrote file [%s]\n",(const char*)tableFileName);
        }   
        else if( answer=="save show file" ) 
        {
            if( nameOfShowFile == "" )
                nameOfShowFile = "genEigs.show";

      // -- save a show file ---
      // Ogshow show( nameOfShowFile );                      // create a show file
            bool useStreamMode=true;  // show files will be saved compressed
            Ogshow show( nameOfShowFile,".",useStreamMode );   
            show.saveGeneralComment("Results from genEigs");    // save a general comment in the show file

            show.setFlushFrequency( flushFrequency ); // save this many solutions per sub-showFile

      // ListOfShowFileParameters showFileParams;
      // showFileParams.push_back(ShowFileParameter("u1Component",u1c));
      // showFileParams.push_back(ShowFileParameter("u2Component",u2c));
      // show.saveGeneralParameters(showFileParams);

            show.startFrame(); 

      // -- save eigenvalues --
            HDF_DataBase *dbp=NULL;

            dbp = show.getFrame();
            assert( dbp!=NULL );

      // save parameters that go in this frame
            HDF_DataBase & db = *dbp;
            db.put(eig,"eig");   // computed eigenvalues 

      // ** we should always compute multiplicities here -- see cgWave 
            if( orthogonalize )
            {
        // --- eigenvectors were orthogonalized ----
                db.put(eigMultiplicity,"eigMultiplicity");
                db.put(eigStartIndex,"eigStartIndex");

        // save integration weights so we do not need to recompute
                Integrate & integrate = ogev.dbase.get<Integrate>("integrate");
                integrate.put(db,"integrate" );
            }

            if( eigenValuesAreKnown )
            {
        // -- OLD WAY --
                db.put(multiplicityTrue,"multiplicityTrue");
                db.put(multIndexTrue,"multIndexTrue");
            }


                                                        // start a new frame
            aString buff;
      // show.saveComment(0,sPrintF(buff,"Eigenvalues computed by genEigs with SLEPc, %s: order=%d",(const char*)problem,orderOfAccuracy));  
      // show.saveComment(1,sPrintF(buffer,"  t=%e ",t));            // comment 1 (shown on plot)

            bool saveEigenvectorsAsTimeSteps=true;
            if( saveEigenvectorsAsTimeSteps )
            {
        // -- save eigenvectors as separate solutions --

                realCompositeGridFunction q(cg,all,all,all);
                q.setName("phi",0);
                Index I1,I2,I3;
                for( int ie=0; ie<numEigenVectors; ie++ )
                {
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        OV_GET_SERIAL_ARRAY(real,q[grid],qLocal);
                        OV_GET_SERIAL_ARRAY(Real,u[grid],uLocal);
                        getIndex(cg[grid].dimension(),I1,I2,I3);
                        qLocal(I1,I2,I3)=uLocal(I1,I2,I3,ie);
                    }

                    if( ie>0 )
                        show.startFrame();                  // start a new frame

                    show.saveComment(0,sPrintF("genEigs: FD%i eig=%d lam=[%.4g,%.4g]",orderOfAccuracy,ie,eig(0,ie),eig(1,ie)));
                    show.saveSolution( q );              // save to show file
                    show.endFrame();

                }

            }
            else
            { // **old way** save EV's as components in one big grid function 

                show.saveComment(0,sPrintF(buff,"genEigs: %s: order=%d",(const char*)problem,orderOfAccuracy));  
                show.saveSolution( u );                                        // save the current grid function
                show.endFrame(); 
            }  

            show.close(); 
            printF("Wrote show file=[%s]\n",(const char*)nameOfShowFile);


        }

        else if( answer=="save matlab file" )
        {
            printF("Saving eigenvalues to a matlab file.\n");
  
            aString fileName;
            fileName = matlabFileName + ".m";
            FILE *matlabFile = fopen((const char*)fileName,"w" );     // Save some tex output here 
            fPrintF(matlabFile,"%% Eigenvalues for grid=[%s]\n",(const char*)nameOfGridFile);
            fPrintF(matlabFile,"%% File created eig/genEigs\n");
            fPrintF(matlabFile,"lambdav=[...\n");
            const int numPerLine=10;
            for( int i=0; i<numEigenVectors; i++ )
            {
                fPrintF(matlabFile,"(%16.10e %+16.10e *1i)",i,eig(0,i),eig(1,i));
                if( i<numEigenVectors-1 ) fPrintF(matlabFile,",");
                if( (i % numPerLine )==numPerLine-1 ) fPrintF(matlabFile,"...\n");
            }
            fPrintF(matlabFile,"\n];\n");

            fPrintF(matlabFile,"eigErr=[...\n");
            for( int i=0; i<numEigenVectors; i++ )
            {
                fPrintF(matlabFile,"%9.2e ",i,eigErr(i));
                if( i<numEigenVectors-1 ) fPrintF(matlabFile,",");
                if( (i % numPerLine )==numPerLine-1 ) fPrintF(matlabFile,"...\n");
            }
            fPrintF(matlabFile,"\n];\n");
            fPrintF(matlabFile,"evectErr=[...\n");
            for( int i=0; i<numEigenVectors; i++ )
            {
                fPrintF(matlabFile,"%9.2e ",i,evectErr(i));
                if( i<numEigenVectors-1 ) fPrintF(matlabFile,",");
                if( (i % numPerLine )==numPerLine-1 ) fPrintF(matlabFile,"...\n");
            }
            fPrintF(matlabFile,"\n];\n");      

            fclose(matlabFile);
            printF("Wrote file [%s]\n",(const char*)fileName);

        }
        else if( answer=="plot eigenvectors" ||
                          answer=="plot true eigenvectors" )
        {
            const int numberOfEigenvectors = numEigenVectors;
            const bool plotTrue = answer=="plot true eigenvectors";

            int vector=0; // plot this eigenvector 


            GUIState dialog;

            dialog.setWindowTitle("Eigenvectors");
            dialog.setExitCommand("done", "done");

            aString cmds[] = {"next",
                                                "previous",
                                                "first",
                                                "contour",
                                                ""};

            int numberOfPushButtons=4;  // number of entries in cmds
            int numRows=(numberOfPushButtons+1)/2;
            dialog.setPushButtons( cmds, cmds, numRows ); 

            const int numberOfTextStrings=8;
            aString textLabels[numberOfTextStrings];
            aString textStrings[numberOfTextStrings];


            int nt=0;
            textLabels[nt] = "vector";  sPrintF(textStrings[nt], "%d",vector);  nt++; 
      
      // null strings terminal list
            textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
            dialog.setTextBoxes(textLabels, textLabels, textStrings);
            int numberOfTextBoxes=nt;

            gi.pushGUI(dialog);      

            aString answer2; 
            Real t=0.; 
            for( int it=0; ; it++ )
            {
                if( it==0 )
                    answer2="plot";
                else
                    gi.getAnswer(answer2,"");

                if( answer2=="done" || answer2=="exit" || answer2=="finish" )
                {
                    break;
                }
                else if( dialog.getTextValue(answer2,"vector","%d",vector) )
                {
                    vector = max(0,min(numberOfEigenvectors-1,vector));
                    printF("Setting vector=%d\n",vector);
                }       
                else if( answer2=="first" )
                {
                    vector=0;
                }
                else if( answer2=="next" )
                {
                    vector = ( vector +1 ) % numberOfEigenvectors;
                }
                else if( answer2=="previous" )
                {
                    vector = ( vector - 1 + numberOfEigenvectors) % numberOfEigenvectors;

                }        
                gi.erase();
                if( abs(eig(1,vector)) < 1.e-10*abs(eig(0,vector)) )
                { // Real eigenvalue
                    psp.set(GI_TOP_LABEL,sPrintF(buff,"Eigenvector %d, eig=%.5g",vector,eig(0,vector)));
                }
                else
                { // complex eigenvalue 
                    psp.set(GI_TOP_LABEL,sPrintF(buff,"Eigenvector %d, eig=%.5g + (%.5g) I",vector,eig(0,vector),eig(1,vector)));
                }

                psp.set(GI_COMPONENT_FOR_CONTOURS,vector);
                if( answer2 == "contour" )
                    psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
                else
                    psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

                if( !plotTrue )
                {
                    PlotIt::contour(gi,u,psp);
                }
                else
                {
                    PlotIt::contour(gi,evTrue,psp);
                }

            }

            gi.popGUI(); // restore the previous GUI
        }
        else if( answer=="contour" ) // old way 
        {
            gi.erase();
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);      
            psp.set(GI_TOP_LABEL,sPrintF("Eigenvector"));  // set title
            PlotIt::contour(gi, u,psp );
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);      
        }    
        else if( dialog.getTextValue(answer,"time:","%g",time) ){}// 
        else if( dialog.getTextValue(answer,"numPoints:","%i",numPoints) ){}// 
        else if( dialog.getTextValue(answer,"show file:","%s",nameOfShowFile) ){}// 
        else if( dialog.getTextValue(answer,"flush frequency:","%i",flushFrequency) ){}// 
        else if( answer=="erase" )
        {
            gi.erase();
        }
        else if( answer=="check" && problem=="ile" )
        {  

      // --- check the ILE solution ----
            printF("--- check the incompressible elasticity solution : eigc=%d---\n",eigc);

      // eigenvalue: 
            const Real lambda = eig(0,eigc);  

            const int psim = 0 + numberOfComponents*(eigc); // index for psi
            const int pm   = 1 + numberOfComponents*(eigc); // index for p             


            Index Iv[3],  &I1  =Iv[0], &I2  =Iv[1], &I3=  Iv[2];
            Index Ibv[3], &Ib1=Ibv[0], &Ib2=Ibv[1], &Ib3=Ibv[2];
            Index Igv[3], &Ig1=Igv[0], &Ig2=Igv[1], &Ig3=Igv[2];
            BoundaryConditionParameters extrapParams;
            extrapParams.orderOfExtrapolation=orderOfAccuracy+1;


      // ---- Check residuals in psi equations ----
            ogev.checkResidualsInPsi( eigc, eig, u[0], cg, cgop, bc, numberOfComponents, orderOfAccuracy, 
                                                                useWideStencils, mu, includePressure ); 

      // ======== COMPUTE u1 and u2 from psi ==============
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                realMappedGridFunction & uvg = uv[grid];
                realMappedGridFunction & pg  = p[grid];
                OV_GET_SERIAL_ARRAY(Real,uvg,uLocal);

                const bool isRectangular = mg.isRectangular();
                Real dx[3]={1.,1.,1.};
                if( isRectangular )
                    mg.getDeltaX(dx);

                getIndex(mg.dimension(),I1,I2,I3);

                uv[grid](I1,I2,I3,0) =  u[grid].y(I1,I2,I3,psim)(I1,I2,I3,psim);   // u1 =  psi.y 
                uv[grid](I1,I2,I3,1) = -u[grid].x(I1,I2,I3,psim)(I1,I2,I3,psim);   // u2 = -psi.x 

        // First extrap all ghost 
                for( int ghost=1; ghost<=numGhost; ghost++ )
                {
                    extrapParams.ghostLineToAssign=ghost;
                    uvg.applyBoundaryCondition(Range(0,1),BCTypes::extrapolate,BCTypes::allBoundaries,0.,0.,extrapParams);
                }  

                if( orderOfAccuracy==2 )
                {
          // apply BC's to u1 and u2 -- maybe not needed if use use CBCs when solving for the eigenfunction
                    ForBoundary(side,axis)
                    {
                        getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                        getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);

                        const int is = 1-2*side; 
                        if( bc(side,axis)==displacement )
                        {
                            if( axis==0 )
                            {
                // u1.x = 0 
                                uvg.applyBoundaryCondition(u1c,BCTypes::neumann,    BCTypes::BCTypes::boundary(side,axis),0.);
                // uvg.applyBoundaryCondition(u2c,BCTypes::extrapolate,BCTypes::BCTypes::boundary(side,axis),0.);

                // mu*u2.xx = p.y - lambda*u2  **CHECK ME** TROUBLE: p is NOT KNOWN YET !!
                // uLocal(Ig1,Ig2,Ig3,u2c) = 2.*uLocal(Ib1,Ib2,Ib3,u2c) - uLocal(Ib1+is,Ib2,Ib3,u2c)
                //   + (SQR(dx[0])/mu)*( pg.y(Ib1,Ib2,Ib3)(Ib1,Ib2,Ib3) - lambda*uLocal(Ib1,Ib2,Ib3,u2c) );
                            }
                            else
                            {
                // u2.y = 0 
                                uvg.applyBoundaryCondition(u2c,BCTypes::neumann,    BCTypes::BCTypes::boundary(side,axis),0.);
                // uvg.applyBoundaryCondition(u1c,BCTypes::extrapolate,BCTypes::BCTypes::boundary(side,axis),0.);
                // mu*u1.yy = p.x - lambda*u1  ** CHECK ME **
                // uLocal(Ig1,Ig2,Ig3,u1c) = 2.*uLocal(Ib1,Ib2,Ib3,u1c) - uLocal(Ib1,Ib2+is,Ib3,u1c)
                //   + (SQR(dx[1])/mu)*( pg.x(Ib1,Ib2,Ib3)(Ib1,Ib2,Ib3) - lambda*uLocal(Ib1,Ib2,Ib3,u1c) ) ;             
                            }
                        }
                        else if( bc(side,axis)==traction )
                        {
              // Set v.x + u.y = 0 
              //     u.x + v.y = 0 
              // printF("Set BCs for traction BC -- finish me\n");
                            if( axis==0 )
                            {
                                uLocal(Ib1-is,Ib2,Ib3,u1c) = uLocal(Ib1+is,Ib2,Ib3,u1c) + 
                                                          (is*dx[0]/dx[1])*( uLocal(Ib1,Ib2+1,Ib3,u2c) - uLocal(Ib1,Ib2-1,Ib3,u2c) );
                                uLocal(Ib1-is,Ib2,Ib3,u2c) = uLocal(Ib1+is,Ib2,Ib3,u2c) + 
                                                          (is*dx[0]/dx[1])*( uLocal(Ib1,Ib2+1,Ib3,u1c) - uLocal(Ib1,Ib2-1,Ib3,u1c) );
                            }
                            else
                            {
                                uLocal(Ib1,Ib2-is,Ib3,u1c) = uLocal(Ib1,Ib2+is,Ib3,u1c) + 
                                                          (is*dx[1]/dx[0])*( uLocal(Ib1+1,Ib2,Ib3,u2c) - uLocal(Ib1-1,Ib2,Ib3,u2c) );
                                uLocal(Ib1,Ib2-is,Ib3,u2c) = uLocal(Ib1,Ib2+is,Ib3,u2c) + 
                                                          (is*dx[1]/dx[0])*( uLocal(Ib1+1,Ib2,Ib3,u1c) - uLocal(Ib1-1,Ib2,Ib3,u1c) );

                            }
                        }

                    }
          // ::display(uvg,"uvg after apply BC","%6.3f ");
                }
                else if( orderOfAccuracy==4 )
                {
          // apply BC's to u1 and u2 -- maybe not needed if use use CBCs when solving for the eigenfunction


          // u.x = 0 
          //   u(-2) -8*u(-1) + 8*u(1) -u(2) = 0 
          //   u(-2) = 5*u(-1) -10*u(0) +10*u(1) -5*u(2) + u(3)
          //   -3*u(1) -10*u(0) +18*u(1) -6*u(2) + u(3) = 0
                    const Real cex41= (-10./3.), cex42=6., cex43=-2., cex44=1./3.; 

                    ForBoundary(side,axis)
                    {
                        getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
            // getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);

                        const int is = 1-2*side; 
                        if( bc(side,axis)==displacement )
                        {
                            if( axis==0 )
                            {
                // u1.x = 0 
                // uvg.applyBoundaryCondition(u1c,BCTypes::neumann,    BCTypes::BCTypes::boundary(side,axis),0.);
                // uvg.applyBoundaryCondition(u2c,BCTypes::extrapolate,BCTypes::BCTypes::boundary(side,axis),0.);

                                uLocal(Ib1-is,Ib2,Ib3,u1c) = cex41*uLocal(Ib1     ,Ib2,Ib3,u1c) +
                                                                                          cex42*uLocal(Ib1+  is,Ib2,Ib3,u1c) +
                                                                                          cex43*uLocal(Ib1+2*is,Ib2,Ib3,u1c) +
                                                                                          cex44*uLocal(Ib1+3*is,Ib2,Ib3,u1c);
                                uLocal(Ib1-2*is,Ib2,Ib3,u1c) =  5.*uLocal(Ib1-1*is,Ib2,Ib3,u1c) +
                                                                                            -10.*uLocal(Ib1     ,Ib2,Ib3,u1c) +
                                                                                            +10.*uLocal(Ib1+1*is,Ib2,Ib3,u1c) +
                                                                                              -5.*uLocal(Ib1+2*is,Ib2,Ib3,u1c) +                                       
                                                                                              +1.*uLocal(Ib1+3*is,Ib2,Ib3,u1c);                                        

                // // mu*u2.xx = p.y - lambda*u2  **CHECK ME**
                // uLocal(Ig1,Ig2,Ig3,u2c) = 2.*uLocal(Ib1,Ib2,Ib3,u2c) - uLocal(Ib1+is,Ib2,Ib3,u2c)
                //   + (SQR(dx[0])/mu)*( pg.y(Ib1,Ib2,Ib3)(Ib1,Ib2,Ib3) - lambda*uLocal(Ib1,Ib2,Ib3,u2c) );

                            }
                            else
                            {
                // u2.y = 0 
                // uvg.applyBoundaryCondition(u2c,BCTypes::neumann,    BCTypes::BCTypes::boundary(side,axis),0.);
                                uLocal(Ib1,Ib2-is,Ib3,u2c) = cex41*uLocal(Ib1,Ib2     ,Ib3,u2c) +
                                                                                          cex42*uLocal(Ib1,Ib2+  is,Ib3,u2c) +
                                                                                          cex43*uLocal(Ib1,Ib2+2*is,Ib3,u2c) +
                                                                                          cex44*uLocal(Ib1,Ib2+3*is,Ib3,u2c);
                                uLocal(Ib1,Ib2-2*is,Ib3,u2c) =  5.*uLocal(Ib1,Ib2-1*is,Ib3,u2c) +
                                                                                            -10.*uLocal(Ib1,Ib2     ,Ib3,u2c) +
                                                                                            +10.*uLocal(Ib1,Ib2+1*is,Ib3,u2c) +
                                                                                              -5.*uLocal(Ib1,Ib2+2*is,Ib3,u2c) +                                       
                                                                                              +1.*uLocal(Ib1,Ib2+3*is,Ib3,u2c);   

      
                // uvg.applyBoundaryCondition(u1c,BCTypes::extrapolate,BCTypes::BCTypes::boundary(side,axis),0.);
                // mu*u1.yy = p.x - lambda*u1  ** CHECK ME **
                // uLocal(Ig1,Ig2,Ig3,u1c) = 2.*uLocal(Ib1,Ib2,Ib3,u1c) - uLocal(Ib1,Ib2+is,Ib3,u1c)
                //   + (SQR(dx[1])/mu)*( pg.x(Ib1,Ib2,Ib3)(Ib1,Ib2,Ib3) - lambda*uLocal(Ib1,Ib2,Ib3,u1c) ) ;             
                            }
                        }
                        else if( bc(side,axis)==traction )
                        {
              // Set v.x + u.y = 0 
              //     u.x + v.y = 0 

                            printF("Set BCs for traction BC -- finish me\n");
                        }            
            // for( int ghost=2; ghost<=numGhost; ghost++ )
            // {
            //   extrapParams.ghostLineToAssign=ghost;
            //   uvg.applyBoundaryCondition(Range(0,1),BCTypes::extrapolate,BCTypes::allBoundaries,0.,0.,extrapParams);
            // }            

                    }
          // ::display(uvg,"uvg after apply BC","%6.3f ");
                }
            }      

      // for now extrapolate ghost -- could do better
      // if( orderOfAccuracy>4 )
      // {
      //   for( int ghost=1; ghost<=numGhost; ghost++ )
      //   {
      //     extrapParams.orderOfExtrapolation=orderOfAccuracy+1; 
      //     extrapParams.ghostLineToAssign=ghost;
      //     uv.applyBoundaryCondition(Range(0,1),BCTypes::extrapolate,BCTypes::allBoundaries,0.,0.,extrapParams);
      //   }
      // }
            uv.finishBoundaryConditions();


      // Normalize eigen-solution (u1,u2) to have max norm of 1
            const Real scaleFactor = eigSign/max( maxNorm(uv,0), maxNorm(uv,1) );
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                uv[grid](all,all,all,Range(0,1)) *= scaleFactor;
            }  

      // Normalize eigen-solution psi to have max norm of 1 -- OR SHOULD WE JUST USE scaleFactor ??
            const Real psiScaleFactor = maxNorm(u,eigc);
            u(all,all,all,eigc) *= eigSign/psiScaleFactor;

            if( !includePressure )
            {
        // Solve for the pressure if it is not computed when computing psi
                ogev.getPressureFromDisplacement( uv, p, bc, orderOfAccuracy,mu );
            }
            else
            {
                if( max(abs(bc-displacement)) ==0 )
                {
          // displacement BCs all around : pressure equation is singular -- adjust mean value of p
                    Real sum=0., numPoints=0.; 
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        MappedGrid & mg = cg[grid];
                        OV_GET_SERIAL_ARRAY(Real,u[grid],uLocal);
                        getIndex(mg.gridIndexRange(),I1,I2,I3);
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            sum += uLocal(i1,i2,i3,pm);
                            numPoints += 1.;
                        }
                    }    
                    const Real mean = sum/numPoints;
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        MappedGrid & mg = cg[grid];
                        OV_GET_SERIAL_ARRAY(Real,u[grid],uLocal);
                        getIndex(mg.dimension(),I1,I2,I3);    
                        uLocal(I1,I2,I3,pm) -= mean;
                    }  
                    printF("Adjusting presssure mean to zero (mean was %8.2e)\n",mean);    
                }
            }


            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                getIndex(mg.dimension(),I1,I2,I3);
                if( includePressure )
                {
                    Real pSign = +1.; // what should this be ?   ************************************************ WHY ???
                    p[grid] = (pSign*scaleFactor)*u[grid](I1,I2,I3,pm);  // pressure computed with psi
                }

                uv[grid](I1,I2,I3,pc) = p[grid](I1,I2,I3);

                uv[grid](I1,I2,I3,psic) = u[grid](I1,I2,I3,psim);  // psi : chosen eigenvector 
            }



            if( nameOfShowFile != "" )
            {
        // -- save a show file ---
                Ogshow show( nameOfShowFile );                               // create a show file
          

                show.saveGeneralComment("Results from genEigs");    // save a general comment in the show file


                ListOfShowFileParameters showFileParams;
                showFileParams.push_back(ShowFileParameter("u1Component",u1c));
                showFileParams.push_back(ShowFileParameter("u2Component",u2c));
                show.saveGeneralParameters(showFileParams);

                show.startFrame();                                             // start a new frame
                aString buff;
                show.saveComment(0,sPrintF(buff,"%s: order=%d, nx=%d, ny=%d",(const char*)problem,orderOfAccuracy,nx,ny));  
        // show.saveComment(1,sPrintF(buffer,"  t=%e ",t));               // comment 1 (shown on plot)
                show.saveSolution( uv );                                        // save the current grid function
                show.endFrame(); 
                show.close();   
                printF("Wrote show file=[%s]\n",(const char*)nameOfShowFile);
            }   


      // ---- check the residual ----
            realCompositeGridFunction res(cg,all,all,all,2);
            res.setName("u1res",0);
            res.setName("u2res",1);
            res = 0.;


            Real resMax; 
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                const IntegerArray & gid = mg.gridIndexRange();
                realMappedGridFunction & uvg = uv[grid];
                realMappedGridFunction & pg = p[grid];

        // int extra=-1;
        // getIndex(mg.gridIndexRange(),I1,I2,I3,extra);
        // ---- Fill equation for psi points ----
        //   bc = displacement : skip boundary for psi
        //      = traction     : include boundary
                getIndex(gid,I1,I2,I3);
                for( int axis=0; axis<numberOfDimensions; axis++ )
                {
                    int ia = bc(0,axis)==displacement ? gid(0,axis)+1 : gid(0,axis); 
                    int ib = bc(1,axis)==displacement ? gid(1,axis)-1 : gid(1,axis); 
                    Iv[axis] = Range(ia,ib); 
                }          

                OV_GET_SERIAL_ARRAY(Real,res[grid],resLocal);

                resLocal(I1,I2,I3,0) = mu*( uvg.xx(I1,I2,I3,0)(I1,I2,I3,0) + uvg.yy(I1,I2,I3,0)(I1,I2,I3,0) ) 
                                                                        - pg.x(I1,I2,I3,0)(I1,I2,I3,0) + lambda*uvg(I1,I2,I3,0);
                resMax = max(fabs(resLocal(I1,I2,I3,0)));
                printF("grid=%d: lambda=%9.3e, max-residual in mu*Delta(u1) -p.x + lambda*u1  =%8.2e\n",grid,lambda,resMax);

                resLocal(I1,I2,I3,1) = mu*( uvg.xx(I1,I2,I3,1)(I1,I2,I3,1) + uvg.yy(I1,I2,I3,1)(I1,I2,I3,1) ) 
                                                                        - pg.y(I1,I2,I3,0)(I1,I2,I3,0) + lambda*uvg(I1,I2,I3,1);
                resMax = max(fabs(resLocal(I1,I2,I3,1)));
                printF("grid=%d: lambda=%9.3e, max-residual in mu*Delta(u2) -p.y + lambda*u2  =%8.2e\n",grid,lambda,resMax);  


                ForBoundary(side,axis)
                {
                    if( bc(side,axis)==traction )
                    {
                        int extra=-1; // skip ends
                        getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3,extra);
                        getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,1,extra);

                        resLocal(Ig1,Ig2,Ig3,0) = uvg.x(Ib1,Ib2,Ib3,u1c)(Ib1,Ib2,Ib3,u1c) + uvg.y(Ib1,Ib2,Ib3,u2c)(Ib1,Ib2,Ib3,u2c);
                        resMax = max(fabs(resLocal(Ig1,Ig2,Ig3,0)));
                        printF("grid=%d: side=%d axis=%d residual(u1.x+u2.y) =%9.3e (saved in ghost)\n",grid,side,axis,resMax);            

                        resLocal(Ig1,Ig2,Ig3,1) = uvg.x(Ib1,Ib2,Ib3,u2c)(Ib1,Ib2,Ib3,u2c) + uvg.y(Ib1,Ib2,Ib3,u1c)(Ib1,Ib2,Ib3,u1c);
                        resMax = max(fabs(resLocal(Ig1,Ig2,Ig3,1)));
                        printF("grid=%d: side=%d axis=%d residual(u2.x+u1.y) =%9.3e (saved in ghost)\n",grid,side,axis,resMax); 

                    }
                }      

            }     
            Real resL2u1 = l2Norm(res,0);
            Real resL2u2 = l2Norm(res,1);
            printF(" L2-norm residuals = [%8.2e,%8.2e]\n",resL2u1,resL2u2);

      // res.display("res");

            gi.erase();      

            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);      
            psp.set(GI_TOP_LABEL,sPrintF("Residual"));  // set title
            PlotIt::contour(gi, res,psp );
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true); 
            gi.erase();      

            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);      
            psp.set(GI_TOP_LABEL,sPrintF("Incompressible elasticity"));  // set title
            PlotIt::contour(gi, uv,psp );
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);  


                    

        }  
    // else if( answer=="movie" )
    // {
    //   // --- plot a movie of the eye motion ---

    //   // fix the plot bounds as the eye moves 
    //   RealArray xBound(2,3);
    //   xBound(0,0)=-1.2; xBound(1,0)= 1.2;
    //   xBound(0,1)=-1.2; xBound(1,1)= 1.2;
    //   xBound(0,2)=-1.;  xBound(1,2)= 1.;

    //   Real dt=.005*Pi, tFinal=10*Pi;
    //   int nStep=int( tFinal/dt + .5 );
    //   dt = tFinal/(nStep);

    //   Real yMax=-1.e10; // keep track of the largest y value of the eye-lid 

    //   RealArray x;
    //   for( int step=0; step<nStep; step++ )
    //   {
    //     Real t=step*dt;
                
    //     eyeCurves.getEyeCurve( x,t,numPoints );

    //     Range R = x.dimension(0);
    //     Real yTop = max(x(R,1));
    //     if( yTop>yMax )
    //     {
    //       yMax=yTop;
    //     }
    //     else if( yTop>yMax*(.99999) )
    //     {
    //       printF("Eye reaches yMax=%9.3e at t=%9.3e t/(2*pi)=%9.3e\n",yMax,t,t/twoPi);
    //     }
                
    //     gi.erase();
    //     psp.set(GI_TOP_LABEL,sPrintF(buff,"genEigs: t=%9.2e, yMax=%8.2e",t,yMax));
    //     plotCurve( x, gi,psp );
    //     gi.setGlobalBound(xBound); // set plot bounds         

    //     gi.redraw(true);
    //   }

    // }

    // else if( answer=="save file" )
    // {
    //   aString fileName = "eyeCurveDataPoints.dat";
    //   eyeCurves.saveEyeCurve( time, numPoints, fileName );
    //   printF("Eye coordinates written to file=[%s]\n",(const char*)fileName);
    // }
        
        else
        {
            printF("Unknown response=[%s]\n",(const char*)answer);
        }
    }
    

    gi.unAppendTheDefaultPrompt();
    gi.popGUI(); // restore the previous GUI


    int ierr = SlepcFinalize();

    Overture::finish();  

    return 0;
}

