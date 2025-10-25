#ifndef OGEV_H
#define OGEV_H "Ogev.h"

//===========================================================================
//
//  Overlapping Grid EigenValue/EigenVector solver
//
//===========================================================================

// Follow the order of include files used in PETScEquationSolver.h


#include "mpi.h"
#include "Overture.h"

// #include <slepceps.h>

// forward declarations
typedef struct _p_Mat *Mat;
class CompositeGridOperators;

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

class Ogev 
{

public:

  enum BoundaryConditionsEnum
  {
     periodic      =-1, 
     interpolation = 0, 
     dirichlet     = 1, 
     neumann       = 2,
     characteristic= 3, // radiation BC
     displacement  = 1, // for elasticity
     traction      = 2, // for elasticity
  };

  Ogev();
  ~Ogev();

  aString bcName( int bc );

  int checkResidualsInPsi( int eigc, RealArray & eig, realMappedGridFunction & u, CompositeGrid & cg,  
                           CompositeGridOperators & cgop, IntegerArray & bc, 
                           int numberOfComponents, int orderOfAccuracy, int useWideStencils, Real mu, int includePressure );  

  bool compositeGridsMatch( CompositeGrid & cg, CompositeGrid & cgsf );

  int computeEigenvalues( const aString & problem, const int numberOfComponents,
                          int orderOfAccuracy, int & numEigenValues, int & numEigenVectors, 
                          RealArray & eig, CompositeGrid & cg, realCompositeGridFunction & ucg, 
                          CompositeGridOperators & cgop, 
                          Real tol, int eigOption, int maxIterations,
                          const int setTargetEigenvalue, const Real targetEigenvalue,
                          IntegerArray & bc, int numGhost, int saveMatlab, int useWideStencils,
                          int maximumProjectedDimension=-1  );

  int 
  fillInterpolationCoefficients( Mat & A, realCompositeGridFunction & uu );

  int 
  fillMatrixIncompressibleElasticity( const int numberOfComponents, int orderOfAccuracy, MappedGrid & mg, Mat & A, Mat & B, 
                                      int numGhost, bool useNew, Real tol, int eigOption, IntegerArray & bc, 
                                      int saveMatlab, int useWideStencils );

  int 
  fillMatrixLaplacian( int orderOfAccuracy, realCompositeGridFunction & ucg, CompositeGridOperators & cgop, 
                       Mat & A, Mat & B, int numGhost, bool useNew, 
                       Real tol, int eigOption, IntegerArray & bc, int saveMatlab, Real lambdaShift );

  // Complex case: quadratic eigenvalue problem:
  int 
  fillMatrixLaplacianComplex( int orderOfAccuracy, realCompositeGridFunction & ucg, CompositeGridOperators & cgop, 
                              Mat & A, Mat & B, int numGhost, bool useNew, 
                              Real tol, int eigOption, IntegerArray & bc, int saveMatlab, Real lambdaShift );  

  Real getEigenPairResidual( Real lambdar, Real lambdai, realCompositeGridFunction & v,
                             realCompositeGridFunction & res,  
                             CompositeGridOperators & operators, 
                             RealArray & resbc,
                             int component =0  );
  int 
  getEigenvaluesBox( int numEigs, RealArray & eigs, CompositeGrid & cg,
                     Real lx =1.0 , Real ly =1.0, Real lz =1.0,
                     RealCompositeGridFunction *eigenvector = NULL,
                     const bool discreteEigenvalues =false  );

  int 
  getEigenvaluesCylinder( int numEigs, RealArray & eigs, CompositeGrid & cg, 
                          Real ra =0.5, Real rb =1.0, Real za = 0.0, Real zb = 1.0,
                          RealCompositeGridFunction *eigenvector = NULL  );
  int 
  getEigenvaluesSphere( int numEigs, RealArray & eigs, CompositeGrid & cg, 
                        Real ra =0.0, Real rb =1.0, 
                        RealCompositeGridFunction *eigenvector = NULL  );

  int 
  getTrueEigenValues( const aString & eigCase, CompositeGrid & cg, int numEigsTrue, RealArray & eigsTrue,
                      bool & eigenValuesAreKnown, bool & eigenVectorsAreKnown, 
                      RealCompositeGridFunction *evTrue=NULL, int discreteEigenvalues=0 );

  int 
  getPressureFromDisplacement( realCompositeGridFunction & uv, realCompositeGridFunction & p, 
                               IntegerArray & bc, int orderOfAccuracy, Real mu  );                          

  Real getRayleighQuotient( realCompositeGridFunction & v, int component, CompositeGrid & cg, CompositeGridOperators & cgop );

  // Return true if this is a complex valued problem (e.g. radiation BCs)
  bool isComplexProblem() const;

  // normalize eigenvectors
  int normalizeEigenvectors( const aString & problem, const int numberOfComponents,
                                 int orderOfAccuracy, int & numEigenValues, int & numEigenVectors, 
                                 RealArray & eig, realCompositeGridFunction & ucg );

  // Force the solution to be from the complex problem
  int setIsComplexProblem( bool isComplex=true );

  // read eigenvectors from a file
  int readEigenvectors(  const aString & evFile, CompositeGrid & cg, const int orderOfAccuracy, 
                         int & numEigenValues, int & numEigenVectors, 
                         RealArray & eig, realCompositeGridFunction & uev, 
                         IntegerArray & eigMultiplicity, IntegerArray & eigStartIndex, 
                         bool & eigenVectorsAreOrthogonal );

  // save eigenvectors to a file 
  int saveEigenvectors(  const aString & evFile, const int & numEigenValues, const int & numEigenVectors, 
                         const RealArray & eig, const realCompositeGridFunction & uev, 
                         const IntegerArray & eigMultiplicity, const IntegerArray & eigStartIndex );

  // // normalize eigenvectors, fixup eigenvectors for multiple eigenvalues.
  // int normalizeEigenvector( int eigNumber, RealArray & eig, IntegerArray & eigMultiplicity, realCompositeGridFunction & uev  );

  // count multiplicities, orthogonalize and normalize eigenvectors
  int orthogonalizeEigenvectors( const aString & problem, const int numberOfComponents,
                                 int orderOfAccuracy, int & numEigenValues, int & numEigenVectors, 
                                 RealArray & eig, realCompositeGridFunction & ucg, 
                                 IntegerArray & eigMultiplicity, IntegerArray & eigStartIndex  );

  // convert coarse grid eigenpairs to fine grid eigenpairs
  int coarseToFine( const aString & eigCase, const int numberOfComponents,
                                 int orderOfAccuracy, int & numEigenValues, int & numEigenVectors, 
                                 RealArray & eig, realCompositeGridFunction & ucg, 
                                 IntegerArray & eigMultiplicity, IntegerArray & eigStartIndex, CompositeGrid & cg, CompositeGridOperators & cgop  );

  // Set interpolation stencil width for coarse-to-fine transfers
  int setInterpolationWidth( int interpolationWidth );

  // Use accurate inner product for the Rayleigh quotient
  int setUseAccurateInnerProduct( int useAccurateInnerProduct );

  // Set shift for eigenvalue problem to avoid a singular A
  int setShift( Real shift );

  Real getShift( ) const;

protected:

  int getAdjustedBoundaryIndex( MappedGrid & mg, int side, int axis, Index & Ib1, Index & Ib2, Index & Ib3 );

  Real getDiscreteSymbol( const Real modeNumber, const Real dx ) const;

  int
  buildGlobalIndexing(CompositeGrid & cg );
  
  int 
  getGlobalIndex( int n, int *iv, int grid, int p );

  int 
  getGlobalIndex( int n, int i1, int i2, int i3, int grid, int p );

  int debug;

  int numberOfComponents;
  int numberOfProcessors;
  int numberOfGridPoints;
  int numberOfGridPointsThisProcessor;

  // --- arrays for global indexing -----
  int *pnab, *pnoffset;

public: // do this for now
   // Here is the place to store parameters:
  DataBase dbase;

};


#endif