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

#include <slepceps.h>



#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

class Ogev 
{

public:
  Ogev();
  ~Ogev();

  aString bcName( int bc );

  int checkResidualsInPsi( int eigc, RealArray & eig, realMappedGridFunction & u, CompositeGrid & cg,  
                           CompositeGridOperators & cgop, IntegerArray & bc, 
                           int numberOfComponents, int orderOfAccuracy, int useWideStencils, Real mu, int includePressure );  

  int computeEigenvalues( const aString & problem, const int numberOfComponents,
                          int orderOfAccuracy, int & numEigenValues, int & numEigenVectors, 
                          RealArray & eig, CompositeGrid & cg, realCompositeGridFunction & ucg, 
                          CompositeGridOperators & cgop, 
                          Real tol, int eigOption,
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
  getPressureFromDisplacement( realCompositeGridFunction & uv, realCompositeGridFunction & p, 
                               IntegerArray & bc, int orderOfAccuracy, Real mu  );                          

protected:

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

   // Here is the place to store parameters:
  DataBase dbase;

};


#endif