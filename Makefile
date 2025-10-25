# --------------------------------------
#     Various routines using SLEPc
# --------------------------------------
#
# Main programs:
#   genEigs      -- Compute eigevalues and eigenvectors of the Laplacian on overset grids 
#   getnEigsILDE -- Compute eigevalues and eigenvectors of incompressible elasticity
#
# NOTE: To compile optimized:
#   setenv COMPILE [opt|dbg]
#

all = genEigs
all: $(all);


include ${Overture}/make.options

usePETSc := on
# usePETSc := off
ifeq ($(usePETSc),on)
	usePETSc   = on
	petscSolver = obj/solvePETSc.o

	OGES_PETSC = buildEquationSolvers.o PETScEquationSolver.o
	PETSC_INCLUDE = -I$(PETSC_DIR)/include  -I$(PETSC_DIR)/$(PETSC_ARCH)/include -DOVERTURE_USE_PETSC -I$(PETSC_LIB)/include -I$(PETSC_DIR)/include/petsc/mpiuni
	# PETSC_INCLUDE = -I$(PETSC_DIR)/include  -I$(PETSC_DIR)/$(PETSC_ARCH)/include -DOVERTURE_USE_PETSC -I$(PETSC_LIB)/include -I$(PETSC_DIR)/include/mpiuni
	PETSC_LIBS = -Wl,-rpath,$(PETSC_LIB) -L$(PETSC_LIB) -lpetsc

  SLEPC_INCLUDE = -I$(SLEPC_DIR) -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(SLEPC_DIR)/include 

  SLEPC_LIBS = $(LIB_ARPACK) -Wl,-rpath,$(SLEPC_DIR)/$(PETSC_ARCH)/lib -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc

else
	usePETSc   = off
	petscSolver = obj/solvePETScNull.o

	OGES_PETSC = 
	PETSC_INCLUDE = 
	PETSC_LIBS = 
endif


# optimization flag:
#  -Ofast = -O3 + disregards strict standard compliance 
# OPTFLAG = -O
OPTFLAG = -O3
# OPTFLAG = -Ofast


CCFLAGS = -I$(SLEPC_DIR) $(SLEPC_INCLUDE)
CCFLAGS += $(PETSC_INCLUDE)

# for macos 
CCFLAGS += -std=c++11 

CCFLAGS += $(OV_CXX_FLAGS) -I. -I$(Overture)/include -I$(APlusPlus)/include -I$(OpenGL)/include $(USE_PPP_FLAG)
# CCFLAGS += -I/home/henshw/software/slepc-3.4.4 -I/home/henshw/software/slepc-3.4.4/linux-gnu-opt/include -I/home/henshw/software/slepc-3.4.4/include 

CFLAGS = $(OV_CC_FLAGS) -I. -I$(Overture)/include -I$(APlusPlus)/include

# Fortran flags, not optimzed unless COMPILE set to opt
FFLAGS  = $(OV_FORTRAN_FLAGS) $(OV_AUTO_DOUBLE_FLAGS) 

# Fortran flags, optimized by default 
FFLAGSO = $(OV_FORTRAN_FLAGS) $(OV_AUTO_DOUBLE_FLAGS) $(OPTFLAG)
#
# Fortran flags for debug always
FFLAGSG = $(OV_FORTRAN_FLAGS) $(OV_AUTO_DOUBLE_FLAGS) -g 

# no auto-double:
FFLAGSSO = $(OV_FORTRAN_FLAGS)

ifeq ($(COMPILE),opt)
  CCFLAGS += $(OPTFLAG)
  FFLAGS  += $(OPTFLAG)
  FFLAGSSO += $(OPTFLAG)
else
	ifeq ($(COMPILE),dbg)
    # debug:
    CCFLAGS += -w -g -finit-real=snan
    FFLAGS  += -g
    FFLAGSO += -g
    FFLAGSSO += g   
  else
    # default case:
    CCFLAGS += -w -g 
    FFLAGS  += -g 
    FFLAGSO += $(OPTFLAG) 
    FFLAGSG += -g
    FFLAGSSO += $(OPTFLAG)   
  endif
endif

# ifeq ($(COMPILE),opt)
# 	CCFLAGS += -w -O 
# 	FFLAGS  += -O
# 	FFLAGSSO += -O3
# else
# 	ifeq ($(COMPILE),dbg)
# 		# debug:
# 		CCFLAGS += -w -g
# 		FFLAGS  += -g
# 		FFLAGSO += -g
# 		FFLAGSSO += g   
# 	else
# 		# default case:
# 		CCFLAGS += -w -g 
# 		FFLAGS  += -g 
# 		FFLAGSO += -O 
# 		FFLAGSG += -g
# 		FFLAGSSO += O   
# 	endif
# endif

LIB_ARPACK = -Wl,-rpath,/home/henshw/software/arpack-ng/lib -L/home/henshw/software/arpack-ng/lib -larpack 

# LAPACK_LIBS =  -Wl,-rpath,$(PGI)/linux86-64/7.0-4/lib -L$(PGI)/linux86-64/7.0-4/lib -llapack -lblas
LAPACK_LIBS =  -L$(LAPACK) -llapack -lblas
# LIBS = $(OVERTURE_LIBRARIES) $(APP_LIBRARIES) $(OV_HDF_LIBRARIES) $(OV_COMPILER_LIBS) 
# LIBS +=  $(OV_FORTRAN_LIBRARIES) $(OV_PERL_LIBRARIES) $(OV_OPENGL_LIBRARIES) $(OV_MOTIF_LIBRARIES) $(OV_X_LIBRARIES) $(LAPACK_LIBS) 

# Ubuntu -- move mesa libs forward
LIBS = $(OVERTURE_LIBRARIES) $(OV_OPENGL_LIBRARIES) $(APP_LIBRARIES) $(OV_HDF_LIBRARIES) $(OV_COMPILER_LIBS) 

# SLEPc
# LIBS += -Wl,-rpath,/home/henshw/software/slepc-3.4.4/linux-gnu-opt/lib -L/home/henshw/software/slepc-3.4.4/linux-gnu-opt/lib -lslepc
LIBS += -Wl,-rpath,$(SLEPC_DIR)/$(PETSC_ARCH)/lib -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc

# LIBS +=  $(PETSC_LIBS) $(OV_FORTRAN_LIBRARIES) $(OV_PERL_LIBRARIES) $(OV_MOTIF_LIBRARIES) $(OV_X_LIBRARIES) $(LAPACK_LIBS) 
LIBS += $(PETSC_LIBS) $(OV_FORTRAN_LIBRARIES) $(OV_PERL_LIBRARIES) $(OV_MOTIF_LIBRARIES) $(OV_X_LIBRARIES) $(LAPACK_LIBS) 


.SUFFIXES:
.SUFFIXES:.C .o .f .o .F .o .bf .f .c .o .cc .o
.C.o:; $(CXX) $(CCFLAGS) -c $<
.cc.o:; $(CXX) $(CCFLAGS) -c $<
.c.o:; $(CC) $(CFLAGS) -c $<
.f.o:; $(FC) $(FFLAGSO) -c $<
.F.o:; $(FC) $(FFLAGSO) -c $<
.bf.f:; bpp $*.bf
.bC.C:; bpp $*.bC
# .bf.o: $*.f ; 
.bC.o: $*.C ; 

# compile some fortran files optimized by default:
$(OBJO) : obj/%.o : %.f
	$(FC) $(FFLAGSO) -o $@ -c $<


# %.o : %.C ; $(CXX) $(CCFLAGS) -c $*.C

# .C: $(LIB_DEPENDENCIES)
#	 $(CXX) $(CCFLAGS) -o $@ $< $(CLIBS) $(FLIBS)  $(GLIBS)

BPP = bpp
%.C : %.bC
	$(BPP) -quiet -clean  $<
%.f : %.bf
	$(BPP) -quiet -clean $<


current = .
mapping = $(current)/../mapping
ogshow = $(current)/../ogshow

VPATH = src:obj



Oges = $(Overture)/Oges
linkFiles:
	ln -sf $(Oges)/PETScEquationSolver.C src/
	ln -sf $(Oges)/PETScSolver.C src/
	ln -sf $(Oges)/buildEquationSolvers.C src/

# This is needed with PETSc
buildEquationSolvers.o : $(Oges)/buildEquationSolvers.C; $(CXX) $(CCFLAGS) -DOVERTURE_USE_PETSC -c $(Oges)/buildEquationSolvers.C
PETScEquationSolver.o : $(Oges)/PETScEquationSolver.C; $(CXX) $(CCFLAGS) -DOVERTURE_USE_PETSC -c $(Oges)/PETScEquationSolver.C

OBJC = obj/genEigs.o obj/Ogev.o obj/computeEigenvalues.o obj/fillMatrixLaplacian.o \
			 obj/fillMatrixIncompressibleElasticity.o obj/fillInterpolationCoefficients.o obj/orthogonalize.o obj/residual.o \
			 obj/coarseToFine.o obj/fillMatrixLaplacianComplex.o  \
			 $(OGES_PETSC)

# Fortran 90 (FN) object files: 
FNOBJO = obj/sumEigenvectors.o 
# Fortran obj:
FOBJO =  obj/rjbesl.o obj/rybesl.o

# OBJO = obj/advWave.o\
#         obj/advWave2dOrder2r.o obj/advWave2dOrder2c.o obj/advWave2dOrder4r.o obj/advWave2dOrder4c.o 


# For regression tests: (Note: the master version of checkCheckFiles is in cg/common/src --> could move to Overture)
checkCheckFiles = obj/checkCheckFiles.o 
checkCheckFiles: $(checkCheckFiles); $(CXX) $(CCFLAGS) -o check/checkCheckFiles $(checkCheckFiles) $(LIBS)

# --- genEigs solver ---
genEigsFiles = $(OBJC) $(OBJO) $(FNOBJO) $(FOBJO)
genEigs: $(genEigsFiles) ; $(CXX) $(CCFLAGS) -o bin/genEigs $(genEigsFiles) $(LIBS)


# --- genEigsILE :  solver for incompressible linear elasticity ILE ---
genEigsILEFiles = $(OBJC) $(OBJO) $(FNOBJO)  $(FOBJO)
genEigsILE: $(genEigsILEFiles) ; $(CXX) $(CCFLAGS) -o bin/genEigsILE $(genEigsILEFiles) $(LIBS)

genEigsOld: genEigs.C
	gcc -o genEigs.o -c -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O   -I/home/henshw/software/petsc-3.4.5-serial/include -I/home/henshw/software/petsc-3.4.5-serial/linux-gnu-opt/include -I/home/henshw/software/petsc-3.4.5-serial/include/mpiuni    -D__INSDIR__=src/eps/examples/tutorials/ -I/home/henshw/software/slepc-3.4.4 -I/home/henshw/software/slepc-3.4.4/linux-gnu-opt/include -I/home/henshw/software/slepc-3.4.4/include genEigs.C
	gcc -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O  -o genEigs genEigs.o -Wl,-rpath,$(SLEPC_DIR)/linux-gnu-opt/lib -L$(SLEPC_DIR)/linux-gnu-opt/lib -lslepc       -Wl,-rpath,/home/henshw/software/petsc-3.4.5-serial/linux-gnu-opt/lib -L/home/henshw/software/petsc-3.4.5-serial/linux-gnu-opt/lib  -lpetsc -llapack -lblas -lX11 -lpthread -lm -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lstdc++ -ldl -lgcc_s -ldl


eveSolverFiles = obj/eveSolver.o $(FNOBJO) $(FOBJO)
eveSolver: $(eveSolverFiles); $(CXX) $(CCFLAGS) -o bin/eveSolver $(eveSolverFiles) $(LIBS)

ex2: ex2.c
	gcc -o ex2.o -c -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O   -I/home/henshw/software/petsc-3.4.5-serial/include -I/home/henshw/software/petsc-3.4.5-serial/linux-gnu-opt/include -I/home/henshw/software/petsc-3.4.5-serial/include/mpiuni    -D__INSDIR__=src/eps/examples/tutorials/ -I/home/henshw/software/slepc-3.4.4 -I/home/henshw/software/slepc-3.4.4/linux-gnu-opt/include -I/home/henshw/software/slepc-3.4.4/include ex2.c
	gcc -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O  -o ex2 ex2.o -Wl,-rpath,$(SLEPC_DIR)/linux-gnu-opt/lib -L$(SLEPC_DIR)/linux-gnu-opt/lib -lslepc       -Wl,-rpath,/home/henshw/software/petsc-3.4.5-serial/linux-gnu-opt/lib -L/home/henshw/software/petsc-3.4.5-serial/linux-gnu-opt/lib  -lpetsc -llapack -lblas -lX11 -lpthread -lm -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lstdc++ -ldl -lgcc_s -ldl 





# SLEPc example - matrix free
ex3: ex3.C
	gcc -o ex3.o -c -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O   -I$(PETSC_DIR)/include -I$(PETSC_DIR)/linux-gnu-opt/include -I$(PETSC_DIR)/include/petsc/mpiuni    -D__INSDIR__=src/eps/examples/tutorials/ -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/linux-gnu-opt/include ex3.C
	gcc -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O  -o ex3 ex3.o $(LIB_ARPACK) -Wl,-rpath,$(SLEPC_DIR)/linux-gnu-opt/lib -L$(SLEPC_DIR)/linux-gnu-opt/lib -lslepc   -Wl,-rpath,$(PETSC_DIR)/linux-gnu-opt/lib -L$(PETSC_LIB) -lpetsc -llapack -lblas -lX11 -lpthread -lm -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lstdc++ -ldl -lgcc_s -ldl 





cpuTest: cpuTest.C;  
	gcc -o cpuTest.o -c cpuTest.C
	gcc -o cpuTest cpuTest.o 	-lm

# cpuTest: cpuTest.C
#   gcc -o cpuTest.o -c cpuTest.C
#   gcc -fPIC -Wall -O  -o cpuTest cpuTest.o 


# # ----- test of coeff matricies with wide stencils -----
# # tcmWideStencilFiles = obj/tcmWideStencil.o $(OBJC) $(OBJO) $(FNOBJO)
# tcmWideStencilFiles = obj/tcmWideStencil.o obj/cgesl1234.o $(OBJC) $(OBJO) $(petscSolver) $(FNOBJO) 
# tcmWideStencil: $(tcmWideStencilFiles); $(CXX) $(CCFLAGS) -o bin/tcmWideStencil $(tcmWideStencilFiles) $(LIBS)

# obj/cgesl1234.o : src/cgesl1234.F; $(FC) $(FFLAGSO) -I.  -DOV_USE_DOUBLE  -o $*.o -c $<
# # cgesl1234.o:
# #   gfortran -O  -fPIC  -fdefault-real-8 -fdefault-double-8   -I/home/henshw/Overture.g/include -I.   -DOV_USE_DOUBLE -c cgesl1234.F  

# obj/tcmWideStencil.o : src/tcmWideStencil.C; $(CXX) $(CCFLAGS) -o $*.o -c $<
# obj/solvePETScNull.o : src/solvePETScNull.C; $(CXX) $(CCFLAGS) -o $*.o -c $<


# --------- CgWaveHoltz ----------

# test: 
# 	-@echo "usePETSc=[$(usePETSc)]"
# 	-@echo "petscSolver=$(petscSolver)"

# # cgwh = obj/cgwh.o obj/CgWaveHoltz.o $(petscSolver) $(OBJC) $(OBJO)
# cgwh = obj/cgwh.o obj/CgWaveHoltz.o obj/solveHelmholtz.o $(petscSolver) $(OGES_PETSC) $(OBJC) $(OBJO) $(FNOBJO)
# cgwh: $(cgwh) 
# 	$(CXX) $(CCFLAGS) -o bin/cgwh $(cgwh) $(PETSC_LIBS) $(LIBS)


# bpp files: 
src/genEigs.C: src/genEigs.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include genEigs.bC

src/computeEigenvalues.C: src/computeEigenvalues.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include computeEigenvalues.bC  
src/fillMatrixLaplacian.C: src/fillMatrixLaplacian.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include fillMatrixLaplacian.bC  
src/fillMatrixLaplacianComplex.C: src/fillMatrixLaplacianComplex.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include fillMatrixLaplacianComplex.bC  
src/fillMatrixIncompressibleElasticity.C: src/fillMatrixIncompressibleElasticity.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include fillMatrixIncompressibleElasticity.bC  

src/fillInterpolationCoefficients.C: src/fillInterpolationCoefficients.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include fillInterpolationCoefficients.bC  

src/genEigsILE.C: src/genEigsILE.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include genEigsILE.bC

src/eveSolver.C: src/eveSolver.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include eveSolver.bC  

src/Ogev.C: src/Ogev.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include Ogev.bC  

src/orthogonalize.C: src/orthogonalize.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include orthogonalize.bC 
src/coarseToFine.C: src/coarseToFine.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include coarseToFine.bC 
src/residual.C: src/residual.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include residual.bC  


#  -- optimized routines
src/sumEigenvectors.f90: src/sumEigenvectors.bf90 
	      @cd src; $(BPP) -clean -quiet -I$(Overture)/include sumEigenvectors.bf90



# src/sumEigenvectors.f90 : src/sumEigenvectors.f90



# # -- optimized BC routine
# src/bcOptWave.f90: src/bcOptWave.bf90
# 	      @cd src; $(BPP) -clean -quiet -I$(Overture)/include bcOptWave.bf90	

# # dependencies
# obj/getDt.o : src/getDt.C; $(CXX) $(CCFLAGS) -o $*.o -c $<
# obj/cgwh.o : src/cgwh.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
# obj/CgWaveHoltz.o : src/CgWaveHoltz.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
# obj/solvePETSc.o : src/solvePETSc.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
# obj/solvePETScNull.o : src/solvePETScNull.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<	
# obj/solveHelmholtz.o : src/solveHelmholtz.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<	

# # --- cgWave ---
# obj/cgWaveMain.o : src/cgWaveMain.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<

obj/genEigs.o : src/genEigs.C; $(CXX) $(CCFLAGS) -o $*.o -c $<

obj/Ogev.o : src/Ogev.C src/Ogev.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/computeEigenvalues.o : src/computeEigenvalues.C src/Ogev.h; $(CXX) $(CCFLAGS) -o $*.o -c $<

obj/fillMatrixLaplacian.o : src/fillMatrixLaplacian.C; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/fillMatrixLaplacianComplex.o : src/fillMatrixLaplacianComplex.C; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/fillMatrixIncompressibleElasticity.o : src/fillMatrixIncompressibleElasticity.C; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/fillInterpolationCoefficients.o : src/fillInterpolationCoefficients.C; $(CXX) $(CCFLAGS) -o $*.o -c $<

obj/genEigsILE.o : src/genEigsILE.C; $(CXX) $(CCFLAGS) -o $*.o -c $<

obj/orthogonalize.o : src/orthogonalize.C; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/coarseToFine.o : src/coarseToFine.C; $(CXX) $(CCFLAGS) -o $*.o -c $<	
obj/residual.o : src/residual.C; $(CXX) $(CCFLAGS) -o $*.o -c $<

obj/eveSolver.o : src/eveSolver.C; $(CXX) $(CCFLAGS) -o $*.o -c $<


# done below: 
# obj/sumEigenvectors.o : src/sumEigenvectors.f90; $(FC) $(FFLAGSO) -ffree-line-length-none -o $*.o -c $<	

obj/rjbesl.o : src/rjbesl.f; $(FC) $(FFLAGSO) -o $@ -c $<  
obj/rybesl.o : src/rybesl.f; $(FC) $(FFLAGSO) -o $@ -c $<  


# compile f90 files optimized by default:
$(FNOBJO) : obj/%.o : %.f90
	$(FC) $(FFLAGSO) -ffree-line-length-none -o $@ -c $<	

# # compile f90 files optimized by default:
# $(FNOBJO) : obj/%.o : %.f90
# 	$(FC) $(FFLAGSO) -ffree-line-length-none -finit-real=snan -o $@ -c $<	

clean:  
	rm -f obj/*.o bin/genEigs bin/eveSolver


.PRECIOUS: 





