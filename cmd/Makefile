
genEigs = $(HOME)/research/eigSolvers/bin/genEigs
WAHO= $(HOME)/Dropbox/AMP/waho


# ============= Files for cgWave checks in checkEigen/grids
checkEigen:
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=square -g=square32.order4.ng3 -numEigenValues=64  -tol=1.0e-14 -bc1=d -show=square32O4EigsEv64.show -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=square -g=square64.order2     -numEigenValues=64  -tol=1.0e-14 -bc1=d -show=square64O2EigsEv64.show -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk   -g=sice4.order2        -numEigenValues=128 -tol=1.0e-12 -bc1=d -show=diskG4O2Ev128.show      -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=square -g=box1.order2         -numEigenValues=64 - tol=1.0e-12 -bc1=d -show=boxG1O2EigsEv64.show   -orthogonalize=1  -go=go

# ============= DISK -- runs to find grid refinement convergence rates =============

diskOrder2:
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=sice4.order2  -matlab=diskG4O2Ev32  -table=diskG4O2Ev32Table  -numEigenValues=32 -tol=1.0e-14 -bc1=d -show=junk.show -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=sice8.order2  -matlab=diskG8O2Ev32  -table=diskG8O2Ev32Table  -numEigenValues=32 -tol=1.0e-14 -bc1=d -show=junk.show -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=sice16.order2 -matlab=diskG16O2Ev32 -table=diskG16O2Ev32Table -numEigenValues=32 -tol=1.0e-14 -bc1=d -show=junk.show -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=sice32.order2 -matlab=diskG32O2Ev32 -table=diskG32O2Ev32Table -numEigenValues=32 -tol=1.0e-14 -bc1=d -show=junk.show -orthogonalize=1 -go=go

diskOrder4:
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=sice4.order4.ng3  -matlab=diskG4O4Ev32  -table=diskG4O4Ev32Table  -numEigenValues=32 -tol=1.0e-14 -bc1=d -show=junk.show -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=sice8.order4.ng3  -matlab=diskG8O4Ev32  -table=diskG8O4Ev32Table  -numEigenValues=32 -tol=1.0e-14 -bc1=d -show=junk.show -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=sice16.order4.ng3 -matlab=diskG16O4Ev32 -table=diskG16O4Ev32Table -numEigenValues=32 -tol=1.0e-14 -bc1=d -show=junk.show -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=sice32.order4.ng3 -matlab=diskG32O4Ev32 -table=diskG32O4Ev32Table -numEigenValues=32 -tol=1.0e-14 -bc1=d -show=junk.show -orthogonalize=1 -go=go


# ============= PIPE -- runs to find grid refinement convergence rates =============

pipeOrder2:
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=pipeze2.order2 -matlab=pipeG2O2Ev32 -table=pipeG2O2Ev32Table -numEigenValues=32 -tol=1.0e-12 -bc1=d -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -show=junk.show -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=pipeze4.order2 -matlab=pipeG4O2Ev32 -table=pipeG4O2Ev32Table -numEigenValues=32 -tol=1.0e-12 -bc1=d -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -show=junk.show -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=pipeze8.order2 -matlab=pipeG8O2Ev32 -table=pipeG8O2Ev32Table -numEigenValues=32 -tol=1.0e-12 -bc1=d -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -show=junk.show -orthogonalize=1 -go=go


pipeOrder4:
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=pipeze2.order4.ng3 -matlab=pipeG2O4Ev32 -table=pipeG2O4Ev32Table -numEigenValues=32 -tol=1.0e-12 -bc1=d -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -show=junk.show -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=pipeze4.order4.ng3 -matlab=pipeG4O4Ev32 -table=pipeG4O4Ev32Table -numEigenValues=32 -tol=1.0e-12 -bc1=d -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -show=junk.show -orthogonalize=1 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -eigCase=disk -g=pipeze8.order4.ng3 -matlab=pipeG8O4Ev32 -table=pipeG8O4Ev32Table -numEigenValues=32 -tol=1.0e-12 -bc1=d -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -show=junk.show -orthogonalize=1 -go=go


sphereOrder2:
	$(genEigs) -noplot eigs.cmd -problem=laplace -g=spheree2.order2.hdf -matlab=sphereG2O2Ev32 -table=sphereG2O2Ev32Table -orthogonalize=1 -eigCase=sphere -numEigenValues=32 -tol=1.0e-12 -show=junk -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -g=spheree4.order2.hdf -matlab=sphereG4O2Ev32 -table=sphereG4O2Ev32Table -orthogonalize=1 -eigCase=sphere -numEigenValues=32 -tol=1.0e-12 -show=junk -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -go=go
	$(genEigs) -noplot eigs.cmd -problem=laplace -g=spheree8.order2.hdf -matlab=sphereG8O2Ev32 -table=sphereG8O2Ev32Table -orthogonalize=1 -eigCase=sphere -numEigenValues=32 -tol=1.0e-12 -show=junk -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -go=go

sphereOrder2a:
	$(genEigs) -noplot eigs.cmd -problem=laplace -g=spheree3.order2.hdf -matlab=sphereG3O2Ev32 -table=sphereG3O2Ev32Table -orthogonalize=1 -eigCase=sphere -numEigenValues=32 -tol=1.0e-12 -show=junk -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -go=go > sphereG3O4Ev32
	$(genEigs) -noplot eigs.cmd -problem=laplace -g=spheree6.order2.hdf -matlab=sphereG6O2Ev32 -table=sphereG6O2Ev32Table -orthogonalize=1 -eigCase=sphere -numEigenValues=32 -tol=1.0e-12 -show=junk -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -go=go > sphereG6O4Ev32

sphereOrder4:
	setsid $(genEigs) -noplot eigs.cmd -problem=laplace -g=spheree2.order4.ng3.hdf -matlab=sphereG2O4Ev32 -table=sphereG2O4Ev32Table -orthogonalize=1 -eigCase=sphere -numEigenValues=32 -tol=1.0e-12 -show=junk -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -go=go
	setsid $(genEigs) -noplot eigs.cmd -problem=laplace -g=spheree4.order4.ng3.hdf -matlab=sphereG4O4Ev32 -table=sphereG4O4Ev32Table -orthogonalize=1 -eigCase=sphere -numEigenValues=32 -tol=1.0e-12 -show=junk -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -go=go
	setsid $(genEigs) -noplot eigs.cmd -problem=laplace -g=spheree8.order4.ng3.hdf -matlab=sphereG8O4Ev32 -table=sphereG8O4Ev32Table -orthogonalize=1 -eigCase=sphere -numEigenValues=32 -tol=1.0e-12 -show=junk -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -go=go

sphereOrder4a:
	setsid $(genEigs) -noplot eigs.cmd -problem=laplace -g=spheree3.order4.ng3.hdf -matlab=sphereG3O4Ev32 -table=sphereG3O4Ev32Table -orthogonalize=1 -eigCase=sphere -numEigenValues=32 -tol=1.0e-12 -show=junk -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -go=go > sphereG3O4Ev32.out
	setsid $(genEigs) -noplot eigs.cmd -problem=laplace -g=spheree6.order4.ng3.hdf -matlab=sphereG6O4Ev32 -table=sphereG6O4Ev32Table -orthogonalize=1 -eigCase=sphere -numEigenValues=32 -tol=1.0e-12 -show=junk -st_ksp_type gmres -st_pc_type ilu -st_pc_factor_levels 2 -st_ksp_rtol 1e-12 -go=go > sphereG6O4Ev32.out
