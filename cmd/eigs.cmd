#
#  Command file for genEigs
# Usage: 
#     genEigs [-noplot] eigs.cmd -show=<s> -matlab=[0|1] -eigenSolver=[KrylovSchur|ARPACK|LAPACK] -go=[go|read] <command line options for genEigs> 
#
#    NOTE: Most command line arguments are handled in genEigs.C
#
$show=""; $go="og"; $matlab=1; $coarseToFine=0; $matlabFileName="genEigs"; 
$eigenSolver="KrylovSchur";
GetOptions( "show=s"=>\$show,"coarseToFine=i"=>\$coarseToFine,"matlab=i"=>\$matlab,"matlabFileName=s"=>\$matlabFileName,\
            "eigenSolver=s"=>\$eigenSolver,"go=s"=>\$go );
#
$eigenSolver
#
if( $go eq "read" ){ $cmd = "read eigenvectors from a file\n orthogonalize"; }else{ $cmd = "compute"; }
if( $coarseToFine ){ $cmd = "read eigenvectors from a file\n coarse to fine\n orthogonalize"; }
$cmd 
# compute
# open graphics
if( $show ne "" ){ $cmd="save show file"; }else{ $cmd="#"; }
$cmd
#
if( $matlab eq "1" ){ $cmd="save matlab file"; }else{ $cmd="#"; }
$cmd
$cmd="#"; 
if( $go eq "og" ){ $cmd = "open graphics"; }
if( $go eq "go" ){ $cmd = "exit"; }
$cmd

check

exit
exit
exit
open graphics

