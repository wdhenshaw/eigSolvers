#
#  Compute eigenvalues -- incompressible elasticity
#
#   comp compEigs.cmd -logFile=squareEigsILEorder2.log
#   comp compEigs.cmd -logFile=squareEigsILEorder4.log -order=4
# 
#
$solution=1; $name="compEigs"; $logFile="compEigs.log"; $order=2; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"solution=i"=>\$solution,"logFile=s"=>\$logFile,"order=i"=>\$order  );
# 
specify files
  if( $order eq 2 ){ $files="eigILE40Order2.show\n eigILE80Order2.show\n  eigILE160Order2.show"; }\
              else{  $files="eigILE40Order4.show\n eigILE80Order4.show\n  eigILE160Order4.show"; }
  $files 
exit
#
output file name: $logFile
#
choose a solution
  $solution
#
compute errors