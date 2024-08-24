#
#  Command file for genEigs
#    Most command line arguments are handled in genEigs.C
#
$show=""; $go="og"; 
GetOptions( "show=s"=>\$show,"go=s"=>\$go );
#
compute
# open graphics
if( $show ne "" ){ $cmd="save show file"; }else{ $cmd="#"; }
$cmd
#
save matlab file
$cmd="#"; 
if( $go eq "og" ){ $cmd = "open graphics"; }
if( $go eq "go" ){ $cmd = "exit"; }
$cmd

check

exit
exit
exit
open graphics

