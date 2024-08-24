Getopt::Long::Configure("prefix_pattern=(--square|--|-)"); 
GetOptions("N=i"=>\$N, "bc=s"=>\$bc, "saveGrid=i"=>\$saveGrid,"grid_order=i"=>\$grid_order); 
if ( !$N ) { $N=10;}; 
if ( !$bc ) { $bc="notPeriodic";}; 
if (!$grid_order) { $grid_order =2;}; 
$order_str = $grid_order==2 ? "second" : ($grid_order==4 ? "fourth" : ( $grid_order==6 ? "sixth" : ( $grid_order==8 ?  "eighth" : "unknown"))); 
$order_str = $order_str." order"; 
$Np = $N; 
$name = "square$N".($bc eq "periodic" ? "p" : "").".hdf"; 
$bc_l = $bc eq "periodic" ? "periodicity\n1 1\n" : "show parameters"; 
$finish = "exit"; 
$Ng = int($grid_order/2); 
$Nge= $Ng+1; 
if ( $saveGrid ) { $finish = "save an overlapping grid\n$name\nsquare\nexit"; } 
* 
* make a simple square (periodic BCs) 
* 
create mappings 
rectangle 
lines 
$Np $Np 
$bc_l 
exit 
