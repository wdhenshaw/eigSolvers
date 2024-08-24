#
# plotStuff plotEig.cmd -show=diskEigsEve32.show -name=disk -numToSave=16
# plotStuff plotEig.cmd -show=ellipseEve512.show -name=ellipse -numToSave=16
# plotStuff plotEig.cmd -show=shapes64.show -name=shapes -numToSave=16
# plotStuff plotEig.cmd -show=cicEigs128.show -name=cic -numToSave=16
# plotStuff plotEig.cmd -show=rpiEigs128.show -name=rpi -numToSave=16
# plotStuff plotEig.cmd -show=rpiEigs1024.show -name=rpiSelected -start=64 -stride=64 -numToSave=1024 -lines=0
#
$show="diskEigsEve32.show"; $dsf=.1; $cf=1; 
$numToSave=4; $start=0; $stride=1; $lines=1; 
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"dsf=f"=>\$dsf,"cf=i"=>\$cf,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"start=i"=>\$start,"stride=i"=>\$stride,\
      "field=s"=>\$field,"vmin=f"=>\$vmin,"vmax=f"=>\$vmax,"lines=i"=>\$lines );
#
$show
contour
  coarsening factor 1
  if( $lines eq "0" ){ $cmd="plot contour lines (toggle)"; }else{ $cmd="#"; }
  $cmd 
exit
# 
DISPLAY AXES:0 0
x-:0
plot:psi0
pause
#
# -- save plots ---
#
$cmd = "";
for( $i=$start; $i < $numToSave; $i=$i+$stride ){\
  $cmd .= "plot:psi$i\n"; \
  $plotName = $name . "EigenVector$i.ps"; \
  $cmd .= "hardcopy file name:0 $plotName\n"; \
  $cmd .= "hardcopy save:0\n"; \
}
$cmd .="#";
# execute commands 
$cmd


#
DISPLAY AXES:0 0
# DISPLAY LABELS:0 0
# DISPLAY COLOUR BAR:0 0
DISPLAY SQUARES:0 0
#  
#
$plotName = $name . "Displacement.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0  
#
pause
#
erase

plot:u1
$plotName = $name . "u1.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0 
#
plot:u2
$plotName = $name . "u2.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0  
#
plot:p
$plotName = $name . "p.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0 
#
plot:psi
$plotName = $name . "psi.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0  
