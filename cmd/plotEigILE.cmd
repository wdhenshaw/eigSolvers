#
# plotStuff plotEig.cmd -show=eigILE40Order2.show -name=squareEigILE40
# plotStuff plotEig.cmd -show=eigILE40Eig3Order4.show -name=squareEigILE40Eig3 -dsf=.075 
# 
# plotStuff plotEig.cmd -show=square41TDDD.show -name=square41TDDD
#
# plotStuff plotEig.cmd -show=eigILESquare250BCtddd.show -name=eigILESquare250BCtddd -cf=5
#
$show="eigILE40Order2.show"; $dsf=.1; $cf=1; 
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"dsf=f"=>\$dsf,"cf=i"=>\$cf,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field,"vmin=f"=>\$vmin,"vmax=f"=>\$vmax );
#
$show
#
displacement
  displacement scale factor $dsf
  coarsening factor $cf
  exit this menu
displacement
  exit this menu
#
pause
x-
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
contour
  coarsening factor 1
  exit
# 
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
