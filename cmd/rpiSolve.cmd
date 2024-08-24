read eigenvectors
process eigenvectors
fit function
plot fit
  plot contour lines (toggle)
  plot:u0err
  plot:u1Fit
  plot:u1err
  exit
solve wave
set time 1
set time 1.5
contour
  min max -.2 1.3
  plot contour lines (toggle)
  set view:0 0.0241692 0.111782 0 2.16458 1 0 0 0 1 0 0 0 1
  reset min max
  reset:0
  plot contour lines (toggle)
  exit
set time 2
done
beta, x0, y0, z0, k0: 20, 0, -2, 0, 2
fit function
plot fit
  x-r:0
  y-r:0
  y-r:0
  y-r:0
  y-r:0
  y-r:0
  y-r:0
  reset:0
  y-r:0
  y-r:0
  y-r:0
  x-r:0
  x-r:0
  x-r:0
  x+r:0
  exit
solve wave
set time 1.5
reset:0
x-r:0
x-r:0
x-r:0
x-r:0
x-r:0
x-r:0
x-r:0
y-r:0
y-r:0
y-r:0
x+r:0
x+r:0
contour
  vertical scale factor 0.25
  reset:0
  exit
done
beta, x0, y0, z0, k0: 20, 0, -2, 0, 4
fit function
solve wave
set time 1.5
contour
  min max -1.4 2.
  min max -1.4 1.5
  min max -1. 1.5
  exit
times to plot 0.1
next
next
previous
previous
contour
  min max -1 1.4
  min max -.8 1.4
  min max -0.6 1.4
  min max -0.5 1.4
  min max -0.7 1.4
  min max -0.7 1.2
  min max -0.7 1.
  min max -0.7 1.5
  min max -0.7 2
  min max -0.7 1.8
  min max -0.9 1.8
  min max -1 1.8
  min max -1 1.4
  hardcopy file name:0 rpiScat4096t1p5.ps
  hardcopy save:0
  exit
set time 2
hardcopy file name:0 rpiScat4096t2p0.ps
hardcopy save:0
contour
  min max -1 1.
  min max -1 1.1
