set terminal postscript eps enhanced color size font 'Times-Bold,24'
set output 'test.eps'
set pm3d
set view map
unset colorbox
set xlabel "Distance [m]" offset 0.,1.5
set xrange [-20:220]
set ylabel "Depth [m]" offset 2.,0.
set yrange [220:-20]

set size square
#set size 3,3
set ydim=5
set xdim=5
    
set xtics nomirror out offset 0.,1.0
set ytics nomirror out offset 1.0,0.0


set palette grey 
splot 'snapz00801' binary format="%float" origin=(-20,-20,0.) dx=0.5 dy=0.5 array=(481,481) with pm3d notitle