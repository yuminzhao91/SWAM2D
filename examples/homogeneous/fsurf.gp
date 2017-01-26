set terminal postscript portrait enhanced color font 'Times-Roman,11'
set output 'validation_fsurf.ps'
unset colorbox

set xlabel "Distance [m]" offset 0.,0.75
set xrange [-20:220]
set ylabel "Depth [m]" offset 2.5,0.
set yrange [220:-20]

set cbrange [-5.e-11:5.0e-11]
set size ratio -1

set xtics 40 out nomirror offset 0.,0.5
set ytics 40 out nomirror offset 0.75,0.

set palette grey

fb(x) = (x >= 0. && x <= 200.) ? 200. : 1/0
ft(x) = (x >= 0. && x <= 200.) ? 0. : 1/0
fr(z) = (z >= 0. && z <= 200.) ? 200. : 1/0
fl(z) = (z >= 0. && z <= 200.) ? 0. : 1/0

# >> SETTING MULTIPLOT (3,2)
set multiplot layout 3,2 spacing 0.15,0.0
set linestyle 1 lt 0 lw 2

# >> FIGURE (1,1)
set title "vaccum  approach" offset 0.0,-0.5 font "Times-Bold,11"
set xlabel " "
set label 1 "t=0.2s" tc rgb "black" font "Times-Bold,11" front at 10.,190.
plot 'validation_fsurf0/snapz00801' binary format="%float" dx=0.5 dy=0.5 center=(100.,100.) array=481x481 with image notitle, \
     'validation_fsurf0/pml.txt' using 1:2 with line ls 1 notitle 

# >> FIGURE (1,2)
set title "image theory approach" offset 0.0,-0.5 font "Times-Bold,11"
set xlabel " "
set label 1 "t=0.2s" tc rgb "black" font "Times-Bold,11" front at 10.,190.
plot 'validation_fsurf1/snapz00801' binary format="%float" dx=0.5 dy=0.5 center=(100.,100.) array=481x481 with image notitle, \
     'validation_fsurf0/pml.txt' using 1:2 with line ls 1 notitle 

# >> FIGURE (2,1)
set title " " offset -16.7,-0.5
set xlabel " "
set label 1 "t=0.3s" tc rgb "black" font "Times-Bold,11" front at 10.,190.
plot 'validation_fsurf0/snapz01201' binary format="%float" dx=0.5 dy=0.5 center=(100.,100.) array=481x481 with image notitle, \
     'validation_fsurf0/pml.txt' using 1:2 with line ls 1 notitle 

# >> FIGURE (2,2)
set title " " offset -16.7,-0.5
set xlabel " "
set label 1 "t=0.3s" tc rgb "black" font "Times-Bold,11" front at 10.,190.
plot 'validation_fsurf1/snapz01201' binary format="%float" dx=0.5 dy=0.5 center=(100.,100.) array=481x481 with image notitle, \
     'validation_fsurf0/pml.txt' using 1:2 with line ls 1 notitle 

# >> FIGURE (3,1)
set title " " offset -16.7,-0.5
set xlabel "Distance [m]" offset 0.,0.75
set label 1 "t=0.4s" tc rgb "black" font "Times-Bold,11" front at 10.,190.
plot 'validation_fsurf0/snapz01601' binary format="%float" dx=0.5 dy=0.5 center=(100.,100.) array=481x481 with image notitle, \
     'validation_fsurf0/pml.txt' using 1:2 with line ls 1 notitle 

# >> FIGURE (3,2)
set title " " offset -16.7,-0.5
set xlabel "Distance [m]" offset 0.,0.75
set label 1 "t=0.4s" tc rgb "black" font "Times-Bold,11" front at 10.,190.
plot 'validation_fsurf1/snapz01601' binary format="%float" dx=0.5 dy=0.5 center=(100.,100.) array=481x481 with image notitle, \
     'validation_fsurf0/pml.txt' using 1:2 with line ls 1 notitle 