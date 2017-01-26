set terminal postscript portrait enhanced color font 'Times-Roman,11'
set output 'validation_courant.ps'
unset colorbox

set xlabel "Distance [m]" offset 0.,0.75
set xrange [-20:220]
set ylabel "Depth [m]" offset 2.5,0.
set yrange [220:-20]

set cbrange [-5.e-11:5.0e-11]
set size ratio -1

set xtics 40 out nomirror offset 0.,0.5
set ytics out nomirror offset 0.75,0.

set palette grey

# >> SETTING MULTIPLOT (3,2)
set multiplot layout 3,2 spacing 0.15,0.0


# >> FIGURE (1,1)
set title "grid points per wavelength > 4" offset 0.0,-0.5 font "Times-Bold,11"
set xlabel " "
set label 1 "t=0.2s" tc rgb "black" font "Times-Bold,11" front at 10.,190.
plot 'validation_courant1/snapz00801' binary format="%float" dx=0.5 dy=0.5 center=(100.,100.) array=481x481 with image notitle

# >> FIGURE (1,2)
set title "grid points per wavelength < 4" offset 0.0,-0.5 font "Times-Bold,11"
set xlabel " "
set label 1 "t=0.2s" tc rgb "black" font "Times-Bold,11" front at 10.,190.
plot 'validation_courant0/snapz00801' binary format="%float" dx=2.0 dy=2.0 center=(100.,100.) array=121x121 with image notitle

# >> FIGURE (2,1)
set title " " offset -16.7,-0.5
set xlabel " "
set label 1 "t=0.3s" tc rgb "black" font "Times-Bold,11" front at 10.,190.
plot 'validation_courant1/snapz01201' binary format="%float" dx=0.5 dy=0.5 center=(100.,100.) array=481x481 with image notitle

# >> FIGURE (2,2)
set title " " offset -16.7,-0.5
set xlabel " "
set label 1 "t=0.3s" tc rgb "black" font "Times-Bold,11" front at 10.,190.
plot 'validation_courant0/snapz01201' binary format="%float" dx=2.0 dy=2.0 center=(100.,100.) array=121x121 with image notitle

# >> FIGURE (3,1)
set title " " offset -16.7,-0.5
set xlabel "Distance [m]" offset 0.,0.75
set label 1 "t=0.4s" tc rgb "black" font "Times-Bold,11" front at 10.,190.
plot 'validation_courant1/snapz01601' binary format="%float" dx=0.5 dy=0.5 center=(100.,100.) array=481x481 with image notitle

# >> FIGURE (3,2)
set title " " offset -16.7,-0.5
set xlabel "Distance [m]" offset 0.,0.75
set label 1 "t=0.4s" tc rgb "black" font "Times-Bold,11" front at 10.,190.
plot 'validation_courant0/snapz01601' binary format="%float" dx=2.0 dy=2.0 center=(100.,100.) array=121x121 with image notitle