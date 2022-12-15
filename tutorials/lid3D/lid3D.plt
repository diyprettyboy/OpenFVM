set title 'Convergence'
set xlabel 'Iteration'
set ylabel 'Residual'
#set logscale xy
set logscale y
plot 'lid3D.res' using 1:2 title 'u' with lines, 'lid3D.res' using 1:3 title 'v' with lines, 'lid3D.res' using 1:5 title 'p' with lines
#plot 'lid3D.res' using 1:2 title 'u' with lines, 'lid3D.res' using 1:3 title 'v' with lines, 'lid3D.res' using 1:4 title 'w' with lines, 'lid3D.res' using 1:5 title 'p' with lines, 'lid3D.res' using 1:6 title 'T' with lines, 'lid3D.res' using 1:7 title 's' with lines
pause -1
