set title 'Convergence'
set xlabel 'Iteration'
set ylabel 'Residual'
#set logscale xy
set logscale y
plot 'testheat.res' using 1:6 title 'T' with lines
pause -1
