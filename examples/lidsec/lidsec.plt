set title 'Convergence'
set xlabel 'Iteration'
set ylabel 'Residual'
#set logscale xy
set logscale y
plot 'lidsec.res' using 1:2 title 'u' with lines, 'lidsec.res' using 1:3 title 'v' with lines, 'lidsec.res' using 1:4 title 'w' with lines, 'lidsec.res' using 1:5 title 'p' with lines, 'lidsec.res' using 1:6 title 'T' with lines, 'lidsec.res' using 1:7 title 's' with lines
pause -1
