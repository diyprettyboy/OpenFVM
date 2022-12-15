set title 'Convergence'
set xlabel 'Iteration'
set ylabel 'Residual'
#set logscale xy
set logscale y
plot 'lidteta.res' using 1:2 title 'u' with lines, 'lidteta.res' using 1:3 title 'v' with lines, 'lidteta.res' using 1:5 title 'p' with lines
#plot 'lidteta.res' using 1:6 title 'T' with lines
pause -1
#pause 10
#reread 
