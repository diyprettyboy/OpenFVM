set title 'Convergence'
set xlabel 'Iteration'
set ylabel 'Residual'
#set logscale xy
set logscale y
plot 'still.res' using 1:2 title 'u' with lines, 'still.res' using 1:3 title 'v' with lines, 'still.res' using 1:5 title 'p' with lines, 'still.res' using 1:7 title 's' with lines
pause -1
#pause 10
#reread 
