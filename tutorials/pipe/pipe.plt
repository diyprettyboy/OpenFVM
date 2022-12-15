set terminal postscript eps            # set output to eps file
set output 'pipe.eps' 
set title 'Convergence'
set xlabel 'Iteration'
set ylabel 'Residual'
#set logscale xy
set logscale y
plot 'pipe.res' using 1:6 title 'Temperature' with lines
pause -1

