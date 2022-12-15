#set terminal postscript eps            # set output to eps file
#set output 'tubo.eps' 
set title 'Convergence'
set xlabel 'Iteration'
set ylabel 'Residual'
#set logscale xy
set logscale y
plot 'tubo.res' using 1:6 title 'Temperature' with lines
pause -1

