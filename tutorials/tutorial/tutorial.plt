#set terminal postscript eps enhanced color            # set output to eps file
#set output 'cavity.eps' 
#set title 'Convergence'
set xlabel 'Iteration'
set ylabel 'Residual'
#set logscale xy
set logscale y
plot 'tutorial.res' using 1:2 title 'u' with lines, \
     'tutorial.res' using 1:3 title 'v' with lines, \
     'tutorial.res' using 1:5 title 'p' with lines
pause -1
