set title 'Convergence'
set xlabel 'Iteration'
set ylabel 'Residual'
set logscale xy
plot 'backstep.res' using 1:2 title 'u' with lines, 'backstep.res' using 1:3 title 'v' with lines,  'backstep.res' using 1:5 title 'p' with lines
pause 10
reread
