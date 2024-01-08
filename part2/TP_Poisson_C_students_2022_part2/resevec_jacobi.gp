set terminal png size 800, 600 
set output "plot/resvec_jacobi.png"

set xlabel "Residual norm"
set ylabel "Number of iteration"
set style data line
plot  "data/Jacobi/jac_01.dat" linecolor rgb "#FF96C5"
