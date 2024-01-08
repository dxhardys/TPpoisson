set terminal png size 800, 600 
set output "plot/resvec_gs.png"

set xlabel "Residual norm"
set ylabel "Number of iteration"
set style data line
plot  "data/GS/gs_01.dat" linecolor rgb "#FF96C5"
