set terminal png size 800, 600 
set output "plot/resvec_richardson.png"

set ylabel "Residual norm"
set xlabel "Number of iteration"
set style data line
plot  "data/Richardson/rich_01.dat" linecolor rgb "#FF96C5"
