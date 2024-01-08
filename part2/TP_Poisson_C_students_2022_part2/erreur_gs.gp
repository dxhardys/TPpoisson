set terminal png size 800, 600 
set output "plot/erreur_gs.png"

set xlabel "Data size"
set ylabel "Forward error"
set style data line
plot    "data/GS/erreur_avant_gs.dat" using 1 linecolor rgb "#00A5E3" title "GS",\

