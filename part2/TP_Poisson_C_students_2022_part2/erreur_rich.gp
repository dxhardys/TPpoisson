set terminal png size 800, 600 
set output "plot/erreur_rich.png"

set xlabel "Data size"
set ylabel "Forward error"
set style data line
plot    "data/Richardson/erreur_avant_rich.dat" using 1 linecolor rgb "#00A5E3" title "Richardson",\

