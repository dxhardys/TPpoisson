set terminal png size 800, 600 
set output "plot/erreur_jac.png"

set xlabel "Data size"
set ylabel "Forward error"
set style data line
plot    "data/Jacobi/erreur_avant_jac.dat" using 1 linecolor rgb "#00A5E3" title "Jacobi",\

