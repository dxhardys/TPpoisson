set terminal png enhanced font "arial,15" fontscale 1.0 size 1000, 1000 
set output "graphs/direct/direct_graphd.png"

set xlabel "Taille des données"
set ylabel "Temps en secondes"
set style data line
plot    "complx/direct_graphd.dat" using 1 title 'dgbtrftridiag + dgbtrs' linecolor rgb "#00A5E3",\
        "complx/direct_graphd.dat" using 2 title 'dgbtrf + dgbtrs' linecolor rgb "#8DD7BF",\
        "complx/direct_graphd.dat" using 3 title 'dgbsv' linecolor rgb "#FF96C5"
