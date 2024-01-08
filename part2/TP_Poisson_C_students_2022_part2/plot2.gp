set terminal png size 1000, 1000 
set output "plot/dgbtrf.png"

set xlabel "Data size"
set ylabel "Time (s)"
set style data line
plot    "complx/dgbtrf.dat" using 1 title 'dgbtrftridiag' linecolor rgb "#00A5E3",\
        "complx/dgbtrf.dat" using 2 title 'dgbtrf' linecolor rgb "#8DD7BF"
