set terminal png size 1000, 1000 
set output "plot/direct.png"

set xlabel "Data size"
set ylabel "Time (s)"
set style data line
plot    "complx/direct.dat" using 1 title 'dgbtrftridiag + dgbtrs' linecolor rgb "#00A5E3",\
        "complx/direct.dat" using 2 title 'dgbtrf + dgbtrs' linecolor rgb "#8DD7BF",\
        "complx/direct.dat" using 3 title 'dgbsv' linecolor rgb "#FF96C5"
