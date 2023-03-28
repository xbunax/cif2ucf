#!/bin/bash
cat $1 | tail -n +$2 > text
cut -f $3,$4 text  > data
rm text
gnuplot -persist << EOF
set term png
set output "data.png"
plot "data" with lines
replot
set output
EOF
rm data 
