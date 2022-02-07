#!/usr/bin/gnuplot

bin(x,w) = w*floor(x/w)

set style fill solid 0.3
set terminal png
set size square
unset key

filename = ARG1

titulo = filename
set title titulo
figname = sprintf("%s.png", filename)
set output figname
plot filename u (bin($3, 1)):(1.0) smooth freq w boxes lt 22
