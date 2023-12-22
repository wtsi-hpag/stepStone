#!/bin/tcsh

function plotcmd
{
printf "set terminal svg\n"
printf "set style line 1 lc rgb 'blue' pt 6\n"
printf "set style line 2 lc rgb 'red' pt 5\n"
printf "set xlabel \"Copy Numbers\"\n"
printf "set ylabel \"Standard Error - SDV\"\n"
printf "plot [ 1.5 to  2.5 ] [ 0.001 to 0.6 ]  \"pca-deveil-1.dat\" title \"Devil Cancer DFT1 \" with points ls 1,\"pca-deveil-2.dat\" title \"Devil Cancer DFT2 \" with points ls 2"
}
plotcmd | gnuplot > data-chr6.svg
inkscape -z --export-text-to-path --export-pdf data-chr6.pdf data-chr6.svg
gs -r600 -dNOPAUSE -dBATCH -sDEVICE=png256 -sOutputFile=human_chr.png data-chr6.pdf


