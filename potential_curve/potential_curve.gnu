#!/usr/bin/env gnuplot


# scope

  @png

  set output "curve.png"

  curve = "curve.dat"

# default settings

  set loadpath "~/.gnuplot.d/config" "~/.gnuplot.d/palettes"
  load "all.cfg"
  load "spectral.pal"

# variables / functions

# fit

# stats

  stats curve nooutput
  n = STATS_columns

# settings

  set xrange [0:*]
  set xlabel "radial separation"

  set yrange [*:*]
  set ylabel "HF energy"

  set key top right

# plots

  set title "HF Energy of Diatomic Helium (++) against Radial Separation"
  plot for [i = 2:n] curve u 1:i ls i lw 2 t sprintf("%i", i - 1)
