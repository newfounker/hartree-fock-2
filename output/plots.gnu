#!/usr/bin/env gnuplot


# scope

  @png
  set output "plots.png"

  wavefunctions = "core_plots.dat"
  coulomb = "core_coulomb.dat"
  potential = "core_potential.dat"

# default settings

  set loadpath "~/.gnuplot.d/config" "~/.gnuplot.d/palettes"
  load "all.cfg"
  load "spectral.pal"

# variables / functions

  J(r) = 1 / r

# fit

# stats

  stats wavefunctions nooutput
  n = STATS_columns

# settings

  set xrange [0:5]
  set xlabel "r"

# plots

  set multiplot layout 2, 2 title sprintf("Partial Waves for {/Symbol l} = 0, ..., %i", n-2)

  set title "{{/Symbol Y}_{/Symbol l}}(r)"
  set key top right
  set yrange [*:*]
  plot \
    for [i = 2:n] wavefunctions u 1:i ls i lw 2 t sprintf("%i", i-2)

  set title "{|{/Symbol Y}_{/Symbol l}|}^{2}(r)"
  set key top right
  set yrange [*:*]
  plot \
    for [i = 2:n] wavefunctions u 1:(column(i) ** 2) ls i lw 2 t sprintf("%i", i-2)

  set title "{J_{/Symbol l}}(r)"
  set key top right
  set yrange [0:3]
  plot \
    for [i = 2:n] coulomb u 1:i ls i lw 2 t sprintf("%i", i-2), \
    coulomb u 1:(2*J(column(1))) lc black lw 2 t sprintf("2/r")

  set title "{V_{/Symbol l}}(r)"
  set key bottom right
  set yrange [-3:0]
  plot \
    for [i = 2:n] potential u 1:i ls i lw 2 t sprintf("%i", i-2), \
    coulomb u 1:(-2*J(column(1))) lc black lw 2 t sprintf("-2/r")

  unset multiplot