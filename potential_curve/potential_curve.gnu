#!/usr/bin/env gnuplot


# scope

  @png

  set output "potential_curve.png"

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
  plot for [i = 2:n] curve u 1:i ls i lw 2 t sprintf("%i", i-2)

  # set multiplot layout 2, 2 title sprintf("Partial Waves for {/Symbol l} = 0, ..., %i", n-2)

  # set title "{{/Symbol Y}_{/Symbol l}}(r)"
  # set yrange [*:*]
  # plot \
  #   for [i = 2:n] wavefunctions u 1:i ls i lw 2 t sprintf("%i", i-2)

  # set title "{|{/Symbol Y}_{/Symbol l}|}^{2}(r)"
  # set yrange [*:*]
  # plot \
  #   for [i = 2:n] wavefunctions u 1:(column(i) ** 2) ls i lw 2 t sprintf("%i", i-2)

  # set title "{J_{/Symbol l}}(r)"
  # set yrange [0:3]
  # plot \
  #   for [i = 2:n] coulomb u 1:i ls i lw 2 t sprintf("%i", i-2), \
  #   coulomb u 1:(2*J(column(1))) lc black lw 2 t sprintf("2/r")

  # set title "{V_{/Symbol l}}(r)"
  # set yrange [-3:0]
  # plot \
  #   for [i = 2:n] potential u 1:i ls i lw 2 t sprintf("%i", i-2), \
  #   coulomb u 1:(-2*J(column(1))) lc black lw 2 t sprintf("-2/r")

  # unset multiplot