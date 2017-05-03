#!/usr/bin/env gnuplot


# scope

# replacing: @png
  set term pngcairo size 1600, 875 dashed enhanced font "Fira mono, 9"
# done

  set output "plots.png"

  wavefunctions = "core_plots.dat"
  coulomb = "core_coulomb.dat"
  potential = "core_potential.dat"

# default settings

# replacing:
  # set loadpath "~/.gnuplot.d/config" "~/.gnuplot.d/palettes"
  # load "all.cfg"

# border
  set border 31 front lc rgb '#000000' lt 1 lw 0.5

# grid
  set grid front lc rgb '#000000' lt 0 lw 1

# key
  set key top right box lc rgb '#000000' lt 1

# tics
  set format "%g"
  set tics nomirror out scale 0.75

# style
  set style textbox
  set style data lines
  set style function lines
#done

# replacing: load "spectral.pal"
# line styles
  set style line 1 lt 1 lc rgb '#D53E4F' # red
  set style line 2 lt 1 lc rgb '#F46D43' # orange
  set style line 3 lt 1 lc rgb '#FDAE61' # pale orange
  set style line 4 lt 1 lc rgb '#FEE08B' # pale yellow-orange
  set style line 5 lt 1 lc rgb '#E6F598' # pale yellow-green
  set style line 6 lt 1 lc rgb '#ABDDA4' # pale green
  set style line 7 lt 1 lc rgb '#66C2A5' # green
  set style line 8 lt 1 lc rgb '#3288BD' # blue

# palette
  set palette defined ( 0 '#D53E4F',\
      	    	      1 '#F46D43',\
  		      2 '#FDAE61',\
  		      3 '#FEE08B',\
  		      4 '#E6F598',\
  		      5 '#ABDDA4',\
  		      6 '#66C2A5',\
  		      7 '#3288BD' )
# done

# variables / functions

  J(r) = 1 / r

# fit

# stats

# replacing:
  # stats wavefunctions nooutput
  # n = STATS_columns
  n = 20
#done

# settings

  set xrange [0:5]
  set xlabel "r"
  set key center right

# plots

  set multiplot layout 2, 2 title sprintf("Partial Waves for {/Symbol l} = 0, ..., %i", n-2)

  set title "{{/Symbol Y}_{/Symbol l}}(r)"
  set yrange [*:*]
  plot \
    for [i = 2:n] wavefunctions u 1:i ls i lw 2 t sprintf("%i", i-2)

  set title "{|{/Symbol Y}_{/Symbol l}|}^{2}(r)"
  set yrange [*:*]
  plot \
    for [i = 2:n] wavefunctions u 1:(column(i) ** 2) ls i lw 2 t sprintf("%i", i-2)

  set title "{J_{/Symbol l}}(r)"
  set yrange [0:3]
  plot \
    for [i = 2:n] coulomb u 1:i ls i lw 2 t sprintf("%i", i-2), \
    coulomb u 1:(2*J(column(1))) lc black lw 2 t sprintf("2/r")

  set title "{V_{/Symbol l}}(r)"
  set yrange [-3:0]
  plot \
    for [i = 2:n] potential u 1:i ls i lw 2 t sprintf("%i", i-2), \
    coulomb u 1:(-2*J(column(1))) lc black lw 2 t sprintf("-2/r")

  unset multiplot