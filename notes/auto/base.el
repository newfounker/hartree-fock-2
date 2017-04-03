(TeX-add-style-hook
 "base"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("nag" "l2tabu" "orthodox") ("babel" "english") ("inputenc" "utf8x") ("caption" "font={small, it}") ("todo" "marginpar") ("enumitem" "shortlabels")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "nag"
    "babel"
    "inputenc"
    "microtype"
    "xifthen"
    "etoolbox"
    "amsmath"
    "amssymb"
    "amsfonts"
    "mathtools"
    "physics"
    "siunitx"
    "tensor"
    "geometry"
    "caption"
    "pgfplots"
    "tikz"
    "graphicx"
    "titling"
    "booktabs"
    "fancyhdr"
    "hyperref"
    "todo"
    "framed"
    "enumitem"
    "verbatim")
   (TeX-add-symbols
    "integer"
    "rational"
    "complex")
   (LaTeX-add-mathtools-DeclarePairedDelimiters
    '("lr" "")
    '("lrsq" "")
    '("lrset" "")
    '("lrang" "")
    '("lrabs" "")
    '("lrnorm" "")))
 :latex)

