(TeX-add-style-hook
 "hartree_fock_notes"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "draft")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "base")
   (TeX-add-symbols
    '("irrep" ["argument"] 0)
    '("harmonic" ["argument"] 0)
    '("laguerre" ["argument"] 0)
    '("sturmian" ["argument"] 0)
    '("electron" ["argument"] 0)
    '("orbitalpw" ["argument"] 0)
    '("orbital" ["argument"] 0)
    '("orbitalspin" ["argument"] 0)
    '("molecule" ["argument"] 0)
    "spin"
    "parity"
    "proj")
   (LaTeX-add-labels
    "sec:frozen-core-hartree"
    "sec:coulomb"
    "sec:exchange"))
 :latex)

