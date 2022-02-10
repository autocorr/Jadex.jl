# Jadex

[![Docs]()](https://jadex.readthedocs.io)
[![License](https://img.shields.io/badge/License-GPL-blue?style=flat)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-replace-blue?style=flat)](https://doi.org/)

Port of the Fortran 77 code [RADEX](https://personal.sron.nl/~vdtak/radex/index.shtml) (van der Tak et al. [2007, A&A 468, 627](https://www.aanda.org/articles/aa/abs/2007/23/aa6820-06/aa6820-06.html)) to the [Julia](https://julialang.org/) programming language. A Python wrapper is provided using `PyJulia`. Distinguishing features of this implementation include:

  * Improved performance (~150x).
  * (*in progress*) Parameter grid evaluation with parallel execution.
  * User definable escape probability and background radiation field.
  * Test and validation suite.
  * Extensible design.

For cases where the same input parameters are used, results from Jadex are expected to match RADEX within five significant figures. These differences arise in-part from the use of higher precision mathematical constants and general numerical instability for levels with very small populations.


## Installation
To install Jadex, open an [interactive Julia session](https://docs.julialang.org/en/v1/manual/getting-started/), press the `]` key to enter the package management mode, and execute the command `add Jadex`. To execute the test suite, run `test Jadex` from package mode.

To use the Python wrapper, first install Jadex per the above instruction and then follow the [PyJulia installation instructions](https://pyjulia.readthedocs.io/en/latest/installation.html). Jadex can then be imported from Python by calling `from julia import Jadex`.

For validation purposes, optional compilation instructions are included in `src/wrap_slatec.jl` for compiling and linking the `slatec.f` Fortran file from RADEX into a shared library. The resulting `libslatec.so` is then wrapped and can be called to factor the rate matrix and solve for the level populations.


## Documentation
Please refer to the online https://jadex.readthedocs.io for the [Quickstart]() guide, [User Manual](), and [API]() reference. The documentation source files are also supplied in the `docs/` folder distributed with Jadex.


## Citing this work
If you use Jadex in an academic work, we ask that you cite the following references, including the original publication for RADEX (van der Tak et al. 2007):

```latex
@ARTICLE{2007A&A...468..627V,
       author = {{van der Tak}, F.~F.~S. and {Black}, J.~H. and {Sch{\"o}ier}, F.~L. and {Jansen}, D.~J. and {van Dishoeck}, E.~F.},
        title = "{A computer program for fast non-LTE analysis of interstellar line spectra. With diagnostic plots to interpret observed line intensity ratios}",
      journal = {\aap},
     keywords = {radiative transfer, methods: numerical, radio lines: ISM, infrared: ISM, submillimeter, Astrophysics},
         year = 2007,
        month = jun,
       volume = {468},
       number = {2},
        pages = {627-635},
          doi = {10.1051/0004-6361:20066820},
archivePrefix = {arXiv},
       eprint = {0704.0155},
 primaryClass = {astro-ph},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2007A&A...468..627V},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

If one uses the collision rate files from the [LAMDA](https://home.strw.leidenuniv.nl/~moldata/) database, the following citation should be included in addition to the source references listed on the page for the specie(s) used.

```latex
@ARTICLE{2005A&A...432..369S,
       author = {{Sch{\"o}ier}, F.~L. and {van der Tak}, F.~F.~S. and {van Dishoeck}, E.~F. and {Black}, J.~H.},
        title = "{An atomic and molecular database for analysis of submillimetre line observations}",
      journal = {\aap},
     keywords = {astronomical data bases: miscellaneous, atomic data, molecular data, radiative transfer, ISM: atoms, ISM: molecules, Astrophysics},
         year = 2005,
        month = mar,
       volume = {432},
       number = {1},
        pages = {369-379},
          doi = {10.1051/0004-6361:20041729},
archivePrefix = {arXiv},
       eprint = {astro-ph/0411110},
 primaryClass = {astro-ph},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2005A&A...432..369S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```


## License and Acknowledgements
Copyright Brian Svoboda (2021) and distributed under the terms of the GPL v3 software license. RADEX is authored by [Floris van der Tak](https://personal.sron.nl/~vdtak/) and contributors.
