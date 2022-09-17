```@meta
Author = "Brian Svoboda"
```

# Jadex

High-performance and extensible re-implementation of the Fortran 77 code [RADEX](https://personal.sron.nl/~vdtak/radex/index.shtml) (van der Tak et al. [2007, A&A 468, 627](https://www.aanda.org/articles/aa/abs/2007/23/aa6820-06/aa6820-06.html)) in the [Julia](https://julialang.org/) programming language. A Python wrapper is provided using `PyJulia`.

Jadex & RADEX are non-LTE radiative transfer codes for calculating atomic and molecular spectral line intensities. They assume a uniform medium (i.e., single zone) with energy level populations in statistical equillibrium. Optical depth is calculated using an escape probability defined for different cloud geometries (e.g., LVG, uniform sphere, slab).

For cases where the same input parameters and constants are used, results from Jadex should match RADEX to within a tolerance of five significant figures. When using an un-patched copy of RADEX values are expected to match within four significant figures due to higher precision mathematical constants in Jadex. Jadex has been validated against the RADEX wrapper [SpectralRadex](https://github.com/uclchem/SpectralRadex) for a suite of species and physical conditions (see `test/validation.jl`).


## Package features

* Improved performance (~110x).
* Improved convergence rate through Ng-acceleration
* User definable escape probability and background radiation field.
* Multi-threaded parameter grid calculations with interpolation.
* Cross-platform: tested on Linux, MacOS, and Windows.
* Python interface provided using PyJulia.
* Test and validation suite.
* Extensible design.


## Installing

To install Jadex, open an [interactive Julia session](https://docs.julialang.org/en/v1/manual/getting-started/), press the `]` key to enter the package management mode, and execute the command `add Jadex`. For further instructions, see the [Installation](@ref) page.


## Quickstart guide

First download the spectroscopic and collision rate data `.dat` files from [LAMDA](https://home.strw.leidenuniv.nl/~moldata/) for the specie(s) of interest. Import the functions and types from Jadex into the namespace by executing:

```julia
using Jadex
```

Now read and parse the LAMDA file into a `Specie` object. The identifier passed to `Specie` should be the same as the file name except with the `.dat` extension removed. Thus for the `hco+@xpol.dat` file, use the identifier `hco+@xpol`. After downloading the file, it is worth occasionally checking LAMDA for updates.

```julia
mol = Specie("hco+@xpol", datadir="path/to/LAMDA/files")
```

Note that Julia is a JIT-compiled language, and so the first execution will have significant run-time overhead due to code compilation. In Jadex, the majority of the overhead occurs when calling `Specie` because the CSV.jl and DataFrames.jl libraries used to parse the LAMDA rate files are large. Subsequent executions will be much faster.

Next define the physical properties of the run. The default escape probability geometry is the expanding (or collapsing) sphere using the "large velocity gradient" (LVG) approximation. The default background radiation field is a blackbody with the CMB temperature.

```julia
rdf = RunDef(mol, density=Dict("h2" => 1e4), tkin=20.0, cdmol=1e13, deltav=1.0)
```

Now solve for the level populations, optical depths, and excitation temperatures for all transitions. Emergent radiation temperatures and integrated intensities are then calculated for all transitions within the selected frequency range.

```julia
df = get_results(rdf, freq_max=500)  # GHz
```

Results are returned in a `DataFrame` object. Please refer to the DataFrames.jl [documentation](https://dataframes.juliadata.org/stable/) for further information on working with `DataFrame`s.

Please refer to the [User Guide](@ref) for additional information regarding running and using Jadex.


## Citing this work
If you use Jadex in an academic work, we request that you cite the following references, including the original publication for RADEX (van der Tak et al. 2007):

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


## Similar tools

Other similar tools for calculating spectral line intensities include:

* [RADEX](https://personal.sron.nl/~vdtak/radex/index.shtml): Single-zone, non-LTE radiative transfer code written in Fortran 77.
* [pyradex](https://github.com/keflavich/pyradex): Python interface to RADEX with Fortran bindings auto-generated by `f2py`.
* [SpectralRadex](https://github.com/uclchem/SpectralRadex): Python interface to RADEX with Fortran bindings auto-generated by `f2py`, grid evaluation, and creation of model spectra.
* [myradex](https://github.com/fjdu/myRadex): Single-zone, non-LTE radiative transfer code written in Fortran 90. Level populations are computed using an ODE solver.
* [ndradex](https://github.com/astropenguin/ndradex): Command-line wrapper for the RADEX executable designed for evaluating parameter grids.
* [DESPOTIC](https://bitbucket.org/krumholz/despotic/): Python package to perform calculations related to line emission and thermal behavior in cold interstellar clouds.
* [MOLPOP-CEP](https://github.com/aasensio/molpop-cep): Multi-zone model using the coupled escape probability (CEP) method.


## License and Acknowledgements
Copyright Brian Svoboda (2021) and distributed under the terms of the GPL v3 software license. RADEX is authored by [Floris van der Tak](https://personal.sron.nl/~vdtak/) and contributors.

