```@meta
Author = "Brian Svoboda"
```

# Installation
To install Jadex, open a Julia interactive session (see the Julia documentation's [Getting Started](https://docs.julialang.org/en/v1/manual/getting-started/) page for further instructions) and then enter the package management mode by entering `]` and then `add Jadex`:

```
julia> ]

pkg> add Jadex
```

and, optionally, to run the test suite call:

```
pkg> test Jadex
```

Alternatively, to install from the main branch hosted on GitHub, call:

```
pkg> add https://github.com/autocorr/autocorr/Jadex.jl
```


## Installing in Python
To call Jadex from Python we can use the [PyJulia](https://pyjulia.readthedocs.io/en/stable/) Python package to call-out to our Julia environment. Following the PyJulia [installation instructions](https://pyjulia.readthedocs.io/en/stable/installation.html):

__(1)__ Ensure that Julia is installed and that the Julia executable is in one's `PATH` (a symbolic link is okay but not an alias).

__(2)__ Install Jadex within Julia (see the previous section above).

### Python installation without Anaconda
Now we need to install PyJulia, which can be imported as the Python `julia` module. If your system uses Anaconda, then please see the alternate instructions in the following section.

__(3a)__ Install PyJulia in Python by calling `python3 -m pip install --user julia` from the command line.

__(4a)__ Install the Julia dependencies of PyJulia (i.e., PyCall) from within Python by calling:

```python
import julia
julia.install()
```

To import Jadex and call it, run:

```python
from julia import Jadex
mol = Jadex.Specie("hco+@xpol", datadir="path/to/data")
```


### Python installation with Anaconda
If using Anaconda, we need to mitigate a known limitation of PyJulia/PyCall for distributions of Python that use static linking, which is the case for Anaconda (see PyJulia's [Troubleshooting](https://pyjulia.readthedocs.io/en/stable/troubleshooting.html) page for more details).

__(3b)__ Install PyJulia in Python by calling `pip install julia` from the conda environment you wish to use (i.e., after calling `conda activate myenv` or similar).

__(4b)__ Before using PyJulia (i.e., `import julia`), turn off the Julia module compilation cache and then install the Julia dependencies of PyJulia:

```python
from julia.api import Julia
jl = Julia(compiled_modules=False)
import julia
julia.install()
```

!!! note
    Unfortunately turning off module compilation and cacheing adds some initial
    startup overhead because pre-compiled modules can't be re-used. This is not
    so bad for Jadex, however, but would be more problematic if one was calling
    Julia modules for plotting.

