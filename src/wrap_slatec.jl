@doc raw"""
# Wrap SLATEC
Wrap the `slatec.f` Fortran file distributed as part of RADEX. The
`slatec.f` library is used to solve the system of equations for statistical
equillibrium and return the populations for each level. The rate matrix is
factored into lower and upper diagonal components (LU decomposition) and solved
using partial-pivoting. The subroutines `ludcmp` (LU decompose) and `lubksb`
(LU back-substition) format the inputs for the main `SGEIR` subroutine, which
itself calls the LINPACK subroutines `SGEFA` and `SGESL`.

The Julia standard library includes LINPACK and uses LU-decomposition with
partial-pivoting as the default algorithm for matrix inversion when calling
`A \ b`, thus this module is not necessary but included for the purposes of
verifying the output against the exact output of RADEX.

# Usage
First compile the `slatec.f` file distributed in the `src` directory of RADEX
into a shared library `libslatec.so`:

```bash
    gfortran -shared -O2 slatec.f -o libslatec.so -fPIC
```

The Fortran standard does not specify how symbol names are expressed in shared
libraries and are mangled in different ways by different compilers and even
versions of the same compiler. To check, use the `nm` program to list symbol
names from the object file:

```bash
    nm libslatec.so
```

Look for entries in the output of `nm` that are similar to "ludcmp" and
"lubksb". The gfortran compiler typically appends an underscore to the symbol
name, thus `ludcmp_` and `lubksb_`. If these are not the names, then this
file will need to edited to reflect those names, as Julia's `ccall` must be
given literals for the symbol, library name, output type, and input types
(although this could be avoided with a macro that replaces these values at
compile time).

Finally, ensure that the directory containing the library file is included
in the shell environment variable `LD_LIBRARY_PATH`:

```bash
    # Run before starting Julia or include in shell configuration.
    export LD_LIBRARY_PATH=<path-to-lib>:$LD_LIBRARY_PATH
    # Or include when starting Julia.
    LD_LIBRARY_PATH=<path-to-lib> julia
```

To use, call the `call_slatec!` function to modify the `rhs` level population
work array in-place.
"""
module WrapSlateC

export call_slatec!, solve_rate_eq_slatec!


@doc raw"""
    call_slatec!(yrate, rhs)

Solve the system of equations in `yrate` using LU-decomposition and
partial-pivoting to perform the matrix inversion operation `yrate \ rhs`. The
results are stored by modifying the vector `rhs` in-place. The Fortran routines
in RADEX are called natively and included for verification purposes.
"""
function call_slatec!(yrate::Matrix{Float64}, rhs::Vector{Float64})
    @assert size(yrate, 1) == size(yrate, 2) == length(rhs)
    maxlev = length(rhs)
    nplus = maxlev + 1
    indx = 0
    dsign = 0.0
    # matrix.f: call ludcmp(yrate,nplus,maxlev,indx,dsign)
    ccall((:ludcmp_, "libslatec.so"),
          Cvoid,
          (Ptr{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}),
          yrate, nplus, maxlev, indx, dsign)
    # matrix.f: call lubksb(yrate,nplus,maxlev,indx,rhs)
    ccall((:lubksb_, "libslatec.so"),
          Cvoid,
          (Ptr{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{Float64}),
          yrate, nplus, maxlev, indx, rhs)
    return nothing
end


"For testing results with LU factorization as used by RADEX."
function solve_rate_eq_slatec!(xpop, yrate, rhs)
    call_slatec!(yrate, rhs)
    xpop .= rhs
    F = eltype(rhs)
    rhs .= zero(F)
    rhs[end] .= one(F)
    xpop
end
function solve_rate_eq_slatec!(sol)
    sole_rate_eq_slatec!(sol.xpop, sol.yrate, sol.rhs)
    sol
end


end  # module
