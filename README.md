@note this repository is a copy of [Camelio (2022), "hydro-bulk-1D", zenodo, https://doi.org/10.5281/zenodo.6478023](https://zenodo.org/record/6478023)

# hydro-bulk-1D

1-dimensional general relativistic hydrodynamic code with bulk viscosity.

`hydro-bulk-1D` tests the correspondence between the bulk viscous formulation
of the multi-component fluid and that of the Israel-Stewart theories
(see Gavassino, Antonelli, and Haskell, 2021, CQG 38:075001)
in the context of neutron stars.

The theory and the code are described in the following papers:

* [1] Camelio, Gavassino, Antonelli, Bernuzzi, and Haskell, "Simulating bulk viscosity in neutron stars I: formalism", submitted to PRD (2022).
* [2] Camelio, Gavassino, Antonelli, Bernuzzi, and Haskell, "Simulating bulk viscosity in neutron stars II: evolution in spherical symmetry", submitted to PRD (2022).

`hydro-bulk-1D` is released under the MIT License.

If you use `hydro-bulk-1D` in your research, please cite Refs. [1] and [2].

## Content

* `README.md` -- this file
* `hydro-bulk-1D.c` -- the hydrodynamic code
* `parameters.h` -- parameter file included in `hydro-bulk-1D.c`
* `par/par-*.h` -- parameter files used in [2]
* `Doxyfile` -- Doxygen configuration file
* `LICENSE.txt` -- the MIT License

## Usage

The program is the file `hydro-bulk-1D.c`, which directly include the parameter file `parameters.h`.

The parameter files for the simulations and tests presented in Ref. [2] are `par/par-*.h`.

For example, to simulate the shocktube test using the `gcc` compiler, execute on a Linux shell:

```
$ cp par/par-shocktube.h parameters.h
$ gcc hydro-bulk-1D.c -O2 -lm -o shocktube.x
$ ./shocktube.x
```

In order to generate the automatic documentation using Doxygen, run:

```
$ doxygen
```

You can read the automatic documentation openening `html/index.html` in a browser.
Note that I have set `HAVE_DOT = YES` in the Doxygen configuration file;
this option requires the `dot` program from the `graphviz` package.

## Output

`hydro-bulk-1D` generates in output a profile file and a log file.

The profile file name can be set with the `OUTPUT_FILE` macro.
It is a `gnuplot`-friendly text file, which means that each time snapshot is separated by two blank lines.
The format of the output is the following:

* `r|x` -- radial or Cartesian coordinate,
* `ρ` -- rest mass density,
* `Wv` --  Lorentz factor times velocity,
* `u` -- internal specific (per unit mass) energy density,
* `[Ye|Π]` -- electron fraction or bulk stress (optional),
* `[Yμ]` -- muon fraction (optional),
* `[mg]` -- enclosed gravitational mass (optional),
* `cs` -- speed of sound.

The log file name can be set with the `LOG_FILE` macro.
It is a text file, and its format is the following:

* `t` -- time,
* `ρ` -- rest mass density in the first non-ghost cell,
* `Wv` -- Lorentz factor times velocity in the first non-ghost cell,
* `u` -- internal specific (per unit mass) energy density in the first non-ghost cell,
* `[Ye|Π]` -- electron fraction or bulk stress in the first non-ghost cell (optional),
* `[Yμ]` -- muon fraction in the first non-ghost cell (optional),
* `cs` -- speed of sound in the first non-ghost cell,
* `[Mg]` -- total gravitational mass (optional),
* `[Mb]` -- total baryon mass (optional),
* `ρend` -- rest mass density in the last non-ghost cell.

Unless otherwise specified, all quantities in the code, in the parameter files, and in the output files, are in code units: c = G = M☉ = kB = 1.

## Acknowledgments

This work was supported by the Polish National Science Centre (NCN) grant number OPUS 2019/33/B/ST9/00942.

I am grateful to Sebastiano Bernuzzi, Lorenzo Gavassino, Marco Antonelli, and Brynmor Haskell for useful discussions.
