# Changelog

## v0.1.11
- Add SizeMax to meshing resolution settings
- Add examples for coplanar waveguides and crosstalk
- Improve plotting function

## v0.1.10
- Add backward compatability for "Fix Typo: Culomb -> Coulomb" until the end of the year

## v0.1.9
- Improvements on the Docs
- Fix Typo: Culomb -> Coulomb

## v0.1.8

- Fix build-backend in pyproject.toml

## v0.1.7

- Add __init__.py to femwell/maxwell

## v0.1.6

- Add build-backend to pyproject.toml

## v0.1.5

- Add function to calculate intensities
- Add function to plot intensities

## v0.1.4

- Fix meshing if a LineString is cut by the geometry

## v0.1.3

- Update minimum version of scikit-fem

## v0.1.2

- Fix n_eff guess to make scipy solver work

## v0.1.1

- Add mesh refinement
- Add scipy benchmarks
- Fix n_eff-guess for compute_modes (->seems to fix mode calculations with scipy)

## v0.1.0

- Rewrite maxwell waveguide mode solver - now completely object oriented!
- Start implementing schrödinger's equation
- Start docs for schrödinger's equation
- Add potential well

## v0.0.18

- Add function to sort the modes by power within a given area

## v0.0.17

- Allow to disable normalizing the modes in compute_modes
- Add cache to compute_modes
- Add automatic changelog generation
- Add changelog to docs
- Fixed meshing
