# Changelog

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
