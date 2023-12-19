# ace-pac-utils (APU)
Tool to process and analyze the output of the chemistry codes ACE and PAC

## Introduction

Welcome the acepacutils (APU) repository. APU is software package that can be used to analyse and post-process the output of the chemical kinetics code ([Agúndez et al., 2014](https://www.aanda.org/articles/aa/full_html/2014/04/aa22895-13/aa22895-13.html)). The code consists of a 1D chemical equilibrium solver (ACE) and a 1D and pseudo-2D chemical disequilibrium solver (PAC). With APU you can easily
- plot the convergence evolution, steady-state abundances, pressure-temperature structure, input stellar spectrum, and eddy diffusion (Kzz) profile of any ACE/PAC simulation,  
- read-in the simulation parameters that were used (such as R_pl, M_pl, ...) in the directory of the ACE/PAC simulation
- post-process the output of a chemistry simulation into the radiative transfer code [petitRADTRANS](https://gitlab.com/mauricemolli/petitRADTRANS) to generate a transmission spectrum.


## Installation
APU is currently only installable by cloning this repository.

## Requirements 
- petitRADTRANS
- Xarray

## Documentation
t.b.d.

## Contact

If you have questions regarding the code, or you need further help, feel free to contact me ([Thomas Konings](thomaskonings.github.io)).
