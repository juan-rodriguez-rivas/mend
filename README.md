
# MEND
**MEND** is a correction for DCA coevolutionary scores meant to improve the prediction of inter-protein physical contacts. It takes as input a paired alignment codified in two FASTA files and returns the MEND (corrected) coevolutionary scores for each pair inter-protein positions. It is based on the observation that the maximum APC scores obtained by randomized the sequences pairings provide an effective case-specific background estimation. Given its robustness regarding the number of sequences, it is particularly useful to detect inter-protein contact predictions when the paired alignment contain few sequences.

MEND uses the following software:
- [fmpl](https://github.com/simomarsili/fmpl) implementation of plmDCA developed by Simone Marsili (@github/simomarsili) to compute the DCA coevolutionary model.

## Installation

The following instructions are given for its installation in Ubuntu. The installation in any linux distribution should be very similar, with the use of a different package manager.

MEND requires python3 and a fortran compiler:
> sudo apt-get install python3 gfortran

The following python packages are needed
> pip install bio numpy

Compile fmpl program
> cd src/mpl/
> make

To recompile, use ``make clean; make``


## Test
The included test should take around 5 minutes
> ./run_test.py

The test will quantify if there are relevant differences between the computed contact prediction against a stored one.

