# Boltzmann machine learning for Potts models of biological data

## Authors

Anna Paola Muntoni and Francesco Zamponi

## Description

This is an implementation of the Boltzmann machine learning to infer several maximum-entropy statistical models of Potts or Ising variables given a set of observables. More precisly, it infers the couplings and the fields of a set of generalized Direct Coupling Analysis (DCA) models given a Multiple Sequence Alignment (MSA) of protein or RNA sequences. It is also possible to infer an Ising model from a set of binary configurations. The learning is performed via a gradient ascent of the likelihood of the data in which the model observables are computed via a Markov Chain Monte Carlo (MCMC) sampling.

`adabmDCA` has been used in:
 - *Aligning biological sequences by exploiting residue conservation and coevolution* - A. Muntoni, A. Pagnani, M. Weigt, F. Zamponi (on [Phys. Rev. E](https://link.aps.org/doi/10.1103/PhysRevE.102.062409)) for learning the Potts model and pseudo Hidden Markov model (see option `-H` in `Advanded options/Available maximum entropy models`) of the studied proteins and RNA sequences;
 - *Sparse generative modeling of protein-sequence families* - P. Barrat-Charlaix, A. Muntoni, K. Shimagaki, M. Weigt, F. Zamponi (on [arXiv](https://arxiv.org/abs/2011.11259)) for the learning of the dense model and to perform the decimation procedure (see `Advanced options/Pruning the coupling`)


## Installation

This code in written in C/C++ language. To properly install `adabmDCA`, run
```
make
```
on `adabmDCA/src` folder. It suffices a `g++` compiler.

## Usage

All the possible routines (model choices, decimation procedure, input and output files) implemented in `adabmDCA` can be shown typing
```
./adabmDCA -h
```

For a detailed description see the [documentation](https://github.com/anna-pa-m/adabmDCA/blob/master/docs/documentation.md) file.

