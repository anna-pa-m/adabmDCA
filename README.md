# Adaptive Boltzmann machine learning for Potts models of biological data

## Description

This is an implementation of the Boltzmann machine learning to infer several maximum-entropy statistical models of Potts or Ising variables given a set of observables. More precisely, it infers the couplings and the fields of a set of generalized Direct Coupling Analysis (DCA) models given a Multiple Sequence Alignment (MSA) of protein or RNA sequences. It is also possible to infer an Ising model from a set of spin configurations. The learning is performed via a gradient ascent of the likelihood of the data in which the model observables are computed via a Markov Chain Monte Carlo (MCMC) sampling.

More details are available in [adabmDCA: Adaptive Boltzmann machine learning for biological sequences](https://arxiv.org/abs/2109.04105). Please cite this paper if you use (even partially) this code.

`adabmDCA` has been used in:
 - *Aligning biological sequences by exploiting residue conservation and coevolution* - A. Muntoni, A. Pagnani, M. Weigt, F. Zamponi (on [Phys. Rev. E](https://link.aps.org/doi/10.1103/PhysRevE.102.062409)) to learn the Potts model and pseudo Hidden Markov model (see the option `-H` in [Advanced options/Available maximum entropy models](https://github.com/anna-pa-m/adabmDCA/blob/master/docs/documentation.md)) of the studied protein and RNA seed alignment;
 - *Sparse generative modeling of protein-sequence families* - P. Barrat-Charlaix, A. Muntoni, K. Shimagaki, M. Weigt, F. Zamponi (on [Phys. Rev. E](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.104.024407)) to learn the dense model and to perform the information-based pruning of the coupling parameters (see [Advanced options/Pruning the coupling](https://github.com/anna-pa-m/adabmDCA/blob/master/docs/documentation.md)).


## Installation

This code is written in C/C++ language. To properly install `adabmDCA`, run
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

