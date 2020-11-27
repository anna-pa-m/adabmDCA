# Boltzmann machine learning for Potts models of biological data

This algorithm infers several maximum-entropy statistical models of Potts variables given a set of observables (usually a set of biological sequences). More precisly, it is able to learn the couplings and the fields parameters of a set of generalized Direct Coupling Analysis (DCA) models given a Multiple Sequence Alignment (MSA) of protein or RNA domains, binary sequences (for Ising variables) in FASTA format. The learning is performed via a gradient ascent of the likelihood of the data; the moel observables are computed via a Markov Chain Monte Carlo (MCMC) sampling.

## Installation

This code in written in C++ language. To properly install ```adabmDCA```, run
```
make
```
on adabmDCA/src folder. It suffices a ```g++``` compiler.

## Usage

All the possible features (model choices, input and output files) implemented in ```adabmDCA``` can be shown typing
```
./adabmDCA -h
```
We show here the most important features of the algorithm.

## Basic run
Suppose that we would like to learn a DCA model associated with a given MSA, i.e. that in `test/PF00018.fasta`. The command
```
./adabmDCA -f ../test/PF00018.fasta -z 1 -m 500 -c 1e-2 
```
will run a Boltzmann learning, using the default MCMC, printing to files every 500 iterations the parameters and the statistics of MSA against those of the model. The run will stop when either we reach the maximum number of iterations (default: 2000) or we reach convergence, i.e. the maximum fitting error among all the connected covariances is 1e-2 (as specified in the command line using `-c`).
The output will produce a list of information concerning the input data and the adopted learning strategy
```
***** Boltzmann machine for DCA model ******
****** Initializing data structures ******
Using alphabet: -ACDEFGHIKLMNPQRSTVWY
Reading MSA from../test/PF00018.fasta
Reading alignment completed.
M = 55784 L = 48 q = 21
Computing weights...done
Computing empirical statistics...Meff: 7e+03
No three-points correlations indices specified
****** Initializing model ******
Performing Metropolis-Hastings MC
Adaptive sampling time. Initial sampling time: 10 initial equilibration time: 20
Using 1000 seeds and tot. number of points 50000
MC chains are randomly initialized at the beginning of each iteration
Learning strategy: 
Using standard gradient descent with constant learning rate (for J 0.05 for h 0.05 )
Using pseudo-count: 0.0001
Zero-parameters initialization...done
Using full Potts model
Sparsity after initialization: 0
****** Starting learning loop ******
Printing output every 1 iterations - parameters every 500

```
At every step (the frequency can be specified by the `-z` flag)  it writes to the `stdout` an update of the learning, for instance
```
it: 4 el_time: 48 N: 50000 Teq: 12 Twait: 6 merr_fm: 8.8e-01 merr_sm: 8.7e-01 averr_fm: 5.4e-02 averr_sm: 3.4e-03 cov_err: 2.0e-01 corr: 0.0057 sp: 0 lrav: 0.05
```
where:

  - `it:` gives the number of the current iteration
  - `el_time:` is the running time (in seconds)
  - `N:` is the number of samples used for computing the model observables
  - `Teq:` is the equilibration time (in MC sweep) of each MC chain
  - `Twait:` is the de-correlation time (in MC sweep) between two sampled configurations, and within each chain
  - `merr_fm:` maximum value of the fitting error for the first moments
  - `merr_sm` maximum value of the fitting error for the second moments
  - `averr_fm:` mean value of the fitting error for the first moments
  - `averr_sm:` mean value of the fitting error for the second moments
  - `cov_err:` maximum value of the fitting error for the connected covariances (used for checking the convergence)
  - `corr:` Pearson correlation coefficients between the model and MSA connected covariances
  - `sp:` sparsity of the couplings parameters, i.e. the number of couplings fixed to zero divided by q<sup>2</sup>L(L-1)/2
  - `lrav:` average learning rate
  
## Input/Output files

## License
Apache?

## Citations
The performances of the adaptive Boltzmann Machine are shown in MISSING LINK. Please cite this article if you use this software or any part of it.

```adabmDCA``` has been used in:
 - https://arxiv.org/abs/2005.08500 for learning the Potts model and pseudo Hidden Markov model (see option ```-H``` in Advanded options/Available maximum entropy models) of the studied proteins and RNA sequences;
 - MISSING LINK for the learning of the dense model and the decimation procedure (see Advanced options/Pruning the coupling)

## Advanced options

### Tuning the Monte Carlo Markov Chain

### Compute non-fitted third order statistics

### Available maximum entropy models

### Pruning the couplings

### Learning strategies

## To do list

Several things:

  - implement ```adam``` learning strategy
  - parallelization of the MCMC




