# Boltzmann machine learning for Potts models of biological data

## Authors

Anna Paola Muntoni and Francesco Zamponi

## Description

This is an implementation of the Boltzmann machine learning to infer several maximum-entropy statistical models of Potts variables given a set of observables (usually a set of biological sequences). More precisly, it is able to learn the couplings and the fields of a set of generalized Direct Coupling Analysis (DCA) models given a Multiple Sequence Alignment (MSA) of protein or RNA domains, binary sequences (for Ising variables) in FASTA format. The learning is performed via a gradient ascent of the likelihood of the data in which the model observables are computed via a Markov Chain Monte Carlo (MCMC) sampling.

`adabmDCA` has been used in:
 - *Aligning biological sequences by exploiting residue conservation and coevolution* - A. Muntoni, A. Pagnani, M. Weigt, F. Zamponi (on [arXiv](https://arxiv.org/abs/2005.08500)) for learning the Potts model and pseudo Hidden Markov model (see option `-H` in Advanded options/Available maximum entropy models) of the studied proteins and RNA sequences;
 - *Sparse generative modeling of protein-sequence families* - P. Barrat-Charlaix, A. Muntoni, K. Shimagaki, M. Weigt, F. Zamponi (on [arXiv](https://arxiv.org/abs/2011.11259)) for the learning of the dense model and to perform the decimation procedure (see Advanced options/Pruning the coupling)


## Installation

This code in written in C/C++ language. To properly install `adabmDCA`, run
```
make
```
on `adabmDCA/src` folder. It suffices a `g++` compiler.

## Usage

All the possible features (model choices, decimation procedure, input and output files) implemented in `adabmDCA` can be shown typing
```
./adabmDCA -h
```
In the following, we present in details the main features of the algorithm.

## Basic run
Suppose that we would like to learn a DCA model associated with a given MSA, i.e. that in `test/PF00018.fasta`. The command
```
./adabmDCA -f ../test/PF00018.fasta -z 1 -m 500 -c 1e-2 
```
will run a Boltzmann learning, using the default MCMC, printing to files every 500 iterations (as specified in the option `-m`) the parameters and the statistics of MSA against those of the learned model. The run will stop when either we reach the maximum number of iterations (default: 2000) or we reach convergence, i.e. the maximum fitting error among all the connected covariances is 1e-2 (as specified in the command line using `-c`).
The program will output a list of information concerning the input data and the adopted learning strategy
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
As mentioned in the output, the algorithm performs a reweighting of the statistical significance of every sequence as explained [here](https://www.pnas.org/content/108/49/E1293), the sequence similarity threshold can be set using `-l` input flag. Then, the algorithm computes the data statistics (first and second moments) and eventually corrects the empirical moments introducing a pseudo-count (read [here](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.012132) for details) whose value is set to `Meff` and can be modified using `-d`. 
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

`
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




