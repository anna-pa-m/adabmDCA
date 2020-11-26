# Boltzmann machine learning for Potts models of biological data

## Installation

This code in written in C++ language. To properly install ```adabmDCA```, run
```
make
```
on adabmDCA/src folder. It suffices a ```g++``` compiler.

## Usage
This algorithm is able to infer several Max-Ent statistical models of Potts variables given a set of sequences. More precisly, it is able to learn the couplings and the fields parameters of a set of Direct Coupling Analysis (DCA) models given a Multiple Sequence Alignment (MSA) of protein, RNA, binary sequences (for Ising variables) in FASTA format. The learning is performed via a gradient ascent of the likelihood of the data, computed via a Monte Carlo sampling, either using Metropolis or Gibbs steps.
All the possible features (model choices, input and output files) implemented in ```adabmDCA``` can be shown typing
```
./adabmDCA -h
```
## Basic run
Suppose we would like to learn a DCA model associated with a given MSA, i.e. in ```test/PF00018.fasta```. The command
```
./adabmDCA -f ../test/PF00018.fasta -z 1 -m 500 -c 1e-2 
```
will run a Boltzmann learning, using the default MCMC, printing to files every 500 iterations the parameters and the statistics of MSA against thoses of the model. The run will stop when either we reach the maximum number of iterations (default: 2000) or we reach convergence, i.e. the maximum fitting error among all the connected covariances is 1e-2 (tunable using the flag ```-c```).
The output will produce a list of information concerning the adopted learning strategy
```
****** Boltzmann machine for DCA model ******
****** Initializing data structures ******
Using alphabet: -ACDEFGHIKLMNPQRSTVWY
Reading MSA from ../test/PF00018.fasta
Reading alignment completed.
M = 17370 L = 48 q = 21 
Computing weights...done
Computing empirical statistics...Meff: 3742.451636
No three-points correlations indices specified
****** Initializing model ******
Performing Metropolis-Hastings MC
Adaptive sampling time. Initial sampling time: 10, initial equilibration time: 20
Using 1000 seeds and tot. number of points 50000
MC chains are randomly initialized at the beginning of each iteration
Learning strategy: Using standard gradient descent with constant learning rate (for J 5.000e-02, for h 5.000e-02)
Using pseudo-count: 3e-04
Zero-parameters initialization...done
Using full Potts model
Sparsity after initialization: 0.000000
```
At every step (the frequency can be specified by the ```-z``` flag)  it writes to ```stdout``` an update of the learning
```
it: 0 el_time: 7 N: 50000 Teq: 18 Twait: 9 merr_fm: 9.2e-01 merr_sm: 8.6e-01 averr_fm: 5.7e-02 averr_sm: 3.5e-03 cov_err: 2.0e-01 corr: -0.00 sp: 0.0e+00 lrav: 5.0e-02
```
where:

  - ```it:``` gives the current iteration
  - ```el_time:``` is the running time (in seconds)
  - ```N:``` is the number of samples used for computing the model averages
  - ```Teq:``` is the equilibration time (in sweep) of each MC chain
  - ```Twait:``` is the de-correlation time (in sweep) between two sampled configurations within each chain
  - ```merr_fm:``` maximum value for the fitting error for the first moments
  - ```merr_sm:``` maximum value for the fitting error for the second moments
  - ```averr_fm:``` mean value for the fitting error for the first moments
  - ```averr_sm:``` mean value for the fitting error for the second moments
  - ```cov_err:``` maximum value for the fitting error for the connected covariances
  - ```corr:``` Pearson correlation coefficients between the model and MSA connected covariances
  - ```sp:``` sparsity of the couplings parameters, i.e. the number of couplings fixed to zero divided by q<sup>2</sup>L(L-1)/2
  - ```lrav:``` average learning rate
  
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




