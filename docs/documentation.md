

## Documentation 

### Basic run
Let us infer a DCA model, i.e. a Potts model of 21 colors on a fully connected topology, associated with the MSA in `test/PF00018.fasta`. The command
```
./adabmDCA -f ../test/PF00018.fasta -z 1 -m 500 -c 1e-2 
```
will run a Boltzmann learning, using the default MCMC, printing to files every 500 iterations (as specified in the option `-m`) the temporarily parameters of the model, and both data and model statistics. The run will stop when either we reach the maximum number of iterations (default: 2000) or we reach convergence, i.e. the maximum fitting error among all the connected covariances is 1e-2 (as specified in the command line using `-c`).
The program will output a list of information concerning the input data and the adopted learning strategy:
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
In the early stage of the running, the algorithm performs a reweighting of the original sequences and modifies their statistical significance as explained [here](https://www.pnas.org/content/108/49/E1293); the sequence similarity threshold can be set using `-l` input flag. Besides, it is also possible to give the set of weigths in a file using the flag `-w`. Then, the algorithm computes the data statistics (one and two-site frequencies) and eventually corrects the empirical moments introducing a pseudo-count (read [here](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.012132) for details) whose value is set to `1/Meff` and can be modified using `-d`.
Finally, the algorithm itervatively updates the parameters `(J,h)` of the Potts model by standard gradient ascent (some adaptive learning strategies are described in `Advanced options/Learning strategies`). The model observables are computed from 1000 independent MC chains which, at equilibrium, sample 50 configurations each. The number of MC chains and the number of sampled configurations can be modified by `-s` and `-n` respectively. `adabmDCA` allows for an advanced tuning of the sampling, see `Advanced options/Tuning the Markov Chain Monte Carlo`.
At every step (the frequency can be specified by the `-z` flag)  it writes to the `stdout` an update of the learning, like the following:
```
it: 4 el_time: 48 N: 50000 Teq: 12 Twait: 6 merr_fm: 8.8e-01 merr_sm: 8.7e-01 averr_fm: 5.4e-02 averr_sm: 3.4e-03 cov_err: 2.0e-01 corr: 0.0057 sp: 0 lrav: 0.05
```
where:

  - `it:` gives the number of the current iteration
  - `el_time:` is the running time (in seconds)
  - `N:` is the number of samples used for computing the model observables
  - `Teq:` is the equilibration time (in MC sweep) of each MC chain
  - `Twait:` is the de-correlation time (in MC sweep) between two sampled configurations, for all the chains
  - `merr_fm:` maximum value of the fitting error for the first moments
  - `merr_sm` maximum value of the fitting error for the second moments
  - `averr_fm:` mean value of the fitting error for the first moments
  - `averr_sm:` mean value of the fitting error for the second moments
  - `cov_err:` maximum value of the fitting error for the connected covariances (used for checking the convergence)
  - `corr:` Pearson correlation coefficient between the model and the data connected covariances
  - `sp:` sparsity of the couplings parameters, i.e. the number of couplings fixed to zero divided by q<sup>2</sup>L(L-1)/2
  - `lrav:` average learning rate
  
### Input/Output files

#### Raw data

`adabmDCA` takes as input a FASTA file whose name mus follow the flag `-f` as well as a set of given empirical statistics collected in a file (use `-q` here). In the latter case, the input file must be formatted as
```
s i j symbol_i symbol_j value
...
m i symbol_i value
```
where the `s` rows contain the two-site frequencies of the `i j` sites for colors `symbol_i symbol_j` and the `m` lines the empirical one-site frequencies of site `i` for color `symbol_i`. 

### Advanced options

#### Tuning the Markov Chain Monte Carlo

##### Equilibration test

The standard run of this implementation of the Boltzmann learning ensures that the model statistics is estimated using a MCMC sampling performed at equilibrium. To do this the number of MC sweeps between each pair of sampled configurations, `Twait`, is tuned at each iteration, and the number of MC sweeps before the first collected configuration, `Teq`, is set equal to `2Twait`. This last choice is likely to provide a first equilibrium configuration, considering how we fix `Twait`. Let us call the configuration of chain `i` sampled after `n` steps of the MCMC as s<sup>i</sup><sub>n</sub>(Teq + n Twait). For each iteration we compute the average value (among both chains and n) and deviations from them,  of the following quantities:

  - Q<sup>exp</sup> = δ<sub> s<sup>i</sup><sub>n</sub> s<sup>k</sup><sub>n</sub> </sub> 
  - Q<sup>int,1</sup> = δ<sub> s<sup>i</sup><sub>n</sub> s<sup>i</sup><sub>n+1</sub> </sub>
  - Q<sup>int,2</sup> = δ<sub> s<sup>i</sup><sub>n</sub> s<sup>i</sup><sub>n+2</sub> </sub>
  
  If the overlap between independent chains is not similar to the overlap of two samples in the same chain, distant `2Twait`, i.e. Q<sup>exp</sup> < Q<sup>int,2</sup>, then increase Twait: Twait <- Twait + L
  If the overlap between independent chains is similar to the intra-chain overlap between two samples at distance `Twait`, i.e. Q<sup>exp</sup> ~ Q<sup>int,1</sup>, then decrease Twait: Twait <- Twait -L
  
  If this way, two consecutive samples of the same chain can be slightly correlated but for sure two configurations at distance `2Twait` are reasonably de-correlated. One may assume that even starting from an arbitrary sample and waiting `2Twait`, the final sample would be at equilibrium: this is not proven but very likely the case. Besides, if one uses persistent chains (see below), it is fair to consider the time to equilibrate equals to the de-correlation time.
  
  ##### Advanced settings
  
  It is possible to avoid the equilibration test and to sample at fixed Teq and Twait using the flags
  ```
  -L -e Teq -t Twait
  ```
  It is also possible to perform a persistent sampling (flag `-P`), meaning that the chains are initialized as uniformly random configurations only at the first iteration, and then kept persistent. The initial configurations can be extracted from data points using the flag `-Q`.
  
  The standard implementation of `adabmDCA` uses as default the Metropolis update rule, but a Gibbs sampling can be used by adding the `-G` flag.
  

#### Compute non-fitted third order statistics

#### Available maximum entropy models

#### Pruning the couplings

#### Learning strategies

### To do list

Several things:

  - implement ```adam``` learning strategy
  - parallelization of the MCMC




