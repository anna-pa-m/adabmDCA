

# Documentation 

## Basic run
Let us infer a [DCA model](https://en.wikipedia.org/wiki/Direct_coupling_analysis), i.e. a Potts model of 21 colors on a fully connected topology, associated with the multiple sequence alignment in `test/PF00018.fasta`. This file has been downloaded from [Pfam](http://pfam.xfam.org/) database. The command
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
In the early stage of the running, the algorithm performs a reweighting of the original sequences and modifies their statistical significance as explained [here](https://www.pnas.org/content/108/49/E1293); the sequence similarity threshold can be set using `-l` input flag. Besides, it is also possible to give the set of weights in a file using the flag `-w`. Then, the algorithm computes the data statistics (one and two-site frequencies) and eventually corrects the empirical moments introducing a pseudo-count (read [here](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.012132) for details) whose value is set to `1/Meff` and can be modified using `-d`.
Finally, the algorithm iteratively updates the parameters `(J,h)` of the Potts model by standard gradient ascent using the default learning rates, or the given values `-u lrate(J) -v lrate(h)`. The model observables are computed from 1000 independent MC chains which, at equilibrium, sample 50 configurations each. The number of MC chains and the number of sampled configurations can be modified by `-s` and `-n` respectively. `adabmDCA` allows for advanced tuning of the sampling, see `Advanced options/Tuning the Markov Chain Monte Carlo`.
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
  
## Input/Output files

`adabmDCA` takes as input the set of configuration in [FASTA](https://en.wikipedia.org/wiki/FASTA_format). By default, the program assumes to read protein sequences; alternative alphabets can be set by using the `-b` flag followed by the letter `n` for RNA sequences or `i` for Ising variables. 
Notice that the spin configurations must be reported as `0/1` sequences and the `0` character is interpreted as `-1` by the program.

Every X iterations (X is specified by `-m X`), and at convergence, `adabmDCA` prints to file the one-site and two-site statistics of the data/model as well as the set of parameters of the Potts model and the (average corrected) Frobenius norms associated with them. The files are written respectively as `First_mom_label.dat`, `Sec_mom_label.dat`, and `Parameters_label.dat`, `Scores_label.dat` where `label` is a string that can be modified by the option `-k label`. 
The file `First_mom_label.dat` contains the list of one-site frequencies using the format:
```
i a m_MSA(i,a) m_model(i,a)
```
where `i` and `a` are the site and color indices, `m_X(i,a)` is the one-site frequency computed using `X`. For the file `Sec_mom_label.dat` we use:
```
i j a b s_MSA(i,j,a,b) s_model(i,j,a,b) c_MSA(i,j,a,b) c_model(i,j,a,b)
```
where `i j` run over the site indices and `a b` over all pairs of colors; `s_X(i,j,a,b)` and `c_X(i,j,a,b)` are the two-site frequency and the second connected moment, accordin to `X`, computed for positions `i j` and color `a_i = a`,`a_j = b`. The parameters file is structured as:
```
J i j a b value
...
h i a value
```
while the scores file is organized as:
```
i j score
```
The gap state `-` is always mapped to the `0` color.

### Additional input - MSA statistics

Alternatively to the set of configurations, `adabmDCA` can directly read a set of given empirical statistics collected in a file (use `-q` here). In the latter case, the input file must be formatted as
```
s i j symbol_i symbol_j value
...
m i symbol_i value
```
where the `s` rows contain the two-site frequencies of the `i j` sites for colors `symbol_i symbol_j` and the `m` lines the empirical one-site frequencies of site `i` for color `symbol_i`. 

## Initialization of the parameters

By default, the parameters of the DCA model are all initialized to 0. However, it is possible to initialize the Boltzmann machine to the set of parameters of the profile model (read [here](https://iopscience.iop.org/article/10.1088/1361-6633/aa9965/meta)) using the flag `-I` or to a set of given parameters stored in a file, using the flag `-p file`. The `file` must be formatted as:
```
J i j a b value
...
h i a value

```

## Advanced options

Here is a list of auxiliary functions that can be performed by `adabmDCA`.

### Tuning the Markov Chain Monte Carlo

#### Equilibration test

The standard run of this implementation of the Boltzmann machine learning ensures that the model statistics is estimated using a MCMC sampling performed at equilibrium. To do this, the number of MC sweeps between each pair of sampled configurations, `Twait`, is tuned at each iteration, and the number of MC sweeps before the first collected configuration, `Teq`, is set equal to `2Twait`. This last choice is likely to provide a first well equilibrated configuration. 
Let us call the configuration of chain `i` sampled after `n` steps of the MCMC as s<sup>i</sup><sub>n</sub>(Teq + n Twait). For each iteration we compute the average value (among both chains and `n`) and deviations from them,  of the following quantities, the so-called overlaps:

  - Q<sup>ext</sup> = δ<sub> s<sup>i</sup><sub>n</sub> s<sup>k</sup><sub>n</sub> </sub> 
  - Q<sup>int,1</sup> = δ<sub> s<sup>i</sup><sub>n</sub> s<sup>i</sup><sub>n+1</sub> </sub>
  - Q<sup>int,2</sup> = δ<sub> s<sup>i</sup><sub>n</sub> s<sup>i</sup><sub>n+2</sub> </sub>
  
  If the overlap between independent chains is not similar to the overlap of two samples in the same chain, distant `2Twait`, i.e. Q<sup>ext</sup> < Q<sup>int,2</sup>, then, we increase Twait: `Twait <- Twait + 1`.
  If the overlap between independent chains is similar to the intra-chain overlap between two samples at distance `Twait`, i.e. Q<sup>ext</sup> ~ Q<sup>int,1</sup>, then, we decrease Twait: `Twait <- Twait - 1`.
  
  In this way, two consecutive samples of the same chain can be slightly correlated but for sure two configurations at distance `2Twait` are reasonably de-correlated. One may assume that even starting from an arbitrary sample and waiting `2Twait` sweeps, the final sample would be at equilibrium: this is not proven but very likely the case. Besides, if one uses persistent chains (see below), it is fair to consider the time to equilibrate equals to the de-correlation time.
  
#### Advanced settings
  
  It is possible to avoid the equilibration test and to sample at fixed `Teq` and `Twait` using the flags
  ```
  -L -e Teq -t Twait
  ```
  It is also possible to perform a persistent sampling (flag `-P`), meaning that the chains are initialized as uniformly random configurations only at the first iteration, and then kept persistent, i.e. for the next iterations each chain is initialized to the last configuration of the previous iteration. The initial configurations can be extracted from empirical samples using the flag `-Q`.
  
  The standard implementation of `adabmDCA` uses as default the Metropolis update rule, but a Gibbs sampling can be used by adding the `-G` flag.
  
### Sampling

`adabmDCA` can be easily used to sample a given model. To do this, one has to give the parameters file using the flag `-p` and set the maximum number of iterations used for the training as 0, that is `-i 0`. A FASTA or a frequency input file must be specified also for this procedure. The sampled configurations, and the energies associated with them, are stored in a FASTA file whose names must be given using the following flag. Energies are written in the sequences name.
```
-S sample_file 
```
It is also possible to multiply by a constant (an inverse-temperature in the statistical physics jargon) the set of parameters; use the flag `-J` to this purpose.

### Learning on a fixed topology

The Boltzmann machine learning assumes to learn all the possible couplings associated with a fully connected interaction graph. Alternatively, if the interaction graph is known, one may give to the programm the list of component-wise couplings to be larned. This list must be stored in a `file`, and passed to the programm using `-C file`. The file must be formatted as
```
i j a b
```
or 
```
i j a b value
```
where value is some kind of measurement associated with the tuple `(i,j,a,b)`. This latter format encompasses the possibility of a priori determined a fixed topology from a statistical measure (i.e. connected correlations) and to directly give the list to `adabmDCA` (the `value` column is not used by the Boltzmann machine learning and can be omitted).

### Compute non-fitted third order statistics

At convergence, the model returned by `adabmDCA` fits, up to the convergence error, the one-site and two-site data statistics. To facilitate the testing on the generative properties of the learned model, it is possible to give to the program a list of tuple `(i,j,k,a,b,c)` on which the data and model third order connected moments are computed and printed to an output file `Third_mom_label.dat`. The list must satisfy the format
```
i j k a b c
```
and the `file` must be given using `-T file`. The output file will contain a list of the type
```
i j k a b c third_MSA(i,j,k,a,b,c) third_model(i,j,k,a,b,c)
```
where `third_X(i,j,k,a,b,c)` is the third connected moment computed using X of the sites `i ,j ,k` for colors `a_i = a, a_j = b, a_k = c`.

### Available maximum entropy models

`adabmDCA` gives the possibility of fitting a partial set of observables associated with the `-` symbols; one possibility is the so-called `pseudo Hidden Markov Model` (see [here](https://link.aps.org/doi/10.1103/PhysRevE.102.062409))

<img src="https://latex.codecogs.com/gif.latex?\mathcal{H}_{phmm}(\boldsymbol{S}&space;|&space;\boldsymbol{J},&space;\boldsymbol{h})&space;=&space;-&space;\sum_{i,i&plus;1}J_{i,i&plus;1}(-,-)&space;\delta_{S_i,-}&space;\delta_{S_{i&plus;1},-}&space;-&space;\sum_{i}h_{i}(S_{i})" title="\mathcal{H}_{phmm}(\boldsymbol{S} | \boldsymbol{J}, \boldsymbol{h}) = - \sum_{i,i+1}J_{i,i+1}(-,-) \delta_{S_i,-} \delta_{S_{i+1},-} - \sum_{i}h_{i}(S_{i})" />

to be set using `-H` flag. A second one is the `nearest-neighbors gaps` model defined as

<img src="https://latex.codecogs.com/gif.latex?\mathcal{H}_{nng}(\boldsymbol{S}&space;|&space;\boldsymbol{J},&space;\boldsymbol{h})&space;=&space;-&space;\sum_{i}&space;h_{i}&space;(S_{i})&space;-&space;\sum_{i<j}J_{ij}(S_i,S_j)&space;(1&space;-&space;\delta_{S_{i},-})&space;(1&space;-&space;\delta_{S_{j},-})&space;-&space;\sum_{i,i&plus;1}&space;J_{i,i&plus;1}&space;(-,&space;-)&space;\delta_{S_{i},-}&space;\delta_{S_{i&plus;1},-}" title="\mathcal{H}_{nng}(\boldsymbol{S} | \boldsymbol{J}, \boldsymbol{h}) = - \sum_{i} h_{i} (S_{i}) - \sum_{i<j}J_{ij}(S_i,S_j) (1 - \delta_{S_{i},-}) (1 - \delta_{S_{j},-}) - \sum_{i,i+1} J_{i,i+1} (-, -) \delta_{S_{i},-} \delta_{S_{i+1},-}" />

whose corresponding input flag is `-N`.

### Regularizations

By default, the Boltzmann machine learning does not assume any regularization of the learned parameters. Still, it is possible to add an L1 or L2 regularizations, of strength `lambda`, using the flags `-g lambda` or `-r lambda` respectively.

### Pruning the couplings

In `adabmDCA` it is possible to run several procedures to prune the couplings of the fully connected Potts model up to a desired sparsity, specified using `-x value`. At convergence of the fully connected model, or alternatively every X iteration (use `-X value` here), `adabmDCA` computes a score associated with each non-zero coupling and prune the 1% of these parameters showing the lowest scores. 
There are three choices for the score: 
  - the empirical frequency `s_MSA(i,j,a,b)`associated with `J(i,j,a,b)` (`-U` flag)
  - the absolute value of the coupling `J(i,j,a,b)` (use `-V` flag) 
  - the symmetric Kullbak-Leibler distance between the models with or without a coupling. 
  
The recommended procedure is the last one which corresponds to an information-based and component-wise pruning, achievable using the `-W` flag. A block-wise pruning can be set using `-B` (not yet implemented). For more details about the procedure and the performances of the pruned models, see [here](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.104.024407). 

#### Remove gauge invariance

Each pruning strategy described above breaks the gauge invariance of the target Potts model. To start the learning in a fixed gauge, we propose to initially set to 0 a subset of the couplings determined as follows. For each `q x q` coupling matrix `J(i,j)`, we sort the associated connected second moments `c(i,j,a,b)`. Then for the `2q-1` couples of colors `(a.b)` having the smallest connected moments, we set `J(i,j,a,b) = 0`.





