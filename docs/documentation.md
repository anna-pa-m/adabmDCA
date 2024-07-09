

# Documentation 

## Basic run
Let us infer a [DCA model](https://en.wikipedia.org/wiki/Direct_coupling_analysis), i.e. a Potts model of 21 colors on a fully connected topology, associated with the multiple sequence alignment in `test/PF00018.fasta`. This file has been downloaded from [Pfam](http://pfam.xfam.org/) database. The command
```
./adabmDCA -f ../test/PF00018.fasta -z 1 -m 500 -c 1e-2 -k PF00018 
```
will run a Boltzmann learning, using the default MCMC, printing to files every 500 iterations (as specified in the option `-m`) the temporary parameters of the model, and both data and model statistics. Output files are saved in the folder `PF00018` as specified in the `-k` flag. The run will stop when either we reach the maximum number of iterations (default: 2000) or we reach convergence, i.e. the maximum fitting error among all the connected covariances is 1e-2 (as specified in the command line using `-c`).
The program will output a list of information concerning the input data and the adopted learning strategy:
```
***** Boltzmann machine for DCA model ******
Running from user USERNAME on host HOSTNAME
Output files in PF00018 folder
Seed of random number generator: 1720109177
****** Initializing data structures ******
Using alphabet: -ACDEFGHIKLMNPQRSTVWY
Reading MSA from ../test/PF00018.fasta
Reading alignment completed.
M = 55784 L = 48 q = 21
Computing weights...done
Computing empirical statistics...Meff: 7170.02
No three-point correlation indices specified
****** Initializing model ******
Performing Metropolis-Hastings MC
Adaptive sampling time. Initial sampling time: 10 initial equilibration time: 20
Using 1000 seeds and tot. number of points 50000
MC chains are randomly initialized at the beginning of each iteration
Learning strategy: 
Using standard gradient descent with constant learning rate (for J 0.05 for h 0.05)
Using pseudo-count: 0.00
Zero-parameters initialization...done
Using full Potts model
Sparsity after initialization: 0.00
****** Starting learning loop ******
Printing output every 1 iteration - summary every 500

```
In the early stage of the running, the algorithm performs a reweighting of the original sequences and modifies their statistical significance as explained [here](https://www.pnas.org/content/108/49/E1293); the sequence similarity threshold can be set using `-l` input flag. Besides, it is also possible to give the set of weights in a file using the flag `-w`. Then, the algorithm computes the data statistics (one and two-site frequencies) and eventually corrects the empirical moments introducing a pseudo-count (read [here](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.012132) for details) whose value is set to `1/Meff` and can be modified using `-d`.
Finally, the algorithm iteratively updates the parameters `(J,h)` of the Potts model by standard gradient ascent using the default learning rates, or the given values `-u lrate(J) -v lrate(h)`. The model observables are computed from 1000 independent MC chains which, at equilibrium, sample 50 configurations each. The number of MC chains and the number of sampled configurations can be modified by `-s` and `-n` respectively. `adabmDCA` allows for advanced tuning of the sampling, see `Advanced options/Tuning the Markov Chain Monte Carlo`.
At every step (the frequency can be specified by the `-z` flag)  it writes to the `stdout` an update of the learning, like the following:
```
it: 1 el_time: 8 Teq: 20 Twait: 10  lrav: 5.00e-02 sp: 0.00e+00 cov_err: 1.99e-01 pears_act: 3.13e-03 pears_all: 3.13e-03
```
where:

  - `it:` gives the number of the current iteration
  - `el_time:` is the running time (in seconds)
  - `Teq:` is the equilibration time (in MC sweep) of each MC chain
  - `Twait:` is the de-correlation time (in MC sweep) between two sampled configurations, for all the chains  
  - `lrav:` average learning rate
  - `sp:` sparsity of the couplings parameters, i.e. the number of couplings fixed to zero divided by q<sup>2</sup>L(L-1)/2
  - `cov_err:` maximum value of the fitting error for the connected covariances (used for checking the convergence)
  - `pears_act:` Pearson correlation coefficient between the model and the data connected covariances associated with the active (or not pruned) couplings
  - `pears_all:` Pearson correlation coefficient between the model and the data connected covariances 


## Input/Output files

`adabmDCA` takes as input the set of configurations in [FASTA](https://en.wikipedia.org/wiki/FASTA_format). By default, the program assumes to read protein sequences; alternative alphabets can be set by using the `-b` flag followed by the letter `n` for RNA sequences or `i` for Ising variables. 
Notice that the spin configurations must be reported as `0/1` sequences and the `0` character is interpreted as `-1` by the program.

Every X iteration (X is specified by `-m X`), and at convergence, `adabmDCA` prints to file the one-site and two-site statistics of the data/model as well as the set of parameters of the Potts model and the (average corrected) Frobenius norms associated with them. To restore a stopped run from the last temporarily saved files, use the flag `--restore` as
```
./adabmDCA -f ../test/PF00018.fasta -k PF00018 --restore
```
The output files take the names `First_mom_$label.dat`, `Sec_mom_$label.dat`, and `Parameters_label.dat`, `Scores_label.dat` where `label` is a string that can be modified by the option `-k label`. 
The file `First_mom_$label.dat` contains the list of one-site frequencies using the format:
```
i a m_MSA(i,a) m_model(i,a)
```
where `i` and `a` are the site and color indices, `m_X(i,a)` is the one-site frequency computed using `X`. For the file `Sec_mom_label.dat` we use:
```
i j a b s_MSA(i,j,a,b) s_model(i,j,a,b) c_MSA(i,j,a,b) c_model(i,j,a,b)
```
where `i j` runs over the site indices and `a b` over all pairs of colors; `s_X(i,j,a,b)` and `c_X(i,j,a,b)` are the two-site frequency and the second connected moment, according to `X`, computed for positions `i j` and color `a_i = a`,`a_j = b`. The parameters file is structured as:
```
J i j a b value
...
h i a value
```
while the scores file is organized as:
```
i j score
```
The gap state `-` is always mapped to the `0` color. Inactive (or pruned) couplings are not printed in the output. 

### Additional input - MSA statistics

Alternatively to the set of configurations, `adabmDCA` can directly read a set of given empirical statistics collected in a file (use `-q` here). In the latter case, the input file must be formatted as
```
s i j symbol_i symbol_j value
...
m i symbol_i value
```
where the `s` rows contain the two-site frequencies of the `i j` sites for colors `symbol_i symbol_j` and the `m` lines the empirical one-site frequencies of site `i` for color `symbol_i`. 

## Initialization of the parameters

By default, the parameters of the DCA model are all initialized to 0. However, using the flag `I`, it is possible to initialize the Boltzmann machine to the set of parameters of the profile model (read [here](https://iopscience.iop.org/article/10.1088/1361-6633/aa9965/meta)) for the fields, and the (rescaled, by argument) covariance matrix for the couplings. Alternatively, `adabmDCA` can consider a set of parameters stored in a file, using the flag `-p file`. The `file` must be formatted as:
```
J i j a b value
...
h i a value

```
`adabmDCA` assumes zero (inactive or pruned) couplings for the missing rows in the input parameter files, and these are not updated during the training.

## Advanced options

Here is a list of auxiliary functions that can be performed by `adabmDCA`.

### Tuning the Markov Chain Monte Carlo

#### Equilibration test

The standard run of this implementation of the Boltzmann machine learning ensures that the model statistics are estimated using an MCMC sampling performed at equilibrium. To do this, the number of MC sweeps between each pair of sampled configurations, `Twait`, is tuned at each iteration, and the number of MC sweeps before the first collected configuration, `Teq`, is set equal to `2Twait`. This last choice is likely to provide a first well-equilibrated configuration. 
Let us call the configuration of chain `i` sampled after `n` steps of the MCMC as s<sup>i</sup><sub>n</sub>(Teq + n Twait). For each iteration, we compute the average value (among both chains and `n`) and deviations from them,  of the following quantities, the so-called overlaps:

  - Q<sup>ext</sup> = δ<sub> s<sup>i</sup><sub>n</sub> s<sup>k</sup><sub>n</sub> </sub> 
  - Q<sup>int,1</sup> = δ<sub> s<sup>i</sup><sub>n</sub> s<sup>i</sup><sub>n+1</sub> </sub>
  - Q<sup>int,2</sup> = δ<sub> s<sup>i</sup><sub>n</sub> s<sup>i</sup><sub>n+2</sub> </sub>
  
  If the overlap between independent chains is not similar to the overlap of two samples in the same chain, distant `2Twait`, i.e. Q<sup>ext</sup> < Q<sup>int,2</sup>, then, we increase Twait: `Twait <- 2Twait`.
  If the overlap between independent chains is similar to the intra-chain overlap between two samples at a distance `Twait`, i.e. Q<sup>ext</sup> ~ Q<sup>int,1</sup>, then, we decrease Twait by computing the average between the current value of the waiting time and the value of `Twait` before the last increasing step. This guarantees to keep the waiting time bounded in the correct interval of values within the learning process.
  
  Finally, two consecutive samples of the same chain can be slightly correlated but for sure two configurations at distance `2Twait` are reasonably de-correlated. One may assume that even starting from an arbitrary sample and waiting `2Twait` sweeps, the final sample would be at equilibrium: this is not proven but very likely the case. Besides, if one uses persistent chains (see below), it is fair to consider the time to equilibrate equal to the de-correlation time.
  
#### Advanced settings
  
  It is possible to avoid the equilibration test and to sample at fixed `Teq` and `Twait` using the flags
  ```
  -L -e Teq -t Twait
  ```
  It is also possible to perform a persistent sampling (flag `-P`), meaning that the chains are initialized as uniformly random configurations only at the first iteration, and then kept persistent, i.e. for the next iterations each chain is initialized to the last configuration of the previous iteration. The initial configurations can be extracted from empirical samples using the flag `-Q`.
  In a restored run, the first configurations are read from the `Last_configurations_$label.dat` file of the results folder.
  
  The standard implementation of `adabmDCA` uses as default the Metropolis update rule, but a Gibbs sampling can be used by adding the `-G` flag.
  
### Sampling

`adabmDCA` can be easily used to sample a given model. To do this, one has to give the parameters file using the flag `-p` and set the maximum number of iterations used for the training as 0, that is `-i 0`. A FASTA or a frequency input file must be specified also for this procedure. The sampled configurations are stored in a FASTA file if the flag `-S` has been given in input. Energies are written in the sequence name.

It is also possible to multiply by a constant (an inverse-temperature in the statistical physics jargon) the set of parameters; use the flag `-J` and `-H` for this purpose for rescaling the couplings or the fields respectively.

### Learning on a fixed topology

The Boltzmann machine learning assumes to learn all the possible couplings associated with a fully connected interaction graph. Alternatively, if the interaction graph is known, one may give to the program the list of component-wise couplings to be learned. This list must be stored in a `file`, and passed to the program using `-C file`. The file must be formatted as
```
i j a b
```
or 
```
i j a b value
```
where value is some kind of measurement associated with the tuple `(i,j,a,b)`. This latter format encompasses the possibility of a priori determined a fixed topology from a statistical measure (i.e. connected correlations) and directly gives the list to `adabmDCA` (the `value` column is not used by the Boltzmann machine learning and can be omitted).
If the topology is known, and furthermore a set of initial coupling matrices is given, use the `-p` option.

### Compute non-fitted third-order statistics

At convergence, the model returned by `adabmDCA` fits, up to the convergence error, the one-site and two-site data statistics. To facilitate the testing on the generative properties of the learned model, it is possible to give the program a list of tuples `(i,j,k,a,b,c)` on which the data and model third-order connected moments are computed and printed to an output file `Third_mom_label.dat`. The list must satisfy the format
```
i j k a b c
```
and the `file` must be given using `-T file`. The output file will contain a list of the type
```
i j k a b c third_MSA(i,j,k,a,b,c) third_model(i,j,k,a,b,c)
```
where `third_X(i,j,k,a,b,c)` is the third connected moment computed using X of the sites `i ,j ,k` for colors `a_i = a, a_j = b, a_k = c`.

### Available maximum entropy models

`adabmDCA` gives the possibility of fitting a partial set of observables associated with the `-` symbols; one possibility is the so-called `pseudo-Hidden Markov Model` (see [here](https://link.aps.org/doi/10.1103/PhysRevE.102.062409))

```math
\mathcal{H}_{phmm}(\boldsymbol{S} | \boldsymbol{J}, \boldsymbol{h}) = - \sum_{i,i+1}J_{i,i+1}(-,-) \delta_{S_i,-} \delta_{S_{i+1},-} - \sum_{i}h_{i}(S_{i})
```

to be set using `-D` flag. A second one is the `nearest-neighbors gaps` model defined as

```math
\mathcal{H}_{nng}(\boldsymbol{S} | \boldsymbol{J}, \boldsymbol{h}) = - \sum_{i} h_{i} (S_{i}) - \sum_{i < j} J_{ij} (S_i, S_j ) (1 - \delta_{S_{i}, -}) (1 - \delta_{S_{j}, -} ) - \sum_{i,i+1} J_{i,i+1} (-, -) \delta_{S_{i}, -} \delta_{S_{i+1}, -}
```
whose corresponding input flag is `-N`.

#### Data-driven model
If a set of configurations and an experimental score associated with them are known, one may use the training procedure described [here](https://www.nature.com/articles/srep37812). The experimental data can be given as input by using the flag `-E $file`: `adabmDCA` will look for a sequences file in FASTA format named `$file.fasta` and a text file containing the corresponding fitnesses `$file.fit`.

### Regularizations

By default, the Boltzmann machine learning does not assume any regularization of the learned parameters. Still, it is possible to add $L_1$ or $L_2$ regularizations, of strength `lambda`, using the flags `-g lambda` or `-r lambda` respectively.

### Pruning/activating the couplings

In `adabmDCA` it is possible to run several procedures to prune the couplings of the fully connected Potts model up to a desired sparsity, or, alternatively, to activate initially-set zero couplings from a profile model. The required sparsity can be specified using `-x value`. For a sparsity close to the target value, `adabmDCA` uses both pruning and activation schemes to reach convergence.

#### Pruning 

At the convergence of the fully connected model, or alternatively every X iteration (use `-X value` here), `adabmDCA` computes a score associated with each non-zero coupling and prunes the 1% of these parameters showing the lowest scores. 
There are three choices for the score: 
  - the empirical frequency `s_MSA(i,j,a,b)`associated with `J(i,j,a,b)` (`-U` flag)
  - the absolute value of the coupling `J(i,j,a,b)` (use `-V` flag) 
  - the symmetric Kullbak-Leibler distance between the models with or without a coupling. 
  
The recommended procedure is the last one which corresponds to an information-based and component-wise pruning, achievable using the `-W` flag. For more details about the procedure and the performances of the pruned models, see [here](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.104.024407).

#### Remove gauge invariance

Each pruning strategy described above breaks the gauge invariance of the target Potts model. To start the learning in a fixed gauge, we propose to initially set to 0 a subset of the couplings determined as follows. For each `q x q` coupling matrix `J(i,j)`, we sort the associated connected second moments `c(i,j,a,b)`. Then for the `2q-1` couples of colors `(a.b)` having the smallest connected moments, we set `J(i,j,a,b) = 0`.

### Activation

From a profile model or an initial sparse model, `adabmDCA` activates at each iteration 0.1% of the inactive parameters up to the sparsity required through the `-x` flag. These couplings are chosen as the ones associated with the largest gap between the model and the data two-site statistics. For more details, see [here](https://doi.org/10.1093/nar/gkae289).







