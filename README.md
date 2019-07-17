# adabmDCA 
This is an (adaptive) Boltzmann machine for Direct Coupling Analysis model. It is written in C language.

# Installation

To properly install adabmDCA, run
```
make
```
on adabmDCA/src folder. It suffices a gcc compiler.

# Usage
This algorithm is able to infer several Max-Ent statistical models of Potts variables given a set of observations. More precisly, it is able to learn the couplings and the fields parameters of a set of Direct Coupling Analysis (DCA) models given a Multiple Sequence Alignment (MSA) of protein, RNA, binary sequences. The learning is performed via a gradient ascent of the likelihood of the data, computed via a Monte Carlo sampling, either using Metropolis or Gibbs sampling.
All the possible features (model choices, input and output files) implemented in adabmDCA can be shown typing
```
./adabmDCA -h
```
## Basic input/output files


