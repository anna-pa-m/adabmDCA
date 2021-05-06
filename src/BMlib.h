// BOLTZMANN MACHINE CODE - version of March, 2021
// Authors: Anna Paola Muntoni and Francesco Zamponi
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include <iostream>
#include "BMaux.h"

using namespace std;

#ifndef BMLIB_H
#define BMLIB_H


class Params {
 public:
  char * file_msa, * file_freq, * file_w , * file_params, init, * label, ctype, * file_3points, *file_cc;
  bool print_samples, Metropolis, Gibbs, nprinteq, rmgauge, dgap, gapnn, phmm, blockwise, compwise, persistent, initdata, overwrite, adapt, dec_sdkl, dec_f, dec_J;
  double sparsity, rho, w_th,  regJ1, regJ2, lrateJ, lrateh, conv, pseudocount, beta;
  int tau, seed, learn_strat, nprint, nprintfile, Teq, Nmc_starts, Nmc_config, Twait, maxiter, dec_steps, num_threads;
  
  Params();

  int read_params (int & argc, char ** argv);
  void print_learning_strategy();
  void construct_filenames(int iter, bool conv, char * par, char * par_zsum, char * ene, char * score, char * first, char * sec, char * third);
};

class Data {
 public:
  int q, L, M;
  double Meff;
  vector< vector<unsigned char> > msa;
  vector<float> w;
  vector<float> fm;
  vector< vector<float> > sm;
  vector< vector<float> > cov;
  vector<float> tm; // Pay attention. File contains connected 3rd order correlations
  vector< vector<int> > tm_index;
  Params * params;


  Data(Params * _params);

  void read_msa(); 
  void compute_w(); 
  void alloc_structures();
  void compute_empirical_statistics(); 
  void read_freq(); 

  /******************** METHODS FOR 3RD ORDER STATISTICS ***********************************************************/

  void load_third_order_indices();
  /******************** METHODS FOR OUTPUT ***********************************************************/
  
  void print_msa(char *filename);
  int print_statistics(char *file_sm, char *file_fm, char *file_tm, vector<float> & fm_s, vector< vector<float> > & sm_s, vector<float> & tm_s); 
   
};





#endif
