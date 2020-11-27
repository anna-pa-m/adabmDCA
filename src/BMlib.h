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
  char * file_msa, * file_freq, * file_w , * file_params, init, * label, * ctype, * file_3points, *file_cc, *file_samples, *file_en;
  bool Metropolis, Gibbs, nprinteq, rmgauge, dgap, gapnn, phmm, blockwise, compwise, persistent, initdata, overwrite, adapt, dec_sdkl, dec_f, dec_J;
  double sparsity, rho, w_th,  regJ1, regJ2, lrateJ, lrateh, conv, pseudocount, beta;
  int tau, seed, learn_strat, nprint, nprintfile, Teq, Nmc_starts, Nmc_config, Twait, maxiter, dec_steps;

  Params();

  int read_params (int & argc, char ** argv);
  void print_learning_strategy();
  void construct_filenames(int iter, bool conv, char * par, char * par_zsum, char * score, char * first, char * sec, char * third);
};


class Data {
 public:
  int q, L, M;
  double Meff;
  vector< vector<int> > msa;
  vector<double> w;
  vector<double> fm;
  vector< vector<double> > sm;
  vector< vector<double> > cov;
  vector<double> tm; // Pay attention. File contains connected 3rd order correlations
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
  int print_statistics(char *file_sm, char *file_fm, char *file_tm, vector<double> & fm_s, vector< vector<double> > & sm_s, vector<double> & tm_s); 
   
};





#endif