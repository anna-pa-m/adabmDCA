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
#include <fstream>
#include <vector>
#include <valarray>
#include "BMlib.h"

using namespace std;

#ifndef BMMC_H
#define BMMC_H

struct Errs {
  double errnorm;
  double averrh;
  double averrJ; 
  double merrh;
  double merrJ;
};

class Model {
  public:
  int q, L;
  vector< vector<int> > curr_state;
  vector<float> h, Gh;
  vector< vector<float> > J, GJ;
  vector< vector<unsigned char > > decJ;
  vector<float> fm_s;
  vector< vector<float> > sm_s;
  vector<float> tm_s;
  vector< vector<int> > * tm_index;
  bool Gibbs;
  Params * params;
  double alpha, acc;  // for FIRE
  int counter;   // for FIRE
  double model_sp;
  vector< vector<int> > idx;
  vector<int> tmp_idx;
  vector<float> sorted_struct;

  Model(int _q, int _L, Params * _params, vector< vector<int> > & msa, int _ntm, vector< vector<int> > * _tm_index);

  void init_current_state(vector< vector<int> > & msa);
  int remove_gauge_freedom(vector< vector<float> > & cov);
  int initialize_parameters(vector<float> & fm);
  void initial_decimation(vector< vector<float> > & cov);
  int print_model(char *filename);
  double prof_energy(vector<int> & seq);
  double DCA_energy(vector<int> & seq);
  double energy(vector<int> & seq);
  void metropolis_step(vector<int> & x);
  void gibbs_step(vector<int> & x);
  void MC_sweep(vector<int> & x);
  valarray<int> mc_chain(vector<int> & x1, vector<int> & x2, valarray<double> & corr);
  bool sample(vector< vector<int> > & msa);
  void init_statistics();
  void update_statistics(vector<int> & x, FILE * fp, FILE * fe);
  void update_tm_statistics(vector<int> & x);
  void compute_third_order_correlations();
  int compute_errors(vector<float> & fm, vector< vector<float> > & sm, vector< vector<float> > & cov, Errs & errs);
  double pearson(vector< vector<float> > & cov);
  double update_parameters(vector<float> & fm, vector< vector<float> > & sm, int iter);
  int n_links();
  void init_decimation_variables();
  int decimate_compwise(int c, int iter);


};


#endif





