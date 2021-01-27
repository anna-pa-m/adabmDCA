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
  vector<double> h, Gh;
  vector< vector<double> > J, decJ, GJ;
  vector<double> fm_s;
  vector< vector<double> > sm_s;
  vector<double> tm_s;
  vector< vector<int> > * tm_index;
  bool Gibbs;
  Params * params;
  double alpha, acc;  // for FIRE
  int counter;   // for FIRE
  double model_sp;
  vector< vector<int> > idx;
  vector<int> tmp_idx;
  vector<double> sorted_struct;

  Model(int _q, int _L, Params * _params, vector< vector<int> > & msa, int _ntm, vector< vector<int> > * _tm_index);

  void init_current_state(vector< vector<int> > & msa);
  int remove_gauge_freedom(vector< vector<double> > & cov);
  int initialize_parameters(vector<double> & fm);
  void initial_decimation(vector< vector<double> > & cov);
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
  int compute_errors(vector<double> & fm, vector< vector<double> > & sm, vector< vector<double> > & cov, Errs & errs);
  double pearson(vector< vector<double> > & cov);
  double update_parameters(vector<double> & fm, vector< vector<double> > & sm, int iter);
  int n_links();
  void init_decimation_variables();
  int decimate_compwise(int c, int iter);


};

#endif





