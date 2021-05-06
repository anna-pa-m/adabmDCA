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

struct Stats {

  vector<float> fm_s;
  vector< vector<float> > sm_s;
  vector<float> tm_s;
  vector<int> qs;
  vector < vector <int> > qs_t;
  vector < vector <unsigned char> > old_state1;
  vector < vector <unsigned char> > old_state2;
  vector < vector <unsigned char> > oldold_state1;
  vector < vector <unsigned char> > oldold_state2;
  vector < vector <unsigned char> > x1i;
  vector < vector <unsigned char> > x2i;
  valarray<float> corr;
  vector< vector < vector <unsigned char> > > synth_msa;
  vector< vector < vector <unsigned char> > > curr_state;
  omp_lock_t lock;

};

class Model {
  public:
  int q, L;
  vector<float> h, Gh;
  vector< vector<float> > J, GJ;
  vector< vector<unsigned char > > decJ;
  vector< vector<int> > * tm_index;
  bool Gibbs;
  Params * params;
  Stats * mstat;
  double alpha, acc;  // for FIRE
  int counter;   // for FIRE
  double model_sp;
  vector< vector<int> > idx;
  vector<int> tmp_idx;
  vector<float> sorted_struct;

  Model(int _q, int _L, Params * _params, Stats * _mstat, vector< vector<unsigned char> > & msa, int _ntm, vector< vector<int> > * _tm_index);


  void update_synth_msa(vector<unsigned char> & x1, vector<unsigned char> & x2);
  void init_model_stat(int ntm);
  void init_current_state(vector< vector<unsigned char> > & msa);
  void init_current_state_ising(vector< vector<unsigned char> > & msa);
  int remove_gauge_freedom(vector< vector<float> > & cov);
  int initialize_parameters(vector<float> & fm);
  void initial_decimation(vector< vector<float> > & cov);
  int print_model(char *filename);
  double prof_energy(vector<unsigned char> & seq);
  double DCA_energy(vector<unsigned char> & seq);
  double energy(vector<unsigned char> & seq);
  void metropolis_step(vector<unsigned char> & x);
  void metropolis_step_ising(vector<unsigned char> & x);
  void gibbs_step(vector<unsigned char> & x);
  void gibbs_step_ising(vector<unsigned char> & x);
  void MC_sweep(vector<unsigned char> & x);
  void MC_sweep_ising(vector<unsigned char> & x);
  void mc_chain_ising(vector<unsigned char> & x1, vector<unsigned char> & x2, int s);
  void mc_chain(vector<unsigned char> & x1, vector<unsigned char> & x2, int s);
  void update_overlap(vector <int> & qs);
  void update_corr(int i, int value);
  bool sample(vector< vector<unsigned char> > & msa);
  bool sample_ising(vector< vector<unsigned char> > & msa);
  void init_statistics();
  void update_statistics();
  void update_statistics_ising();
  void update_statistics_lock(vector<unsigned char> & x, FILE * fp, FILE * fe);
  void update_tm_statistics(vector<unsigned char> & x);
  void compute_third_order_correlations();
  int compute_errors(vector<float> & fm, vector< vector<float> > & sm, vector< vector<float> > & cov, Errs & errs);
  double pearson(vector< vector<float> > & cov);
  double update_parameters(vector<float> & fm, vector< vector<float> > & sm, int iter);
  int n_links();
  void init_decimation_variables();
  int decimate_compwise(int c, int iter);
  void print_samples(char * filename);
  void print_samples_ising(char * filename);


};


#endif





