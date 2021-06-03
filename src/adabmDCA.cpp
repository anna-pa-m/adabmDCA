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
#include <iomanip>
#include "BMaux.h"
#include "BMlib.h"
#include "BMmc.h"

using namespace std;

int main(int argc, char **argv)
{

  cout << setprecision(1);
  /* START INITIALIZATION */
  cout << "****** Boltzmann machine for DCA model ******" << endl;
  Params params;
  params.read_params(argc, argv);
  srand(params.seed ? params.seed : time(NULL));
  Data data(&params);
  Stats mstat;
  params.print_learning_strategy();
  Model model(data.q, data.L, &params, &mstat, data.msa, data.tm.size(), &data.tm_index);
  model.initialize_parameters(data.fm);
  model.initial_decimation(data.cov);
  /* END INITIALIZATION */

  /* START ITERATION LOOP */
  char par[1000];
  char par_zsum[1000];
  char ene[1000];
  char score[1000];
  char first[1000];
  char sec[1000];
  char third[1000];
  char corr[1000];
  int iter = 0;
  int in_time = time(NULL);
  bool conv = (params.maxiter > 0) ? false : true;
  Errs errs;
  double lrav = params.lrateJ;
  if (!conv) {
    cout << "****** Starting learning loop ******" << endl;
    cout << "Printing output every " << params.nprint << " iterations - parameters every " << params.nprintfile << endl;
    fflush(stdout);
  }
  while (!conv) {
    if(params.ctype == 'i')
      model.sample_ising(data.msa);
    else
      model.sample(data.msa);
    model.compute_errors(data.fm, data.sm, data.cov, errs);
    if (iter % params.nprint == 0) {
      cout << setprecision(1);
      cout.unsetf(ios::fixed | ios::scientific);
      cout << "it: " << iter << " el_time: "  << int(time(NULL) - in_time) << " N: " << params.Nmc_config * params.Nmc_starts * params.num_threads << " Teq: " << params.Teq << " Twait: " << params.Twait;
      cout.setf(ios::scientific);
      cout << " merr_fm: " << errs.merrh << " merr_sm: " << errs.merrJ << " averr_fm: " << errs.averrh << " averr_sm: " << errs.averrJ << " cov_err: " << errs.errnorm;
      cout.unsetf(ios::scientific);
      cout << setprecision(2);
      cout << " corr: " << model.pearson(data.cov) << " sp: " << model.model_sp << " lrav: " << lrav << endl;
      fflush(stdout);
    }
    if (iter > 0 && (iter % params.nprintfile == 0)) {
      params.construct_filenames(iter, conv, par, par_zsum, ene, corr, score, first, sec, third);
      print_frobenius_norms(model.h, model.J, model.L, model.q, score, par_zsum);
      model.print_model(par);
      if(params.print_samples) 
        model.print_samples(ene);
      if (data.tm.size() > 0)
        model.compute_third_order_correlations();
      data.print_statistics(sec, first, third, corr, model.mstat->corr, model.mstat->fm_s, model.mstat->sm_s, model.mstat->tm_s);
    }
    lrav = model.update_parameters(data.fm, data.sm, iter);
    if (iter > 0 && (params.compwise || params.blockwise) && (errs.errnorm < params.conv || iter % params.dec_steps == 0)) {
      // Print converged parameters before decimation
      params.construct_filenames(iter, conv, par, par_zsum, ene, corr, score, first, sec, third);
      print_frobenius_norms(model.h, model.J, model.L, model.q, score, par_zsum);
      model.print_model(par);
      if (data.tm.size() > 0)
        model.compute_third_order_correlations();
      data.print_statistics(sec, first, third, corr, model.mstat->corr, model.mstat->fm_s, model.mstat->sm_s, model.mstat->tm_s);
      // Then decimate
      int aux = ceil(model.n_links() / 100);
      if(params.compwise) {
        if(params.ctype == 'i')
          model.decimate_ising(aux, iter);
        else
          model.decimate_compwise(aux, iter);
      } else if(params.blockwise) {
          if(params.ctype == 'i')
            model.decimate_ising(aux,iter);
          else
            model.decimate_blockwise(iter);
      }
    }
    if (errs.errnorm < params.conv && !params.compwise && !params.blockwise) {
      conv = true;
      cout << "Reached convergence of error, end of learning" << endl;
    }
    if (model.model_sp >= params.sparsity && params.sparsity > 0 && errs.errnorm < params.conv) {
      conv = true;
      cout << "Reached convergence of error and desired sparsity, end of learning" << endl;
    }
    if (iter >= params.maxiter) {
      cout << "Reached maximum number of iterations, end of learning" << endl;
      conv = true;
    }
    iter++;
  }
  /* END ITERATION LOOP */

  /* FINAL OPERATIONS */
  cout << "****** Final sampling ******" << endl;
  fflush(stdout);
  if (!params.maxiter) {
    bool eqmc = false;
    while (!eqmc) {
      if(params.ctype == 'i'){
        eqmc = model.sample_ising(data.msa);
      } else
        eqmc = model.sample(data.msa);
      model.compute_errors(data.fm, data.sm, data.cov, errs);
      cout << "N: " << params.Nmc_config * params.Nmc_starts * params.num_threads << " Teq: " << params.Teq << " Twait: " << params.Twait;
      cout.setf(ios::scientific);
      cout << " merr_fm: " << errs.merrh << " merr_sm: " << errs.merrJ << " averr_fm: " << errs.averrh << " averr_sm: " << errs.averrJ << " cov_err: " << errs.errnorm;
      cout.unsetf(ios::scientific);
      cout << " corr: " << model.pearson(data.cov) << endl;
      fflush(stdout);
    }
  }
  else {
    if(params.ctype == 'i'){
      model.sample_ising(data.msa);
    } else
      model.sample(data.msa);
  }

  params.construct_filenames(iter, conv, par, par_zsum, ene,  corr, score, first, sec, third);
  print_frobenius_norms(model.h, model.J, model.L, model.q, score, par_zsum);
  model.print_model(par);
  if (data.tm.size() > 0)
    model.compute_third_order_correlations();
  data.print_statistics(sec, first, third, corr, model.mstat->corr, model.mstat->fm_s, model.mstat->sm_s, model.mstat->tm_s);
  if(params.print_samples) {
    if(params.ctype == 'i')
      model.print_samples_ising(ene);
    else
      model.print_samples(ene);
  }
  cout << "****** Execution completed ******" << endl;

  fflush(stdout);
  /* FINAL OPERATIONS */

  return 0;
}
