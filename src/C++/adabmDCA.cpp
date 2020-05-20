// BOLTZMANN MACHINE CODE - version of March 13, 2020


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
#include "BMlib.h"
#include "BMmc.h"

using namespace std;


int main(int argc, char ** argv) {

  /* START INITIALIZATION */
  cout << "****** Boltzmann machine for DCA model ******" << endl;
  Params params;
  params.read_params(argc,argv);
  srand(params.seed ? params.seed : time(NULL));
  Data data(&params);
  params.print_learning_strategy();
  Model model(data.q,data.L,&params,data.msa,data.tm.size(),&data.tm_index);
  model.initialize_parameters(data.fm);
  model.initial_decimation(data.cov);
  /* END INITIALIZATION */

  /* START ITERATION LOOP */
  char par[1000];
  char par_zsum[1000];
  char score[1000];
  char first[1000];
  char sec[1000];
  char third[1000];
  int iter = 0;
  long in_time = time(NULL);
  bool conv = (params.maxiter > 0) ?  false : true;
  Errs errs;
  double lrav=params.lrateJ;
  if (!conv) {
    cout << "****** Starting learning loop ******" << endl;
    cout << "Printing output every " << params.nprint << " iterations - parameters every " << params.nprintfile << endl;
    fflush(stdout);
  }
  while(!conv) { 
    model.sample(data.msa);
    model.compute_errors(data.fm,data.sm,data.cov,errs);
    if(iter % params.nprint == 0) {
      fprintf(stdout, "it: %i el_time: %li N: %i Teq: %i Twait: %i merr_fm: %.1e merr_sm: %.1e averr_fm: %.1e averr_sm: %.1e cov_err: %.1e corr: %.2f sp: %.1e lrav: %.1e\n", iter, time(NULL)-in_time, params.Nmc_config * params.Nmc_starts, params.Teq, params.Twait, errs.merrh, errs.merrJ, errs.averrh, errs.averrJ, errs.errnorm, model.pearson(data.cov), model.model_sp, lrav);
      fflush(stdout);
    }
    if(iter > 0 && (iter % params.nprintfile == 0)) {
      params.construct_filenames(iter, conv, par, par_zsum, score, first, sec, third);
      print_frobenius_norms(model.h,model.J,model.L,model.q,score,par_zsum);
      model.print_model(par);
      if(data.tm.size()>0)
	model.compute_third_order_correlations();
      data.print_statistics(sec, first, third, model.fm_s, model.sm_s, model.tm_s);
    }
    lrav=model.update_parameters(data.fm,data.sm,iter);
    if(iter > 0 && params.compwise && (errs.errnorm<params.conv || iter % params.dec_steps == 0)) {
      // Print converged parameters before decimation
      params.construct_filenames(iter, conv, par, par_zsum, score, first, sec, third);
      print_frobenius_norms(model.h,model.J,model.L,model.q,score,par_zsum);
      model.print_model(par);
      if(data.tm.size()>0)
	model.compute_third_order_correlations();
      data.print_statistics(sec, first, third, model.fm_s, model.sm_s, model.tm_s);
      // Then decimate
      int aux = ceil( model.n_links() / 100);
      model.decimate_compwise(aux,iter);
    }
    if(errs.errnorm<params.conv && !params.compwise) {
      conv = true;
      fprintf(stdout,"Reached convergence of error, end of learning\n");
    }
    if(model.model_sp >= params.sparsity && params.sparsity > 0 && errs.errnorm<params.conv) {
      conv = true;
      fprintf(stdout,"Reached convergence of error and desired sparsity, end of learning\n");
    }
    if (iter >= params.maxiter) {
      fprintf(stdout,"Reached maximum number of iterations, end of learning\n");
      conv = true;
    }
    iter++;
  }
  /* END ITERATION LOOP */

  /* FINAL OPERATIONS */
  fprintf(stdout, "****** Final sampling ******\n");
  fflush(stdout);
  if(!params.maxiter) {
    bool eqmc = false;
    while(!eqmc) {
      eqmc = model.sample(data.msa);
      model.compute_errors(data.fm,data.sm,data.cov,errs);
      fprintf(stdout, "N: %i Teq: %i Twait: %i merr_fm: %.1e merr_sm: %.1e averr_fm: %.1e averr_sm: %.1e cov_err: %.1e corr: %.2f \n", params.Nmc_config * params.Nmc_starts, params.Teq, params.Twait, errs.merrh, errs.merrJ, errs.averrh, errs.averrJ, errs.errnorm, model.pearson(data.cov));
      fflush(stdout);
    }
  } else {
    model.sample(data.msa);
  }

  params.construct_filenames(iter, conv, par, par_zsum, score, first, sec, third);
  print_frobenius_norms(model.h,model.J,model.L,model.q,score, par_zsum);
  model.print_model(par);
  if(data.tm.size()>0)
    model.compute_third_order_correlations();
  data.print_statistics(sec, first, third, model.fm_s, model.sm_s, model.tm_s);
  fprintf(stdout, "****** Execution completed ******\n");
  fflush(stdout);
  /* FINAL OPERATIONS */

  return 0;

}














