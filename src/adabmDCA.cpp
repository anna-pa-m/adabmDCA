// BOLTZMANN MACHINE CODE - version of October, 2023
// Authors: Anna Paola Muntoni and Francesco Zamponi

#include <sys/stat.h>
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
#include <limits.h>
#include <unistd.h>

using namespace std;
#ifndef HOST_NAME_MAX
#define HOST_NAME_MAX 255
#endif

#ifndef LOGIN_NAME_MAX
#define LOGIN_NAME_MAX 255
#endif

int main(int argc, char **argv)
{

  cout << setprecision(1);
  /* START INITIALIZATION */
  cout << "****** Boltzmann machine for DCA model ******" << endl;
  char hostname[HOST_NAME_MAX];
  char username[LOGIN_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  getlogin_r(username, LOGIN_NAME_MAX);
  cout << "Running from user " << username << " on host " << hostname << endl;
  Params params;
  params.read_params(argc, argv);
  long myseed = params.seed ? params.seed : time(NULL);
  srand(myseed);
  cout << "Output files in " << params.label << " folder" << endl;
  if (mkdir(params.label, 0777) == -1)
    cerr << "Is the folder " << params.label << " already present?" << endl;
  else
    cout << "Directory created" << endl;
  fflush(stdout);

  cout << "Seed of random number generator: " << myseed << endl;
  Data data(&params);
  Stats mstat;
  params.print_learning_strategy();
  Model model(data.q, data.L, &params, &mstat, data.msa, data.tm.size(), &data.tm_index);
  model.initialize_parameters(data.fm, data.cov);
  if(!params.restore_flag && !params.file_params)
    model.initial_decimation(data.cov);
  Data_e data_e(&params, data.q, data.L, data.Meff);
  char ene_fit[1000];
  if (params.file_msa_e)
  {
    if (data.L != data_e.L)
    {
      cout << "Error: the two alignments have different length" << endl;
      exit(EXIT_FAILURE);
    }
    data_e.energy = model.energy(data_e.msa);
    data_e.set_tg_energy();
    sprintf(ene_fit, "%s/Fitness_init_%s.dat", params.label, params.label);
    data_e.print_energy(ene_fit);
    cout << "Spearman coefficient for energy fitting: " << spearman(data_e.fitness, data_e.energy) << endl;
  }

  data_e.compute_grad();

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
  char lchain[1000];
  int iter = 0;
  int in_time = time(NULL);
  bool conv = (params.maxiter > 0) ? false : true;
  Errs errs;
  double lrav = params.lrateJ;

  if (!conv)
  {

    cout << "****** Starting learning loop ******" << endl;
    cout << "Printing output every " << params.nprint << " iterations - summary every " << params.nprintfile << endl;
  }
  while (!conv)
  {

    if (strcmp(params.ctype, "i") == 0)
      model.sample_ising(data.msa);
    else
      model.sample(data.msa);
    model.compute_errors(data.fm, data.sm, data.cov, errs);
    if (iter % params.nprint == 0)
    {
      double spe(0.), stde(0.), peae(0.);
      if (params.file_msa_e)
      {
        data_e.energy = model.energy(data_e.msa);
        spe = spearman(data_e.tg_energy, data_e.energy);
        peae = pearson(data_e.tg_energy, data_e.energy);
        stde = data_e.destd();
      }
      fprintf(stdout, "it: %d el_time: %d Teq: %d Twait: %d  ", iter, int(time(NULL)-in_time), params.Teq, params.Twait);
      fprintf(stdout, "lrav: %.2e sp: %.2e ", lrav, model.model_sp);
      fprintf(stdout, "cov_err: %.2e pears_act: %.2e pears_all: %.2e", errs.errnorm, model.pearson(data.cov, false), model.pearson(data.cov, true));
      //cout << " merr_fm: " << setprecision(2) << errs.merrh << " merr_sm: " << setprecision(2) << errs.merrJ << " cov_err: " << setprecision(2) << errs.errnorm;
      // cout << " averr_fm: " << errs.averrh << " averr_sm: " << errs.averrJ;
      if (params.file_msa_e)
        fprintf(stdout, "sp_e: %.1e pears_e: %.1e stde: %.1e \n", spe, peae, stde);
        //cout << " sp_e: " << spe << " pears_e: " << peae << " stde: " << stde << endl;
      else
        //cout << endl;
        fprintf(stdout, "\n");
      fflush(stdout);
    }
    if (iter > 0 && (iter % params.nprintfile == 0))
    {
      params.construct_filenames(iter, conv, par, par_zsum, ene, corr, score, first, sec, third, lchain);
      print_frobenius_norms(model.h, model.J, model.L, model.q, score, par_zsum);
      model.print_model(par);
      model.print_last_chain(lchain);
      if (params.print_samples)
        model.print_samples(ene);
      if (params.file_msa_e)
      {
        data_e.energy = model.energy(data_e.msa);
        data_e.print_energy(ene_fit);
      }
      if (data.tm.size() > 0)
        model.compute_third_order_correlations();
      data.print_statistics(sec, first, third, corr, model.mstat->corr, model.mstat->fm_s, model.mstat->sm_s, model.mstat->tm_s);
    }
    lrav = model.update_parameters(data.fm, data.sm, iter, data_e);
    if (iter > 0 && params.sparsity && (errs.errnorm < params.conv || iter % params.dec_steps == 0))
    {
      // Print converged parameters before decimation
      params.construct_filenames(iter, conv, par, par_zsum, ene, corr, score, first, sec, third, lchain);
      print_frobenius_norms(model.h, model.J, model.L, model.q, score, par_zsum);
      model.print_model(par);
      model.print_last_chain(lchain);
      if (params.print_samples)
        model.print_samples(ene);
      if (params.file_msa_e)
      {
        data_e.energy = model.energy(data_e.msa);
        data_e.print_energy(ene_fit);
      }
      if (data.tm.size() > 0)
        model.compute_third_order_correlations();
      data.print_statistics(sec, first, third, corr, model.mstat->corr, model.mstat->fm_s, model.mstat->sm_s, model.mstat->tm_s);
      // Then decimate or activate
      if (model.model_sp < params.sparsity)
      {
        int aux = fabs(model.n_active() - model.n_total() * (1. - params.sparsity)) > model.n_active() / 100. ? ceil(model.n_active() / 100.) : ceil(fabs(model.n_active() - model.n_total() * (1. - params.sparsity)));
        model.decimate_compwise(aux, iter);
      }
      else
      {
        int aux = fabs(model.n_active() - model.n_total() * (1. - params.sparsity)) > (model.n_total() - model.n_active()) / 1000. ? ceil((model.n_total() - model.n_active()) / 1000.) : ceil(fabs(model.n_active() - model.n_total() * (1. - params.sparsity)));
        model.activate_compwise(aux, iter, data.sm);
      }
    }
    if (errs.errnorm < params.conv && params.sparsity == 0)
    {
      conv = true;
      cout << "Reached convergence of error, end of learning" << endl;
    }
    if (fabs(model.model_sp - params.sparsity) <= 0.01 * params.sparsity && params.sparsity > 0 && errs.errnorm < params.conv)
    {
      conv = true;
      cout << "Reached convergence of error and desired sparsity, end of learning" << endl;
    }
    if (iter >= params.maxiter)
    {
      cout << "Reached maximum number of iterations, end of learning" << endl;
      conv = true;
    }
    iter++;
  }
  /* END ITERATION LOOP */

  /* FINAL OPERATIONS */
  cout << "****** Final sampling ******" << endl;
  fflush(stdout);
  if (!params.maxiter)
  {
    bool eqmc = false;
    while (!eqmc)
    {
      if (strcmp(params.ctype, "i") == 0)
      {
        eqmc = model.sample_ising(data.msa);
      }
      else
        eqmc = model.sample(data.msa);
      model.compute_errors(data.fm, data.sm, data.cov, errs);
      fprintf(stdout, "it: %d el_time: %d Teq: %d Twait: %d  ", iter, int(time(NULL)-in_time), params.Teq, params.Twait);
      fprintf(stdout, "lrav: %.2e sp: %.2e ", lrav, model.model_sp);
      fprintf(stdout, "cov_err: %.2e pears_act: %.2e pears_all: %.2e\n", errs.errnorm, model.pearson(data.cov, false), model.pearson(data.cov, true));
      fflush(stdout);
    }
  }
  else
  {
    if (strcmp(params.ctype, "i") == 0)
    {
      model.sample_ising(data.msa);
    }
    else
      model.sample(data.msa);
  }

  params.construct_filenames(iter, conv, par, par_zsum, ene, corr, score, first, sec, third, lchain);
  print_frobenius_norms(model.h, model.J, model.L, model.q, score, par_zsum);
  model.print_model(par);
  model.print_last_chain(lchain);

  if (params.print_samples)

    if (data.tm.size() > 0)
      model.compute_third_order_correlations();
  data.print_statistics(sec, first, third, corr, model.mstat->corr, model.mstat->fm_s, model.mstat->sm_s, model.mstat->tm_s);
  if (params.print_samples)
  {
    if (!strcmp(params.ctype, "i"))
      model.print_samples_ising(ene);
    else
    {
      model.print_samples(ene);
      sprintf(ene, "%s/Sample_natural_conv_%s.dat", params.label, params.label);
      model.print_natural_samples(ene, data.msa);
      if (params.file_msa_e)
      {
        data_e.energy = model.energy(data_e.msa);
        data_e.print_energy(ene_fit);
      }
    }
  }
  cout << "****** Execution completed ******" << endl;

  fflush(stdout);
  /* FINAL OPERATIONS */

  return 0;
}
