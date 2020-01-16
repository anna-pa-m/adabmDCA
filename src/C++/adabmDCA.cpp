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

// structures and important parameters
Params params;
Model model(0,0,0,&params);
int ** msa;
double * w;
int ** idx;
int * tmp_idx;
double * sorted_struct;
double * fm;
double ** sm;
double ** cov;
double * tm;
int ** tm_index;



int q = 21, L, M;
double Meff;
double maxsdkl;
double model_sp;

// list of functions
int alloc_structures();
int decimate(int c);
int decimate_compwise(int c, int iter);


int main(int argc, char ** argv) {

  /* START INITIALIZATION */
  params.read_params(argc,argv);
  fprintf(stdout, "Boltzmann machine for DCA model\n");
  srand(params.seed ? params.seed : time(NULL));
  q=print_alphabet(params.ctype);
  msa=read_msa(params.file_msa,params.ctype,M,L,q);
  if(!params.Twait)
    params.Twait = 10*L;
  if(!params.Teq)
    params.Teq = 20*L;
  if(params.Metropolis)
    fprintf(stdout, "Performing Metropolis-Hastings MC.\nInitial sampling time: %d\nInitial equilibration time: %d\nUsing %d seeds and tot. number of points %d\n", params.Twait, params.Teq, params.Nmc_starts, params.Nmc_starts * params.Nmc_config);
  else if(params.Gibbs) {	
    fprintf(stdout, "Performing Gibbs sampling.\nInitial sampling time: %d\nInitial equilibration time: %d\nUsing %d seeds and tot. number of points %d\n", params.Twait, params.Teq, params.Nmc_starts, params.Nmc_starts * params.Nmc_config);
  }
  fprintf(stdout, "Learning strategy: ");
  switch(params.learn_strat) {
  case 0:
    fprintf(stdout, "Using standard gradient descent with constant learning rate (for J %.3e, for h %.3e)\n", params.lrateJ, params.lrateh);
    break;
  case 1:
    fprintf(stdout, "Using adagrad\n");
    break;
  case 2:
    fprintf(stdout, "Using RMSprop with reinforcement %.2f\n", params.rho);
    break;
  case 3:
    fprintf(stdout, "Using search and converge with decay time %d and learning rate (for J %.3e, for h %.3e)\n", params.tau, params.lrateJ, params.lrateh);
    break;
  case 4:
    fprintf(stdout, "Using adam : NOT YET IMPLEMENTED, EXIT\n");
    exit(EXIT_FAILURE);
    break;
  }
  if(params.sparsity > 0.0) {
    params.regJ = 0;
    fprintf(stdout, "Sparsity %.3f, using direct decimation on couplings\n", params.sparsity);
  }
  if(params.regJ > 0)
    fprintf(stdout, "L1 regularization on couplings: lambda %.1e\n", params.regJ);
  if(params.pseudocount)
    fprintf(stdout, "Using pseudo-count: %1.e\n", params.pseudocount);
  w=compute_w(params.file_w,params.label,params.w_th,msa,M,L);
  alloc_structures();
  Meff=compute_empirical_statistics(fm,sm,cov,params.pseudocount,msa,w,M,L,q);
  model.resize(q,L,M);
  model.initialize_parameters(fm);
  int n = model.initialize_model(cov);
  model_sp = 0.0;
  load_third_order_indices();
  /* END INITIALIZATION */

  /* START ITERATION LOOP */
  char score[1000];
  char sec[1000];
  char first[1000];
  char third[1000];
  char par[1000];
  char sc;

  fprintf(stdout, "Printing every %d iteration(s)\n", params.nprint);
  int iter = 1;
  long in_time = time(NULL);
  bool conv = false;
  Errs errs;
  while(!conv && iter < params.maxiter) {
    bool print_aux = false;
    model.init_statistics();
    model.sample();
    double lrav=model.update_parameters(fm,sm,iter);
    model.compute_errors(fm,sm,cov,errs);
    if(model_sp < params.sparsity && iter % 10 == 0) {
      fprintf(stdout, "Decimating..");
      decimate(ceil((n*10)/(params.maxiter*0.5)));
    } else if(params.compwise && errs.errnorm<params.conv) {
      fprintf(stdout,"Decimating..");
      int aux = ceil( (1.0 - model_sp) * n / 100);
      decimate_compwise(aux,iter);
      print_aux = true;
    }
    if(errs.errnorm<params.conv && params.sparsity == 0 && !params.compwise)
      conv = true;
    if(model_sp >= params.sparsity && params.sparsity > 0 && errs.errnorm<params.conv)
      conv = true;
    if(iter % params.nprint == 0) {
      fprintf(stdout, "it: %i el_time: %li N: %i Teq: %i Twait: %i merr_fm: %.1e merr_sm: %.1e averr_fm: %.1e averr_sm: %.1e cov_err: %.1e corr: %.2f sp: %.1e lrav: %.1e\n", iter, time(NULL)-in_time, params.Nmc_config * params.Nmc_starts, params.Teq, params.Twait, errs.merrh, errs.merrJ, errs.averrh, errs.averrJ, errs.errnorm, model.pearson(cov), model_sp, lrav);
      fflush(stdout);
    }
    if(iter % params.nprintfile == 0 || print_aux) {
      sc = (params.Gibbs == 0) ? 'M' : 'G';
      sprintf(par, "Parameters_tmp_%d_zerosum_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", iter, params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
      sprintf(score, "Score_tmp_%d_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", iter, params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
      print_frobenius_norms(model.h,model.J,L,q,score, par);
      sprintf(par, "Parameters_tmp_%d_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", iter, params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
      model.print_model(par);
      sprintf(sec, "Sec_mom_tmp_%d_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", iter, params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
      sprintf(first, "First_mom_tmp_%d_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", iter, params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
      sprintf(third, "Third_mom_tmp_%d_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", iter, params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
      if(compute_tm)
	compute_third_order_correlations();
      print_statistics(sec, first, third);
    }
    iter++;
  }
  /* END ITERATION LOOP */

  /* FINAL OPERATIONS */
  model.init_statistics();
  model.sample(); // compute 3rd order moments through sampling and possibly print the sequences
  if(compute_tm)
    compute_third_order_correlations();
  sc = (params.Gibbs == 0) ? 'M' : 'G';
  sprintf(par, "Parameters_conv_zerosum_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
  sprintf(score, "Score_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat",params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
  print_frobenius_norms(model.h,model.J,L,q,score, par);
  sprintf(sec, "Sec_mom_conv_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
  sprintf(first, "First_mom_conv_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
  sprintf(third, "Third_order_connected_corr_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
  print_statistics(sec, first, third);
  sprintf(par, "Parameters_conv_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
  model.print_model(par);
  return 0;

}







///////FZ: BLOCK OF ROUTINES THAT SHOULD BE MOVED SOMEWHERE //////////////////////////////////////////////////////////////////////

int alloc_structures() {
  // one-point frequencies and fields
  fm = (double *)calloc(L*q, sizeof(double));
  // two-point frequencies and fields
  sm = (double **)calloc(L*q, sizeof(double *));
  cov = (double **)calloc(L*q, sizeof(double *));
  for(int i = 0; i < L*q; i++) {
    sm[i] = (double *)calloc(L*q, sizeof(double));
    cov[i] = (double *)calloc(L*q, sizeof(double));
  }
  // variables for sorting
  if(params.sparsity > 0 || params.compwise || params.blockwise) {
    int n = L*(L-1)*q*q/2;
    idx = (int **)calloc(n, sizeof(int *));
    sorted_struct = (double *)calloc(n, sizeof(double));
    tmp_idx = (int *)calloc(n, sizeof(int));
    for(int i = 0; i < n;i++)
      idx[i] = (int *)calloc(4, sizeof(int));
    int k = 0;
    for(int i = 0; i < L; i++) {
      for(int j = i+1; j < L; j++) {
	for(int a = 0; a < q; a++) {
	  for(int b = 0; b < q; b++) {
	    idx[k][0] = i;
	    idx[k][1] = j;
	    idx[k][2] = a;
	    idx[k][3] = b;
	    sorted_struct[k] = 0.0;
	    k += 1;
	  }
	}
      }
    }
  }
  return 0;
}

int load_third_order_indices() {
  if(params.file_3points) {
    fprintf(stdout, "Reading three points correlations indices...");
    FILE *file3;
    if(!(file3 = fopen(params.file_3points, "r"))) {
      fprintf(stderr, "File %s not found \n", params.file_3points);
      return EXIT_FAILURE;
    } else {
      int i, j, k, a, b, c, m;
      int ind;
      double value;
      int line = 0;
      char buffer[1000];
      char ch;
      ntm = 0;
      compute_tm = true;
      while(!feof(file3)) {
	ch = fgetc(file3);
	if(ch == '\n')
	  ntm++;
      }
      rewind(file3);
      fprintf(stdout, "Number of indices %d\n", ntm);
      tm_index = (int **)calloc(ntm, sizeof(int *));
      for(ind = 0; ind < ntm; ind++)
	tm_index[ind] = (int *)calloc(6, sizeof(int));
      tm = (double *)calloc(ntm, sizeof(double)); // Pay attention. File contains connected 3rd order correlations
      tm_s = (double *)calloc(ntm, sizeof(double));
      while (!feof(file3) && fgets(buffer, 1000, file3)) {
	if(sscanf(buffer, "%d %d %d %d %d %d %lf \n", &i, &j, &k, &a, &b, &c, &value) == 7) {
	  tm_index[line][0] = i;
	  tm_index[line][1] = j;
	  tm_index[line][2] = k;
	  tm_index[line][3] = a;
	  tm_index[line][4] = b;
	  tm_index[line][5] = c;
	  tm[line] = 0.0;	
	  line++;
	}
      }
      // compute 3rd order moments of msa, slow!
      for(ind =0; ind < ntm; ind++) {
	i = tm_index[ind][0];
	j = tm_index[ind][1];
	k = tm_index[ind][2];
	a = tm_index[ind][3];
	b = tm_index[ind][4];
	c = tm_index[ind][5];
	for(m = 0; m < M; m++) 
	  if(msa[m][i] == a && msa[m][j] == b && msa[m][k] == c)
	    tm[ind] += w[m]/Meff;
      }
    }
    return 0;
  } else {
    fprintf(stdout, "No three-points correlations indices specified\n");
    return 0;
  }  
}

int update_tm_statistics(vector<int> & x)
{
	int i, j, k;
	int Ns = params.Nmc_starts * params.Nmc_config;
	int a, b, c;
	int ind;
	for(ind = 0; ind < ntm; ind ++) {
		i = tm_index[ind][0];
		j = tm_index[ind][1];
		k = tm_index[ind][2];
		a = tm_index[ind][3];
		b = tm_index[ind][4];
		c = tm_index[ind][5];
		if(x[i] == a && x[j] == b && x[k] == c)
			tm_s[ind] += 1.0/Ns;
	}
	return 0;
}

int compute_third_order_correlations() {
  int ind, i, j, k, a, b, c;
  for(ind = 0; ind < ntm; ind++) {
    i = tm_index[ind][0];
    j = tm_index[ind][1];
    k = tm_index[ind][2];
    a = tm_index[ind][3];
    b = tm_index[ind][4];
    c = tm_index[ind][5];
    //tm[ind] = tm[ind] - sm[i*q+a][j*q+b]*fm[k*q+c] - sm[i*q+a][k*q+c]*fm[j*q+b] - sm[j*q+b][k*q+c]*fm[i*q+a] + 2*fm[i*q+a]*fm[j*q+b]*fm[k*q+c];
    tm_s[ind] = tm_s[ind] - sm_s[i*q+a][j*q+b]*fm_s[k*q+c] - sm_s[i*q+a][k*q+c]*fm_s[j*q+b] - sm_s[j*q+b][k*q+c]*fm_s[i*q+a] + 2*fm_s[i*q+a]*fm_s[j*q+b]*fm_s[k*q+c]; 
  }
  return 0;
}

int print_statistics(char *file_sm, char *file_fm, char *file_tm) {
  int i, j , a, b;
  FILE *fs, *ff, *ft;
  fs = fopen(file_sm, "w");
  ff = fopen(file_fm, "w");
  for(i = 0; i < L; i++) {
    for(j = i+1; j < L; j++) {
      for(a = 0; a < q; a++) {
	for(b = 0; b < q; b++)
	  fprintf(fs, "%d %d %d %d %.5f %.5f\n",i, j,a,b, sm[i*q+a][j*q+b], sm_s[i*q+a][j*q+b]);
      }
    }
  }
  for(i = 0; i < L; i++) {
    for(a = 0; a < q; a++)
      fprintf(ff, "%i %i %.5f %.5f\n", i, a, fm[i*q+a], fm_s[i*q+a]);
  }
  if(params.file_3points) {
    ft = fopen(file_tm, "w");
    int k, c, ind;
    double aux;
    for(ind = 0; ind < ntm; ind++) {
      i = tm_index[ind][0];
      j = tm_index[ind][1];
      k = tm_index[ind][2];
      a = tm_index[ind][3];
      b = tm_index[ind][4];
      c = tm_index[ind][5]; 
      aux = tm[ind] - sm[i*q+a][j*q+b]*fm[k*q+c] - sm[i*q+a][k*q+c]*fm[j*q+b] - sm[j*q+b][k*q+c]*fm[i*q+a] + 2*fm[i*q+a]*fm[j*q+b]*fm[k*q+c];
      fprintf(ft, "%d %d %d %d %d %d %.5f %.5f\n", i, j, k, a, b,c, aux, tm_s[ind]);
    }
    fflush(ft);
    fclose(ft);
  }
  fflush(fs);
  fflush(ff);
  fclose(fs);
  fclose(ff);
  return 0;
}

///////FZ: CHECKED UP TO THIS POINT //////////////////////////////////////////////////////////////////////












int decimate_compwise(int c, int iter)
{
	int k, i, j, a, b;
	int index;
	int n, m = 0;
	int neff;
	FILE *fileout;
	double auxsm;
	char filename_aux[1000];
	sprintf(filename_aux, "sDKL_couplings_iter_%i.dat", iter);
	fileout = fopen(filename_aux, "w");
	maxsdkl = -1e50;
	printf("n terms: %d\n",c);
	if(params.dgap)
		n = (L*(L-1)*(q-1)*(q-1))/2;
	else if(params.gapnn)
		n = (L*(L-1)*(q-1)*(q-1))/2 + L - 1;
	else if(params.phmm)
		n = L - 1;
	else {
		neff = (L*(L-1)*q*q)/2 - (L*(L-1)*(2*q-1))/2 + L;
		n = (L*(L-1)*q*q)/2;
	}

	for(k = 0; k < n; k++) {
		tmp_idx[k] = k;
		sorted_struct[k] = 0.0;
		i = idx[k][0];
		a = idx[k][2];
		j = idx[k][1];
		b = idx[k][3];
		if(model.decJ[i*q + a][j*q + b] > 0) {
			m += 1;
			auxsm = sm_s[i*q+a][j*q+b] == 0 ? params.pseudocount : sm_s[i*q+a][j*q+b];
			sorted_struct[k] = model.J[i*q+a][j*q+b] * auxsm
					- (model.J[i*q+a][j*q+b] * exp(-model.J[i*q+a][j*q+b]) * auxsm) / (exp(-model.J[i*q+a][j*q+b]) *auxsm + 1 - auxsm);
			maxsdkl = max(maxsdkl, sorted_struct[k]);
			fprintf(fileout, "%i %i %i %i %.2e %f %f %f\n", i, j, a, b, sorted_struct[k], model.J[i*q+a][j*q+b], model.decJ[i*q+a][j*q+b], sm_s[i*q+a][j*q+b]);
			//sorted_struct[k] += 1e-4 * rand01();
		} else {
			double f = rand01();
			sorted_struct[k] =  n * f; // to be optimized: elements should be removed instead of putting large numbers
			//fprintf(stderr, "n %d fake %f rnd %f\n", n, sorted_struct[k], f);
		}
	}
	fprintf(stdout, "Non-zeros parameters %d / %d \n",m, neff);
	quicksort(sorted_struct, tmp_idx, 0, n-1);
	for(k = 0; k < c; k++) {
		index = tmp_idx[k];
		i = idx[index][0];
		a = idx[index][2];
		j = idx[index][1];
		b = idx[index][3];
		//fprintf(stderr, "%i %f\n", index, J[idx[index][0]*q+idx[index][2]][idx[index][1]*q + idx[index][3]]);
		model.J[i*q + a][j*q + b] = 0.0;
		model.J[j*q + b][i*q + a] = 0.0;
	       	model.decJ[i*q+a][j*q + b] = 0.0;
		model.decJ[j*q+b][i*q + a] = 0.0;
	}
	for(k = 0; k<=c; k++) {
		index = tmp_idx[k];
		i = idx[index][0];
		a = idx[index][2];
		j = idx[index][1];
		b = idx[index][3];
		fprintf(fileout, "%i %i %i %i %.2e %f * \n", i, j, a, b, sorted_struct[k], model.J[i*q+a][j*q+b]);
	}
	index = tmp_idx[c];
	fflush(fileout);
	fclose(fileout);
	fprintf(stdout, "Smallest sDKL associated with the first kept coupling is %.2e (i: %i j: %i a: %i b: %i)\n", sorted_struct[c], idx[index][0], idx[index][1], idx[index][2], idx[index][3]);
	model_sp += 1.0*c/n;
	return 0;
}

int decimate(int c)
{
	int k, i, j, a, b;
	int index;
	int n, m = 0;
	printf("c: %d\n",c);
	if(params.dgap)
		n = (L*(L-1)*(q-1)*(q-1))/2;
	else if(params.gapnn)
		n = (L*(L-1)*(q-1)*(q-1))/2 + L - 1;
	else if(params.phmm)
		n = 2*L - 1;
	else
		n = (L*(L-1)*q*q)/2;
	for(k = 0; k < n; k++) {
		i = idx[k][0];
		a = idx[k][2];
		j = idx[k][1];
		b = idx[k][3];
		if(model.decJ[i*q + a][j*q + b] > 0) {
			m += 1;
			sorted_struct[k] = fabs(model.J[i*q+a][j*q+b]);
		} else
			sorted_struct[k] = n * rand01(); // to be optimized: elements should be removed instead of putting large numbers
	}
	fprintf(stdout, "Non-zeros parameters %d / %d \n",m,n);
	for(k = 0; k < n; k++)
		tmp_idx[k] = k;
	quicksort(sorted_struct, tmp_idx, 0, n-1);
	for(k = 0; k < c; k++) {
		index = tmp_idx[k];
		i = idx[index][0];
		a = idx[index][2];
		j = idx[index][1];
		b = idx[index][3];
		//fprintf(stderr, "%i %f\n", index, J[idx[index][0]*q+idx[index][2]][idx[index][1]*q + idx[index][3]]);
		model.J[i*q + a][j*q + b] = 0.0;
		model.J[j*q + b][i*q + a] = 0.0;
	       	model.decJ[i*q+a][j*q + b] = 0.0;
		model.decJ[j*q+b][i*q + a] = 0.0;
	}
	model_sp += 1.0*c/n;
	return 0;
}







