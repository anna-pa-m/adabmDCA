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
#include <vector>
#include <valarray>
#include "BMaux.h"
#include "BMlib.h"

using namespace std;

#ifndef BMmc
#define BMmc

/////////////////////////////////////INCAPSULATE ALL THIS///////////////////////
vector<double> fm_s;
vector< vector<double> > sm_s;
// BEGIN - THIRD MOMENT
bool compute_tm = false;
int ntm = 0;
double * tm_s;
int load_third_order_indices();
int update_tm_statistics(vector<int> & x);
int compute_third_order_correlations();
int print_statistics(char * file_sm, char *file_fm, char *file_tm);
// END - THIRD MOMENT
/////////////////////////////////////INCAPSULATE ALL THIS///////////////////////


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
  vector<double> h, Gh;
  vector< vector<double> > J, decJ, GJ;
  bool Gibbs;
  Params * params;
  double alpha, acc;  // for FIRE
  int counter;   // for FIRE

 Model(int _q, int _L, Params * _params):
  q(_q),L(_L),h(L*q,0),J(L*q,h),decJ(L*q,h),Gibbs(false),params(_params),alpha(0.1),acc(1),counter(0) {
    if (params->learn_strat == 1 || params->learn_strat == 2 || params->learn_strat == 5) {
      Gh.clear();
      Gh.resize(L*q,0);
      GJ.clear();
      GJ.resize(L*q,Gh);
    }
  }

  /******************** METHODS FOR INIT AND OUTPUT ***********************************************************/

  void resize(int _q, int _L) {
    q=_q;
    L=_L;
    h.clear();
    h.resize(L*q,0);
    J.clear();
    J.resize(L*q,h);
    decJ.clear();
    decJ.resize(L*q,h);
    if (params->learn_strat == 1 || params->learn_strat == 2 || params->learn_strat == 5) {
      Gh.clear();
      Gh.resize(L*q,0);
      GJ.clear();
      GJ.resize(L*q,Gh);
    }
    fm_s.clear();
    fm_s.resize(L*q,0);
    sm_s.clear();
    sm_s.resize(L*q,fm_s);
  }

  int remove_gauge_freedom(double pseudocount, vector< vector<double> > & cov) {
    double sorted_matrix[q*q];
    int mapping[q*q];
    int idx_aux[q*q][2];
    int neff = 0;
    for(int i = 0; i < L; i++) {
      for(int j = i+1; j < L;j++) {
	int k;
	for(k = 0; k < q*q; k++)
	  sorted_matrix[k] = 0;
	k = 0;
	for(int a = 0; a < q; a++) {
	  for(int b = 0; b < q; b++){
	    mapping[k] = k;
	    idx_aux[k][0] = a;
	    idx_aux[k][1] = b;
	    sorted_matrix[k] = fabs(cov[i*q+a][j*q+b]) + pseudocount * 0.001 * rand01();
	    k += 1;
	  }
	}
	quicksort(sorted_matrix, mapping, 0, q*q - 1);
	for(k = 0; k < 2*q-1; k++) {
	  static int a = idx_aux[mapping[k]][0];
	  static int b = idx_aux[mapping[k]][1];
	  J[i*q+a][j*q+b] = 0.0;
	  decJ[i*q+a][j*q+b] = 0.0;
	  J[j*q+b][i*q+a] = 0.0;
	  decJ[j*q+b][i*q+a] = 0.0;
	}
	neff += 2*q -1;
      }
    }
    fprintf(stdout, "%d couplings have been removed\n", neff);
    return 0;
  }
  
  
  int initialize_parameters(vector<double> & fm) {
    if (params->Metropolis && !params->Gibbs) {
      Gibbs=false;
    } else if (!params->Metropolis && params->Gibbs) {
      Gibbs=true;
    } else {
	fprintf(stderr, "Conflict in Gibbs-Metropolis initialization\n");
	exit(EXIT_FAILURE);
    }
    /* {READ J,H FROM FILE} OR {SET J,H=0 OR H=IND.SITE.MODEL} */
    if(params->file_params) {
      fprintf(stdout, "Reading input parameters from file...");
      FILE *filep;
      if(!(filep = fopen(params->file_params, "r"))) {
	fprintf(stderr, "File %s not found\n", params->file_params);
	exit(EXIT_FAILURE);
      } else {
	int  i,j,a, b;
	char buffer[100];
	char c;
	double tmp;
	while(!feof(filep) && fgets(buffer,100, filep) && sscanf(buffer, "%c ", &c) == 1) {
	  switch (c) {
	  case 'J':
	    sscanf(buffer, "J %d %d %d %d %lf \n", &i, &j, &a, &b, &tmp);
	    J[i*q + a][j*q + b] = tmp;
	    J[j*q + b][i*q + a] = tmp;
	    break;
	  case 'j':
	    sscanf(buffer, "j %d %d %d %d %lf \n", &i, &j, &a, &b, &tmp);
	    J[i*q + a][j*q + b] = tmp;
	    J[j*q + b][i*q + a] = tmp;
	    break;
	  case 'H':
	    sscanf(buffer, "H %d %d %lf \n", &i, &a, &tmp);
	    h[i*q + a] = tmp;
	    break;
	  case 'h':
	    sscanf(buffer, "h %d %d %lf \n", &i, &a, &tmp);
	    h[i*q + a] = tmp;
	    break;
	  }
	}
	fclose(filep);
      }
      fprintf(stdout, "done\n");
    } else {	  
      for(int i = 0; i < L*q; i++) {
	h[i] = 0.0;
	for(int j = i; j < L*q; j++) {
	  J[i][j] = 0.0;
	  J[j][i] = 0.0;
	}
      }
      switch(params->init) {
      case 'R':
	fprintf(stdout, "Zero-parameters initialization...done\n");
	fflush(stdout);
	break;
      case 'I':
	fprintf(stdout, "Initializing parameters using independent sites approximation...");
	double mean_all;
	for(int i = 0; i < L; i++) {
	  mean_all = 0;
	  for(int a = 0; a < q; a++) {
	    h[i*q +a] = log(fm[i*q+a]);
	    mean_all += h[i*q+a];
	  }
	  mean_all /= q;
	  for(int a = 0; a < q; a++)
	    h[i*q+a] -= mean_all;
	}
	fprintf(stdout, "done\n");
	fflush(stdout);
      }
    }
    return 0;
  }
  
  int initialize_model(vector< vector<double> > & cov) {
    int n=0;
    if(!params->file_cc) {
      if(params->dgap) {
	fprintf(stdout, "Using DGap model\n");
	n = (L*(L-1)*(q-1)*(q-1))/2;
      } else if(params->gapnn) {
	fprintf(stdout, "Using GapNN model\n");
	n = (L*(L-1)*(q-1)*(q-1))/2 + L - 1;
      } else if(params->phmm) {
	fprintf(stdout, "Using Hmmer-like model\n");
	n = L - 1;
      } else if(params->rmgauge) {
	fprintf(stdout, "Using Potts model with Gauge fixing via cc-decimation\n");
	n = (L*(L-1)*(q-1)*(q-1))/2;
      } else {
	fprintf(stdout, "Using full Potts model\n");
	n = (L*(L-1)*q*q)/2;
      }
      for(int i = 0; i < L; i++) {
	for(int j = 0; j < L; j++) {
	  for(int a = 0; a < q; a++) {
	    for(int b = 0; b < q; b++) {
	      decJ[i*q + a][j*q + b] = 1.0;	
	      if(i == j) {
		decJ[i*q+a][j*q+b] = 0.0;
		J[i*q+a][j*q+b] = 0.0;
	      } else if(params->dgap && (a == 0 || b == 0)) {
		decJ[i*q + a][j*q + b] = 0.0;
		J[i*q + a][j*q + b] = 0.0;
	      } else if(params->gapnn) {
		if(abs(i - j) > 1 && (a == 0 || b == 0)) {
		  decJ[i*q + a][j*q + b] = 0.0;
		  J[i*q + a][j*q + b] = 0.0;
		}
		if(abs(i-j)==1 && ((a == 0 && b > 0) || (b == 0 && a > 0))) {
		  decJ[i*q + a][j*q + b] = 0.0;
		  J[i*q + a][j*q + b] = 0.0;
		} 
	      } else if(params->phmm && !(a == 0 && b == 0 && abs(i-j) == 1)) {
		decJ[i*q + a][j*q + b] = 0.0;
		J[i*q + a][j*q + b] = 0.0;
	      } 
	    }
	  }
	}
      }
      if(params->rmgauge) remove_gauge_freedom(params->pseudocount,cov);
    } else {
      fprintf(stdout, "Reading interaction graph from input file...");
      FILE *filep;
      if(!(filep = fopen(params->file_cc, "r"))) {
	fprintf(stderr, "File %s not found\n", params->file_params);
	return EXIT_FAILURE;
      } else {
	int  i,j, a, b;
	char buffer[100];
	double tmp;
	for(i = 0; i < L*q; i++) {
	  for(j = 0; j < L*q; j++)
	    decJ[i][j] = 0.0;
	}
	while(!feof(filep) && fgets(buffer,100, filep) && 
	      ( sscanf(buffer, "%d %d %d %d %lf \n", &i, &j, &a, &b, &tmp) == 5 || sscanf(buffer, "%d %d %d %d \n", &i, &j, &a, &b) == 4) ) {
	  decJ[i*q+a][j*q+b] = 1.0;
	  decJ[j*q+b][i*q+a] = 1.0;
	  n++;
	}
	fclose(filep);
	for(i = 0; i < q*L; i++) {
	  for(j = 0; j < q*L; j++) {
	    if(decJ[i][j] == 0.0) {
	      J[i][j] = 0.0;
	    }
	  }
	}
      }
      fprintf(stdout, "done\nNumber of links %d\n", n);
    }
    return n;	
  }

  int print_model(char *filename) {
    FILE *fp;
    fp = fopen(filename, "w");
    for(int i = 0; i < L; i++) {
      for(int j = 0; j < L; j++) if(i < j) {
	  for(int a = 0; a < q; a++) {
	    for(int b = 0; b < q; b++)
	      fprintf(fp, "J %d %d %d %d %f\n",i, j,a,b, J[i*q+a][j*q+b]);
	  }
	}
    }
    for(int i = 0; i < L; i++) {
      for(int a = 0; a < q; a++)
	fprintf(fp, "h %i %i %f\n", i, a, h[i*q+a]);
    }
    fflush(fp);
    fclose(fp);
    return 0;
  }


  /******************** METHODS FOR MONTE CARLO ***********************************************************/

  double energy(vector<int> & seq) {
    double en = 0;  
    for(int i = 0; i < L; i++) {
      en += -h[i*q + seq[i]];
      for(int j = i+1; j < L; j++) 
	en += -J[i*q + seq[i]][j*q +seq[j]];
    }
    return en;
  }

  int metropolis_step(vector<int> & curr_state) {
    int i = (int)rand() % L;
    int a = (int)rand() % q;
    while(a == curr_state[i])
      a = (int)rand() %q;
    double deltaE = -h[i*q + a] + h[i*q + curr_state[i]];
    for(int j = 0; j < L; j++) if(j != i) {
	deltaE += - J[i*q + a][j*q + curr_state[j]] + J[i*q + curr_state[i]][j*q + curr_state[j]];
      }
    double p = rand01();
    if (exp(-deltaE) > p) {
      curr_state[i] = a;
    }
    return 0;
  }

  int gibbs_step(vector<int> & curr_state) {
    double H, cum[q];
    int i = (int)rand() % L;
    for(int a = 0; a < q; a++) {
      H = -h[i*q + a];
      for(int j = 0; j < L; j++) if(j != i) {
	  H += -J[i*q + a][j*q + curr_state[j]];
	}
      if(a==0) cum[a]=exp(-H);
      else cum[a] = cum[a-1] + exp(-H);
    }
    double r = cum[q-1]*rand01();
    int a=0;
    while(r>cum[a]) {a++;}
    curr_state[i] = a;
    return 0;
  }

  int MC_step(vector<int> & curr_state) {
    if (!Gibbs) {return metropolis_step(curr_state);} else {return gibbs_step(curr_state);}
  }

  valarray<int> mc_chain(vector<int> & curr_state1, vector<int> & curr_state2) {
    FILE * fp = 0, * fe = 0;
    if(params->file_samples)
      fp  = fopen(params->file_samples, "a");
    if(params->file_en)
      fe = fopen(params->file_en, "a");
    for(int t=0; t < params->Teq * L; t++) {
      MC_step(curr_state1);
      MC_step(curr_state2);
    }
    valarray<int> qs(6);
    vector<int> old_state1, old_state2, oldold_state1, oldold_state2;
    update_statistics(curr_state1);        
    update_statistics(curr_state2);
    if(compute_tm) {
      update_tm_statistics(curr_state1);
      update_tm_statistics(curr_state2);
    }
    double o12=overlap(curr_state1,curr_state2);
    qs[0]+=o12;
    qs[1]+=(o12*o12);
    for(int n = 0; n < params->Nmc_config - 1; n++) { 
      oldold_state1=old_state1;
      oldold_state2=old_state2;
      old_state1=curr_state1;
      old_state2=curr_state2;
      for(int t=0; t < params->Twait * L; t++) {
	MC_step(curr_state1);
	MC_step(curr_state2);
      }
      update_statistics(curr_state1);
      update_statistics(curr_state2);
      if(compute_tm) {
	update_tm_statistics(curr_state1);
	update_tm_statistics(curr_state2);
      }
      o12=overlap(curr_state1,curr_state2);
      qs[0]+=o12;
      qs[1]+=(o12*o12);
      if (n>0) {
	double o1=overlap(old_state1,curr_state1);
	double o2=overlap(old_state2,curr_state2);
	qs[2]+=(o1+o2);
	qs[3]+=(o1*o1+o2*o2);
      }
      if (n>1) {
	double oo1=overlap(oldold_state1,curr_state1);
	double oo2=overlap(oldold_state2,curr_state2);
	qs[4]+=(oo1+oo2);
	qs[5]+=(oo1*oo1+oo2*oo2);
      }
      if(params->file_samples) {
	for(int i = 0; i < L; i++)
	  fprintf(fp, "%d ", curr_state1[i]);
	fprintf(fp, "\n");
	for(int i = 0; i < L; i++)
	  fprintf(fp, "%d ", curr_state2[i]);
	fprintf(fp, "\n");
	fflush(fp);
      }
      if(params->file_en) {
	fprintf(fe, "%lf\n", energy(curr_state1));
	fprintf(fe, "%lf\n", energy(curr_state2));
      }
    }
    if(params->file_samples)
      fclose(fp);
    if(params->file_en)
      fclose(fe);
    return qs; 
  }

  int sample() {
    init_statistics();
    vector<int> curr_state1(L);
    vector<int> curr_state2(L);
    valarray<int> qs(6);
    cout<<"Start sampling..."<<flush;
    for(int s = 0; s < params->Nmc_starts/2; s++) {
      for(int i = 0; i < L; i++) {
	curr_state1[i] = (int)rand() % q;
	curr_state2[i] = (int)rand() % q;
      }
      qs+=mc_chain(curr_state1,curr_state2);
    } 
    double nse=params->Nmc_config*(params->Nmc_starts/2);
    double qext=qs[0]/nse;
    double dqext=sqrt(qs[1]/(nse-1)-qs[0]*qs[0]/nse/(nse-1))/sqrt(nse);
    double nsi1=(params->Nmc_config-1)*params->Nmc_starts;
    double qin1=nsi1>0 ? qs[2]/nsi1 : 0;
    double dqin1=nsi1>1 ? sqrt(qs[3]/(nsi1-1)-qs[2]*qs[2]/nsi1/(nsi1-1))/sqrt(nsi1) : 0;
    double nsi2=(params->Nmc_config-2)*params->Nmc_starts;
    double qin2=nsi2>0 ? qs[4]/nsi2 : 0;
    double dqin2=nsi2>1 ? sqrt(qs[5]/(nsi2-1)-qs[4]*qs[4]/nsi2/(nsi2-1))/sqrt(nsi2) : 0;
    int test1=(abs(qext-qin1)<3*sqrt(dqext*dqext+dqin1*dqin1) ? 1 : 0);
    int test2=(abs(qext-qin2)<3*sqrt(dqext*dqext+dqin2*dqin2) ? 1 : 0);
    cout<<"done, q_ext: "<<qext<<" +- "<<dqext<<" q_int_1: "<<qin1<<" +- "<<dqin1<<" q_int_2: "<<qin2<<" +- "<<dqin2<<" Test_eq1: "<<test1<<" Test_eq2: "<<test2<<endl;
    if (test1) {
      if (params->Twait > 1) params->Twait-=1;
      params->Teq = 2*params->Twait;
      cout<<"Reduced equilibration time, Teq="<<params->Teq<<" Twait="<<params->Twait<<endl;
    } else if (!test2) {
      params->Twait+=1;
      params->Teq = 2*params->Twait;
      cout<<"Increased equilibration time, Teq="<<params->Teq<<" Twait="<<params->Twait<<endl;
    }
    return 0;
  }



  /******************** METHODS FOR STATISTICS ***********************************************************/

  int init_statistics() {
    for(int i = 0; i < L*q; i++) {
      fm_s[i] = 0;
      for(int j = 0; j < L*q; j++) {
	sm_s[i][j] = 0;
      }
    }
    for(int ind = 0; ind < ntm; ind++) 
      tm_s[ind] = 0;  
    return 0;
  }
  
  int update_statistics(vector<int> & x) {
    int Ns = params->Nmc_starts * params->Nmc_config;
    for(int i = 0; i < L; i++) {
      fm_s[i*q + x[i]] += 1.0/Ns;
      sm_s[i*q + x[i]][i*q + x[i]] += 1.0/Ns;
      for(int j = i+1; j < L; j++) {
	sm_s[i*q + x[i]][j*q + x[j]] += 1.0/Ns;
	sm_s[j*q + x[j]][i*q + x[i]] += 1.0/Ns;
      }
    }
    return 0;
  }

  int compute_errors(vector<double> & fm, vector< vector<double> > & sm, vector< vector<double> > & cov, Errs & errs) {
    errs.errnorm = 0;
    errs.merrh = 0;
    errs.merrJ = 0;
    errs.averrh = 0;
    errs.averrJ = 0;
    for(int i = 0; i < L; i++) {
      for(int a = 0; a < q; a++) {
	errs.averrh += fabs(fm_s[i*q+a] - fm[i*q+a]);
	errs.merrh =  max(errs.merrh, fabs(fm_s[i*q+a] - fm[i*q+a]));	
	for(int j = i+1; j < L; j++) {
	  for(int b = 0; b < q; b++) {
	    errs.errnorm = max(errs.errnorm, decJ[i*q+a][j*q+b] * fabs(cov[i*q+a][j*q+b] - sm_s[i*q+a][j*q+b] + fm_s[i*q+a]*fm_s[j*q+b]));
	    errs.averrJ += fabs(sm_s[i*q+a][j*q+b] - sm[i*q+a][j*q+b]);
	    errs.merrJ =  max(errs.merrJ, fabs(sm_s[i*q+a][j*q+b] - sm[i*q+a][j*q+b]));
	  }
	}
      }
    }
    errs.averrh /= L*q;
    errs.averrJ /= (L*(L-1)/2)*q*q;
    return 0;
  }

  double pearson(vector< vector<double> > & cov) {
    double mean_cov_s = 0.0;
    double mean_cov = 0.0;
    double mean_prod = 0.0;
    double mean_x2 = 0.0;
    double mean_y2 = 0.0;
    int n = 0;
    for(int i = 0; i < L; i++) {
      for(int j = i + 1; j < L; j++) {
	for(int a = 0; a < q; a++) {
	  for(int b = 0; b < q; b++) {
	    if(decJ[i*q + a][j*q +b]) {
	      n += 1;
	      double cov_s = sm_s[i*q+a][j*q+b] - fm_s[i*q+a]*fm_s[j*q+b];
	      mean_cov_s += cov_s;
	      mean_cov += cov[i*q+a][j*q+b];
	      mean_prod += cov[i*q+a][j*q+b] * cov_s;
	      mean_x2 += cov[i*q+a][j*q+b] * cov[i*q+a][j*q+b];
	      mean_y2 += cov_s * cov_s;
	    }
	  }
	}
      }
    }
    mean_cov_s /= n;
    mean_cov /= n;
    mean_prod /= n;
    mean_y2 /= n;
    mean_x2 /=n;
    double covxy = mean_prod - mean_cov_s * mean_cov;
    double std_cov_s = sqrt(mean_y2 - mean_cov_s * mean_cov_s);
    double std_cov = sqrt(mean_x2 - mean_cov * mean_cov);
    double rho = covxy / (std_cov_s * std_cov);
    return rho;
  }
  

  /******************** METHODS FOR STATISTICS ***********************************************************/

  double update_parameters(vector<double> & fm, vector< vector<double> > & sm, int iter) {
    if (params->learn_strat != 5) {
      double lrav=0;
      int n=0;
      for(int i = 0; i < L; i++) {
	for(int a = 0; a < q; a++) {
	  double gradh = fm[i*q+a] - fm_s[i*q+a];	
	  double lrh=params->lrateh;
	  if (params->learn_strat == 1) {
	    Gh[i*q+a]+=gradh*gradh;
	    lrh /= sqrt(Gh[i*q+a] + 1e-12);
	  } else if (params->learn_strat == 2) {
	    Gh[i*q+a] = params->rho * Gh[i*q+a] + (1.0 - params->rho)*gradh*gradh;
	    lrh /= sqrt(Gh[i*q+a] + 1e-12);
	  } else if(params->learn_strat == 3) {
	    lrh /= (1.0 + (double)iter/params->tau);
	  }
	  h[i*q+a] += lrh * gradh;
	  lrav+=lrh;
	  n++;
	  for(int j = i+1; j < L; j++) {
	    for(int b = 0; b <q; b++) {
	      double gradJ=sm[i*q+a][j*q+b] - sm_s[i*q+a][j*q+b] - params->regJ * ( (J[i*q +a][j*q + b] > 0)  - (J[i*q +a][j*q+b] < 0) );
	      double lrJ=params->lrateJ;
	      if (params->learn_strat == 1) {
		GJ[i*q+a][j*q+b]+=gradJ*gradJ;
		lrJ /= sqrt(GJ[i*q+a][j*q+b] + 1e-12);
	      } else if (params->learn_strat == 2) {
		GJ[i*q+a][j*q+b] = params->rho * GJ[i*q+a][j*q+b] + (1.0 - params->rho)*gradJ*gradJ;
		lrJ /= sqrt(GJ[i*q+a][j*q+b] + 1e-12);
	      } else if(params->learn_strat == 3) {
		lrJ /= (1.0 + (double)iter/params->tau);
	      }
	      J[i*q + a][j*q + b] += lrJ * decJ[i*q + a][j*q + b] * gradJ;
	      J[j*q + b][i*q + a] = J[i*q + a][j*q + b];
	      lrav+=lrJ;
	      n++;
	    }
	  }
	}
      } 
      return lrav/n;
    } else {
      double P=0;
      double modF=0;
      double modv=0;
      for(int i = 0; i < L; i++) {
	for(int a = 0; a < q; a++) {
	  double gradh = fm[i*q+a] - fm_s[i*q+a];
	  h[i*q+a] += params->lrateh * acc * Gh[i*q+a];
	  Gh[i*q+a] += params->lrateh * acc * gradh;
	  P+=gradh*Gh[i*q+a];
	  modF+=gradh*gradh;
	  modv+=Gh[i*q+a]*Gh[i*q+a];
	  for(int j = i+1; j < L; j++) {
	    for(int b = 0; b <q; b++) {
	      double gradJ=sm[i*q+a][j*q+b] - sm_s[i*q+a][j*q+b];
	      J[i*q + a][j*q + b] += params->lrateJ * acc * decJ[i*q + a][j*q + b] * GJ[i*q + a][j*q + b];
	      J[j*q + b][i*q + a] = J[i*q + a][j*q + b];
	      GJ[i*q + a][j*q + b] += params->lrateJ * acc * decJ[i*q + a][j*q + b] * gradJ;
	      GJ[j*q + b][i*q + a] = GJ[i*q + a][j*q + b];
	      P+=gradJ*GJ[i*q + a][j*q + b];
	      modF+=gradJ*gradJ;
	      modv+=GJ[i*q+a][j*q + b]*GJ[i*q+a][j*q + b];
	    }
	  }
	}      
      }
      modF=sqrt(modF);
      modv=sqrt(modv);
      for(int i = 0; i < L; i++) {
	for(int a = 0; a < q; a++) {
	  Gh[i*q+a] = (1-alpha)*Gh[i*q+a] + alpha*(fm[i*q+a] - fm_s[i*q+a])/modF*modv;
	  for(int j = i+1; j < L; j++) {
	    for(int b = 0; b <q; b++) {
	      GJ[i*q+a][j*q + b] = (1-alpha)*GJ[i*q+a][j*q + b] + alpha*(sm[i*q+a][j*q+b] - sm_s[i*q+a][j*q+b])/modF*modv;
	    }
	  }
	}
      }
      if (P>=0) {
	counter++;
	if (counter>5) {
	  acc = min(acc*1.1,10);
	  alpha*=0.99;
	}
      } else {
	acc=0.5*acc;
	alpha=0.1;
	fill(Gh.begin(), Gh.end(), 0);
	fill(GJ.begin(), GJ.end(), Gh);
	counter=0;
      }
      cout<<"FIRE step - acc: "<<acc<<" alpha: "<<alpha<<" counter: "<<counter<<" P: "<<P/modF/modv<<endl;
      return acc*params->lrateh;
    }
  }


  
};















#endif




