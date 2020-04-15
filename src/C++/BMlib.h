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

#ifndef BMlib
#define BMlib




class Params {
 public:
  char * file_msa, * file_freq, * file_w , * file_params, init, * label, * ctype, * file_3points, *file_cc, *file_samples, *file_en;
  bool Metropolis, Gibbs, rmgauge, dgap, gapnn, phmm, blockwise, compwise, persistent, initdata, overwrite;
  double sparsity, rho, w_th,  regJ, lrateJ, lrateh, conv, pseudocount;
  int tau, seed, learn_strat, nprint, nprintfile, Teq, Nmc_starts, Nmc_config, Twait, maxiter, dec_steps;
  Params() {
    file_msa = 0;
    file_w = 0;
    file_freq = 0;
    file_params =0;
    init = 'R';
    label = (char *) "nolabel";
    ctype = (char *) "a";
    file_3points =0;
    file_cc = 0;
    file_samples = 0;
    file_en=0;
    Metropolis = true;
    Gibbs = false;
    rmgauge = false;
    dgap = false;
    gapnn = false;
    phmm = false;
    blockwise = false;
    compwise = false;
    persistent = false;
    initdata = false;
    overwrite = true;
    sparsity = 0;
    rho = 0.9; // RMSprop reinforcement (learn_strat = 2)
    w_th = 0.2;
    regJ = 0.0;
    lrateJ = 5e-2;
    lrateh = 5e-2;
    conv = 8e-3;
    pseudocount = 0;
    tau = 1000; // tau parameter for search and converge (learn_strat = 3)
    seed = 0;
    learn_strat = 0;
    nprint = 100;
    nprintfile = 500;
    Teq = 20;
    Nmc_starts = 1000;
    Nmc_config = 50;
    Twait = 10;
    maxiter = 2000;
    dec_steps = INT_MAX;
  }

  int read_params (int & argc, char ** argv) {
    int c;
    while ((c = getopt(argc, argv, "y:b:f:w:l:u:v:s:n:m:p:j:t:o:i:a:c:z:g:e:k:x:S:d:T:C:X:MPQGIRAhDNFE:HBWq:")) != -1) {
		switch (c) {
			case 'b':
				ctype = optarg;
				break;
			case 'd':
				pseudocount = atof(optarg);
				break;
			case 'T':
				file_3points = optarg;
				break;
			case 'C':
				file_cc = optarg;
				break;
			case 'B':
				blockwise = true;
				break;
			case 'W':
				compwise = true;
				break;
			case 'S':
				file_samples = optarg;
				break;
			case 'E':
				file_en = optarg;
				break;
			case 'H':
				phmm = true;
				break;
			case 'D':
				dgap = true;
				break;
			case 'N':
				gapnn = true;
				break;
			case 'y':
				seed = atoi(optarg);
				break;
			case 'f':
				file_msa = optarg;
				break;
		        case 'q':
 				file_freq = optarg;
 				break;
			case 'w':
				file_w = optarg;
				break;
			case 'x':
				sparsity = atof(optarg);
				break;
			case 'X':
				dec_steps = atoi(optarg);
				break;
			case 'l':
				w_th = atof(optarg);
				break;
			case 's':
				Nmc_starts = atoi(optarg);
				break;
			case 'n':
				Nmc_config = atoi(optarg);
				break;
			case 'p':
				file_params = optarg;
				break;
			case 'e':
				Teq = atoi(optarg);
				break;
			case 't':
				Twait = atoi(optarg);
				break;
			case 'i':
				maxiter = atoi(optarg);
				break;
			case 'u':
				lrateJ = atof(optarg);
				break;
			case 'v':
				lrateh = atof(optarg);
				break;
			case 'a':
				learn_strat = atoi(optarg);
				break;
			case 'c':
				conv = atof(optarg);
				break;
			case 'z':
				nprint = atoi(optarg);
				break;
			case 'm':
				nprintfile = atoi(optarg);
				break;
			case 'F':
				overwrite = false;
				break;
			case 'g':
				regJ = atof(optarg);
				break;
			case 'A':
				rmgauge = true;
				break;
			case 'M':
			 	Metropolis = true;
				Gibbs = false;
				break;
			case 'G':
				Metropolis = false;
				Gibbs = true;
				break;
			case 'P':
				persistent = true;
				break;
			case 'Q':
				initdata = true;
				break;
			case 'R':
				init = 'R';
				break;
			case 'I':
				init = 'I';
				break;
			case 'k':
				label = optarg;
				break;
			case 'h':
				fprintf(stdout, "Here's the list of the instructions\n");
				fprintf(stdout, "-f : MSA alignment in FASTA format\n");
				fprintf(stdout, "-q : Read frequencies from file - no MSA\n");
				fprintf(stdout, "-d : Pseudo-count, default: 1/M \n");
				fprintf(stdout, "-b : Alphabet. \n  \ta : amino-acids. \n \tn : nucleic acids. \n \ti : present/absent. \n \te : epigenetic data. \n \tDefault: a\n");
				fprintf(stdout, "-w : (optional file) weights file\n");
				fprintf(stdout, "-S : (optional file) file name in which print configurations at convergence\n");
				fprintf(stdout, "-E : (optional file) file name in which print energy's configurations at convergence\n");
				fprintf(stdout, "-T : (optional file) (i j k a b c) indices for third order correlations\n");
				fprintf(stdout, "-C : (optional file) (i j a b) or (i j a b corr) input file for interaction graph \n");
				fprintf(stdout, "-l : Threshold for computing weigts, default: %.1f\n", w_th);
				fprintf(stdout, "-g : L1 Regularization (J parameters), default :%.3e\n", regJ);
				fprintf(stdout, "-k : Label used in output files\n");
				fprintf(stdout, "-x : Required sparsity. Add -B for block-wise decimation or -W for component-wise decimation\n");
				fprintf(stdout, "-X : Decimate every x steps even if not converged (default: infinite)\n");
				fprintf(stdout, "-B : (flag) A block-wise decimation is applied to the couplings using sDKL as a criterion. Default: false\n");
				fprintf(stdout, "-W : (flag) A component-wise decimation is applied to the couplings using sDKL as a criterion. Default: false\n");
				fprintf(stdout, "-R : (flag) Zero J,H initialization\n");
				fprintf(stdout, "-I : (flag) Indipendent model initialization\n");
				fprintf(stdout, "-A : (flag) Remove gauge invariance\n");
				fprintf(stdout, "-D : (flag) DGap model: Only first moments of gap statistics are fitted\n");
				fprintf(stdout, "-N : (flag) GapNN model: Fit first moments and nearest-neighbors second moments gap statistics\n");
				fprintf(stdout, "-H : (flag) Hmmer model: profile + couplings(gap, gap) for nearest-neighbours\n");
				fprintf(stdout, "-G : (flag) Using Gibbs sampling\n");
				fprintf(stdout, "-M : (flag) Using Metropolis-Hastings sampling\n");
				fprintf(stdout, "-P : (flag) Use persistent MC chains\n");
				fprintf(stdout, "-Q : (flag) Initialize MC chains in data points\n");
				fprintf(stdout, "-y : Seed of random number generator, default: %d\n", seed);
				fprintf(stdout, "-s : Metropolis chains, default: %d\n", Nmc_starts);
				fprintf(stdout, "-n : Number of MC configurations per chain, default: %d\n", Nmc_config);
				fprintf(stdout, "-p : (optional file) Initial parameters J, h\n");
				fprintf(stdout, "-c : Convergence tolerance, default: %.3e\n", conv);
				fprintf(stdout, "-e : Initial MC equilibration time (in MCsweeps), default: 20 \n");
				fprintf(stdout, "-t : Initial sampling time of MC algorithm (in MCsweeps), default: 10 \n");
				fprintf(stdout, "-i : Maximum number of iterations, default: %d\n", maxiter);
				fprintf(stdout, "-z : Print output every x iterations, default: %d\n", nprint);
				fprintf(stdout, "-m : Print Frobenius norms and parameters every x iterations, default: %d\n", nprintfile);
				fprintf(stdout, "-F : (flag) Do not overwrite temporary output\n");
				fprintf(stdout, "-u : Learning rate for couplings, default: %.e\n", lrateJ);
				fprintf(stdout, "-v : Learning rate for fields, default: %.e\n", lrateh);
				fprintf(stdout, "-a : Learning strategy.\n \t0: standard gradient descent\n \t1: adagrad\n \t2. RMSprop\n \t3. search then converge\n \t4. adam (currently not implemented)\n \t5. FIRE\n \tDefault: %d\n", learn_strat);
				exit(0);
			default:
			        fprintf(stdout, "Incorrect input. Run with -h to get a list of inputs.\n");
				exit(1);
		}
	}

    return 0;
  }
  
  void print_learning_strategy() {
    fprintf(stdout, "****** Initializing model ******\n");
    if(Nmc_starts % 2 == 1 || Nmc_starts <=0) {
      fprintf(stderr, "You need an even number of MC chains to check equilibration\n");
      exit(EXIT_FAILURE);
    }
    if(Metropolis)
      fprintf(stdout, "Performing Metropolis-Hastings MC.\nInitial sampling time: %d\nInitial equilibration time: %d\nUsing %d seeds and tot. number of points %d\n", Twait, Teq, Nmc_starts, Nmc_starts * Nmc_config);
    else if(Gibbs) {	
      fprintf(stdout, "Performing Gibbs sampling.\nInitial sampling time: %d\nInitial equilibration time: %d\nUsing %d seeds and tot. number of points %d\n", Twait, Teq, Nmc_starts, Nmc_starts * Nmc_config);
    }
    if (initdata) {
      fprintf(stdout, "MC chains are initialized using MSA sequences");
    } else {
      fprintf(stdout, "MC chains are randomly initialized");
    }
    if (persistent) {
      fprintf(stdout, " only at the beginning, and are then persistent\n");
    } else {
      fprintf(stdout, " at the beginning of each iteration\n");
    }

    fprintf(stdout, "Learning strategy: ");
    switch(learn_strat) {
    case 0:
      fprintf(stdout, "Using standard gradient descent with constant learning rate (for J %.3e, for h %.3e)\n", lrateJ, lrateh);
      break;
    case 1:
      fprintf(stdout, "Using adagrad\n");
      break;
    case 2:
      fprintf(stdout, "Using RMSprop with reinforcement %.2f\n", rho);
      break;
    case 3:
      fprintf(stdout, "Using search and converge with decay time %d and learning rate (for J %.3e, for h %.3e)\n", tau, lrateJ, lrateh);
      break;
    case 4:
      fprintf(stdout, "Using adam : NOT YET IMPLEMENTED, EXIT\n");
      exit(EXIT_FAILURE);
      break;
    case 5:
      fprintf(stdout, "Using FIRE\n");
      break;
    }
    if(sparsity > 0.0) {
      if (regJ != 0.0) {
	fprintf(stdout, "L1 regularization is not compatible with sparsification\n");
	exit(1);
      }
      fprintf(stdout, "Required sparsity %.3f", sparsity);
      if (!compwise && !blockwise) {
	fprintf(stdout, "\n Please use either -W or -B flags to specify decimation type! EXIT\n");
	exit(EXIT_FAILURE);
      } else if(compwise && !blockwise) {
	fprintf(stdout, " using (component-wise) Kullback-Leibler based decimation");
      } else if (!compwise && blockwise) {
	fprintf(stdout, " using (block-wise) Kullback-Leibler based decimation -- NOT YET IMPLEMENTED!\n");
	exit(EXIT_FAILURE);
      } else if (compwise && blockwise) {
	fprintf(stdout, "\n Please do not use -W and -B flags together! EXIT\n");
	exit(EXIT_FAILURE);
      }
      if (dec_steps==INT_MAX)
	fprintf(stdout, " at convergence\n");
      else
	fprintf(stdout, " every at most %d steps\n", dec_steps);
    }
    if(regJ > 0)
      fprintf(stdout, "L1 regularization on couplings: lambda %.1e\n", regJ);
    if(pseudocount)
      fprintf(stdout, "Using pseudo-count: %1.e\n", pseudocount);

  }

  void construct_filenames(int iter, bool conv, char * par, char * par_zsum, char * score, char * first, char * sec, char * third) {
      char sc = (Gibbs == 0) ? 'M' : 'G';
      if (!conv) {
	if (overwrite) {
	  sprintf(par, "Parameters_tmp_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", label, sc, init, lrateJ, lrateh, learn_strat);
	  sprintf(par_zsum, "Parameters_tmp_zerosum_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", label, sc, init, lrateJ, lrateh, learn_strat);
	  sprintf(score, "Score_tmp_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", label, sc, init, lrateJ, lrateh, learn_strat);
	  sprintf(first, "First_mom_tmp_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", label, sc, init, lrateJ, lrateh, learn_strat);
	  sprintf(sec, "Sec_mom_tmp_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", label, sc, init, lrateJ, lrateh, learn_strat);
	  sprintf(third, "Third_mom_tmp_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", label, sc, init, lrateJ, lrateh, learn_strat);
	} else {
	  sprintf(par, "Parameters_tmp_%d_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", iter, label, sc, init, lrateJ, lrateh, learn_strat);
	  sprintf(par_zsum, "Parameters_tmp_%d_zerosum_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", iter, label, sc, init, lrateJ, lrateh, learn_strat);
	  sprintf(score, "Score_tmp_%d_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", iter, label, sc, init, lrateJ, lrateh, learn_strat);
	  sprintf(first, "First_mom_tmp_%d_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", iter, label, sc, init, lrateJ, lrateh, learn_strat);
	  sprintf(sec, "Sec_mom_tmp_%d_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", iter, label, sc, init, lrateJ, lrateh, learn_strat);
	  sprintf(third, "Third_mom_tmp_%d_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", iter, label, sc, init, lrateJ, lrateh, learn_strat);
	}
      } else {
	sprintf(par, "Parameters_conv_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", label, sc, init, lrateJ, lrateh, learn_strat);
	sprintf(par_zsum, "Parameters_conv_zerosum_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", label, sc, init, lrateJ, lrateh, learn_strat);
	sprintf(score, "Score_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat",label, sc, init, lrateJ, lrateh, learn_strat);
	sprintf(first, "First_mom_conv_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", label, sc, init, lrateJ, lrateh, learn_strat);
	sprintf(sec, "Sec_mom_conv_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", label, sc, init, lrateJ, lrateh, learn_strat);
	sprintf(third, "Third_order_connected_corr_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", label, sc, init, lrateJ, lrateh, learn_strat);
      }
  }
  
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

 Data(Params * _params):
  params(_params) {
    fprintf(stdout, "****** Initializing data structures ******\n");
    q=print_alphabet(params->ctype);
    if(params->file_msa) {
      read_msa();
      compute_w();
      alloc_structures();
      compute_empirical_statistics();
    } else if(params->file_freq) {
      read_freq();
    }
    load_third_order_indices();
  }

  void read_msa() {
    char * filename=params->file_msa;
    fprintf(stdout, "Reading MSA from %s\n", filename);
    FILE * filemsa;
    char ch;
    int readseq = 0, newseq = 0;
    if(!filename || !(filemsa = fopen(filename, "r"))) {
      fprintf(stderr, "I couldn't open %s\n", filename);
      exit(EXIT_FAILURE);
    }
    vector<int> auxseq;
    M = 0;
    L = 0;
    msa.clear();
    while((ch = fgetc(filemsa)) != EOF) {
      if (ch == '>') {
	newseq = 1;
	readseq = 0;
      } else if (ch == '\n' && newseq == 1) {
	readseq = 1;
	newseq = 0;
	auxseq.clear();
      } else if (ch != '\n' && newseq == 0 && readseq == 1) {
	if(!strcmp(params->ctype, "a"))
	  auxseq.push_back(convert_char_amino(ch));
	else if(!strcmp(params->ctype, "n"))
	  auxseq.push_back(convert_char_nbase(ch));
	else if(!strcmp(params->ctype, "i"))
	  auxseq.push_back(convert_char_ising(ch));
	else if(!strcmp(params->ctype, "e"))
	  auxseq.push_back(convert_char_epi(ch));
      } else if (ch == '\n' && newseq == 0 && readseq == 1) {
	if (L == 0) {
	  L = auxseq.size();
	} else if (L!=auxseq.size()) {
	  cout<<"MSA reading error!"<<endl;
	  exit(1);
	}
	readseq = 0;
	msa.push_back(auxseq);
      }
    }
    M = msa.size();
    fprintf(stdout, "Reading alignment completed.\nM = %i L = %i q = %i \n", M, L, q);
    fclose(filemsa);
  }

  void compute_w() {
    w.clear();
    w.resize(M,0);
    char * filename = params->file_w;
    char * label = params->label;
    if(filename) {
      fprintf(stdout, "Reading weights from file...");
      FILE *filew;
      if(!(filew = fopen(filename, "r"))) {
	fprintf(stderr, "File %s not found\n", filename);
	exit(EXIT_FAILURE);
      } else {
	int i = 0;
	double tmp;
	while((fscanf(filew, "%lf", &tmp)) != EOF) {
	  w[i] = tmp;
	  i += 1;
	}
	fclose(filew);
      }
    } else {
      fprintf(stdout, "Computing weights...");
      fflush(stdout);
      FILE *fw;
      char file_w[1000];
      sprintf(file_w,  "Weights_%s.dat", label);
      fw = fopen(file_w, "w");
      int m = 0, n = 0, d = 0, l = 0;
      for (m = 0; m < M; m++)
	w[m] = 1.0;
      double h = L * (1 - params->w_th);
      for(m = 0; m < M-1; m++) {
	for (n = m+1; n < M; n++ ) {
	  d = 0;
	  for (l = 0; l < L; l++) {
	    if(msa[m][l] == msa[n][l])
	      d++;
	  }
	  if(d > h) {
	    w[n] += 1;
	    w[m] += 1;
	  }
	}
      }
      for (m = 0; m < M; m++) {
	w[m] = 1.0 / w[m];
	fprintf(fw,"%f\n", w[m]);
	fflush(fw);
      }
      fclose(fw);
    }
    fprintf(stdout,"done\n");
    fflush(stdout);
  }

  void alloc_structures() {
    fm.clear();
    fm.resize(L*q,0);
    sm.clear();
    sm.resize(L*q,fm);
    cov.clear();
    cov.resize(L*q,fm);
  }

  void compute_empirical_statistics() {
    fprintf(stdout, "Computing empirical statistics...");
    fflush(stdout);
    Meff = 0;
    for(int m = 0; m < M; m++) {
      Meff += w[m];
      for(int i = 0; i < L; i++) {
	fm[i*q + msa[m][i]] += w[m];
	for(int j = i+1; j < L; j++) {
	  sm[i*q + msa[m][i]][j*q + msa[m][j]] += w[m];
	}
      }
    }
    if(!params->pseudocount)
      params->pseudocount = 1.0/(1+Meff);
    for(int i = 0; i < L*q; i++) {
      fm[i] = (1-params->pseudocount)*fm[i]/Meff + params->pseudocount/q;
      sm[i][i] = fm[i];
    }
    for(int i = 0; i < L; i++) {
      for(int j = i+1; j < L; j++) {
	for (int a = 0; a < q; a++) {
	  for(int b = 0; b < q; b++) {
	    sm[i*q+a][j*q+b] = (1-params->pseudocount)*sm[i*q+a][j*q+b]/Meff + params->pseudocount/(q*q);
	    sm[j*q+b][i*q+a] = sm[i*q+a][j*q+b];
	  }
	}
      }
    }
    for(int i = 0; i < L*q; i++) {
      for(int j = 0; j < L*q; j++) {
	cov[i][j] = sm[i][j] - fm[i]*fm[j];
      }
    }
    fprintf(stdout, "Meff: %lf\n", Meff);
    fflush(stdout);
  }
  
  void read_freq() {    
    FILE * filefreq;
    int i, j, a, b;
    char ch, cha,chb, t;
    char tmp[1024];
    double aux;
    if(!params->file_freq || !(filefreq = fopen(params->file_freq, "r"))) {
      fprintf(stderr, "I couldn't open %s\n", params->file_freq);
      exit(EXIT_FAILURE);
    } else {
      fprintf(stdout, "Reading frequencies from %s\n", params->file_freq);
    }
    L = 0;
    while(!feof(filefreq) && fgets(tmp, 1024, filefreq) && sscanf(tmp, "%c ", &t) == 1) {
      switch (t) {
      case 'm':
	sscanf(tmp, "m %d %c %lf \n", &i, &ch, &aux);
	if(i+1 > L)
	  L = i+1;
	break;
      }
    }
    fprintf(stdout, "L = %i, M = %i, q = %i, alphabet = %s\n", L, M, q, params->ctype);
    alloc_structures();
    rewind(filefreq);
    while (!feof(filefreq) && fgets(tmp, 1024, filefreq) && sscanf(tmp, "%c ", &t) == 1) {
      switch(t) {
      case 's':
	sscanf(tmp, "s %d %d %c %c %lf \n", &i, &j, &cha, &chb, &aux);  
	if(!strcmp(params->ctype, "a")) {
	  a = convert_char_amino(cha); 
	  b = convert_char_amino(chb);
	} else if(!strcmp(params->ctype, "n")) {
	  a = convert_char_nbase(cha);
	  b = convert_char_nbase(chb);
	} else if(!strcmp(params->ctype, "i")) {
	  a = convert_char_ising(cha);
	  b = convert_char_ising(chb);
	} else if(!strcmp(params->ctype, "e")) {
	  a = convert_char_epi(cha);
	  b = convert_char_epi(chb);
	}
	if(i != j) {
	  sm[i*q+a][j*q+b] = aux;
	  sm[j*q+b][i*q+a] = aux;
	}
	break;
      case 'm':
	sscanf(tmp, "m %d %c %lf \n", &i, &ch, &aux);
	if(!strcmp(params->ctype, "a"))  
	  a = convert_char_amino(ch);
	else if(!strcmp(params->ctype, "n")) 
	  a = convert_char_nbase(ch);
	else if(!strcmp(params->ctype, "i")) 
	  a = convert_char_ising(ch);
	else if(!strcmp(params->ctype, "e")) 
	  a = convert_char_epi(ch);
	fm[i*q+a] = aux;
	break;
      }
    }    
    for(i = 0; i < q*L; i++) 
      for(j = 0; j < q*L; j++) 
	cov[i][j] = sm[i][j] - fm[i]*fm[j];  
  }


  /******************** METHODS FOR 3RD ORDER STATISTICS ***********************************************************/

  void load_third_order_indices() {
    if(params->file_3points) {
      fprintf(stdout, "Reading three points correlations indices...");
      FILE *file3;
      if(!(file3 = fopen(params->file_3points, "r"))) {
	fprintf(stderr, "File %s not found \n", params->file_3points);
	exit(1);
      } else {
	int i, j, k, a, b, c;
	double value;
	char buffer[1000];
	tm_index.clear();
	tm.clear();
	vector<int> tmp_vec(6,0);
	while (!feof(file3) && fgets(buffer, 1000, file3)) {
	  if(sscanf(buffer, "%d %d %d %d %d %d %lf \n", &i, &j, &k, &a, &b, &c, &value) == 7) {
	    tmp_vec[0] = i;
	    tmp_vec[1] = j;
	    tmp_vec[2] = k;
	    tmp_vec[3] = a;
	    tmp_vec[4] = b;
	    tmp_vec[5] = c;
	    tm_index.push_back(tmp_vec);
	    tm.push_back(0.0);	
	  }
	}
	fprintf(stdout, "Number of indices %d\n", (int)tm_index.size());
	// compute 3rd order moments of msa, slow!
	for(int ind =0; ind < tm_index.size(); ind++) {
	  i = tm_index[ind][0];
	  j = tm_index[ind][1];
	  k = tm_index[ind][2];
	  a = tm_index[ind][3];
	  b = tm_index[ind][4];
	  c = tm_index[ind][5];
	  for(int m = 0; m < M; m++) 
	    if(msa[m][i] == a && msa[m][j] == b && msa[m][k] == c)
	      tm[ind] += w[m]/Meff;
	}
      }
    } else {
      fprintf(stdout, "No three-points correlations indices specified\n");
    }  
  }
  
  /******************** METHODS FOR OUTPUT ***********************************************************/
  
  void print_msa(char *filename) {
    FILE *fp;
    fp = fopen(filename, "w");    
    for(int m = 0; m < M;m++) {
      for(int i = 0; i < L; i++)
	fprintf(fp, "%d ", msa[m][i]);
      fprintf(fp, "\n");
    }
    fflush(fp);
    fclose(fp);
  }

  int print_statistics(char *file_sm, char *file_fm, char *file_tm, vector<double> & fm_s, vector< vector<double> > & sm_s, vector<double> & tm_s) {
    FILE *fs, *ff, *ft;
    fs = fopen(file_sm, "w");
    ff = fopen(file_fm, "w");
    for(int i = 0; i < L; i++) {
      for(int j = i+1; j < L; j++) {
	for(int a = 0; a < q; a++) {
	  for(int b = 0; b < q; b++)
	    fprintf(fs, "%d %d %d %d %.5f %.5f %.5f %.5f\n",i, j,a,b, sm[i*q+a][j*q+b], sm_s[i*q+a][j*q+b], cov[i*q+a][j*q+b], sm_s[i*q+a][j*q+b]-fm_s[i*q+a]*fm_s[j*q+b]);
	}
      }
    }
    for(int i = 0; i < L; i++) {
      for(int a = 0; a < q; a++)
	fprintf(ff, "%i %i %.5f %.5f\n", i, a, fm[i*q+a], fm_s[i*q+a]);
    }
    if(tm_index.size()>0) {
      ft = fopen(file_tm, "w");
      int i, j, k, a, b, c;
      double aux;
      for(int ind = 0; ind < tm_index.size(); ind++) {
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
  
  
  
};





#endif
