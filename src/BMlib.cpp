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
#include <valarray>
#include "BMaux.h"
#include "BMlib.h"
#include "BMmc.h"
using namespace std;


Params::Params() {
	
    file_msa = 0;
    file_w = 0;
    file_freq = 0;
    file_params =0;
    init = 'R';
    label = (char *) "nolabel";
    ctype = 'a';
    file_3points =0;
    file_cc = 0;
    print_samples = false;
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
    adapt = true;
    dec_sdkl = true;
    dec_f = false;
    dec_J = false;
    sparsity = 0;
    num_threads = 1;
    rho = 0.9; // RMSprop reinforcement (learn_strat = 2)
    w_th = 0.2;
    regJ1 = 0.0;
    regJ2 = 0.0;
    lrateJ = 5e-2;
    lrateh = 5e-2;
    conv = 8e-3;
    pseudocount = 0;
    beta = 1.0;
    tau = 1000; // tau parameter for search and converge (learn_strat = 3)
    seed = 0;
    learn_strat = 0;
    nprinteq = false;
    nprint = 100;
    nprintfile = 500;
    Teq = 20;
    Nmc_starts = 1000;
    Nmc_config = 50;
    Twait = 10;
    Twait_last = 10;
    maxiter = 2000;
    dec_steps = INT_MAX;
  }

  int Params::read_params (int & argc, char ** argv) {
    int c;
    while ((c = getopt(argc, argv, "a:b:c:d:e:f:g:hi:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:ABC:DFGHIJ:KLMNPQRST:UVWX:")) != -1) {
		switch (c) {
			case 'b':
				ctype = optarg[0];
				break;
			case 'd':
				pseudocount = atof(optarg);
				break;
			case 'T':
				file_3points = optarg;
				break;
			case 'K':
				nprinteq = true;
				break;
			case 'C':
				file_cc = optarg;
				break;
			case 'B':
				blockwise = true;
				dec_f = false;
				dec_J = false;
				dec_sdkl = true;
				break;
			case 'W':
				compwise = true;
				dec_f = false;
				dec_J = false;
				dec_sdkl = true;
				break;
			case 'U':
				dec_sdkl = false;
				dec_f = true;
				dec_J = false;
				break;
			case 'V':
				dec_sdkl = false;
				dec_J = true;
				dec_f = false;
				break;
			case 'S':
				print_samples = true;
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
    		case 'J':
		        beta = atof(optarg);
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
				regJ1 = atof(optarg);
				break;
			case 'r':
				regJ2 = atof(optarg);
				break;
			case 'A':
				rmgauge = true;
				break;
			case 'L':
				adapt = false;
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
				cout << "Here is the list of all the instructions" << endl;
				cout << "### Basic run ###" << endl;
				cout << "-f : (file) MSA alignment in FASTA format" << endl;
				cout << "-b : (letter) Alphabet. " << endl << "\ta : amino-acids. " << endl << "\tn : nucleic acids. " << endl << "\ti : Ising " << endl << "\te : epigenetic data. " << endl << "\tDefault: a" << endl;
				cout << "-k : (string) Label used in output files" << endl;
				cout << "-c : (number) Convergence tolerance, default: " << conv << endl;
				cout << "-i : (number) Maximum number of iterations, default: " << maxiter << endl;
				cout << "-l : (number) Threshold used within the reweighting process, default: " << w_th << endl;
				cout << "-d : (number) Pseudo-count, default: 1/Meff" << endl;
				cout << "-O : (number) Number of threads" << endl;
				cout << "### Additional I/O files ###" << endl;
				cout << "-q : (file) Read frequencies from file - no MSA" << endl;
				cout << "-w : (file) weights file" << endl;
				cout << "-p : (file) Initial parameters J, h" << endl;
				cout << "-S : (file) file name used to print the MCMC configurations at convergence" << endl;
				cout << "-E : (file) file name used to print energy's configurations at convergence" << endl;
				cout << "-T : (file) (i j k a b c) indices for third order correlations" << endl;
				cout << "-C : (file) (i j a b) or (i j a b corr) input interaction graph" << endl;
				cout << "-m : (number) Print temporary Frobenius norms and parameters every x iterations, default: " << nprintfile << endl;
				cout << "-F : (flag) Do not overwrite temporary output" << endl;
				cout << "### Settings of the learning ###" << endl;
			    cout << "-e : (number) Initial MC equilibration time (in MCsweeps), default: " << Teq << endl;
				cout << "-t : (number) Initial sampling time of MC algorithm (in MCsweeps), default: " << Twait << endl;
				cout << "-s : (number) Number of the MC chains per thread, default: " << Nmc_starts << endl;
				cout << "-n : (number) Number of MC sampled configurations per chain, default: " << Nmc_config << endl;
				cout << "-g : (number) L1 Regularization (J parameters), default: " << regJ1 << endl;
				cout << "-r : (number) L2 Regularization (J parameters), default: " << regJ2 << endl;
				cout << "-R : (flag) Zero J,H initialization" << endl;
				cout << "-I : (flag) Indipendent model initialization" << endl;
				cout << "-A : (flag) Remove gauge-invariance" << endl;
				cout << "-G : (flag) Using Gibbs sampling" << endl;
				cout << "-M : (flag) Using Metropolis-Hastings sampling" << endl;
				cout << "-P : (flag) Use persistent MC chains" << endl;
				cout << "-Q : (flag) Initialize MC chains in data points" << endl;
				cout << "-L : (flag) Do not adapt Teq and Twait to achieve equilibration" << endl;
				cout << "-u : (number) Learning rate for couplings, default: " << lrateJ << endl;
				cout << "-v : (number) Learning rate for fields, default: " << lrateh << endl;
				cout << "-a : (number) Learning strategy." << endl << "\t0: standard gradient descent" << endl << "\t1: adagrad" << endl << "\t2. RMSprop" << endl << "\t3. search then converge" << endl << "\t4. pseudo-Newton" << endl << "\t5. FIRE" << endl << "\tDefault: " << learn_strat << endl;
				cout << "### Sparse Potts model ###" << endl;
				cout << "-x : (number) Required sparsity. Add -B for block-wise decimation or -W for component-wise decimation. Default: 1.0" << endl;
				cout << "-X : (number) Decimate every x steps even if not converged, default: infinite" << endl;
				cout << "-B : (flag) A block-wise decimation is applied to the couplings using the symKL. Default: false" << endl;
				cout << "-W : (flag) A component-wise decimation is applied to the couplings using the symKL. Default: false" << endl << "\tadd -U (flag) to use frequencies instead" << endl << "\tadd -V (flag) to use couplings instead" << endl;
				cout << "### Other max-ent models ###" << endl;
				cout << "-D : (flag) DGap model: Only the first moments of the gap statistics are fitted" << endl;
				cout << "-N : (flag) GapNN model: Fit the first moments and the nearest-neighbors second moments gap statistics" << endl;
				cout << "-H : (flag) Hmmer-like model: profile + couplings(gap, gap) for nearest-neighbours" << endl;
				cout << "### Miscellaneous ###" << endl;
				cout << "-y : (number) Seed of random number generator, default: " << seed << endl;
				cout << "-J : (number) Rescale initial given parameters (J,h) by argument, default: " << beta << endl;
				cout << "-K : (flag) Print information on equilibration, default: " << nprinteq << endl;
				cout << "-z : (number) Print output every x iterations, default: " << nprint << endl;
				exit(0);
			default:
			        cout << "Incorrect input. Run with -h to get a list of inputs." << endl;
				exit(1);
		}
	}

    return 0;
  }
  
  void Params::print_learning_strategy() {
    cout << "****** Initializing model ******" << endl;
    if(Nmc_starts % 2 == 1 || Nmc_starts <=0) {
      cerr <<"You need an even number of MC chains to check equilibration" << endl;
      exit(EXIT_FAILURE);
    }
    if(Metropolis) {
      cout << "Performing Metropolis-Hastings MC" << endl;
    } else if(Gibbs) {
      cout << "Performing Gibbs sampling" << endl;
    }
    if (adapt) {
      cout << "Adaptive sampling time. Initial sampling time: " << Twait << " initial equilibration time: " << Teq << endl;
    } else {
      cout << "Fixed sampling time: " << Twait << " fixed equilibration time: " << Teq << endl;
    }
    cout << "Using " << Nmc_starts << " seeds and tot. number of points " << Nmc_starts * Nmc_config * num_threads << endl;
    if (initdata) {
      cout << "MC chains are initialized using MSA sequences";
    } else {
      cout << "MC chains are randomly initialized";
    }
    if (persistent) {
      cout << " only at the beginning, and are then persistent" << endl;
    } else {
      cout << " at the beginning of each iteration" << endl;
    }

    cout << "Learning strategy: " << endl;
    switch(learn_strat) {
    case 0:
      cout << "Using standard gradient descent with constant learning rate (for J " << lrateJ << " for h " << lrateh << ")" << endl;
      break;
    case 1:
      cout << "Using adagrad" << endl;
      break;
    case 2:
      cout << "Using RMSprop with reinforcement " << rho << endl;
      break;
    case 3:
      cout << "Using search and converge with decay time " << tau << " and learning rate (for J " << lrateJ << ", for h " << lrateh << endl;
      break;
    case 4:
      cout << "Using pseudo-Newton " << endl;
      break;
    case 5:
      cout << "Using FIRE" << endl;
      break;
    }
    if(sparsity > 0.0) {
      if (regJ1 != 0.0 || regJ2 != 0.0) {
		cout << "Regularization is not compatible with sparsification because of gauge choice conflicts" << endl;
		exit(1);
      }
      cout << "Required sparsity " << sparsity;
      if (!compwise && !blockwise) {
		cerr << endl << "Please use either -W or -B flags to specify decimation type. EXIT" << endl;
		exit(EXIT_FAILURE);
      } else if(compwise && !blockwise) {
			cout << " using (component-wise) ";
			if(dec_sdkl)
				cout << "Kullback-Leibler based decimation" ;
			if(dec_f)
				cout << "second moments based decimation" ;
			if(dec_J)
				cout << "couplings based decimation" ;
      } else if (!compwise && blockwise) {
			cout << " using (block-wise) ";
			if(dec_sdkl)
				cout << "Kullback-Leibler based decimation" ;
			if(dec_f)
				cout << "second moments based decimation (this is equivalent to a random choice!)" ;
			if(dec_J)
				cout << "couplings based decimation" ;
      } else if (compwise && blockwise) {
			cout << endl << "Please do not use -W and -B flags together! EXIT" << endl;
			exit(EXIT_FAILURE);
      }
      if (dec_steps==INT_MAX)
		cout << " at convergence" << endl;
      else
		cout << " every at most "  << dec_steps << " steps" << endl;
    }
    if(regJ1 > 0)
      cout << "L1 regularization on couplings: lambda " << regJ1 << endl;
    if(regJ2 > 0)
      cout << "L2 regularization on couplings: lambda " << regJ2 << endl;
    if(pseudocount)
      cout << "Using pseudo-count: " << pseudocount << endl;

  }

  void Params::construct_filenames(int iter, bool conv, char * par, char * par_zsum, char * ene, char * corr, char * score, char * first, char * sec, char * third) {
      if (!conv) {
		if (overwrite) {
			sprintf(corr, "Corr_tmp_%s.dat", label);
	  		sprintf(par, "Parameters_tmp_%s.dat", label);
	  		sprintf(par_zsum, "Parameters_tmp_zerosum_%s.dat", label);
	  		sprintf(ene, "Sample_tmp_%s.dat", label);
	  		sprintf(score, "Score_tmp_%s.dat", label);
	  		sprintf(first, "First_mom_tmp_%s.dat", label);
	  		sprintf(sec, "Sec_mom_tmp_%s.dat", label);
	  		sprintf(third, "Third_mom_tmp_%s.dat", label);
		} else {
			sprintf(corr, "Corr_tmp_%d_%s.dat", iter, label);
	  		sprintf(par, "Parameters_tmp_%d_%s.dat", iter, label);
	  		sprintf(par_zsum, "Parameters_tmp_%d_zerosum_%s.dat", iter, label);
	  		sprintf(ene, "Sample_tmp_%d_%s.dat", iter, label);
	  		sprintf(score, "Score_tmp_%d_%s.dat", iter, label);
	  		sprintf(first, "First_mom_tmp_%d_%s.dat", iter, label);
	  		sprintf(sec, "Sec_mom_tmp_%d_%s.dat", iter, label);
	 		sprintf(third, "Third_mom_tmp_%d_%s.dat", iter, label);
		}
     } else {
        sprintf(corr, "Corr_conv_%s.dat",label);
		sprintf(par, "Parameters_conv_%s.dat", label);
		sprintf(par_zsum, "Parameters_conv_zerosum_%s.dat", label);
		sprintf(ene, "Sample_conv_%s.dat",  label);
		sprintf(score, "Score_%s.dat",label);
		sprintf(first, "First_mom_conv_%s.dat", label);
		sprintf(sec, "Sec_mom_conv_%s.dat", label);
		sprintf(third, "Third_order_connected_corr_%s.dat", label);
      }
  }
  


Data::Data(Params * _params):
  params(_params) {
    cout << "****** Initializing data structures ******" << endl;
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

void Data::read_msa() {
	char * filename=params->file_msa;
    	cout << "Reading MSA from " << params->file_msa << endl;
    	FILE * filemsa;
    	char ch;
    	int readseq = 0, newseq = 0;
    	if(!filename || !(filemsa = fopen(filename, "r"))) {
		cerr << "I couldn't open " << filename << endl;
		exit(EXIT_FAILURE);
	}
	vector<unsigned char> auxseq;
	M = 0;
	L = 0;
	msa.clear();
	while((ch = fgetc(filemsa)) != EOF) {
		if (ch == '>') {
			newseq = 1;
			readseq = 0;
			if (L == 0 && int(auxseq.size()) > 0) {
		 	 	L = int(auxseq.size());
			} else if (L!=int(auxseq.size())) {
		 	 	cout<<"MSA reading error!"<<endl;
		 	 	exit(1);
			}
			if(int(auxseq.size()) > 0)
				msa.push_back(auxseq);
		} else {
			if (ch == '\n' && newseq == 1) {
				readseq = 1;
				newseq = 0;
				auxseq.clear();
			} else if (ch != '\n' && newseq == 0 && readseq == 1) {
				if(params->ctype == 'a')
		 	 		auxseq.push_back(convert_char_amino(ch));
				else if(params->ctype == 'n')
		 	 		auxseq.push_back(convert_char_nbase(ch));
				else if(params->ctype == 'i')
		 	 		auxseq.push_back(convert_char_ising(ch));
				else if(params->ctype == 'e')
		 	 		auxseq.push_back(convert_char_epi(ch));
			}
		}
	}
	if(int(auxseq.size()) > 0)
		msa.push_back(auxseq); // last sequence
	M = int(msa.size());
	cout << "Reading alignment completed." << endl << "M = " << M << " L = " << L << " q = " << q << endl;
	fclose(filemsa);
  }

/*
  void Data::read_msa_old() {
    char * filename=params->file_msa;
    cout << "Reading MSA from " << params->file_msa << endl;
    FILE * filemsa;
    char ch;
    int readseq = 0, newseq = 0;
    if(!filename || !(filemsa = fopen(filename, "r"))) {
      cerr << "I couldn't open " << filename << endl;
      exit(EXIT_FAILURE);
    }
    vector<unsigned char> auxseq;
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
		if(params->ctype == 'a')
	  		auxseq.push_back(convert_char_amino(ch));
		else if(params->ctype == 'n')
	  		auxseq.push_back(convert_char_nbase(ch));
		else if(params->ctype == 'i')
	  		auxseq.push_back(convert_char_ising(ch));
		else if(params->ctype == 'e')
	  		auxseq.push_back(convert_char_epi(ch));
      } else if (ch == '\n' && newseq == 0 && readseq == 1) {
		if (L == 0) {
	  		L = int(auxseq.size());
		} else if (L!=int(auxseq.size())) {
	  		cout<<"MSA reading error!"<<endl;
	  		exit(1);
		}
		readseq = 0;
		msa.push_back(auxseq);
      }
    }
    M = int(msa.size());
    cout << "Reading alignment completed." << endl << "M = " << M << " L = " << L << " q = " << q << endl;
    fclose(filemsa);
  }
  */

  void Data::compute_w() {
    w.clear();
    w.resize(M,0);
    char * filename = params->file_w;
    char * label = params->label;
    if(filename) {
      cout << "Reading weights from " << params->file_w << "...";
      FILE *filew;
      if(!(filew = fopen(filename, "r"))) {
	cerr << "File " << params->file_w << " not found" << endl;
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
      cout << "Computing weights...";
      fflush(stdout);
      ofstream fw;
      char file_w[1000];
      sprintf(file_w,  "Weights_%s.dat", label);
      fw.open(file_w);
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
	fw << w[m] << endl;
      }
      fw.close();
    }
    cout << "done" << endl;
    fflush(stdout);
  }

  void Data::alloc_structures() {
    fm.clear();
    fm.resize(L*q,0);
    sm.clear();
    sm.resize(L*q,fm);
    cov.clear();
    cov.resize(L*q,fm);
  }

  void Data::compute_empirical_statistics() {
    cout << "Computing empirical statistics...";
    Meff = 0;
    for(int m = 0; m < M; m++) {
      Meff += w[m];
      for(int i = 0; i < L; i++) {
	if(params->ctype == 'i') {
	  fm[i] += w[m] * (2.0*msa[m][i] - 1.0);
	} else {
	  fm[i*q + msa[m][i]] += w[m];
	}
	for(int j = i+1; j < L; j++) {
	  if(params->ctype == 'i') {
	    sm[i][j] += w[m] * (2.0*msa[m][i]-1.0) * (2.0*msa[m][j]-1.0);
	  } else {
	    sm[i*q + msa[m][i]][j*q + msa[m][j]] += w[m];
	  }
	}
      }
    }
    if(!params->pseudocount)
      params->pseudocount = 1.0/(1+Meff);
    for(int i = 0; i < L*q; i++) {
      if(params->ctype == 'i') {
	fm[i] = (1-params->pseudocount)*fm[i]/Meff;
	sm[i][i] = fm[i];
      } else {	
	fm[i] = (1-params->pseudocount)*fm[i]/Meff + params->pseudocount/q;
	sm[i][i] = fm[i];
      }
    }
    for(int i = 0; i < L; i++) {
      for(int j = i+1; j < L; j++) {
	for (int a = 0; a < q; a++) {
	  for(int b = 0; b < q; b++) {
	    if(params->ctype == 'i') {
	      sm[i][j] = (1-params->pseudocount)*sm[i][j]/Meff;
	      sm[j][i] = sm[i][j];
	    } else {
	      sm[i*q+a][j*q+b] = (1-params->pseudocount)*sm[i*q+a][j*q+b]/Meff + params->pseudocount/(q*q);
	      sm[j*q+b][i*q+a] = sm[i*q+a][j*q+b];
	    }
	  }
	}
      }
    }
    for(int i = 0; i < L*q; i++) {
      for(int j = 0; j < L*q; j++) {
		cov[i][j] = sm[i][j] - fm[i]*fm[j];
      }
    }
    cout << "Meff: " << int(Meff) << endl;
  }
  
  void Data::read_freq() {    
    FILE * filefreq;
    int i, j, a = -1, b = -1;
    char ch, cha,chb, t;
    char tmp[1024];
    double aux;
    if(!params->file_freq || !(filefreq = fopen(params->file_freq, "r"))) {
      cerr << "I couldn't open " << params->file_freq << endl;
      exit(EXIT_FAILURE);
    } else {
      cout << "Reading frequencies from " << params->file_freq << endl;
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
    cout << "L = " << L << " M = " << M << " q = " << q << " alphabet = " << params->ctype << endl;
    alloc_structures();
    rewind(filefreq);
    while (!feof(filefreq) && fgets(tmp, 1024, filefreq) && sscanf(tmp, "%c ", &t) == 1) {
      switch(t) {
      case 's':
	sscanf(tmp, "s %d %d %c %c %lf \n", &i, &j, &cha, &chb, &aux);  
	if(params->ctype == 'a') {
	  a = convert_char_amino(cha); 
	  b = convert_char_amino(chb);
	} else if(params->ctype == 'n') {
	  a = convert_char_nbase(cha);
	  b = convert_char_nbase(chb);
	} else if(params->ctype == 'i') {
	  a = convert_char_ising(cha);
	  b = convert_char_ising(chb);
	} else if(params->ctype == 'e') {
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
	if(params->ctype == 'a')  
	  a = convert_char_amino(ch);
	else if(params->ctype == 'n') 
	  a = convert_char_nbase(ch);
	else if(params->ctype == 'i') 
	  a = convert_char_ising(ch);
	else if(params->ctype == 'e') 
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

  void Data::load_third_order_indices() {
    if(params->file_3points) {
      cout << "Reading three points correlations indices..." << endl;
      FILE *file3;
      if(!(file3 = fopen(params->file_3points, "r"))) {
	cerr <<"File " << params->file_3points << " not found " << endl;
	exit(1);
      } else {
	int i, j, k, a, b, c;
	double value;
	char buffer[1000];
	tm_index.clear();
	tm.clear();
	while (!feof(file3) && fgets(buffer, 1000, file3)) {
	  if(params->ctype == 'i'){
	    vector<int> tmp_vec(3,0);
	    if(sscanf(buffer, "%d %d %d %lf \n", &i, &j, &k, &value) == 4) {
	      tmp_vec[0] = i;
	      tmp_vec[1] = j;
	      tmp_vec[2] = k;
	      tm_index.push_back(tmp_vec);
	      tm.push_back(0.0);
	    }
	  } else {
	    vector<int> tmp_vec(6,0);
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
	}
	cout << "Number of indices " << (int)tm_index.size() << endl;
	// compute 3rd order moments of msa, slow!
	for(int ind =0; ind < int(tm_index.size()); ind++) {
	  if(params->ctype == 'i') {
	    i = tm_index[ind][0];
	    j = tm_index[ind][1];
	    k = tm_index[ind][2];
	    for(int m = 0; m < M; m++)
	      tm[ind] += w[m]/Meff * (2.0*msa[m][i]-1.0) * (2.0*msa[m][j]-1.0) * (2.0*msa[m][k] -1.0);
	  } else {
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
      }
    } else {
      cout << "No three-points correlations indices specified" << endl;
    }  
  }
  
  /******************** METHODS FOR OUTPUT ***********************************************************/
  
  void Data::print_msa(char *filename) {
    ofstream fp;
    fp.open(filename);    
    for(int m = 0; m < M;m++) {
      for(int i = 0; i < L; i++)
	fp << msa[m][i] << " ";
      fp << endl;
    }
    fp.close();
  }

  int Data::print_statistics(char *file_sm, char *file_fm, char *file_tm, char * file_c, valarray<float> & corr, vector<float> & fm_s, vector< vector<float> > & sm_s, vector<float> & tm_s) {
    ofstream fs;
    ofstream ff;
    ofstream ft;
	ofstream fc;
    fs.open(file_sm);
    ff.open(file_fm);
    fc.open(file_c);
	
    for (int i=0;i<int(corr.size());i++) 
      fc <<  i << " " << corr[i] << endl;
    fc.close();
    if(params->ctype == 'i') {
      for(int i = 0; i <L; i++) {
		for(int j= i+1; j < L; j++)
		  fs << i << " " << j << " " << sm[i][j] << " " << sm_s[i][j] << " " << cov[i][j] << " " << sm_s[i][j] - fm_s[i]*fm_s[j] << endl;
      }
      for(int i = 0; i < L; i++)
		ff << i << " " << fm[i] <<  " " << fm_s[i] << endl;
      if(int(tm_index.size())>0) {
	ft.open(file_tm);
	int i, j, k;
	double aux;
	for(int ind = 0; ind < int(tm_index.size()); ind++) {
	    i = tm_index[ind][0];
	    j = tm_index[ind][1];
	    k = tm_index[ind][2];
	    aux = tm[ind] - sm[i][j]*fm[k] - sm[i][k]*fm[j] -sm[j][k]*fm[i] +2*fm[i]*fm[j]*fm[k];
	    ft << i << " " << j << " " << k << aux << " " << tm_s[ind] << endl;
	}
      }
    } else {
      for(int i = 0; i < L; i++) {
	for(int j = i+1; j < L; j++) {
	  for(int a = 0; a < q; a++) {
	    for(int b = 0; b < q; b++)
	      fs << i << " " << j << " " << a << " " << b << " " << sm[i*q+a][j*q+b] << " " << sm_s[i*q+a][j*q+b] << " " << cov[i*q+a][j*q+b] << " " << sm_s[i*q+a][j*q+b]-fm_s[i*q+a]*fm_s[j*q+b] << endl;
	  }
	}
      }
      for(int i = 0; i < L; i++) {
	for(int a = 0; a < q; a++){
	  ff <<  i << " " << a << " " << fm[i*q+a] << " " << fm_s[i*q+a] << endl;
	}
      }
      if(int(tm_index.size())>0) {
	ft.open(file_tm);
	int i, j, k, a, b, c;
	double aux;
	for(int ind = 0; ind < int(tm_index.size()); ind++) {
	  i = tm_index[ind][0];
	  j = tm_index[ind][1];
	  k = tm_index[ind][2];
	  a = tm_index[ind][3];
	  b = tm_index[ind][4];
	  c = tm_index[ind][5];
	  aux = tm[ind] - sm[i*q+a][j*q+b]*fm[k*q+c] - sm[i*q+a][k*q+c]*fm[j*q+b] - sm[j*q+b][k*q+c]*fm[i*q+a] + 2*fm[i*q+a]*fm[j*q+b]*fm[k*q+c];
	  ft << i << " " << j << " " << k << " " << a << " " << b << " " << c << " " << aux << " " << tm_s[ind] << endl;
	}
	ft.close();
      }
    }
   
  
    fs.close();
    ff.close();
    return 0;
  }
  
  
  


