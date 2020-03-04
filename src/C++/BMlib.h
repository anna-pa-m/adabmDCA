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

///FZ: TO BE ENCAPSULATED IN A CLASS
int q = 21, L, M;
double Meff;
int ** msa;
double * w;
int ** idx;
int * tmp_idx;
double * sorted_struct;
vector<double> fm;
vector< vector<double> > sm;
vector< vector<double> > cov;
double * tm;
int ** tm_index;
//END



class Params {
 public:
  char * file_msa, * file_freq, * file_w , * file_params, init, * label, * ctype, * file_3points, *file_cc, *file_samples, *file_en;
  bool Metropolis, Gibbs, rmgauge, dgap, gapnn, phmm, blockwise, compwise;
  double sparsity, rho, w_th,  regJ, lrateJ, lrateh, conv, pseudocount;
  int tau, seed, learn_strat, nprint, nprintfile, Teq, Nmc_starts, Nmc_config, Twait, maxiter;
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
  }
  int read_params (int & argc, char ** argv) {
    int c;
    while ((c = getopt(argc, argv, "y:b:f:w:l:u:v:s:n:m:p:j:t:o:i:a:c:z:g:e:k:x:S:d:T:C:MGIRAhDNE:HBWq:")) != -1) {
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
			case 'l':
				w_th = atof(optarg);
				break;
			case 's':
				Nmc_starts = atoi(optarg);
				if(Nmc_starts == 1) {
					fprintf(stderr, "You need at least 2 MC chains to check equilibration\n");
					Nmc_starts = 2;
				}
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
				fprintf(stdout, "-B : (flag) A block-wise decimation is applied to the couplings using sDKL as a criterion. Default: false\n");
				fprintf(stdout, "-W : (flag) A component-wise decimation is applied to the couplings using sDKL as a criterion. Default: false\n");
				fprintf(stdout, "-R : (flag) Zero J,H initialization\n");
				fprintf(stdout, "-I : (flag) Indipendent model initialization\n");
				fprintf(stdout, "-A : (flag) Remove gauge invariance\n");
				fprintf(stdout, "-D : (flag) DGap model: Only first moments of gap statistics are fitted\n");
				fprintf(stdout, "-N : (flag) GapNN model: Fit first moments and nearest-neighbors second moments gap statistics\n");
				fprintf(stdout, "-H : (flag) Hmmer model: profile + couplings(gap, gap) for nearest-neighbours\n");
				fprintf(stdout, "-P : (flag) Decimate J(i,j,a,b) looking at min(sec. mom., parameter)\n");
				fprintf(stdout, "-G : (flag) Using Gibbs sampling\n");
				fprintf(stdout, "-M : (flag) Using Metropolis-Hastings sampling\n");
				fprintf(stdout, "-y : Seed of random number generator, default: %d\n", seed);
				fprintf(stdout, "-s : Metropolis chains, default: %d\n", Nmc_starts);
				fprintf(stdout, "-n : Number of MC configurations per chain, default: %d\n", Nmc_config);
				fprintf(stdout, "-p : (optional file) Initial parameters J, h, default: random [-1e-3, 1e-3]\n");
				fprintf(stdout, "-c : Convergence tolerance, default: %.3e\n", conv);
				fprintf(stdout, "-e : Initial MC equilibration time (in MCsweeps), default: 20 \n");
				fprintf(stdout, "-t : Initial sampling time of MC algorithm (in MCsweeps), default: 10 \n");
				fprintf(stdout, "-i : Maximum number of iterations, default: %d\n", maxiter);
				fprintf(stdout, "-z : Print output every x iterations, default: %d\n", nprint);
				fprintf(stdout, "-m : Print Frobenius norms and parameters every x iterations, default: %d\n", nprintfile);
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
  
};

int convert_char_amino(char a) {
	int i;
	switch(a) {
		case '-':
			i = 0;
			break;
		case 'A':
			i = 1;
			break;
		case 'B':
			i = 0;
			break;
		case 'C':
			i = 2;
			break;
		case 'D':
			i = 3;
			break;
		case 'E':
			i = 4;
			break;
		case 'F':
			i = 5;
			break;
		case 'G':
			i = 6;
			break;
		case 'H':
			i = 7;
			break;
		case 'I':
			i = 8;
			break;
		case 'J':
			i = 0;
			break;
		case 'K':
			i = 9;
			break;
		case 'L':
			i = 10;
			break;
		case 'M':
			i = 11;
			break;
		case 'N':
			i = 12;
			break;
		case 'O':
			i = 0;
			break;
		case 'P':
			i = 13;
			break;
		case 'Q':
			i = 14;
			break;
		case 'R':
			i = 15;
			break;
		case 'S':
			i = 16;
			break;
		case 'T':
			i = 17;
			break;
		case 'U':
			i = 0;
			break;
		case 'V':
			i =18;
			break;
		case 'W':
			i = 19;
			break;
		case 'X':
			i = 0;
			break;
		case 'Y':
			i = 20;
			break;
		case 'Z':
			i = 0;
			break;
		default:
			fprintf(stderr, "%c not recognized\n", a);
			return(EXIT_FAILURE);
	}
	return i;
}

int convert_char_nbase(char a) {
	int i;
	switch(a) {
		case '-':
			i = 0;
			break;
		case 'A':
			i = 1;
			break;
		case 'U':
			i = 2;
			break;
		case 'T':
			i = 2;
			break;
		case 'C':
			i = 3;
			break;
		case 'G':
			i = 4;
			break;
		default:
			fprintf(stderr, "%c not recognized\n", a);
			i = 0;
			break;
	}
	return i;

}


int convert_char_epi(char a) {
	int i;
	switch(a) {
		case '-':
			i = 0;
			break;
		case 'A':
			i = 1;
			break;
		case 'F':
			i = 2;
			break;
		case '5':
			i = 3;
			break;
		case 'T':
			i = 4;
			break;
		case 't':
			i = 5;
			break;
		case 'G':
			i = 6;
			break;
		case 'E':
			i = 7;
			break;
		case 'Z':
			i = 8;
			break;
		case 'h':
			i = 9;
			break;
		case 'B':
			i = 10;
			break;
		case 'b':
			i = 11;
			break;
		case 'e':
			i = 12;
			break;
		case 'R':
			i = 13;
			break;
		case 'r':
			i = 14;
			break;
		case 'q':
			i = 15;
			break;
		default:
			fprintf(stderr, "%c not recognized, assuming '-'\n", a);
			i = 0;
			break;
			//return(EXIT_FAILURE);
	}
	return i;
}

int convert_char_ising(char a){
	int i;
	switch(a) {
		case 'A':
			i = 0;
			break;
		case 'P':
			i = 1;
			break;
		default:
			fprintf(stderr, "%c not recognized\n", a);
			return(EXIT_FAILURE);
	}
	return i;
}

int print_alphabet(char * ctype) {
  if(ctype == NULL) {
    fprintf(stdout, "Error in alphabet pointer\n");
    exit(EXIT_FAILURE);
  }  
  int q=0;
  if(!strcmp(ctype, "a")) {
    fprintf(stdout, "Using alphabet: -ACDEFGHIKLMNPQRSTVWY\n");
    q=21;
  } else if(!strcmp(ctype, "n")) {
    fprintf(stdout, "Using alphabet: -AUCG \n");
    q = 5;
  } else if(!strcmp(ctype, "i")) {
    fprintf(stdout, "Using alphabet: AP (0,1) \n");
    q = 2;
  } else if(!strcmp(ctype, "e")) {
    fprintf(stdout, "Using alphabet: -AF5TtGEZhBbeRrq \n");
    q = 16;
  } else {
    fprintf(stderr, "Use 'a' for amino-acids or 'n' for nitrogenous bases\n");
    return EXIT_FAILURE;
  }
  return q;
}



int alloc_structures(Params & params) {
  fm.clear();
  fm.resize(L*q,0);
  sm.clear();
  sm.resize(L*q,fm);
  cov.clear();
  cov.resize(L*q,fm);
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


int ** read_msa(char * filename, char * ctype, int & M, int & L, int & q) {
        fprintf(stdout, "Reading MSA from %s\n", filename);
	FILE * filemsa;
	char ch;
	int readseq = 0, newseq = 0, i;
        int l = 0, m = -1;
	char * auxseq;

	if(!filename || !(filemsa = fopen(filename, "r"))) {
		fprintf(stderr, "I couldn't open %s\n", filename);
		exit(EXIT_FAILURE);
	}

	auxseq = (char *) malloc(1);
	M = 0;
	L = 0;
	while((ch = fgetc(filemsa)) != EOF) {
		if (ch == '>')
		  M+= 1;
	}
	rewind(filemsa);
        int ** msa;
	while((ch = fgetc(filemsa)) != EOF) {
		if (ch == '>') {
			newseq = 1;
			readseq = 0;
			m += 1;
		} else if (ch == '\n' && newseq == 1) {
			readseq = 1;
			newseq = 0;
			l = 0;
		} else if (ch != '\n' && newseq == 0 && readseq == 1) {
			l += 1;
			if (L == 0 && l > 1) {
			  auxseq = (char *) realloc(auxseq, l * sizeof(int));
				memset(auxseq + l, 0, sizeof(int));
			}
			if(!strcmp(ctype, "a"))
				auxseq[l-1] = convert_char_amino(ch);
			else if(!strcmp(ctype, "n"))
				auxseq[l-1] = convert_char_nbase(ch);
			else if(!strcmp(ctype, "i"))
				auxseq[l-1] = convert_char_ising(ch);
			else if(!strcmp(ctype, "e"))
				auxseq[l-1] = convert_char_epi(ch);
		} else if (ch == '\n' && newseq == 0 && readseq == 1) {
			if (L == 0) {
				L = l;
				msa = (int **)calloc(M, sizeof(int *));
				for (i = 0; i < M; i++)
					msa[i] = (int *)calloc(L, sizeof(int));
			}
			readseq = 0;
			for (i = 0; i < L; i++)
				msa[m][i] = auxseq[i];
		}
	}

	fprintf(stdout, "Reading alignment completed.\nM = %i L = %i q = %i \n", M, L, q);
	free(auxseq);
	fclose(filemsa);
	return msa;
}

int read_freq(vector<double> & fm, vector< vector<double> > & sm, vector< vector<double> > & cov, Params & params,int &M, int &L, int &q) {    /// FZ: TO BE CHECKED

  FILE * filefreq;
  int i, j, a, b;
  char ch, cha,chb, t;
  char tmp[1024];
  double aux;
  
  if(!params.file_freq || !(filefreq = fopen(params.file_freq, "r"))) {
    fprintf(stderr, "I couldn't open %s\n", params.file_freq);
    exit(EXIT_FAILURE);
  } else {
    fprintf(stdout, "Reading frequencies from %s\n", params.file_freq);
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
  fprintf(stdout, "L = %i, M = %i, q = %i, alphabet = %s\n", L, M, q, params.ctype);
  alloc_structures(params);
  rewind(filefreq);
  while (!feof(filefreq) && fgets(tmp, 1024, filefreq) && sscanf(tmp, "%c ", &t) == 1) {
    switch(t) {
    case 's':
      sscanf(tmp, "s %d %d %c %c %lf \n", &i, &j, &cha, &chb, &aux);  
      if(!strcmp(params.ctype, "a")) {
	a = convert_char_amino(cha); 
	b = convert_char_amino(chb);
      } else if(!strcmp(params.ctype, "n")) {
	a = convert_char_nbase(cha);
	b = convert_char_nbase(chb);
      } else if(!strcmp(params.ctype, "i")) {
	a = convert_char_ising(cha);
	b = convert_char_ising(chb);
      } else if(!strcmp(params.ctype, "e")) {
	a = convert_char_epi(cha);
	b = convert_char_epi(chb);
      }
      //	      printf("%d %d %d %d %lf\n", i,a,j,b,aux);
      if(i != j) {
	sm[i*q+a][j*q+b] = aux;
	sm[j*q+b][i*q+a] = aux;
      }
      break;
    case 'm':
      sscanf(tmp, "m %d %c %lf \n", &i, &ch, &aux);
      if(!strcmp(params.ctype, "a"))  
	a = convert_char_amino(ch);
      else if(!strcmp(params.ctype, "n")) 
	a = convert_char_nbase(ch);
      else if(!strcmp(params.ctype, "i")) 
	a = convert_char_ising(ch);
      else if(!strcmp(params.ctype, "e")) 
	a = convert_char_epi(ch);
      //		printf("%d %d %lf\n", i,a,aux);
      fm[i*q+a] = aux;
      break;
    }
  }
  
  for(i = 0; i < q*L; i++) 
    for(j = 0; j < q*L; j++) 
      cov[i][j] = sm[i][j] - fm[i]*fm[j];
  
  return 0;
  
}


double * compute_w(char * filename, char * label, double w_th, int ** msa, int & M, int & L) {
        double * w;
	if(filename) {
		fprintf(stdout, "Reading weights from file...");
		FILE *filew;
		if(!(filew = fopen(filename, "r"))) {
			fprintf(stderr, "File %s not found\n", filename);
			exit(EXIT_FAILURE);
		} else {
   		        w = (double *) calloc(M, sizeof(double));
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
		w = (double *) calloc(M, sizeof(double));
		int m = 0, n = 0, d = 0, l = 0;
		for (m = 0; m < M; m++)
			w[m] = 1.0;
		double h = L * (1 - w_th);
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
	return w;
}

int compute_empirical_statistics(vector<double> & fm, vector< vector<double> > & sm, vector< vector<double> > & cov, double & pseudocount, int ** msa, double * w, int M, int L, int q) {
  fprintf(stdout, "Computing empirical statistics...");
  fflush(stdout);
  double Meff = 0;
  for(int m = 0; m < M; m++) {
    Meff += w[m];
    for(int i = 0; i < L; i++) {
      fm[i*q + msa[m][i]] += w[m];
      for(int j = i+1; j < L; j++) {
	sm[i*q + msa[m][i]][j*q + msa[m][j]] += w[m];
      }
    }
  }

  if(!pseudocount)
    pseudocount = 1.0/Meff;
  for(int i = 0; i < L*q; i++) {
    fm[i] = fm[i] / Meff + pseudocount;
    sm[i][i] = fm[i];
  }
  for(int i = 0; i < L; i++) {
    for(int j = i+1; j < L; j++) {
      for (int a = 0; a < q; a++) {
	for(int b = 0; b < q; b++) {
	  sm[i*q+a][j*q+b] = sm[i*q+a][j*q+b] / Meff + pseudocount;
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
  return Meff;
}

int print_msa(char *filename, int ** msa, int M, int L) {
	FILE *fp;
	fp = fopen(filename, "w");

	for(int m = 0; m < M;m++) {
		for(int i = 0; i < L; i++)
			fprintf(fp, "%d ", msa[m][i]);
		fprintf(fp, "\n");
	}
	fflush(fp);
	fclose(fp);
	return 0;
}






#endif
