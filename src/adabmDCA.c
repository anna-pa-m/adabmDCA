#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
//#include <omp.h>


// files and parameters
struct Params {

	char * file_msa, * file_w , * file_params, init, * label, * ctype, * file_3points, *file_cc, *file_samples, *file_en;
	bool Metropolis, Gibbs, dgap, gapnn, empdec, phmm;
	double sparsity, rho, w_th,  regJ, lrateJ, lrateh, conv, pseudocount;
	int tau, seed, learn_strat, nprint, nprintfile, Teq, Nmc_starts, num_threads, Nmc_config, Nmc_config_max, Tcheck, Twait_max, Twait, maxiter;
} params = {
	.phmm = false,
	.ctype = 0,
	.pseudocount = 0,
	.gapnn = false,
	.dgap = false,
	.empdec = false,
	.seed = 0,
	.label = 0,
	.init = 'R',
	.tau = 1000, // tau parameter for search and converge (learn_strat = 3)
	.regJ = 0.0,
	.lrateJ = 5e-2,
	.lrateh = 5e-2,
	.conv = 8e-3,
	.file_msa = 0,
	.file_w = 0,
	.file_cc = 0,
	.file_params =0,
	.file_samples = 0,
	.Nmc_starts = 1000,
	.Nmc_config = 50,
	.w_th = 0.2,
	.num_threads = 1,
	.Twait = 2000,
	.Teq = 5000,
	.maxiter = 2000,
	.nprint = 100,
	.Tcheck = 50,
	.nprintfile = 500,
	.Twait_max = 500,
	.Nmc_config_max = 100,
	.learn_strat = 0, // learning strategy: learning_rate / (1 + iter/tau)
	.rho = 0.98, // adadelta reinforment
	.Metropolis = true,
	.Gibbs = false,
	.sparsity = 0
};

// structures and important parameters
double * w;
int ** idx;
int * tmp_idx;
double * sortedJ;
int ** msa;
double * fm;
double * h;
double * Gh;
double ** sm;
double ** cov;
double * tm;
double ** J;
double ** Gj;
double * fm_s;
double ** sm_s;
double * tm_s;
double ** mtj;
double ** vtj;
double * mth;
double * vth;
double ** dec;
double ** Jzs;
double * hzs;
int ** tm_index;
double * fb;
double * fi;
double * mean_ai;
double * mean_aj;
int * chain_init_1;
int * chain_init_2;
int * chain_fin_1;
int * chain_fin_2;
int * intra_chain;
int * half_chain;
int idx_chain_1, idx_chain_2;
int * curr_state;


int q = 21, L, M;
double gibbs_en = 0.0, averrh, averrJ, merrh, merrJ, errnorm, model_sp;
int iter = 0;
bool compute_tm = false;
bool print_samples = false;
bool print_en = false;
int ntm = 0;
double Meff;

// list of functions
double max(double x, double y);
double min(double x, double y);
int read_msa();
int preprocess_msa();
int compute_w();
int alloc_structures();
int compute_empirical_statistics();
int compute_third_order_moments();
int compute_statistics();
int initialize_parameters();
int print_model(char * filename);
int print_msa(char * filename);
int print_statistics(char * file_sm, char *file_fm, char *file_tm);
int sample();
int mc_chain(bool test_chain_sampl, bool test_chain_half, int test_chain_eq);
int gibbs_step();
int metropolis_step();
void permute(int * vector, int s);
double energy(int *seq);
int update_parameters(double regJ);
double rand01();
double randrange(int xmin, int xmax);
double adaptive_learn(int i, int a, int j, int b, char type, double grad);
int compute_frobenius_norms(char *filename, char *parfile);
double compute_gibbs_z(int i, int *x);
int init_statistics();
int update_statistics(int *x);
int update_tm_statistics(int *x);
bool check_convergence();
double pearson();
int decimate(int c);
int quicksort(double *x, int *tmp_idx, int first, int last);
int load_third_order_indices();
int compute_third_order_correlations();
int equilibration_test();


int main(int argc, char ** argv)
{
	int c, n;
	bool conv = false;
	char score[1000];
	char sec[1000];
	char first[1000];
	char third[1000];
	char par[1000];
	char sc;
	while ((c = getopt(argc, argv, "y:b:f:w:l:u:v:s:n:m:p:j:t:i:a:c:z:g:e:k:x:S:d:T:C:MGIRhDNE:HPo:")) != -1) {
		switch (c) {
			case 'o':
				params.Tcheck = atoi(optarg);
				break;
			case 'b':
				params.ctype = optarg;
				break;
			case 'd':
				params.pseudocount = atof(optarg);
				break;
			case 'T':
				params.file_3points = optarg;
				break;
			case 'C':
				params.file_cc = optarg;
				break;
			case 'S':
				params.file_samples = optarg;
				break;
			case 'E':
				params.file_en = optarg;
				break;
			case 'P':
				params.empdec = true;
				break;
			case 'H':
				params.phmm = true;
				break;
			case 'D':
				params.dgap = true;
				break;
			case 'N':
				params.gapnn = true;
				break;
			case 'y':
				params.seed = atoi(optarg);
				break;
			case 'f':
				params.file_msa = optarg;
				break;
			case 'w':
				params.file_w = optarg;
				break;
			case 'x':
				params.sparsity = atof(optarg);
				break;
			case 'l':
				params.w_th = atof(optarg);
				break;
			case 's':
				params.Nmc_starts = atoi(optarg);
				break;
			case 'n':
				params.Nmc_config = atoi(optarg);
				break;
			case 'p':
				params.file_params = optarg;
				break;
			case 'j':
				params.num_threads = atoi(optarg);
				break;
			case 'e':
				params.Teq = atoi(optarg);
				break;
			case 't':
				params.Twait = atoi(optarg);
				break;
			case 'i':
				params.maxiter = atoi(optarg);
				break;
			case 'u':
				params.lrateJ = atof(optarg);
				break;
			case 'v':
				params.lrateh = atof(optarg);
				break;
			case 'a':
				params.learn_strat = atoi(optarg);
				break;
			case 'c':
				params.conv = atof(optarg);
				break;
			case 'z':
				params.nprint = atoi(optarg);
				break;
			case 'm':
				params.nprintfile = atoi(optarg);
				break;
			case 'g':
				params.regJ = atof(optarg);
				break;
			case 'M':
			 	params.Metropolis = true;
				params.Gibbs = false;
				break;
			case 'G':
				params.Metropolis = false;
				params.Gibbs = true;
				break;
			case 'R':
				params.init = 'R';
				break;
			case 'I':
				params.init = 'I';
				break;
			case 'k':
				params.label = optarg;
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
				fprintf(stdout, "-C : (optional file) (i j a b corr) given graph for correlations compressed\n");
				fprintf(stdout, "-l : Threshold for computing weigts, default: %.1f\n", params.w_th);
				fprintf(stdout, "-g : Regularization (J parameters), default :%.3e\n", params.regJ);
				fprintf(stdout, "-k : Label used in output files\n");
				fprintf(stdout, "-x : Required sparsity. If this value is larger than zero, l1 par = 0 and decimation is applied\n");
				fprintf(stdout, "-R : (flag) Random initialization of couplings and field in the range [-1e-3, 1e-3]\n");
				fprintf(stdout, "-I : (flag) Indipendent model initialization\n");
				fprintf(stdout, "-G : (flag) Using Gibbs sampling\n");
				fprintf(stdout, "-M : (flag) Using Metropolis-Hastings sampling\n");
				fprintf(stdout, "-D : (flag) Only first moments of gap statistics are fitted\n");
				fprintf(stdout, "-N : (flag) Fit first moments and nearest-neighbors second moments gap statistics\n");
				fprintf(stdout, "-P : (flag) Decimate J(i,j,a,b) looking at min(sec. mom., parameter)\n");
				fprintf(stdout, "-H : (flag) Hmmer-like model: profile + couplings(gap, gap) for nearest-neighbours\n");
				fprintf(stdout, "-y : Seed of random number generator, default: %d\n", params.seed);
				fprintf(stdout, "-o : Perform equilibration check every %d iterations\n", params.Tcheck);
				fprintf(stdout, "-s : Metropolis chains, default: %d\n", params.Nmc_starts);
				fprintf(stdout, "-n : Number of MC configurations per chain, default: %d\n", params.Nmc_config);
				fprintf(stdout, "-p : (optional file) Initial parameters J, h, default: random [-1e-3, 1e-3]\n");
				fprintf(stdout, "-c : Convergence tolerance, default: %.3e\n", params.conv);
				fprintf(stdout, "-j : Number of threads, default: %d\n", params.num_threads);
				fprintf(stdout, "-e : MC Equilibration time, default: %d\n", params.Teq);
				fprintf(stdout, "-t : Sampling time of MC algorithm, default: %d\n", params.Twait);
				fprintf(stdout, "-i : Maximum number of iterations, default: %d\n", params.maxiter);
				fprintf(stdout, "-z : Print output every x iterations, default: %d\n", params.nprint);
				fprintf(stdout, "-m : Print Frobenius norms and parameters every x iterations, default: %d\n", params.nprintfile);
				fprintf(stdout, "-u : Learning rate for couplings, default: %.e\n", params.lrateJ);
				fprintf(stdout, "-v : Learning rate for fields, default: %.e\n", params.lrateh);
				fprintf(stdout, "-a : Learning strategy.\n \t0: standard gradient descent\n \t1: adagrad\n \t2. adadelta\n \t3. search then converge\n \t4. adam\n \tDefault: %d\n", params.learn_strat);
				
				return(EXIT_FAILURE);
			default:
				return(EXIT_FAILURE);
		}
	}
	fprintf(stdout, "Boltzmann machine for DCA model\n");
	srand(params.seed ? params.seed : time(NULL));
	if(params.ctype == NULL)
		params.ctype = "a";
	if(!strcmp(params.ctype, "a")) {
		fprintf(stdout, "Using alphabet: -ACDEFGHIKLMNPQRSTVWY\n");
	} else if(!strcmp(params.ctype, "n")) {
		fprintf(stdout, "Using alphabet: -AUCG \n");
		q = 5;
	} else if(!strcmp(params.ctype, "i")) {
		fprintf(stdout, "Using alphabet: AP (0,1) \n");
		q = 2;
	} else if(!strcmp(params.ctype, "e")) {
		fprintf(stdout, "Using alphabet: -AF5TtGEZhBbeRrq \n");
		q = 16;
	} else {
		fprintf(stderr, "Use 'a' for amino-acids or 'n' for nitrogenous bases\n");
		return EXIT_FAILURE;
	}

	read_msa();
	if(params.Metropolis)
		fprintf(stdout, "Performing Metropolis-Hastings MC: sample every %d configurations (eq. time is %d), using %d seeds. \nTot number of points %d\n", params.Twait, params.Teq, params.Nmc_starts, params.Nmc_starts * params.Nmc_config);
	else if(params.Gibbs) {
		fprintf(stdout, "Performing Gibbs sampling: sample every %d configurations (eq. time is %d), using %d seeds.\nTot number of points %d\n", params.Twait, params.Teq, params.Nmc_starts, params.Nmc_starts * params.Nmc_config);
	}
	fprintf(stdout, "Learning strategy:\n");
	switch(params.learn_strat) {
		case 0:
			fprintf(stdout, "Using standard gradient descent with constant learning rate (for J %.3e, for h %.3e)\n", params.lrateJ, params.lrateh);
			break;
		case 1:
			fprintf(stdout, "Using adagrad\n");
			break;
		case 2:
			fprintf(stdout, "Using adadelta with reinforcement %.2f\n", params.rho);
			break;
		case 3:
			fprintf(stdout, "Using search and converge with decay time %d and learning rate (for J %.3e, for h %.3e)\n", params.tau, params.lrateJ, params.lrateh);
			break;
		case 4:
			fprintf(stdout, "Using adam\n");
			break;
	}
	if(params.sparsity > 0.0) {
		params.regJ = 0;
		fprintf(stdout, "Sparsity %.3f, using decimation\n", params.sparsity);
	} else
		fprintf(stdout, "L1 regularization on couplings: lambda %.1e\n", params.regJ);
	if(params.pseudocount)
		fprintf(stdout, "Using pseudo-count: %1.e\n", params.pseudocount);
	//omp_set_num_threads(params.num_threads);
	compute_w();
	alloc_structures();
	compute_empirical_statistics();
	initialize_parameters();
	iter = 1;
	fprintf(stdout, "Printing every %d iterations\n", params.nprint);
	if(params.dgap) {
		n = (L*(L-1)*(q-1)*(q-1))/2;
	} if(params.gapnn) {
	        n = (L*(L-1)*(q-1)*(q-1))/2 + L - 1;
	} if(params.phmm) {
		n = L + L - 1;
	} else {
	        n = (L*(L-1)*q*q)/2;
	}
	while(!conv && iter < params.maxiter) {
		init_statistics();
		sample();
		if(iter % params.Tcheck == 0)
			equilibration_test();
		update_parameters(params.regJ);
		if(model_sp < params.sparsity && iter % 10 == 0) {
			fprintf(stdout, "Decimating..");
			decimate(ceil((n*10)/(params.maxiter*0.5)));
		}
		if(check_convergence() && params.sparsity == 0)
			conv = true;
		if(model_sp >= params.sparsity && params.sparsity > 0 && check_convergence())
			conv = true;
		if(iter % params.nprint == 0) {
			fprintf(stdout, "it: %i N: %i Teq: %i Twait: %i merr_fm: %.3e merr_sm: %.3e averr_fm: %.3e averr_sm: %.3e cov_err: %.3e corr: %.2f sp: %.2e\n", iter, params.Nmc_config * params.Nmc_starts, params.Teq, params.Twait, merrh, merrJ, averrh, averrJ, errnorm, pearson(), model_sp);
			fflush(stdout);
		}
		if(iter % params.nprintfile == 0) {
			sc = (params.Gibbs == 0) ? 'M' : 'G';
			sprintf(par, "Parameters_tmp_zerosum_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
			sprintf(score, "Score_tmp_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
			compute_frobenius_norms(score, par);
			sprintf(par, "Parameters_tmp_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
			print_model(par);
			sprintf(sec, "Sec_mom_tmp_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
			sprintf(first, "First_mom_tmp_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
			sprintf(third, "Third_mom_tmp_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
			print_statistics(sec, first, third);
		}
		iter++;
	}
	load_third_order_indices();
	if(params.file_samples)
		print_samples = true;
	if(params.file_en)
		print_en = true;
	init_statistics();
	sample(); // compute 3rd order moments through sampling and possibly print the sequences
	if(compute_tm)
		compute_third_order_correlations();
	sc = (params.Gibbs == 0) ? 'M' : 'G';
	sprintf(par, "Parameters_conv_zerosum_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
	sprintf(score, "Score_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat",params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
	compute_frobenius_norms(score, par);
	sprintf(sec, "Sec_mom_conv_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
	sprintf(first, "First_mom_conv_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
	sprintf(third, "Third_order_connected_corr_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
	print_statistics(sec, first, third);
	sprintf(par, "Parameters_conv_%s_%c_%c_lJ%.1e_lh%.1e_a%i.dat", params.label, sc, params.init, params.lrateJ, params.lrateh, params.learn_strat);
	print_model(par);
	return 0;
}

int convert_char_amino(char a)
{
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

int convert_char_nbase(char a) 
{
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


int convert_char_epi(char a)
{
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

int convert_char_ising(char a)
{
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

int read_msa()
{
	FILE * filemsa;
	char ch;
	int readseq = 0, newseq = 0, i;
        int l = 0, m = -1;
	char * auxseq;

	if(!params.file_msa || !(filemsa = fopen(params.file_msa, "r"))) {
		fprintf(stderr, "I couldn't open %s\n", params.file_msa);
		return EXIT_FAILURE;
	}

	auxseq = (char *) malloc(1);
	M = 0;
	L = 0;
	while((ch = fgetc(filemsa)) != EOF) {
		if (ch == '>')
			M+= 1;
	}
	rewind(filemsa);
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
				auxseq = realloc(auxseq, l * sizeof(int));
				memset(auxseq + l, 0, sizeof(int));
			}
			if(!strcmp(params.ctype, "a"))
				auxseq[l-1] = convert_char_amino(ch);
			else if(!strcmp(params.ctype, "n"))
				auxseq[l-1] = convert_char_nbase(ch);
			else if(!strcmp(params.ctype, "i"))
				auxseq[l-1] = convert_char_ising(ch);
			else if(!strcmp(params.ctype, "e"))
				auxseq[l-1] = convert_char_epi(ch);
			//auxseq[l-1] = (!strcmp(params.ctype, "a")) ? convert_char_amino(ch) : convert_char_nbase(ch);
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
	return 0;
}

int load_third_order_indices()
{
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

int compute_third_order_correlations() 
{
	int ind, i, j, k, a, b, c;
	for(ind = 0; ind < ntm; ind++) {
		i = tm_index[ind][0];
		j = tm_index[ind][1];
		k = tm_index[ind][2];
		a = tm_index[ind][3];
		b = tm_index[ind][4];
		c = tm_index[ind][5];
		tm[ind] = tm[ind] - sm[i*q+a][j*q+b]*fm[k*q+c] - sm[i*q+a][k*q+c]*fm[j*q+b] - sm[j*q+b][k*q+c]*fm[i*q+a] + 2*fm[i*q+a]*fm[j*q+b]*fm[k*q+c];
		tm_s[ind] = tm_s[ind] - sm_s[i*q+a][j*q+b]*fm_s[k*q+c] - sm_s[i*q+a][k*q+c]*fm_s[j*q+b] - sm_s[j*q+b][k*q+c]*fm_s[i*q+a] + 2*fm_s[i*q+a]*fm_s[j*q+b]*fm_s[k*q+c]; 
	}
	return 0;
}

int compute_w()
{
	if(params.file_w) {
		fprintf(stdout, "Reading weights from file...");
		FILE *filew;
		if(!(filew = fopen(params.file_w, "r"))) {
			fprintf(stderr, "File %s not found\n", params.file_w);
			return EXIT_FAILURE;
		} else {
			w = calloc(M, sizeof(double));
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
		sprintf(file_w,  "Weights_%s.dat", params.label);
		fw = fopen(file_w, "w");
		w = calloc(M, sizeof(double));
		int m = 0, n = 0, d = 0, l = 0;
		for (m = 0; m < M; m++)
			w[m] = 1.0;
		double h = L * (1 - params.w_th);
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
	return 0;
}

int alloc_structures()
{
	int i;
	fm = (double *)calloc(L*q, sizeof(double));
	h = (double *)calloc(L*q, sizeof(double));
	fm_s = (double *)calloc(L*q, sizeof(double));
	sm = (double **)calloc(L*q, sizeof(double *));
	cov = (double **)calloc(L*q, sizeof(double *));
	J = (double **)calloc(L*q, sizeof(double *));
	dec = (double **)calloc(L*q, sizeof(double *));
	if(params.learn_strat == 1 || params.learn_strat == 2) {
		Gj = (double **)calloc(L*q, sizeof(double *));
		Gh = (double *) calloc(L*q, sizeof(double));
	}
	sm_s = (double **)calloc(L*q, sizeof(double *));
	for(i = 0; i < L*q; i++) {
		cov[i] = (double *)calloc(L*q, sizeof(double));
		sm[i] = (double *)calloc(L*q, sizeof(double));
		sm_s[i] = (double *)calloc(L*q, sizeof(double));
		J[i] = (double *)calloc(L*q, sizeof(double));
		dec[i] = (double *)calloc(L*q, sizeof(double));
	}
	if(params.learn_strat == 1 || params.learn_strat == 2) {
		for(i = 0; i < L*q; i++)
			Gj[i] = (double *)calloc(L*q, sizeof(double));
	}
	if(params.learn_strat == 4) {
		mtj = (double **)calloc(L*q, sizeof(double *));
		vtj = (double **)calloc(L*q, sizeof(double *));
		mth = (double *)calloc(L*q, sizeof(double));
		vth = (double *)calloc(L*q, sizeof(double));
		for(i = 0; i < L*q; i++) {
			vtj[i] = (double *)calloc(L*q, sizeof(double));
			mtj[i] = (double *)calloc(L*q, sizeof(double));
		}
	}
	if(params.sparsity > 0) {
		int n = L*(L-1)*q*q/2;
		idx = (int **)calloc(n, sizeof(int *));
		sortedJ = (double *)calloc(n, sizeof(double));
		tmp_idx = (int *)calloc(n, sizeof(int));
		for(i = 0; i < n;i++)
			idx[i] = (int *)calloc(4, sizeof(int));
	}
	Jzs = (double **)calloc(L*q, sizeof(double *));
	for(i = 0; i < L*q; i++)
		Jzs[i] = (double *)calloc(L*q, sizeof(double));
	hzs = (double *)calloc(L*q, sizeof(double));
	fb = (double *)calloc(L*L, sizeof(double));
	fi = (double *)calloc(L, sizeof(double));
	mean_ai = (double *)calloc(q, sizeof(double));
	mean_aj = (double *)calloc(q, sizeof(double));
	chain_init_1 = (int *)calloc(L, sizeof(int));
	chain_init_2 = (int *)calloc(L, sizeof(int));
	half_chain = (int *)calloc(L, sizeof(int));
	chain_fin_1 = (int *)calloc(L, sizeof(int));
	chain_fin_2 = (int *)calloc(L, sizeof(int));
	intra_chain = (int *)calloc(L, sizeof(int));
	curr_state = (int *)calloc(L, sizeof(int));
	return 0;
}

int compute_empirical_statistics()
{
	int m, i, j, a, b;
	Meff = 0;

	for(m = 0; m < M; m++) {
		Meff += w[m];
		for(i = 0; i < L; i++) {
			fm[i*q + msa[m][i]] += w[m];
			sm[i*q + msa[m][i]][i*q + msa[m][i]] += w[m];
			for(j = i+1; j < L; j++) {
				sm[i*q + msa[m][i]][j*q + msa[m][j]] += w[m];
				sm[j*q + msa[m][j]][i*q + msa[m][i]] += w[m];
			}
		}
	}

	if(!params.pseudocount)
		params.pseudocount = 1.0/Meff;
	for(i = 0; i < L; i++) {
		for(j = i+1; j < L; j++) {
			for (a = 0; a < q; a++) {
				for(b = 0; b < q; b++) {
					sm[i*q+a][j*q+b] = (sm[i*q+a][j*q+b] == 0) ? params.pseudocount : sm[i*q+a][j*q+b] / Meff + params.pseudocount;
					sm[j*q+b][i*q+a] = sm[i*q+a][j*q+b];
				}
			}
		}
	}
	for(i = 0; i < L*q; i++) {
		sm[i][i] = (sm[i][i] == 0) ? params.pseudocount : sm[i][i] / Meff + params.pseudocount;
		fm[i] = (fm[i] == 0) ? params.pseudocount : fm[i] / Meff + params.pseudocount;
	}

	for(i = 0; i < L; i++) {
		for(j = 0; j < L; j++) {
			for(a = 0; a < q; a++) {
				for(b = 0; b < q; b++)
					cov[i*q+a][j*q+b] = sm[i*q+a][j*q+b] - fm[i*q+a]*fm[j*q +b];
			}
		}
	}
	fprintf(stdout, "Meff: %lf\n", Meff);
	return 0;
}

int init_statistics()
{
	int i, j, a, b;
	for(i = 0; i < L; i++) {
		for(a = 0; a < q; a++) {
			fm_s[i*q + a] = 0;
			sm_s[i*q + a][i*q + a] = 0;
			for(j = i+1; j < L; j++) {
				for(b = 0; b < q; b++) {
					sm_s[i*q + a][j*q + b] = 0;
					sm_s[j*q + a][i*q + b] = 0;
				}
			}
		}
	}

	return 0;

}

int update_statistics(int *x)
{
	int i, j;
	int Ns = params.Nmc_starts * params.Nmc_config;
	for(i = 0; i < L; i++) {
		fm_s[i*q + x[i]] += 1.0/Ns;
		sm_s[i*q + x[i]][i*q + x[i]] += 1.0/Ns;
		for(j = i+1; j < L; j++) {
			sm_s[i*q + x[i]][j*q + x[j]] += 1.0/Ns;
			sm_s[j*q + x[j]][i*q + x[i]] += 1.0/Ns;
		}
	}
	return 0;
}

int update_tm_statistics(int *x)
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

bool check_convergence()
{
	int i, j, a, b;
	errnorm = 0.0;
	for(i = 0; i < L; i++) {
		for(j = i+1; j < L; j++) {
			for(a = 0; a < q; a++) {
				for(b = 0; b < q; b++)
					errnorm = max(errnorm, dec[i*q + a][j *q + b] * fabs(cov[i*q+a][j*q+b] - sm_s[i*q + a][j*q+b] + fm_s[i*q +a]*fm_s[j*q +b]));
			}
		}
	}

	if(errnorm < params.conv)
		return true;
	else
		return false;
}

double pearson()
{
	double rho;
	double mean_cov_s, mean_cov, mean_prod, covxy, std_cov_s, std_cov;
	int i,j,a,b,n;
	double cov_s;
	double mean_x2, mean_y2;

	mean_cov_s = 0.0;
	mean_cov = 0.0;
	mean_prod = 0.0;
	mean_x2 = 0.0;
	mean_y2 = 0.0;
	n = 0;
	for(i = 0; i < L; i++) {
		for(j = i + 1; j < L; j++) {
			for(a = 0; a < q; a++) {
				for(b = 0; b < q; b++) {
					if(dec[i*q + a][j*q +b]) {
						n += 1;
						cov_s = sm_s[i*q+a][j*q+b] - fm_s[i*q+a]*fm_s[j*q+b];
						mean_cov_s += cov_s;
						mean_cov += cov[i*q+a][j*q +b];
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
	covxy = mean_prod - mean_cov_s * mean_cov;

	std_cov_s = sqrt(mean_y2 - mean_cov_s * mean_cov_s);
	std_cov = sqrt(mean_x2 - mean_cov * mean_cov);

	rho = covxy / (std_cov_s * std_cov);
	return rho;
}

int compute_frobenius_norms(char *filename, char *parfile)
{

	int i, j, a, b;

	double f = 0;
	FILE *fp, *fpp;
	fp = fopen(filename, "w");
	fpp = fopen(parfile, "w");
	double mean_all,mean_h;
	for(i = 0; i < L; i++) {
		for(j = 0; j < L; j++) {
			mean_all = 0.0;
			for(a = 0; a < q; a++) {
				mean_ai[a] = 0.0;
				mean_aj[a] = 0.0;
			}
			for(a = 0; a < q; a++) {
				for(b = 0; b < q; b++) {
					mean_ai[a] += J[i*q +b][j*q + a];
					mean_aj[a] += J[i*q +a][j*q + b];
					mean_all += J[i*q +a][j*q +b];
				}
			}
			for(a = 0; a < q; a++) {
				mean_ai[a] /= q;
				mean_aj[a] /= q;
			}
			mean_all /= q*q;
			for(a =0; a < q;a++) {
				for(b = 0; b <q;b++){
					Jzs[i*q+a][j*q+b] = J[i*q+a][j*q+b] - mean_ai[b] - mean_aj[a] + mean_all;
					Jzs[j*q+b][i*q+a] = Jzs[i*q+a][j*q+b];
				}
			}
		}
		mean_h = 0.0;
		for(a = 0; a < q; a++)
			mean_h += h[i*q+a];
		mean_h /= q;
		for(a = 0; a < q; a++)
			hzs[i*q+a] = h[i*q+a] - mean_h;
		for(j = 0; j < L; j++) {
			mean_all = 0.0;
			for(a = 0; a < q; a++){
				mean_aj[a] = 0;
				for(b = 0; b < q; b++){
					mean_aj[a] += J[i*q+a][j*q+b];
					mean_all += J[i*q+a][j*q+b];
				}
				mean_aj[a] = mean_aj[a] / q;
			}
			mean_all = mean_all / (q*q);
			for(a = 0; a < q; a++)
				hzs[i*q+a] += mean_aj[a] - mean_all;
		}
	}
//	for(i = 0; i < L; i++) {
//		mean_all = 0;
//		for(a = 0; a < q; a ++)
//			mean_all += h[i*q+a];
//		mean_all /= q;
//		for(a = 0; a < q; a++)
//			hzs[i*q+a] = h[i*q+a] - mean_all;
//	}
	for(i = 0; i < L; i++) {
		for(j = i+1; j < L; j++) {
			for(a = 0; a < q; a++) {
				for(b = 0; b < q; b++)
					fprintf(fpp, "J %d %d %d %d %lf\n", i ,j ,a,b,Jzs[i*q+a][j*q+b]);
			}
		}
	}
	for(i = 0; i < L; i++)
		for(a = 0; a < q; a++)
			fprintf(fpp, "h %d %d %f\n", i,a,hzs[i*q+a]);
	int nf = 0;
	for(i = 0; i < L; i++) {
		for(j = 0; j < L; j++) {
			fb[i*L+j] = 0;
			if(i != j) {
				for(a = 1; a < q; a++) {
					for(b = 1; b < q; b++)
						fb[i*L+j] += Jzs[i*q+a][j*q+b] * Jzs[i*q+a][j*q+b];
				}
				fb[i*L+j] = sqrt(fb[i*L+j]);
				fi[i] += fb[i*L+j];
				f += fb[i*L+j];
				nf++;
			}
		}
	}
	for(i = 0; i < L; i++)
		fi[i] /= (L-1);
	//f /= L*L - L;
	f /= nf;
	for(i = 0; i < L-1; i++) {
		for(j = i+1; j < L; j++) {
			fb[i*L+j] -= fi[i]*fi[j]/f;
			fprintf(fp, "%d, %d, %f\n", i+1,j+1,fb[i*L+j]);
		}
	}
	fflush(fp);
	fflush(fpp);
	fclose(fp);
	fclose(fpp);
	return 0;
}

int initialize_parameters()
{
	int i,j, a ,b;
	if(params.file_params) {
		fprintf(stdout, "Reading input parameters from file...");
		FILE *filep;
		if(!(filep = fopen(params.file_params, "r"))) {
			fprintf(stderr, "File %s not found\n", params.file_params);
			return EXIT_FAILURE;
		} else {
			int  a, b;
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
		switch(params.init) {
			case 'R':
				fprintf(stdout, "Zero-parameters initialization...");
				for(i = 0; i < L; i++) {
					for(a = 0; a < q; a++) {
						//h[i*q + a] = 1e-3 * randrange(-1,1);
						J[i*q + a][i*q + a] = 0;
						h[i*q + a] = 0.0;
						for(j = i+1; j < L; j++) {
							for (b = 0; b < q; b ++) {
								J[i*q + a][j*q + b] = 0.0;
								J[j*q + b][i*q + a] = J[i*q + a][j*q + b];
								//J[i*q + a][j*q + b] = 1e-3 * randrange(-1,1);
							}
						}
					}
				}
				fprintf(stdout, "done\n");
				fflush(stdout);
				break;
			case 'I':
				fprintf(stdout, "Initializing parameters using independent sites approximation...");
				double mean_all;
				for(i = 0; i < L; i++) {
					mean_all = 0;
					for(a = 0; a < q; a++) {
						h[i*q +a] = log(fm[i*q+a]);
						mean_all += h[i*q+a];
					}
					mean_all /= q;
					for(a = 0; a < q; a++)
						h[i*q+a] -= mean_all;
				}
				fprintf(stdout, "done\n");
				fflush(stdout);
		}
	}
	if(params.dgap)
		fprintf(stdout, "Using DGap model\n");
	else if(params.gapnn)
		fprintf(stdout, "Using GapNN model\n");
	else if(params.phmm)
		fprintf(stdout, "Using Hmmer-like model\n");
	if(!params.file_cc) {
		for(i = 0; i < L; i++) {
			for(j = i; j < L; j++) {
				for(a = 0; a < q; a++) {
					for(b = 0; b < q; b ++) {
						if(i == j) {
							dec[i*q+a][j*q+b] = 0.0;
							J[i*q+a][j*q+b] = 0.0;
						}
						if(params.dgap && (a == 0 || b == 0)) {
							dec[i*q + a][j*q + b] = 0.0;
							dec[j*q + b][i*q + a] = 0.0;
							J[i*q + a][j*q + b] = 0.0;
							J[j*q + b][i*q + a] = 0.0;
						} else if(params.gapnn) {
							if(abs(i -j) > 1) {
								if(a == 0 || b == 0) {
									dec[i*q + a][j*q + b] = 0.0;
									dec[j*q + b][i*q + a] = 0.0;
									J[i*q + a][j*q + b] = 0.0;
									J[j*q + b][i*q + a] = 0.0;
								} else {
									dec[i*q+a][j*q+b] = 1.0;
									dec[j*q+b][i*q+a] = 1.0;
								}

							} else {
								if((a == 0 && b > 0) || (b == 0 && a > 0)) {
									dec[i*q + a][j*q + b] = 0.0;
									dec[j*q + b][i*q + a] = 0.0;
									J[i*q + a][j*q + b] =0.0;
									J[j*q + b][i*q + a] =0.0;
								}else {
									dec[i*q+a][j*q+b] = 1.0;
									dec[j*a+b][i*q+a] = 1.0;
								}
							}
						} else if(params.phmm) {
							if(a == 0 && b == 0 && abs(i-j) == 1) {
								dec[i*q+a][j*q+b] = 1.0;
								dec[j*q+b][i*q+a] = 1.0;
							} else {
								dec[i*q + a][j*q + b] = 0.0;
								dec[j*q + b][i*q + a] = 0.0;
								J[i*q + a][j*q + b] = 0.0;
								J[j*q + b][i*q + a] = 0.0;
							}
						} else {
							dec[i*q + a][j*q + b] = 1.0;
							dec[j*q + b][i*q + a] = 1.0;
						}
					}
				}
			}
		}
	} else {
		fprintf(stdout, "Reading graph from file for correlation-compressed...");
		FILE *filep;
		int links = 0;
		if(!(filep = fopen(params.file_cc, "r"))) {
			fprintf(stderr, "File %s not found\n", params.file_params);
			return EXIT_FAILURE;
		} else {
			int  i,j, a, b;
			char buffer[100];
			double tmp;
			for(i = 0; i < L*q; i++) {
				for(j = 0; j < L*q; j++)
					dec[i][j] = 0.0;
			}
			while(!feof(filep) && fgets(buffer,100, filep) && 
				 ( sscanf(buffer, "%d %d %d %d %lf \n", &i, &j, &a, &b, &tmp) == 5 || sscanf(buffer, "%d %d %d %d \n", &i, &j, &a, &b) == 4) ) {
				dec[i*q+a][j*q+b] = 1.0;
				links++;

			}
			fclose(filep);
			for(i = 0; i < L; i++) {
				for(j = i+1; j < L; j++) {
					for(a = 0; a < q; a++) {
						for(b = 0; b < q; b ++) {
							if(dec[i*q+a][j*q+b] == 0.0) {
								J[i*q+a][j*q+b] = 0.0;
								J[j*q+b][i*q+a] = 0.0;
							}
						}
					}
				}
			}
		}
		fprintf(stdout, "done\nNumber of links %d\n", links);
	}
	if(params.sparsity > 0) {
		int k = 0;
		for(i = 0; i < L; i++) {
			for(j = i+1; j < L; j++) {
				for(a = 0; a < q; a++) {
					for(b = 0; b < q; b++) {
						idx[k][0] = i;
						idx[k][1] = j;
						idx[k][2] = a;
						idx[k][3] = b;
						sortedJ[k] = J[i*q+a][j*q+b];
						k += 1;
					}
				}
			}
		}
	}
	model_sp = 0.0;
	return 0;

}

int sample()
{
	int t, i, s;
	int test_chain_eq;
	bool test_chain_half;
	bool test_chain_sampl;
	//int Neff = ceil((1.0*params.Nmc_starts) / params.num_threads);
	idx_chain_1 = (int)rand() % params.Nmc_starts;
	idx_chain_2 = (int)rand() % params.Nmc_starts;
	while(idx_chain_1 == idx_chain_2)
		idx_chain_2 = (int)rand() % params.Nmc_starts;
//#pragma omp parallel for private(i,s)
	for(t = 0; t < params.num_threads; t++) {
		for(s = 0; s < params.Nmc_starts; s++) {
			//int x[L];
			//int *aux;
			test_chain_eq = 0;
			for(i = 0; i < L; i++)
				curr_state[i] = (int)rand() % q;
			if(s == idx_chain_1) {
				test_chain_eq = 1;
				test_chain_sampl = true;
				test_chain_half = true;
			} else {
				test_chain_half = false;
				test_chain_sampl = false;
			}
			if(s == idx_chain_2)
				test_chain_eq = 2;
			mc_chain(test_chain_sampl, test_chain_half, test_chain_eq);
			if(s == idx_chain_1) {
				for(i = 0; i < L; i++)
					chain_fin_1[i] = curr_state[i];
			}
			if(s == idx_chain_2) {
				for(i = 0; i < L; i++)
					chain_fin_2[i] = curr_state[i];
			}
		}
	}
	return 0;
}

int equilibration_test() 
{

	int q_int_chain = 0, q_ext_chain = 0;
	int q_half_first = 0, q_first = 0;
	int i;
	for(i = 0; i < L; i++) {
		if(chain_fin_1[i] == chain_fin_2[i])
			q_ext_chain++;
		if(intra_chain[i] == chain_fin_1[i])
			q_int_chain++;
		if(half_chain[i] == chain_init_1[i])
			q_half_first++;
		if(chain_init_1[i] == chain_init_2[i])
			q_first++;
	}
	printf("q_int_chain: %d q_ext_chain: %d\n", q_int_chain, q_ext_chain);
	printf("q_half_fisrt: %d q_first: %d\n", q_half_first, q_first);
	if(1.0*(q_int_chain - q_ext_chain)/L > 0.05) {
		fprintf(stdout, "Increasing sampl. time: q(halftime config.(A),last config.(A)): %d q(last config.(A), last config.(B)): %d\n", q_int_chain, q_ext_chain);
		params.Twait *= 1.05;
	
	}
	if(1.0*(q_half_first - q_first)/L > 0.05) {
		fprintf(stdout,"Increasing eq. time: q(Teq/2(A), Teq(A)): %d q(Teq(A), Teq(B)): %d\n", q_half_first, q_first);
		params.Teq *= 1.05;
	}
	return 0;
}

int metropolis_step() 
{
	int i, a, j;
	double deltaE, p;
	i = (int)rand() % L;
	a = (int)rand() % q;
	while(a == curr_state[i])
		a = (int)rand() %q;
	deltaE = -h[i*q + a] + h[i*q + curr_state[i]];
	for(j = 0; j < L; j++) if(j != i) {
		deltaE += - J[i*q + a][j*q + curr_state[j]]
			  + J[i*q + curr_state[i]][j*q + curr_state[j]];
	}
	p = rand01();
	if (exp(-deltaE) > p) {
		//printf("i: %d, %d -> %d\n", i, curr_state[i], a);
		curr_state[i] = a;
	}
	return 0;
}

int gibbs_step()
{
	int a, i, j;
	double H[q], p[q+1];
	double r, cum;

	i = (int)rand() % L;
	cum = 0.0;
	for(a = 0; a < q; a++) {
		H[a] = gibbs_en;
		H[a] += -h[i*q + a] + h[i*q + curr_state[i]];
		for(j = 0; j < L; j++) if(j != i) {
			H[a] += - J[i*q + a][j*q + curr_state[j]]
			  + J[i*q + curr_state[i]][j*q + curr_state[j]];
		}
		cum += exp(-H[a]);
		p[a+1] = cum;
		//fprintf(fp, "a: %i en: %f H: %f cum: %f, p[a+1]: %f\n", a, gibbs_en,  H[a], cum, p[a+1]);
	}
	p[0] = 0;
	for(a = 0; a < q + 1; a++) {
		p[a] /= p[q];
		//fprintf(fp, "%f \n", p[a]);
	}
	r = rand01();
	for(a = 0; a < q; a ++) {
		if(p[a] < r && p[a+1] > r) {
			//fprintf(fp, "a: %d, %f < %f < %f\n", a, p[a], r, p[a+1]);
			curr_state[i] = a;
			gibbs_en = H[a];
		}
	}
	return 0;
}

int mc_chain(bool test_chain_sampl, bool test_chain_half, int test_chain_eq)
{
	int t = 0,n,i;
	FILE * fp = 0, * fe = 0;
	if(print_samples)
		fp  = fopen(params.file_samples, "a");
	if(print_en)
		fe = fopen(params.file_en, "a");
	if(params.Gibbs)
		gibbs_en = energy(curr_state);

	while(t <= params.Teq) {
		t++;
		if(test_chain_half && t == params.Teq/2)
			for(i = 0; i < L; i++)
				half_chain[i] = curr_state[i];
		if(params.Metropolis)
			metropolis_step(curr_state);
		else
			gibbs_step(curr_state);
	}
	for(n = 0; n < params.Nmc_config; n++) {
		if(n == 0 && test_chain_eq == 1) {
			for(i = 0; i < L; i++)
				chain_init_1[i] = curr_state[i];
		}
		if(n == 0 && test_chain_eq == 2) {
			for(i = 0; i < L; i++)
				chain_init_2[i] = curr_state[i];
		}
		t = 0;
		while(t <= params.Twait) {
			t++;
			if(params.Metropolis)
				metropolis_step(curr_state);
			else
				gibbs_step(curr_state);
		}
		update_statistics(curr_state);
		if(test_chain_sampl == true && n == params.Nmc_config/2) {
			for(i = 0; i < L; i++)
				intra_chain[i] = curr_state[i];
		}
		if(print_samples) {
			for(i = 0; i < L; i++)
				fprintf(fp, "%d ", curr_state[i]);
			fprintf(fp, "\n");
			fflush(fp);
		}
		if(print_en)
			fprintf(fe, "%lf\n", energy(curr_state));
		if(compute_tm)
			update_tm_statistics(curr_state);
	}
	if(print_samples)
		fclose(fp);
	if(print_en)
		fclose(fe);
	return 0;

}

double randrange(int xmin, int xmax)
{
	double scale = rand() / (double) RAND_MAX;
	return xmin + scale * (xmax - xmin);
}

double rand01()
{
	return (double) rand() / (double)((unsigned)RAND_MAX + 1);
}

double max(double x, double y)
{
	return (x > y) ? x : y;
}

double min(double x, double y)
{
	return (x < y) ? x : y;
}

int update_parameters(double regJ)
{
	int i, j, a, b;
	merrh = 0;
	merrJ = 0;
	averrh = 0;
	averrJ = 0;
	for(i = 0; i < L; i++) {
		for(a = 0; a < q; a++) {
			averrh += fabs(fm_s[i*q +a] - fm[i*q +a]);
			merrh =  max(merrh, fabs(fm_s[i*q+a] - fm[i*q+a]));
			h[i*q+a] += params.lrateh * adaptive_learn(i, a, 0, 0, 'h', fm[i*q+a] - fm_s[i*q+a]);
			for(j = i+1; j < L; j++) {
				for(b = 0; b <q; b++) {
					averrJ += fabs(sm_s[i*q +a][j*q +b] - sm[i*q +a][j*q +b]);
					merrJ =  max(merrJ, fabs(sm_s[i*q+a][j*q+b] - sm[i*q+a][j*q+b]));
					J[i*q + a][j*q + b] += params.lrateJ * adaptive_learn(i, a, j, b, 'J', sm[i*q+a][j*q+b] - sm_s[i*q+a][j*q+b])
							- regJ * ( (J[i*q +a][j*q + b] > 0)  - (J[i*q +a][j*q+b] < 0) ) ;
					J[j*q + b][i*q + a] = J[i*q + a][j*q + b];
				}
			}
		}
	}
	averrh /= L*q;
	averrJ /= (L*(L-1)/2)*q*q;

	return 0;

}

double adaptive_learn(int i, int a, int j, int b, char type, double grad)
{
	if(!params.learn_strat) {
		double tmp = 0;
		switch(type) {
			case 'J':
				tmp = grad * dec[i*q+a][j*q+b];
				break;
			case 'h':
				tmp = grad;
				break;
		}
		return tmp;
	} else if(params.learn_strat == 1) {
		double tmp;
		switch (type) {
			case 'J':
				Gj[i*q + a][j*q + b] += pow(sm_s[i*q + a][j*q + b] - sm[i*q + a][j*q + b], 2);
				tmp = dec[i*q+a][j*q+b] / sqrt(Gj[i*q + a][j*q + b] + 1e-12);
				break;
			case 'h':
				Gh[i*q + a] += pow(fm_s[i*q + a] - fm[i*q +a], 2);
				tmp = 1.0 / sqrt(Gh[i*q + a] + 1e-12);
				break;
			default:
				fprintf(stderr, "I don't know any parameter called %c\n", type);
				return(EXIT_FAILURE);
		}
		return tmp * grad;
	} else if(params.learn_strat == 3) {
		double tmp;
		switch(type) {
			case 'J':
				tmp = (dec[i*q+a][j*q+b] * grad) / (1.0 + iter/params.tau);
				break;
			case 'h':
				tmp = grad / (1.0 + iter/params.tau);
				break;
			default:
				fprintf(stderr, "I don't know any parameter called %c\n", type);
				return(EXIT_FAILURE);
		}
		return tmp;
	}else if(params.learn_strat == 2) {
		double tmp;
		switch (type) {
			case 'J':
				Gj[i*q + a][j*q + b] = params.rho * Gj[i*q + a][j*q + b] + (1.0 - params.rho) * pow(sm_s[i*q + a][j*q + b] - sm[i*q + a][j*q + b], 2);
				tmp = dec[i*q+a][j*q+b] / sqrt(Gj[i*q + a][j*q + b] + 1e-12);
				break;
			case 'h':
				Gh[i*q + a] = params.rho * Gh[i*q + a] + (1.0 - params.rho) * pow(fm_s[i*q + a] - fm[i*q +a], 2);
				tmp = 1.0 / sqrt(Gh[i*q + a] + 1e-12);
				break;
			default:
				fprintf(stderr, "I don't know any parameter called %c\n", type);
				return(EXIT_FAILURE);
		}
		return tmp * grad;
	} else if(params.learn_strat == 4) {
		double tmp, vhat, mhat;
		switch (type) {
			case 'J':
				mtj[i*q + a][j*q + b] *= 0.9;
				mtj[i*q + a][j*q + b] += 0.1*grad;
				vtj[i*q + a][j*q + b] *= 0.999;
				vtj[i*q + a][j*q + b] += 0.001 *grad*grad;
				mhat = mtj[i*q+a][j*q+b] / (1.0 - pow(0.9, iter));
				vhat = vtj[i*q+a][j*q+b] / (1.0 - pow(0.999,iter));
				tmp = mhat / (sqrt(vhat) + 1e-8);
				tmp *= dec[i*q+a][j*q+b];
				break;
			case 'h':
				mth[i*q + a] *= 0.9;
				mth[i*q + a] += 0.1*grad;
				vth[i*q + a] *= 0.999;
				vth[i*q + a] += 0.001 * grad *grad;
				mhat = mth[i*q+a] / (1.0 - pow(0.9, iter));
				vhat = vth[i*q+a] / (1.0 - pow(0.999, iter));
				tmp = mhat / (sqrt(vhat) + 1e-8);
				break;
			default:
				fprintf(stderr, "I don't know any parameter called %c\n", type);
				return(EXIT_FAILURE);
		}
		return tmp;
	} else
		return grad;
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
		if(dec[i*q + a][j*q + b] > 0) {
			m += 1;
			sortedJ[k] = params.empdec ? min(fabs(sm[i*q+a][j*q+b]) + rand01() * params.pseudocount, fabs(J[i*q+a][j*q+b])) : fabs(J[i*q+a][j*q+b]);
		} else
			sortedJ[k] = n * rand01(); // to be optimized: elements should be removed instead of putting large numbers
	}
	fprintf(stdout, "Non-zeros parameters %d / %d \n",m,n);
	for(k = 0; k < n; k++)
		tmp_idx[k] = k;
	quicksort(sortedJ, tmp_idx, 0, n-1);
	for(k = 0; k < c; k++) {
		index = tmp_idx[k];
		i = idx[index][0];
		a = idx[index][2];
		j = idx[index][1];
		b = idx[index][3];
		//fprintf(stderr, "%i %f\n", index, J[idx[index][0]*q+idx[index][2]][idx[index][1]*q + idx[index][3]]);
		J[i*q + a][j*q + b] = 0.0;
		J[j*q + b][i*q + a] = 0.0;
		dec[i*q+a][j*q + b] = 0.0;
		dec[j*q+b][i*q + a] = 0.0;
	}
	model_sp += 1.0*c/n;
	return 0;
}

int quicksort(double *x, int *tmp_idx, int first, int last)
{
	int i, j, pivot, index;
	double temp;
	if(first < last) {
		pivot = first;
		i = first;
		j = last;
		while(i < j) {
			while(x[i] <= x[pivot] && i < last)
				i++;
			while(x[j] > x[pivot])
				j--;
			if(i < j) {
				index = tmp_idx[i];
				tmp_idx[i] = tmp_idx[j];
				tmp_idx[j] = index;
				temp = x[i];
				x[i] = x[j];
				x[j] = temp;
			}
		}

		temp = x[pivot];
		x[pivot] = x[j];
		x[j] = temp;
		index = tmp_idx[pivot];
		tmp_idx[pivot] = tmp_idx[j];
		tmp_idx[j] = index;
		quicksort(x, tmp_idx, first, j-1);
		quicksort(x, tmp_idx, j+1,last);
	}
	return 0;
}

double energy(int *seq)
{

	int i, j;
	double en = 0;

	for(i = 0; i < L; i++) {
		en += -h[i*q + seq[i]];
		for(j = i+1; j < L; j++) 
			en += -J[i*q + seq[i]][j*q +seq[j]];
	}

	return en;
}

void permute(int * vector, int s)
{
	int i,j, aux;
	for(i = 0; i < s; i++) {
		j = i + rand() % (s - i);
		aux = vector[i];
		vector[i] = vector[j];
		vector[j] = aux;
	}
}

int print_model(char *filename)
{
	int i, j , a, b;
	FILE *fp;
	fp = fopen(filename, "w");

	for(i = 0; i < L; i++) {
		for(j = 0; j < L; j++) if(i < j) {
			for(a = 0; a < q; a++) {
				for(b = 0; b < q; b++)
					fprintf(fp, "J %d %d %d %d %f\n",i, j,a,b, J[i*q+a][j*q+b]);
			}
		}
	}
	for(i = 0; i < L; i++) {
		for(a = 0; a < q; a++)
			fprintf(fp, "h %i %i %f\n", i, a, h[i*q+a]);
	}
	fflush(fp);
	fclose(fp);
	return 0;
}

int print_statistics(char *file_sm, char *file_fm, char *file_tm)
{
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
		for(ind = 0; ind < ntm; ind++) {
			i = tm_index[ind][0];
			j = tm_index[ind][1];
			k = tm_index[ind][2];
			a = tm_index[ind][3];
			b = tm_index[ind][4];
			c = tm_index[ind][5];
			fprintf(ft, "%d %d %d %d %d %d %.5f %.5f\n", i, j, k, a, b,c, tm[ind], tm_s[ind]);
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

int print_msa(char *filename)
{
	int i,m;
	FILE *fp;
	fp = fopen(filename, "w");

	for(m = 0; m < M;m++) {
		for(i = 0; i < L; i++)
			fprintf(fp, "%d ", msa[m][i]);
		fprintf(fp, "\n");
	}
	fflush(fp);
	fclose(fp);
	return 0;
}

