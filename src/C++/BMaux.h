#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include <vector>

using namespace std;

#ifndef BMaux
#define BMaux

inline double rand01() {
  return (double) rand() / (double)((unsigned)RAND_MAX + 1);
}

inline double randrange(int xmin, int xmax) {
  return xmin + (xmax - xmin) * rand() / (double) RAND_MAX;
}

inline double max(double x, double y) {
  return (x > y) ? x : y;
}

inline double min(double x, double y) {
  return (x < y) ? x : y;
}


int quicksort(vector<double> & x, vector<int> & tmp_idx, int first, int last) {
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

void permute(int * vector, int s) {
  int i,j, aux;
  for(i = 0; i < s; i++) {
    j = i + rand() % (s - i);
    aux = vector[i];
    vector[i] = vector[j];
    vector[j] = aux;
  }
}

int overlap(vector<int> & x1, vector<int> & x2) {
  int q=0;
  for (int i=0; i<x1.size(); i++) {
    if (x1[i]==x2[i]) q++;
  }
  return q;
}


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



int print_frobenius_norms(vector<double> & h, vector< vector<double> > & J, int L, int q, char *normfile, char *parfile) {

  FILE *fp, *fpp;
  fp = fopen(normfile, "w");
  fpp = fopen(parfile, "w");

  vector<double> hzs=h;
  vector< vector<double> > Jzs=J;

  for(int i = 0; i < L; i++) {
    double mean_h = 0.0;
    for(int a = 0; a < q; a++)
      mean_h += h[i*q+a];
    mean_h /= q;
    for(int a = 0; a < q; a++)
      hzs[i*q+a] = h[i*q+a] - mean_h;
    for(int j = 0; j < L; j++) {
      double mean_all = 0.0;
      double mean_a[q];
      double mean_b[q];
      for(int a = 0; a < q; a++) {
	mean_a[a] = 0.0;
	mean_b[a] = 0.0;
      }
      for(int a = 0; a < q; a++) {
	for(int b = 0; b < q; b++) {
	  mean_a[b] += J[i*q+a][j*q+b];
	  mean_b[a] += J[i*q+a][j*q+b];
	  mean_all += J[i*q+a][j*q+b];
	}
      }
      for(int a = 0; a < q; a++) {
	mean_a[a] /= q;
	mean_b[a] /= q;
      }
      mean_all /= q*q;
      for(int a =0; a < q;a++) {
	for(int b = 0; b <q;b++){
	  Jzs[i*q+a][j*q+b] = J[i*q+a][j*q+b] - mean_a[b] - mean_b[a] + mean_all;
	}
      }
      for(int a = 0; a < q; a++)
	hzs[i*q+a] += mean_b[a] - mean_all;
    }
  }
  
  for(int i = 0; i < L; i++) {
    for(int j = i+1; j < L; j++) {
      for(int a = 0; a < q; a++) {
	for(int b = 0; b < q; b++)
	  fprintf(fpp,"J %d %d %d %d %lf\n",i,j,a,b,Jzs[i*q+a][j*q+b]);
      }
    }
  }
  for(int i = 0; i < L; i++)
    for(int a = 0; a < q; a++)
      fprintf(fpp, "h %d %d %f\n",i,a,hzs[i*q+a]);

  int nf = 0;
  double f = 0;
  double fb[L][L];
  double fi[L];
  for(int i = 0; i < L; i++) {
    fi[i]=0;
    for(int j = 0; j < L; j++) {
      fb[i][j] = 0;
      if(i != j) {
	for(int a = 1; a < q; a++) {
	  for(int b = 1; b < q; b++)
	    fb[i][j] += Jzs[i*q+a][j*q+b] * Jzs[i*q+a][j*q+b];
	}
	fb[i][j] = sqrt(fb[i][j]);
	fi[i] += fb[i][j];
	f += fb[i][j];
	nf++;
      }
    }
    fi[i] /= (L-1);
  }
  f /= nf;
  for(int i = 0; i < L; i++) {
    for(int j = i+1; j < L; j++) {
      fb[i][j] -= fi[i]*fi[j]/f;
      fprintf(fp, "%d, %d, %f\n", i+1,j+1,fb[i][j]);
    }
  }
  fflush(fp);
  fflush(fpp);
  fclose(fp);
  fclose(fpp);
  return 0;
}




#endif
