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
#include <climits>
#include <iostream>
#include <fstream>

using namespace std;

#ifndef BMAUX_H
#define BMAUX_H

double rand01();
double randrange(int xmin, int xmax);
double max(double x, double y);
double min(double x, double y);
int quicksort(vector<double> & x, vector<int> & tmp_idx, int first, int last);
void permute(int * vector, int s); 
int overlap(vector<int> & x1, vector<int> & x2); 
int convert_char_amino(char a);
int convert_char_nbase(char a);
int convert_char_epi(char a);
int convert_char_ising(char a);
int print_alphabet(char * ctype); 
int print_frobenius_norms(vector<double> & h, vector< vector<double> > & J, int L, int q, char *normfile, char *parfile); 

#endif
