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
#include <vector>
#include <climits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cctype>

using namespace std;
typedef double MYFLOAT;

#ifndef BMaux_H
#define BMaux_H

double rand01();
double randrange(int xmin, int xmax);
double max(double x, double y);
double min(double x, double y);
void permute(int *vector, int s);
int overlap(vector<unsigned char> &x1, vector<unsigned char> &x2);
int overlap_two(vector<int> &x1, vector<int> &x2);
int convert_char_amino(char a);
int convert_char_nbase(char a);
int convert_char_epi(char a);
int convert_char_ising(char a);
int convert_char_potts(char a);
int convert_char(char ch, char *ctype);
int print_alphabet(char *ctype);
vector<char> alphabet(char *ctype);

/* quicksort from small to large */
template <class T>
int quicksort(vector<T> &x, vector<int> &tmp_idx, int first, int last)
{
	sort(tmp_idx.begin() + first, tmp_idx.begin() + last + 1, [&](int i, int j)
		 { return x[i] < x[j]; });
	vector<T> y((int)(tmp_idx.size()));
	for (int k = 0; k < int(y.size()); k++)
	{
		y[k] = x[tmp_idx[k]];
	}
	x = y;
	return 0;
}

template <class T>
int quicksort_c(vector<T> &x, vector<int> &tmp_idx, int first, int last)
{
	int i, j, pivot, index;
	double temp;
	if (first < last)
	{
		pivot = first;
		i = first;
		j = last;
		while (i < j)
		{
			while (x[i] <= x[pivot] && i < last)
				i++;
			while (x[j] > x[pivot])
				j--;
			if (i < j)
			{
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
		quicksort_c(x, tmp_idx, first, j - 1);
		quicksort_c(x, tmp_idx, j + 1, last);
	}
	return 0;
}

template <class T>
int print_frobenius_norms(vector<T> &h, vector<vector<T>> &J, int L, int q, char *normfile, char *parfile)
{

	ofstream fp;
	ofstream fpp;
	fp.open(normfile);
	fpp.open(parfile);

	vector<T> hzs = h;
	vector<vector<T>> Jzs = J;

	for (int i = 0; i < L; i++)
	{
		double mean_h = 0.0;
		for (int a = 0; a < q; a++)
			mean_h += h[i * q + a];
		mean_h /= q;
		for (int a = 0; a < q; a++)
			hzs[i * q + a] = h[i * q + a] - mean_h;
		for (int j = 0; j < L; j++)
		{
			double mean_all = 0.0;
			double mean_a[q];
			double mean_b[q];
			for (int a = 0; a < q; a++)
			{
				mean_a[a] = 0.0;
				mean_b[a] = 0.0;
			}
			for (int a = 0; a < q; a++)
			{
				for (int b = 0; b < q; b++)
				{
					mean_a[b] += J[i * q + a][j * q + b];
					mean_b[a] += J[i * q + a][j * q + b];
					mean_all += J[i * q + a][j * q + b];
				}
			}
			for (int a = 0; a < q; a++)
			{
				mean_a[a] /= q;
				mean_b[a] /= q;
			}
			mean_all /= q * q;
			for (int a = 0; a < q; a++)
			{
				for (int b = 0; b < q; b++)
				{
					Jzs[i * q + a][j * q + b] = J[i * q + a][j * q + b] - mean_a[b] - mean_b[a] + mean_all;
				}
			}
			for (int a = 0; a < q; a++)
				hzs[i * q + a] += mean_b[a] - mean_all;
		}
	}

	for (int i = 0; i < L; i++)
	{
		for (int j = i + 1; j < L; j++)
		{
			for (int a = 0; a < q; a++)
			{
				for (int b = 0; b < q; b++)
					fpp << "J " << i << " " << j << " " << a << " " << b << " " << std::fixed << setprecision(5) << Jzs[i * q + a][j * q + b] << endl;
			}
		}
	}
	for (int i = 0; i < L; i++)
		for (int a = 0; a < q; a++)
			fpp << "h " << i << " " << a << " " << std::fixed << setprecision(5) << hzs[i * q + a] << endl;

	int nf = 0;
	double f = 0;
	double fb[L][L];
	double fi[L];
	for (int i = 0; i < L; i++)
	{
		fi[i] = 0;
		for (int j = 0; j < L; j++)
		{
			fb[i][j] = 0;
			if (i != j)
			{
				for (int a = 1; a < q; a++)
				{
					for (int b = 1; b < q; b++)
						fb[i][j] += Jzs[i * q + a][j * q + b] * Jzs[i * q + a][j * q + b];
				}
				fb[i][j] = sqrt(fb[i][j]);
				fi[i] += fb[i][j];
				f += fb[i][j];
				nf++;
			}
		}
		fi[i] /= (L - 1);
	}
	f /= nf;
	for (int i = 0; i < L; i++)
	{
		for (int j = i + 1; j < L; j++)
		{
			fb[i][j] -= fi[i] * fi[j] / f;
			fp << i + 1 << " " << j + 1 << " " << std::fixed << setprecision(5) << fb[i][j] << endl;
		}
	}
	fp.close();
	fpp.close();
	return 0;
}

template <class T>
vector<int> ranks(vector<T> x)
{
	vector<int> idx;
	for (int i = 0; i < int(x.size()); i++)
	{
		idx.push_back(i);
	}
	quicksort(x, idx, 0, int(x.size()) - 1);
	vector<int> ris(x.size());
	for (int i = 0; i < int(x.size()); i++)
	{
		ris[idx[i]] = i;
	}
	return ris;
}

template <class T>
double pearson(vector<T> x, vector<T> y)
{
	double mean_x = 0.0;
	double mean_y = 0.0;
	double mean_prod = 0.0;
	double mean_x2 = 0.0;
	double mean_y2 = 0.0;
	int n;
	if (int(x.size()) == int(y.size()))
	{
		n = int(x.size());
	}
	else
	{
		cout << "Error: trying to compute Pearson between vectors of different size" << endl;
		exit(1);
	}
	for (int i = 0; i < n; i++)
	{
		mean_x += x[i];
		mean_y += y[i];
		mean_prod += x[i] * y[i];
		mean_x2 += x[i] * x[i];
		mean_y2 += y[i] * y[i];
	}
	mean_x /= n;
	mean_y /= n;
	mean_prod /= n;
	mean_y2 /= n;
	mean_x2 /= n;
	double covxy = mean_prod - mean_x * mean_y;
	double std_y = sqrt(mean_y2 - mean_y * mean_y);
	double std_x = sqrt(mean_x2 - mean_x * mean_x);
	double rho = covxy / (std_x * std_y);
	return rho;
}

template <class T>
double spearman(vector<T> x, vector<T> y)
{
	return pearson(ranks(x), ranks(y));
}

#endif
