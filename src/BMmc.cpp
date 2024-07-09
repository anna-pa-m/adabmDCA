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
#include <iomanip>
#include <fstream>
#include <vector>
#include <valarray>
#include "BMlib.h"
#include "BMmc.h"
#include "BMaux.h"
#include <limits.h>
#include <unistd.h>
using namespace std;
typedef double MYFLOAT;

Model::Model(int _q, int _L, Params *_params, Stats *_mstat, vector<vector<unsigned char>> &msa, int _ntm, vector<vector<int>> *_tm_index) : q(_q), L(_L), h(L * q, 0), J(L * q, h), decJ(L * q, vector<unsigned char>(L * q, 0)), tm_index(_tm_index),
                                                                                                                                             Gibbs(false), params(_params), mstat(_mstat), alpha(0.1), acc(1), counter(0), model_sp(0)
{

  init_model_stat(_ntm);
  if (strcmp(params->ctype, "i") == 0)
    init_current_state_ising(msa);
  else
  {
    if (params->restore_flag)
      init_last_chain(params->label);
    else
      init_current_state(msa);
  }
  if (params->learn_strat == 1 || params->learn_strat == 2 || params->learn_strat == 5)
  {
    Gh.clear();
    Gh.resize(L * q, 0);
    GJ.clear();
    GJ.resize(L * q, Gh);
  }
  init_decimation_variables();
}

/******************** METHODS FOR INIT AND OUTPUT ***********************************************************/
void Model::init_model_stat(int _ntm)
{

  mstat->fm_s.clear();
  mstat->fm_s.resize(L * q, 0);
  mstat->sm_s.clear();
  mstat->sm_s.resize(L * q, mstat->fm_s);
  mstat->tm_s.clear();
  mstat->tm_s.resize(_ntm, 0);
  mstat->qs.resize(6);
  mstat->corr.resize((params->Nmc_config - 1) * params->Twait);
  mstat->synth_msa.clear();
  mstat->synth_msa.resize(params->num_threads);
  for (int i = 0; i < params->num_threads; i++)
  {
    mstat->synth_msa[i].clear();
    mstat->synth_msa[i].resize(params->Nmc_config * params->Nmc_starts, vector<unsigned char>(L));
  }
  mstat->curr_state.clear();
  mstat->curr_state.resize(params->num_threads);
  for (int i = 0; i < params->num_threads; i++)
  {
    mstat->curr_state[i].clear();
    mstat->curr_state[i].resize(params->Nmc_starts, vector<unsigned char>(L));
  }
  mstat->qs_t.clear();
  mstat->qs_t.resize(params->num_threads, vector<int>(6));
  mstat->old_state1.clear();
  mstat->old_state2.clear();
  mstat->oldold_state1.clear();
  mstat->oldold_state2.clear();
  mstat->x2i.clear();
  mstat->x1i.clear();
  mstat->x2i.resize(params->num_threads, vector<unsigned char>(L));
  mstat->x1i.resize(params->num_threads, vector<unsigned char>(L));
  mstat->old_state1.resize(params->num_threads, vector<unsigned char>(L));
  mstat->old_state2.resize(params->num_threads, vector<unsigned char>(L));
  mstat->oldold_state1.resize(params->num_threads, vector<unsigned char>(L));
  mstat->oldold_state2.resize(params->num_threads, vector<unsigned char>(L));
}

void Model::init_last_chain(char *label)
{
  for (int t = 0; t < params->num_threads; t++)
  {
    string line;
    ifstream file(params->file_last_chain);
    int s = 0;
    while(getline(file, line))
    {
      vector<int> lineData;
      stringstream lineStream(line);
      int value;
      while(lineStream >> value)
        lineData.push_back(value);
      for(int i = 0; i < L; i++) 
        mstat->curr_state[t][s][i] = lineData[i];
      s++;
    }
  }
}

void Model::init_current_state(vector<vector<unsigned char>> &msa)
{

  for (int t = 0; t < params->num_threads; t++)
  {
    if (!params->initdata)
    {
      vector<unsigned char> tmp(L);
      for (int s = 0; s < params->Nmc_starts; s++)
      {
        for (int i = 0; i < L; i++)
          tmp[i] = (unsigned char)(int)rand() % q;
        mstat->curr_state[t][s] = tmp;
      }
    }
    else
    {
      if (int(msa.size()) == 0)
      {
        cerr << "Empty MSA!" << endl;
        exit(EXIT_FAILURE);
      }
      else
      {
        for (int s = 0; s < params->Nmc_starts; s++)
        {
          int i = (int)rand() % int(msa.size());
          mstat->curr_state[t][s] = msa[i];
        }
      }
    }
  }
}

void Model::init_current_state_ising(vector<vector<unsigned char>> &msa)
{

  for (int t = 0; t < params->num_threads; t++)
  {
    if (!params->initdata)
    {
      vector<unsigned char> tmp(L);
      for (int s = 0; s < params->Nmc_starts; s++)
      {
        for (int i = 0; i < L; i++)
        {
          tmp[i] = (unsigned char)(rand01() > 0.5) ? 1 : 0;
        }
        mstat->curr_state[t][s] = tmp;
      }
    }
    else
    {
      if (int(msa.size()) == 0)
      {
        cerr << "Empty MSA!" << endl;
        exit(EXIT_FAILURE);
      }
      else
      {
        for (int s = 0; s < params->Nmc_starts; s++)
        {
          int i = (int)rand() % int(msa.size());
          mstat->curr_state[t][s] = msa[i];
        }
      }
    }
  }
}

int Model::remove_gauge_freedom(vector<vector<MYFLOAT>> &cov)
{
  vector<MYFLOAT> sorted_matrix(q * q, 0);
  vector<int> mapping(q * q);
  double smalln = min(1e-30, params->pseudocount);
  int idx_aux[q * q][2];
  int neff = 0;
  for (int i = 0; i < L; i++)
  {
    for (int j = i + 1; j < L; j++)
    {
      int k = 0;
      for (int a = 0; a < q; a++)
      {
        for (int b = 0; b < q; b++)
        {
          mapping[k] = k;
          idx_aux[k][0] = a;
          idx_aux[k][1] = b;
          sorted_matrix[k] = fabs(cov[i * q + a][j * q + b]) + smalln * rand01();
          k += 1;
        }
      }
      quicksort(sorted_matrix, mapping, 0, q * q - 1);
      for (k = 0; k < 2 * q - 1; k++)
      {
        int a = idx_aux[mapping[k]][0];
        int b = idx_aux[mapping[k]][1];
        if (decJ[i * q + a][j * q + b] > 0)
          neff++;
        J[i * q + a][j * q + b] = 0.0;
        decJ[i * q + a][j * q + b] = 0;
        J[j * q + b][i * q + a] = 0.0;
        decJ[j * q + b][i * q + a] = 0;
      }
    }
  }
  cout << "Gauge fixing: " << neff << " couplings have been removed" << endl;
  return 0;
}

int Model::initialize_parameters(vector<MYFLOAT> &fm, vector<vector<MYFLOAT>> &cov)
{
  int n = 0;
  /* {READ J,H FROM FILE} OR {SET J,H=0 OR H=IND.SITE.MODEL} */
  if (params->Metropolis && !params->Gibbs)
    Gibbs = false;
  else if (!params->Metropolis && params->Gibbs)
    Gibbs = true;
  else
  {
    cerr << "Conflict in Gibbs-Metropolis initialization" << endl;
    exit(EXIT_FAILURE);
  }
  if (params->file_params)
  {
    cout << "Reading input parameters from " << params->file_params << " with betaJ = " << params->betaJ << " and betaH = " << params->betaH << endl;
    FILE *filep;
    if (!(filep = fopen(params->file_params, "r")))
    {
      cerr << "File " << params->file_params << " not found" << params->file_params << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      int i, j, a, b;
      char buffer[100];
      char c;
      double tmp;
      // Only couplings in the file are set to non-zero, and decJ = 1. Otherwise decJ is assumed = 0
      if (strcmp(params->ctype, "i") == 0)
      {
        for (int i = 0; i < L; i++)
        {
          for (int j = 0; j < L; j++)
          {
            decJ[i][j] = 0;
            decJ[j][i] = 0;
          }
        }
      }
      else
      {
        for (int i = 0; i < L; i++)
        {
          for (int j = 0; j < L; j++)
          {
            for (int a = 0; a < q; a++)
            {
              for (int b = 0; b < q; b++)
              {
                decJ[i * q + a][j * q + b] = 0;
              }
            }
          }
        }
      }
      while (!feof(filep) && fgets(buffer, 100, filep) && sscanf(buffer, "%c ", &c) == 1)
      {
        switch (c)
        {
        case 'J':
          if (strcmp(params->ctype, "i") == 0)
          {
            sscanf(buffer, "J %d %d %lf \n", &i, &j, &tmp);
            J[i][j] = params->betaJ * tmp;
            J[j][i] = params->betaJ * tmp;
            decJ[i][j] = 1;
            decJ[j][i] = 1;
            n++;
          }
          else
          {
            sscanf(buffer, "J %d %d %d %d %lf \n", &i, &j, &a, &b, &tmp);
            J[i * q + a][j * q + b] = params->betaJ * tmp;
            J[j * q + b][i * q + a] = params->betaJ * tmp;
            decJ[i * q + a][j * q + b] = 1;
            decJ[j * q + b][i * q + a] = 1;
            n++;
          }
          break;
        case 'j':
          if (strcmp(params->ctype, "i") == 0)
          {
            sscanf(buffer, "j %d %d %lf \n", &i, &i, &tmp);
            J[i][j] = params->betaJ * tmp;
            J[j][i] = params->betaJ * tmp;
            decJ[i][j] = 1;
            decJ[j][i] = 1;
            n++;
          }
          else
          {
            sscanf(buffer, "j %d %d %d %d %lf \n", &i, &j, &a, &b, &tmp);
            J[i * q + a][j * q + b] = params->betaJ * tmp;
            J[j * q + b][i * q + a] = params->betaJ * tmp;
            decJ[i * q + a][j * q + b] = 1;
            decJ[j * q + b][i * q + a] = 1;
            n++;
          }
          break;
        case 'H':
          if (strcmp(params->ctype, "i") == 0)
          {
            sscanf(buffer, "H %d %lf \n", &i, &tmp);
            h[i] = params->betaH * tmp;
          }
          else
          {
            sscanf(buffer, "H %d %d %lf \n", &i, &a, &tmp);
            h[i * q + a] = params->betaH * tmp;
          }
          break;
        case 'h':
          if (strcmp(params->ctype, "i") == 0)
          {
            sscanf(buffer, "h %d %lf \n", &i, &tmp);
            h[i] = params->betaH * tmp;
          }
          else
          {
            sscanf(buffer, "h %d %d %lf \n", &i, &a, &tmp);
            h[i * q + a] = params->betaH * tmp;
          }
          break;
        }
      }
      fclose(filep);
    }
    cout << "done" << endl;
    double nref = (L * (L - 1) * q * q) / 2;
    model_sp = 1. - n / nref;
    cout << "Sparsity after initialization: " << model_sp << endl;
  }
  else
  {
    for (int i = 0; i < L * q; i++)
    {
      h[i] = 0.0;
      for (int j = i; j < L * q; j++)
      {
        J[i][j] = 0.0;
        J[j][i] = 0.0;
      }
    }
    switch (params->init)
    {
    case 'R':
      cout << "Zero-parameters initialization...done" << endl;
      break;
    case 'I':
      cout.unsetf(ios::scientific);
      cout << setprecision(2);
      cout << "Initializing H parameters using independent sites approximation, J = " << params->initst << " times covariance matrix...";
      double ent = 0;
      for (int i = 0; i < L; i++)
      {
        double mean_all = 0;
        for (int a = 0; a < q; a++)
        {
          h[i * q + a] = log(fm[i * q + a]);
          mean_all += h[i * q + a];
        }
        mean_all /= q;
        double num = 0;
        double den = 0;
        for (int a = 0; a < q; a++)
        {
          h[i * q + a] -= mean_all;
          num += h[i * q + a] * exp(h[i * q + a]);
          den += exp(h[i * q + a]);
        }
        ent += (log(den) - num / den);
        for (int j = i; j < L * q; j++)
        {
          J[i][j] = params->initst * cov[i][j];
          J[j][i] = params->initst * cov[i][j];
        }
      }
      cout.setf(ios::fixed, ios::floatfield);
      cout.precision(2);
      cout << "with Sprof = " << ent << " - done" << endl;
    }
  }
  return 0;
}

void Model::initial_decimation(vector<vector<MYFLOAT>> &cov)
{
  int n = 0;
  if (!params->file_cc)
  {
    if (params->dgap)
    {
      cout << "Using DGap model" << endl;
      n = (L * (L - 1) * (q - 1) * (q - 1)) / 2;
    }
    else if (params->gapnn)
    {
      cout << "Using GapNN model" << endl;
      n = (L * (L - 1) * (q - 1) * (q - 1)) / 2 + L - 1;
    }
    else if (params->phmm)
    {
      cout << "Using Hmmer-like model" << endl;
      n = L - 1;
    }
    else if (params->rmgauge)
    {
      // cout << "Using Potts model with Gauge fixing via cc-decimation" << endl;
      n = (L * (L - 1) * (q - 1) * (q - 1)) / 2;
    }
    else
    {
      cout << "Using full Potts model" << endl;
      n = (L * (L - 1) * q * q) / 2;
    }
    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
      {
        for (int a = 0; a < q; a++)
        {
          for (int b = 0; b < q; b++)
          {
            decJ[i * q + a][j * q + b] = 1;
            if (i == j)
            {
              decJ[i * q + a][j * q + b] = 0;
              J[i * q + a][j * q + b] = 0.0;
            }
            else if (params->dgap && (a == 0 || b == 0))
            {
              decJ[i * q + a][j * q + b] = 0;
              J[i * q + a][j * q + b] = 0.0;
            }
            else if (params->gapnn)
            {
              if (abs(i - j) > 1 && (a == 0 || b == 0))
              {
                decJ[i * q + a][j * q + b] = 0;
                J[i * q + a][j * q + b] = 0.0;
              }
              if (abs(i - j) == 1 && ((a == 0 && b > 0) || (b == 0 && a > 0)))
              {
                decJ[i * q + a][j * q + b] = 0;
                J[i * q + a][j * q + b] = 0.0;
              }
            }
            else if (params->phmm && !(a == 0 && b == 0 && abs(i - j) == 1))
            {
              decJ[i * q + a][j * q + b] = 0;
              J[i * q + a][j * q + b] = 0.0;
            }
          }
        }
      }
    }
    if (params->rmgauge)
    {
      cout << "Using Potts model with Gauge fixing via cc-decimation" << endl;
      n = (L * (L - 1) * (q - 1) * (q - 1)) / 2;
      remove_gauge_freedom(cov);
    }
    else if (params->deczero)
    { // Inactivate all couplings
      n = 0;
      for (int i = 0; i < L * q; i++)
      {
        for (int j = 0; j < L * q; j++)
        {
          if (J[i][j] == 0)
            decJ[i][j] = 0;
          else
            n++;
        }
      }
    }
  }
  else
  {
    cout << "Reading interaction graph from " << params->file_cc << "...";
    fflush(stdout);
    FILE *filep;
    if (!(filep = fopen(params->file_cc, "r")))
    {
      cerr << "File " << params->file_cc << "not found" << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      int i, j, a, b;
      char buffer[100];
      double tmp;
      for (i = 0; i < L * q; i++)
      {
        for (j = 0; j < L * q; j++)
          decJ[i][j] = 0;
      }
      while (!feof(filep) && fgets(buffer, 100, filep) &&
             (sscanf(buffer, "%d %d %d %d %lf \n", &i, &j, &a, &b, &tmp) == 5 || sscanf(buffer, "%d %d %d %d \n", &i, &j, &a, &b) == 4))
      {
        decJ[i * q + a][j * q + b] = 1;
        decJ[j * q + b][i * q + a] = 1;
        n++;
      }
      fclose(filep);
      for (i = 0; i < q * L; i++)
      {
        for (j = 0; j < q * L; j++)
        {
          if (decJ[i][j] == 0)
          {
            J[i][j] = 0.0;
          }
        }
      }
    }
    cout << "done " << endl
         << "Number of links " << n << endl;
  }
  double nref = (L * (L - 1) * q * q) / 2;
  model_sp = 1. - n / nref;
  cout << "Sparsity after initialization: " << model_sp << endl;
}


int Model::print_model(char *filename)
{
  ofstream fp;
  fp.open(filename);
  if (strcmp(params->ctype, "i") == 0)
  {
    for (int i = 0; i < L; i++)
    {
      for (int j = i + 1; j < L; j++)
      {
        if (decJ[i][j] == 1)
          fp << "J " << i << " " << j << " " << std::fixed << setprecision(5) << J[i][j] << endl;
      }
    }
    for (int i = 0; i < L; i++)
      fp << "h " << i << " " << std::fixed << setprecision(5) << h[i] << endl;
  }
  else
  {
    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
        if (i < j)
        {
          for (int a = 0; a < q; a++)
          {
            for (int b = 0; b < q; b++)
            {
              if (decJ[i * q + a][j * q + b] == 1)
                fp << "J " << i << " " << j << " " << a << " " << b << " " << std::fixed << setprecision(5) << J[i * q + a][j * q + b] << endl;
            }
          }
        }
    }
    for (int i = 0; i < L; i++)
    {
      for (int a = 0; a < q; a++)
        fp << "h " << i << " " << a << " " << std::fixed << setprecision(5) << h[i * q + a] << endl;
    }
  }

  fp.close();
  return 0;
}

/******************** METHODS FOR MONTE CARLO ***********************************************************/

double Model::prof_energy(vector<unsigned char> &seq)
{
  double en = 0;
  for (int i = 0; i < L; i++)
  {
    if (strcmp(params->ctype, "i") == 0)
      en += -h[i] * (2.0 * seq[i] - 1.0);
    else
      en += -h[i * q + seq[i]];
  }
  return en;
}

double Model::DCA_energy(vector<unsigned char> &seq)
{
  double en = 0;
  for (int i = 0; i < L; i++)
  {
    for (int j = i + 1; j < L; j++)
    {
      if (strcmp(params->ctype, "i") == 0)
        en += -J[i][j] * (2.0 * seq[i] - 1.0) * (2.0 * seq[j] - 1.0);
      else
        en += -J[i * q + seq[i]][j * q + seq[j]];
    }
  }
  return en;
}

double Model::energy(vector<unsigned char> &seq)
{
  double en = 0;
  for (int i = 0; i < L; i++)
  {
    if (strcmp(params->ctype, "i") == 0)
      en += -h[i] * (2.0 * seq[i] - 1.0);
    else
      en += -h[i * q + seq[i]];
    for (int j = i + 1; j < L; j++)
    {
      if (strcmp(params->ctype, "i") == 0)
        en += -J[i][j] * (2.0 * seq[i] - 1.0) * (2.0 * seq[j] - 1.0);
      else
        en += -J[i * q + seq[i]][j * q + seq[j]];
    }
  }
  return en;
}

vector<MYFLOAT> Model::energy(vector<vector<unsigned char>> &msa)
{
  vector<MYFLOAT> ris;
  for (int m = 0; m < int(msa.size()); m++)
  {
    ris.push_back(energy(msa[m]));
  }
  return ris;
}

void Model::metropolis_step_ising(vector<unsigned char> &x)
{
  int i = (int)rand() % L;
  int a;
  double deltaE;
  a = (x[i] == 0) ? 1 : 0;
  deltaE = -h[i] * (2.0 * a - 1.0) + h[i] * (2.0 * x[i] - 1.0);
  for (int j = 0; j < L; j++)
    deltaE += -J[i][j] * (2.0 * a - 1.0) * (2.0 * x[j] - 1.0) + J[i][j] * (2.0 * x[i] - 1.0) * (2.0 * x[j] - 1.0);
  double p = rand01();
  if (exp(-deltaE) > p)
    x[i] = a;
}

void Model::metropolis_step(vector<unsigned char> &x)
{
  int i = (int)rand() % L;
  int a;
  double deltaE;
  a = (int)rand() % q;
  while (a == x[i])
    a = (int)rand() % q;
  deltaE = -h[i * q + a] + h[i * q + x[i]];
  for (int j = 0; j < L; j++)
    if (j != i)
      deltaE += -J[i * q + a][j * q + x[j]] + J[i * q + x[i]][j * q + x[j]];
  double p = rand01();
  if (exp(-deltaE) > p)
    x[i] = a;
}

void Model::gibbs_step_ising(vector<unsigned char> &x)
{
  double H, cum[q];
  int i = (int)rand() % L;
  for (int a = 0; a < q; a++)
  {
    H = -h[i] * (2.0 * a - 1.0);
    for (int j = 0; j < L; j++)
      if (i != j)
        H += -J[i][j] * (2.0 * a - 1.0) * (2.0 * x[j] - 1.0);
    if (a == 0)
      cum[a] = exp(-H);
    else
      cum[a] = cum[a - 1] + exp(-H);
  }
  double r = rand01();
  if (r > cum[0])
    x[i] = 1;
  else
    x[i] = 0;
}
void Model::gibbs_step(vector<unsigned char> &x)
{
  double H, cum[q];
  int i = (int)rand() % L;
  for (int a = 0; a < q; a++)
  {
    H = -h[i * q + a];
    for (int j = 0; j < L; j++)
      if (j != i)
        H += -J[i * q + a][j * q + x[j]];
    if (a == 0)
      cum[a] = exp(-H);
    else
      cum[a] = cum[a - 1] + exp(-H);
  }
  double r = cum[q - 1] * rand01();
  int a = 0;
  while (r > cum[a])
    a++;
  x[i] = a;
}

void Model::MC_sweep(vector<unsigned char> &x)
{
  if (!Gibbs)
  {
    for (int i = 0; i < L; i++)
      metropolis_step(x);
  }
  else
  {
    for (int i = 0; i < L; i++)
      gibbs_step(x);
  }
}

void Model::MC_sweep_ising(vector<unsigned char> &x)
{
  if (!Gibbs)
  {
    for (int i = 0; i < L; i++)
      metropolis_step_ising(x);
  }
  else
  {
    for (int i = 0; i < L; i++)
      gibbs_step_ising(x);
  }
}

void Model::mc_chain_ising(vector<unsigned char> &x1, vector<unsigned char> &x2, int s)
{

  int numt = 0;
  for (int t = 0; t < params->Teq; t++)
  {
    MC_sweep_ising(x1);
    MC_sweep_ising(x2);
  }
  // valarray<int> qs(6);
  // vector<int> old_state1, old_state2, oldold_state1, oldold_state2;
  // update_statistics_lock(x1,fp,fe);
  // update_statistics_lock(x2,fp,fe);
  // update_synth_msa(x1, x2);
  mstat->synth_msa[numt][(2 * s) * params->Nmc_config] = x1;
  mstat->synth_msa[numt][(2 * s + 1) * params->Nmc_config] = x2;

  // if(int(mstat->tm_s.size())>0) {
  //   update_tm_statistics(x1);
  //   update_tm_statistics(x2);
  // }
  double o12 = overlap(x1, x2);
  mstat->qs_t[numt][0] += o12;
  mstat->qs_t[numt][1] += (o12 * o12);
  mstat->x1i[numt] = x1;
  mstat->x2i[numt] = x2;
  for (int n = 0; n < params->Nmc_config - 1; n++)
  {
    mstat->oldold_state1[numt] = mstat->old_state1[numt];
    mstat->oldold_state2[numt] = mstat->old_state2[numt];
    mstat->old_state1[numt] = x1;
    mstat->old_state2[numt] = x2;
    for (int t = 0; t < params->Twait; t++)
    {
      update_corr(n * params->Twait + t, overlap(x1, mstat->x1i[numt]) + overlap(x2, mstat->x2i[numt]));
      MC_sweep_ising(x1);
      MC_sweep_ising(x2);
    }
    // update_statistics_lock(x1,fp,fe);
    // update_statistics_lock(x2,fp,fe);
    // update_synth_msa(x1, x2);

    mstat->synth_msa[numt][(2 * s) * params->Nmc_config + n + 1] = x1;
    mstat->synth_msa[numt][(2 * s + 1) * params->Nmc_config + n + 1] = x2;
    // if(int(mstat->tm_s.size())>0) {
    //   update_tm_statistics(x1);
    //   update_tm_statistics(x2);
    // }
    o12 = overlap(x1, x2);
    mstat->qs_t[numt][0] += o12;
    mstat->qs_t[numt][1] += (o12 * o12);
    double o1 = overlap(mstat->old_state1[numt], x1);
    double o2 = overlap(mstat->old_state2[numt], x2);
    mstat->qs_t[numt][2] += (o1 + o2);
    mstat->qs_t[numt][3] += (o1 * o1 + o2 * o2);
    if (n > 0)
    {
      double oo1 = overlap(mstat->oldold_state1[numt], x1);
      double oo2 = overlap(mstat->oldold_state2[numt], x2);
      mstat->qs_t[numt][4] += (oo1 + oo2);
      mstat->qs_t[numt][5] += (oo1 * oo1 + oo2 * oo2);
    }
  }
}

void Model::mc_chain(vector<unsigned char> &x1, vector<unsigned char> &x2, int s)
{
  int numt = 0;
  for (int t = 0; t < params->Teq; t++)
  {
    MC_sweep(x1);
    MC_sweep(x2);
  }
  // valarray<int> qs(6);
  // vector<int> old_state1, old_state2, oldold_state1, oldold_state2;
  // update_statistics_lock(x1,fp,fe);
  // update_statistics_lock(x2,fp,fe);
  // update_synth_msa(x1, x2);
  mstat->synth_msa[numt][(2 * s) * params->Nmc_config] = x1;
  mstat->synth_msa[numt][(2 * s + 1) * params->Nmc_config] = x2;
  // if(int(mstat->tm_s.size())>0) {
  //   update_tm_statistics(x1);
  //   update_tm_statistics(x2);
  // }
  double o12 = overlap(x1, x2);
  mstat->qs_t[numt][0] += o12;
  mstat->qs_t[numt][1] += (o12 * o12);
  mstat->x1i[numt] = x1;
  mstat->x2i[numt] = x2;
  for (int n = 0; n < params->Nmc_config - 1; n++)
  {
    mstat->oldold_state1[numt] = mstat->old_state1[numt];
    mstat->oldold_state2[numt] = mstat->old_state2[numt];
    mstat->old_state1[numt] = x1;
    mstat->old_state2[numt] = x2;
    for (int t = 0; t < params->Twait; t++)
    {
      update_corr(n * params->Twait + t, overlap(x1, mstat->x1i[numt]) + overlap(x2, mstat->x2i[numt]));
      MC_sweep(x1);
      MC_sweep(x2);
    }
    // update_statistics_lock(x1,fp,fe);
    // update_statistics_lock(x2,fp,fe);
    // update_synth_msa(x1, x2);
    mstat->synth_msa[numt][(2 * s) * params->Nmc_config + n + 1] = x1;
    mstat->synth_msa[numt][(2 * s + 1) * params->Nmc_config + n + 1] = x2;
    // if(int(mstat->tm_s.size())>0) {
    //   update_tm_statistics(x1);
    //   update_tm_statistics(x2);
    // }
    o12 = overlap(x1, x2);
    mstat->qs_t[numt][0] += o12;
    mstat->qs_t[numt][1] += (o12 * o12);
    double o1 = overlap(mstat->old_state1[numt], x1);
    double o2 = overlap(mstat->old_state2[numt], x2);
    mstat->qs_t[numt][2] += (o1 + o2);
    mstat->qs_t[numt][3] += (o1 * o1 + o2 * o2);
    if (n > 0)
    {
      double oo1 = overlap(mstat->oldold_state1[numt], x1);
      double oo2 = overlap(mstat->oldold_state2[numt], x2);
      mstat->qs_t[numt][4] += (oo1 + oo2);
      mstat->qs_t[numt][5] += (oo1 * oo1 + oo2 * oo2);
    }
  }
}

void Model::update_corr(int i, int value)
{

  mstat->corr[i] += value;
}

bool Model::sample_ising(vector<vector<unsigned char>> &msa)
{

  bool eqmc = true;
  init_statistics();
  mstat->corr.resize((params->Nmc_config - 1) * params->Twait);
  for (int i = 0; i < int(mstat->corr.size()); i++)
    mstat->corr[i] = 0.0;
  for (int i = 0; i < int(mstat->qs.size()); i++)
    mstat->qs[i] = 0;
  if (!params->persistent)
    init_current_state_ising(msa);
  int s;
  for (int t = 0; t < params->num_threads; t++)
  {
    for (s = 0; s < params->Nmc_starts / 2; s++)
    {
      mc_chain_ising(mstat->curr_state[t][2 * s], mstat->curr_state[t][2 * s + 1], s);
    }
  }
  update_statistics_ising();
  mstat->corr /= params->Nmc_starts;
  double nse = params->Nmc_config * (params->Nmc_starts / 2);
  double qext = mstat->qs[0] / nse;
  double dqext = sqrt(mstat->qs[1] / (nse - 1) - mstat->qs[0] * mstat->qs[0] / nse / (nse - 1)) / sqrt(nse);
  double nsi1 = (params->Nmc_config - 1) * params->Nmc_starts;
  double qin1 = nsi1 > 0 ? mstat->qs[2] / nsi1 : 0;
  double dqin1 = nsi1 > 1 ? sqrt(mstat->qs[3] / (nsi1 - 1) - mstat->qs[2] * mstat->qs[2] / nsi1 / (nsi1 - 1)) / sqrt(nsi1) : 0;
  double nsi2 = (params->Nmc_config - 2) * params->Nmc_starts;
  double qin2 = nsi2 > 0 ? mstat->qs[4] / nsi2 : 0;
  double dqin2 = nsi2 > 1 ? sqrt(mstat->qs[5] / (nsi2 - 1) - mstat->qs[4] * mstat->qs[4] / nsi2 / (nsi2 - 1)) / sqrt(nsi2) : 0;
  int test1 = (abs(qext - qin1) < 5 * sqrt(dqext * dqext + dqin1 * dqin1) ? 1 : 0);
  int test2 = (abs(qext - qin2) < 5 * sqrt(dqext * dqext + dqin2 * dqin2) ? 1 : 0);
  if (params->adapt)
  {
    if (test1)
    {
      if (params->Twait > 1)
        // params->Twait-=1;
          params->Twait = floor((params->Twait + params->Twait_last) * 0.5);
      params->Teq = 2 * params->Twait;
    }
    else if (!test2)
    {
      eqmc = false;
      params->Twait *= 2;
      params->Teq = 2 * params->Twait;
    }
  }

  if (params->nprinteq)
  {
    cout.unsetf(ios::scientific);
    cout.setf(ios::fixed, ios::floatfield);
    cout << setprecision(2);
    cout << "Sampling info: q_ext: " << qext << " +- " << dqext << " q_int_1: " << qin1 << " +- " << dqin1 << " q_int_2: " << qin2 << " +- " << dqin2 << " Test_eq1: " << test1 << " Test_eq2: " << test2 << endl;
  }
  return eqmc;
}

bool Model::sample(vector<vector<unsigned char>> &msa)
{

  bool eqmc = true;
  init_statistics();
  mstat->corr.resize((params->Nmc_config - 1) * params->Twait);
  for (int i = 0; i < int(mstat->corr.size()); i++)
    mstat->corr[i] = 0.0;
  for (int i = 0; i < int(mstat->qs.size()); i++)
    mstat->qs[i] = 0;
  if (!params->persistent)
    init_current_state(msa);
  int s;
  for (int t = 0; t < params->num_threads; t++)
  {
    for (s = 0; s < params->Nmc_starts / 2; s++)
    {
      mc_chain(mstat->curr_state[t][2 * s], mstat->curr_state[t][2 * s + 1], s);
    }
  }
  update_statistics();
  mstat->corr /= params->Nmc_starts;
  double nse = params->Nmc_config * (params->Nmc_starts / 2);
  double qext = mstat->qs[0] / nse;
  double dqext = sqrt(mstat->qs[1] / (nse - 1) - mstat->qs[0] * mstat->qs[0] / nse / (nse - 1)) / sqrt(nse);
  double nsi1 = (params->Nmc_config - 1) * params->Nmc_starts;
  double qin1 = nsi1 > 0 ? mstat->qs[2] / nsi1 : 0;
  double dqin1 = nsi1 > 1 ? sqrt(mstat->qs[3] / (nsi1 - 1) - mstat->qs[2] * mstat->qs[2] / nsi1 / (nsi1 - 1)) / sqrt(nsi1) : 0;
  double nsi2 = (params->Nmc_config - 2) * params->Nmc_starts;
  double qin2 = nsi2 > 0 ? mstat->qs[4] / nsi2 : 0;
  double dqin2 = nsi2 > 1 ? sqrt(mstat->qs[5] / (nsi2 - 1) - mstat->qs[4] * mstat->qs[4] / nsi2 / (nsi2 - 1)) / sqrt(nsi2) : 0;
  int test1 = (abs(qext - qin1) < 5 * sqrt(dqext * dqext + dqin1 * dqin1) ? 1 : 0);
  int test2 = (abs(qext - qin2) < 5 * sqrt(dqext * dqext + dqin2 * dqin2) ? 1 : 0);
  if (params->adapt)
  {
    if (test1)
    {
      if (params->Twait > 1)
        // params->Twait-=1;
        params->Twait = floor((params->Twait + params->Twait_last) * 0.5);
      params->Teq = 2 * params->Twait;
    }
    else if (!test2)
    {
      eqmc = false;
      // params->Twait+=1;
      params->Twait_last = params->Twait;
      params->Twait *= 2;
      params->Teq = 2 * params->Twait;
    }
  }

  if (params->nprinteq)
  {
    cout.unsetf(ios::scientific);
    cout << setprecision(2);
    cout << "Sampling info: q_ext: " << qext << " +- " << dqext << " q_int_1: " << qin1 << " +- " << dqin1 << " q_int_2: " << qin2 << " +- " << dqin2 << " Test_eq1: " << test1 << " Test_eq2: " << test2 << endl;
  }
  return eqmc;
}

/******************** METHODS FOR STATISTICS ***********************************************************/

void Model::init_statistics()
{
  for (int i = 0; i < L * q; i++)
  {
    mstat->fm_s[i] = 0;
    for (int j = 0; j < L * q; j++)
    {
      mstat->sm_s[i][j] = 0;
    }
  }
  for (int ind = 0; ind < int(mstat->tm_s.size()); ind++)
    mstat->tm_s[ind] = 0;
  for (int i = 0; i < params->num_threads; i++)
  {
    for (int j = 0; j < int(mstat->qs_t[i].size()); j++)
      mstat->qs_t[i][j] = 0;
  }
}

void Model::update_statistics_ising()
{

  int i, j, k;
  int Ns = params->Nmc_starts * params->Nmc_config * params->num_threads;
  vector<unsigned char> x;

  for (int t = 0; t < params->num_threads; t++)
  {
    for (int idx = 0; idx < int(mstat->qs.size()); idx++)
      mstat->qs[idx] += mstat->qs_t[t][idx];
    for (int m = 0; m < int(mstat->synth_msa[t].size()); m++)
    {
      x = mstat->synth_msa[t][m];
      for (i = 0; i < L; i++)
      {
        mstat->fm_s[i] += 1.0 / Ns * (2.0 * x[i] - 1.0);
        mstat->sm_s[i][i] += 1.0 / Ns;
        for (j = i + 1; j < L; j++)
        {
          mstat->sm_s[i][j] += 1.0 / Ns * (2.0 * x[i] - 1.0) * (2.0 * x[j] - 1.0);
          mstat->sm_s[j][i] += 1.0 / Ns * (2.0 * x[j] - 1.0) * (2.0 * x[i] - 1.0);
        }
      }
      for (int ind = 0; ind < int((*tm_index).size()); ind++)
      {
        i = (*tm_index)[ind][0];
        j = (*tm_index)[ind][1];
        k = (*tm_index)[ind][2];
        mstat->tm_s[ind] += 1.0 / Ns * (2.0 * x[i] - 1.0) * (2.0 * x[j] - 1.0) * (2.0 * x[k] - 1.0);
      }
    }
  }
}
// if(params->file_en)
//  fclose(fe);

void Model::print_last_chain(char * filename)
{

  FILE *fp = 0;
  fp = fopen(filename, "w");
  for (int t = 0; t < params->num_threads; t++)
  {
    for (int s = 0; s < int(mstat->curr_state[t].size()); s++) {
      for (int i = 0; i < L; i++)
        fprintf(fp, "%d ", mstat->curr_state[t][s][i]);
    fprintf(fp, "\n");
    }
  }
  fclose(fp);
}


void Model::print_samples(char *filename)
{

  FILE *fp = 0;
  fp = fopen(filename, "w");
  vector<unsigned char> x;
  vector<char> abc = alphabet(params->ctype);
  for (int t = 0; t < params->num_threads; t++)
  {
    for (int m = 0; m < int(mstat->synth_msa[t].size()); m++)
    {
      x = mstat->synth_msa[t][m];
      fprintf(fp, ">THREAD%d_CHAIN%d_SAMPLE%d h_en %lf J_en %lf \n", t,
              m / params->Nmc_config, m % params->Nmc_config, prof_energy(x), DCA_energy(x));
      for (int i = 0; i < L; i++)
        fprintf(fp, "%c", abc[x[i]]);
      fprintf(fp, "\n");
    }
  }

  fclose(fp);
}

void Model::print_samples_ising(char *filename)
{

  FILE *fp = 0;
  fp = fopen(filename, "w");
  vector<unsigned char> x;

  for (int t = 0; t < params->num_threads; t++)
  {
    for (int m = 0; m < int(mstat->synth_msa[t].size()); m++)
    {
      x = mstat->synth_msa[t][m];
      fprintf(fp, ">THREAD%d_CHAIN%d_SAMPLE%d h_en %lf J_en %lf \n", t,
              m / params->Nmc_config, m % params->Nmc_config, prof_energy(x), DCA_energy(x));
      for (int i = 0; i < L; i++)
        fprintf(fp, "%d", x[i]);
      fprintf(fp, "\n");
    }
  }

  fclose(fp);
}

int Model::print_natural_samples(char *filename, vector<vector<unsigned char>> &msa)
{
  FILE *fp;
  fp = fopen(filename, "w");
  vector<char> ab = alphabet(params->ctype);
  for (int m = 0; m < int(msa.size()); m++)
  {
    fprintf(fp, "> CHAIN %d SAMPLE %d h_ENERGY %lf J_ENERGY %lf \n", 0, m, prof_energy(msa[m]), DCA_energy(msa[m]));
    for (int i = 0; i < L; i++)
      fprintf(fp, "%c", ab[msa[m][i]]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

void Model::update_statistics()
{

  int i, j, k, a, b, c;
  int Ns = params->Nmc_starts * params->Nmc_config * params->num_threads;
  vector<unsigned char> x;

  for (int t = 0; t < params->num_threads; t++)
  {
    for (int idx = 0; idx < int(mstat->qs.size()); idx++)
      mstat->qs[idx] += mstat->qs_t[t][idx];
    for (int m = 0; m < int(mstat->synth_msa[t].size()); m++)
    {
      x = mstat->synth_msa[t][m];
      for (i = 0; i < L; i++)
      {
        mstat->fm_s[i * q + x[i]] += 1.0 / Ns;
        mstat->sm_s[i * q + x[i]][i * q + x[i]] += 1.0 / Ns;
        for (j = i + 1; j < L; j++)
        {
          mstat->sm_s[i * q + x[i]][j * q + x[j]] += 1.0 / Ns;
          mstat->sm_s[j * q + x[j]][i * q + x[i]] += 1.0 / Ns;
        }
      }
      for (int ind = 0; ind < int((*tm_index).size()); ind++)
      {
        i = (*tm_index)[ind][0];
        j = (*tm_index)[ind][1];
        k = (*tm_index)[ind][2];
        a = (*tm_index)[ind][3];
        b = (*tm_index)[ind][4];
        c = (*tm_index)[ind][5];
        if (x[i] == a && x[j] == b && x[k] == c)
          mstat->tm_s[ind] += 1.0 / Ns;
      }
    }
  }
}

void Model::update_statistics_lock(vector<unsigned char> &x, FILE *fp, FILE *fe)
{
  int Ns = params->Nmc_starts * params->Nmc_config;
  if (strcmp(params->ctype, "i") == 0)
  {
    for (int i = 0; i < L; i++)
    {
      mstat->fm_s[i] += 1.0 / Ns * (2.0 * x[i] - 1.0);
      mstat->sm_s[i][i] += 1.0 / Ns;
      for (int j = i + 1; j < L; j++)
      {
        mstat->sm_s[i][j] += 1.0 / Ns * (2.0 * x[i] - 1.0) * (2.0 * x[j] - 1.0);
        mstat->sm_s[j][i] += 1.0 / Ns * (2.0 * x[j] - 1.0) * (2.0 * x[i] - 1.0);
      }
    }
  }
  else
  {
    for (int i = 0; i < L; i++)
    {
      mstat->fm_s[i * q + x[i]] += 1.0 / Ns;
      mstat->sm_s[i * q + x[i]][i * q + x[i]] += 1.0 / Ns;
      for (int j = i + 1; j < L; j++)
      {
        mstat->sm_s[i * q + x[i]][j * q + x[j]] += 1.0 / Ns;
        mstat->sm_s[j * q + x[j]][i * q + x[i]] += 1.0 / Ns;
      }
    }
  }
}

void Model::update_tm_statistics(vector<unsigned char> &x)
{
  int i, j, k, a, b, c;
  int Ns = params->Nmc_starts * params->Nmc_config;
  for (int ind = 0; ind < int((*tm_index).size()); ind++)
  {
    if (strcmp(params->ctype, "i") == 0)
    {
      i = (*tm_index)[ind][0];
      j = (*tm_index)[ind][1];
      k = (*tm_index)[ind][2];
      mstat->tm_s[ind] += 1.0 / Ns * (2.0 * x[i] - 1.0) * (2.0 * x[j] - 1.0) * (2.0 * x[k] - 1.0);
    }
    else
    {
      i = (*tm_index)[ind][0];
      j = (*tm_index)[ind][1];
      k = (*tm_index)[ind][2];
      a = (*tm_index)[ind][3];
      b = (*tm_index)[ind][4];
      c = (*tm_index)[ind][5];
      if (x[i] == a && x[j] == b && x[k] == c)
        mstat->tm_s[ind] += 1.0 / Ns;
    }
  }
}

void Model::compute_third_order_correlations()
{
  int ind, i, j, k, a, b, c;
  if (strcmp(params->ctype, "i") == 0)
  {
    for (ind = 0; ind < int((*tm_index).size()); ind++)
    {
      i = (*tm_index)[ind][0];
      j = (*tm_index)[ind][1];
      k = (*tm_index)[ind][2];
      mstat->tm_s[ind] = mstat->tm_s[ind] - mstat->sm_s[i][j] * mstat->fm_s[k] -
                         mstat->sm_s[i][k] * mstat->fm_s[j] - mstat->sm_s[j][k] * mstat->fm_s[i] +
                         2 * mstat->fm_s[i] * mstat->fm_s[j] * mstat->fm_s[k];
    }
  }
  else
  {
    for (ind = 0; ind < int((*tm_index).size()); ind++)
    {
      i = (*tm_index)[ind][0];
      j = (*tm_index)[ind][1];
      k = (*tm_index)[ind][2];
      a = (*tm_index)[ind][3];
      b = (*tm_index)[ind][4];
      c = (*tm_index)[ind][5];
      mstat->tm_s[ind] = mstat->tm_s[ind] - mstat->sm_s[i * q + a][j * q + b] * mstat->fm_s[k * q + c] -
                         mstat->sm_s[i * q + a][k * q + c] * mstat->fm_s[j * q + b] - mstat->sm_s[j * q + b][k * q + c] * mstat->fm_s[i * q + a] +
                         2 * mstat->fm_s[i * q + a] * mstat->fm_s[j * q + b] * mstat->fm_s[k * q + c];
    }
  }
}

int Model::compute_errors(vector<MYFLOAT> &fm, vector<vector<MYFLOAT>> &sm, vector<vector<MYFLOAT>> &cov, Errs &errs)
{
  errs.errnorm = 0;
  errs.merrh = 0;
  errs.merrJ = 0;
  errs.averrh = 0;
  errs.averrJ = 0;
  for (int i = 0; i < L; i++)
  {
    for (int a = 0; a < q; a++)
    {
      errs.averrh += fabs(mstat->fm_s[i * q + a] - fm[i * q + a]);
      errs.merrh = max(errs.merrh, fabs(mstat->fm_s[i * q + a] - fm[i * q + a]));
      for (int j = i + 1; j < L; j++)
      {
        for (int b = 0; b < q; b++)
        {
          errs.errnorm = max(errs.errnorm, decJ[i * q + a][j * q + b] * fabs(cov[i * q + a][j * q + b] - mstat->sm_s[i * q + a][j * q + b] + mstat->fm_s[i * q + a] * mstat->fm_s[j * q + b]));
          errs.averrJ += fabs(mstat->sm_s[i * q + a][j * q + b] - sm[i * q + a][j * q + b]);
          errs.merrJ = max(errs.merrJ, fabs(mstat->sm_s[i * q + a][j * q + b] - sm[i * q + a][j * q + b]));
        }
      }
    }
  }
  errs.averrh /= L * q;
  errs.averrJ /= (L * (L - 1) / 2) * q * q;
  return 0;
}

double Model::pearson(vector<vector<MYFLOAT>> &cov, bool nodec = false)
{
  double mean_cov_s = 0.0;
  double mean_cov = 0.0;
  double mean_prod = 0.0;
  double mean_x2 = 0.0;
  double mean_y2 = 0.0;
  int n = 0;
  for (int i = 0; i < L; i++)
  {
    for (int j = i + 1; j < L; j++)
    {
      for (int a = 0; a < q; a++)
      {
        for (int b = 0; b < q; b++)
        {
          if (nodec || decJ[i * q + a][j * q + b])
          {
            n += 1;
            double cov_s = mstat->sm_s[i * q + a][j * q + b] - mstat->fm_s[i * q + a] * mstat->fm_s[j * q + b];
            mean_cov_s += cov_s;
            mean_cov += cov[i * q + a][j * q + b];
            mean_prod += cov[i * q + a][j * q + b] * cov_s;
            mean_x2 += cov[i * q + a][j * q + b] * cov[i * q + a][j * q + b];
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
  mean_x2 /= n;
  double covxy = mean_prod - mean_cov_s * mean_cov;
  double std_cov_s = sqrt(mean_y2 - mean_cov_s * mean_cov_s);
  double std_cov = sqrt(mean_x2 - mean_cov * mean_cov);
  double rho = covxy / (std_cov_s * std_cov);
  cout.unsetf(ios::scientific);
  cout << setprecision(1);
  return rho;
}

/******************** METHODS FOR LEARNING ***********************************************************/

double Model::update_parameters(vector<MYFLOAT> &fm, vector<vector<MYFLOAT>> &sm, int iter, Data_e &data_e)
{
  if (params->file_msa_e)
  {
    data_e.energy = energy(data_e.msa);
    data_e.compute_grad();
  }
  if (params->learn_strat != 5)
  {
    double lrav = 0;
    int n = 0;
    for (int i = 0; i < L; i++)
    {
      for (int a = 0; a < q; a++)
      {
        double gradh = fm[i * q + a] - mstat->fm_s[i * q + a] + data_e.gradh[i * q + a];
        double lrh = params->lrateh;
        if (params->learn_strat == 1)
        {
          Gh[i * q + a] += gradh * gradh;
          lrh /= sqrt(Gh[i * q + a] + 1e-12);
        }
        else if (params->learn_strat == 2)
        {
          Gh[i * q + a] = params->rho * Gh[i * q + a] + (1.0 - params->rho) * gradh * gradh;
          lrh /= sqrt(Gh[i * q + a] + 1e-12);
        }
        else if (params->learn_strat == 3)
        {
          lrh /= (1.0 + (double)iter / params->tau);
        }
        else if (params->learn_strat == 4)
        {
          lrh /= fm[i * q + a] * (1.0 - fm[i * q + a]);
        }
        h[i * q + a] += lrh * gradh;
        lrav += lrh;
        n++;
        for (int j = i + 1; j < L; j++)
        {
          for (int b = 0; b < q; b++)
          {
            double gradJ = sm[i * q + a][j * q + b] - mstat->sm_s[i * q + a][j * q + b] + data_e.gradJ[i * q + a][j * q + b] -
                           params->regJ1 * ((J[i * q + a][j * q + b] > 0) - (J[i * q + a][j * q + b] < 0)) - params->regJ2 * J[i * q + a][j * q + b];
            double lrJ = params->lrateJ;
            if (params->learn_strat == 1)
            {
              GJ[i * q + a][j * q + b] += gradJ * gradJ;
              lrJ /= sqrt(GJ[i * q + a][j * q + b] + 1e-12);
            }
            else if (params->learn_strat == 2)
            {
              GJ[i * q + a][j * q + b] = params->rho * GJ[i * q + a][j * q + b] + (1.0 - params->rho) * gradJ * gradJ;
              lrJ /= sqrt(GJ[i * q + a][j * q + b] + 1e-12);
            }
            else if (params->learn_strat == 3)
              lrJ /= (1.0 + (double)iter / params->tau);
            else if (params->learn_strat == 4)
            {
              lrJ /= sm[i * q + a][j * q + b] * (1.0 - sm[i * q + a][j * q + b]);
              // lrJ /= (params->maxiter - iter);
            }
            J[i * q + a][j * q + b] += lrJ * decJ[i * q + a][j * q + b] * gradJ;
            J[j * q + b][i * q + a] = J[i * q + a][j * q + b];
            lrav += lrJ;
            n++;
          }
        }
      }
    }
    return lrav / n;
  }
  else
  {
    double P = 0;
    double modF = 0;
    double modv = 0;
    for (int i = 0; i < L; i++)
    {
      for (int a = 0; a < q; a++)
      {
        double gradh = fm[i * q + a] - mstat->fm_s[i * q + a] + data_e.gradh[i * q + a];
        h[i * q + a] += params->lrateh * acc * Gh[i * q + a];
        Gh[i * q + a] += params->lrateh * acc * gradh;
        P += gradh * Gh[i * q + a];
        modF += gradh * gradh;
        modv += Gh[i * q + a] * Gh[i * q + a];
        for (int j = i + 1; j < L; j++)
        {
          for (int b = 0; b < q; b++)
          {
            double gradJ = sm[i * q + a][j * q + b] - mstat->sm_s[i * q + a][j * q + b] + data_e.gradJ[i * q + a][j * q + b] - params->regJ1 * ((J[i * q + a][j * q + b] > 0) - (J[i * q + a][j * q + b] < 0)) - params->regJ2 * J[i * q + a][j * q + b];
            J[i * q + a][j * q + b] += params->lrateJ * acc * decJ[i * q + a][j * q + b] * GJ[i * q + a][j * q + b];
            J[j * q + b][i * q + a] = J[i * q + a][j * q + b];
            GJ[i * q + a][j * q + b] += params->lrateJ * acc * decJ[i * q + a][j * q + b] * gradJ;
            GJ[j * q + b][i * q + a] = GJ[i * q + a][j * q + b];
            P += gradJ * GJ[i * q + a][j * q + b];
            modF += gradJ * gradJ;
            modv += GJ[i * q + a][j * q + b] * GJ[i * q + a][j * q + b];
          }
        }
      }
    }
    modF = sqrt(modF);
    modv = sqrt(modv);
    for (int i = 0; i < L; i++)
    {
      for (int a = 0; a < q; a++)
      {
        Gh[i * q + a] = (1 - alpha) * Gh[i * q + a] + alpha * (fm[i * q + a] - mstat->fm_s[i * q + a]) / modF * modv;
        for (int j = i + 1; j < L; j++)
        {
          for (int b = 0; b < q; b++)
            GJ[i * q + a][j * q + b] = (1 - alpha) * GJ[i * q + a][j * q + b] + alpha * (sm[i * q + a][j * q + b] - mstat->sm_s[i * q + a][j * q + b]) / modF * modv;
        }
      }
    }
    if (P >= 0)
    {
      counter++;
      if (counter > 5)
      {
        acc = min(acc * 1.1, 10);
        alpha *= 0.99;
      }
    }
    else
    {
      acc = 0.5 * acc;
      alpha = 0.1;
      fill(Gh.begin(), Gh.end(), 0);
      fill(GJ.begin(), GJ.end(), Gh);
      counter = 0;
    }
    cout << "FIRE step - acc: " << acc << " alpha: " << alpha << " counter: " << counter << " P: " << P / modF / modv << endl;
    return acc * params->lrateh;
  }
}

/******************** METHODS FOR DECIMATION ***********************************************************/
int Model::n_total()
{
  return (L * (L - 1) / 2) * q * q;
}

int Model::n_active()
{
  return (int)((1.0 - model_sp) * (L * (L - 1) / 2) * q * q);
}

int Model::n_links()
{
  return (int)((1.0 - model_sp) * (L * (L - 1) / 2) * q * q);
}

void Model::init_decimation_variables()
{
  if (params->sparsity > 0 && params->compwise)
  {
    idx.clear();
    vector<int> tmp(4, 0);
    int n = L * (L - 1) * q * q / 2;
    idx.resize(n, tmp);
    sorted_struct.clear();
    sorted_struct.resize(n, 0);
    tmp_idx.clear();
    tmp_idx.resize(n, 0);
    int k = 0;
    for (int i = 0; i < L; i++)
    {
      for (int j = i + 1; j < L; j++)
      {
        for (int a = 0; a < q; a++)
        {
          for (int b = 0; b < q; b++)
          {
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
  else if (params->blockwise && params->sparsity > 0)
  {
    int n = L * (L - 1) / 2;
    vector<int> tmp(2, 0);
    idx.resize(n, tmp);
    sorted_struct.clear();
    sorted_struct.resize(n, 0);
    tmp_idx.clear();
    tmp_idx.resize(n, 0);
    int k = 0;
    for (int i = 0; i < L; i++)
    {
      for (int j = i + 1; j < L; j++)
      {
        idx[k][0] = i;
        idx[k][1] = j;
        sorted_struct[k] = 0.0;
        k += 1;
      }
    }
  }
}

// THIS NEEDS TO BE CHECKED
int Model::decimate_blockwise(int iter)
{
  int i, j, a, b, index, m = 0;
  double smalln = min(1e-30, params->pseudocount * 0.03);
  double maxsdkl = -1e50;
  int sumJ;
  cout << "Decimating one block" << endl;
  for (int k = 0; k < int(tmp_idx.size()); k++)
  {
    tmp_idx[k] = k;
    i = idx[k][0];
    j = idx[k][1];
    sumJ = 0;
    for (a = 0; a < q; a++)
    {
      for (b = 0; b < q; b++)
        sumJ += decJ[i * q + a][j * q + b];
    }
    if (sumJ > 0)
    {
      m += 1;
      if (params->dec_sdkl)
      {
        double auxsm = smalln * rand01() + mstat->sm_s[i * q + a][j * q + b];
        double term0 = 0.0;
        double num = 0.0;
        double den = 0.0;
        for (a = 0; a < q; a++)
        {
          for (b = 0; b < q; b++)
          {
            term0 += J[i * q + a][j * q + b] * auxsm;
            num += J[i * q + a][j * q + b] * auxsm * exp(-J[i * q + a][j * q + b]);
            den += auxsm * exp(-J[i * q + a][j * q + b]);
          }
        }
        sorted_struct[k] = term0 - num / den;
        sorted_struct[k] += rand01() * smalln;
      }
      else if (params->dec_f)
      {
        for (a = 0; a < q; a++)
        {
          for (b = 0; b < q; b++)
          {
            sorted_struct[k] += smalln * rand01() + fabs(mstat->sm_s[i * q + a][j * q + b]);
          }
        }
      }
      else if (params->dec_J)
      {
        for (a = 0; a < q; a++)
        {
          for (b = 0; b < q; b++)
          {
            sorted_struct[k] += smalln * rand01() + fabs(J[i * q + a][j * q + b]);
          }
        }
      }
      maxsdkl = max(maxsdkl, sorted_struct[k]);
    }
    else
    {
      double f = rand01();
      sorted_struct[k] = int(tmp_idx.size()) + f; // to be optimized: elements should be removed instead of putting large numbers
    }
  }
  cout << "Non-zeros blocks before decimation " << m << endl;
  quicksort(sorted_struct, tmp_idx, 0, int(tmp_idx.size()) - 1);

  index = tmp_idx[0]; // take the one with smallest sDKL/sumof J / sumof sm
  i = idx[index][0];
  j = idx[index][1];
  for (a = 0; a < q; a++)
  {
    for (b = 0; b < q; b++)
    {
      J[i * q + a][j * q + b] = 0.0;
      J[j * q + b][i * q + a] = 0.0;
      decJ[i * q + a][j * q + b] = 0;
      decJ[j * q + b][i * q + a] = 0;
    }
  }

  index = tmp_idx[0];
  cout << "The sDKL associated with the removed block is " << sorted_struct[0] << " (i: " << idx[index][0] << " j: " << idx[index][1] << ")" << endl;
  double nref = (L * (L - 1) * q * q) / 2;
  model_sp += (q * q) / nref;
  cout << "Sparsity after decimation is " << model_sp << endl;
  return 0;
}

int Model::decimate_compwise(int c, int iter)
{
  int i, j, a, b, index, m = 0;
  double smalln = min(1e-30, params->pseudocount * 0.03);
  double maxsdkl = -1e50;
  cout << "Decimating " << c << " couplings" << endl;
  for (int k = 0; k < int(tmp_idx.size()); k++)
  {
    tmp_idx[k] = k;
    i = idx[k][0];
    j = idx[k][1];
    a = idx[k][2];
    b = idx[k][3];
    if (decJ[i * q + a][j * q + b] > 0)
    {
      m += 1;
      if (params->dec_sdkl)
      {
        double auxsm = smalln * rand01() + mstat->sm_s[i * q + a][j * q + b];
        sorted_struct[k] = J[i * q + a][j * q + b] * auxsm - (J[i * q + a][j * q + b] * exp(-J[i * q + a][j * q + b]) * auxsm) / (exp(-J[i * q + a][j * q + b]) * auxsm + 1 - auxsm);
        sorted_struct[k] += rand01() * smalln;
      }
      else if (params->dec_f)
      {
        sorted_struct[k] = smalln * rand01() + fabs(mstat->sm_s[i * q + a][j * q + b]);
      }
      else if (params->dec_J)
      {
        sorted_struct[k] = smalln * rand01() + fabs(J[i * q + a][j * q + b]);
      }
      maxsdkl = max(maxsdkl, sorted_struct[k]);
    }
    else
    {
      double f = rand01();
      sorted_struct[k] = int(tmp_idx.size()) + f; // to be optimized: elements should be removed instead of putting large numbers
    }
  }
  cout << "Active couplings before decimation " << m << endl;
  quicksort(sorted_struct, tmp_idx, 0, int(tmp_idx.size()) - 1);
  for (int k = 0; k < c; k++)
  {
    index = tmp_idx[k];
    i = idx[index][0];
    j = idx[index][1];
    a = idx[index][2];
    b = idx[index][3];
    if (decJ[i * q + a][j * q + b] == 0)
      cerr << "Error: coupling " << i << " " << j << " " << a << " " << b << "was already decimated" << endl;
    J[i * q + a][j * q + b] = 0.0;
    J[j * q + b][i * q + a] = 0.0;
    decJ[i * q + a][j * q + b] = 0;
    decJ[j * q + b][i * q + a] = 0;
  }
  index = tmp_idx[c];
  fprintf(stdout, "Smallest sDKL associated with removed couplings is %.2e (i: %i j: %i a: %d b: %d)\n", sorted_struct[c], idx[index][0], idx[index][1], idx[index][2], idx[index][3]);
  double nref = (L * (L - 1) * q * q) / 2;
  model_sp += c / nref;
  fprintf(stdout, "Sparsity after decimation is %.3f\n", model_sp);
  return 0;
}

int Model::decimate_ising(int c, int iter)
{
  int i, j, index, m = 0;
  double num, den;
  double smalln = min(1e-30, params->pseudocount * 0.03);
  double maxsdkl = -1e50;
  cout << "Decimating " << c << " couplings" << endl;
  for (int k = 0; k < int(tmp_idx.size()); k++)
  {
    tmp_idx[k] = k;
    i = idx[k][0];
    j = idx[k][1];
    if (decJ[i][j] > 0)
    {
      m += 1;
      if (params->dec_sdkl)
      {
        double auxsm = smalln * rand01() + mstat->sm_s[i][j];
        num = (1 + auxsm) * J[i][j] * exp(-J[i][j]) - (1 - auxsm) * J[i][j] * exp(J[i][j]);
        den = (1 + auxsm) * exp(-J[i][j]) + (1 - auxsm) * exp(J[i][j]);
        sorted_struct[k] = J[i][j] * auxsm - num / den;
        sorted_struct[k] += rand01() * smalln;
      }
      else if (params->dec_f)
      {
        sorted_struct[k] = smalln * rand01() + fabs(mstat->sm_s[i][j]);
      }
      else if (params->dec_J)
      {
        sorted_struct[k] = smalln * rand01() + fabs(J[i][j]);
      }
      maxsdkl = max(maxsdkl, sorted_struct[k]);
    }
    else
    {
      double f = rand01();
      sorted_struct[k] = int(tmp_idx.size()) + f; // to be optimized: elements should be removed instead of putting large numbers
    }
  }
  cout << "Non-zeros parameters before decimation " << m << endl;
  quicksort(sorted_struct, tmp_idx, 0, int(tmp_idx.size()) - 1);
  for (int k = 0; k < c; k++)
  {
    index = tmp_idx[k];
    i = idx[index][0];
    j = idx[index][1];
    if (decJ[i][j] == 0)
      cerr << "Error: coupling " << i << " " << j << " was already decimated" << endl;
    J[i][j] = 0.0;
    J[j][i] = 0.0;
    decJ[i][j] = 0;
    decJ[j][i] = 0;
  }
  index = tmp_idx[c];
  fprintf(stdout, "Smallest sDKL associated with removed couplings is %.2e (i: %i j: %i)\n", sorted_struct[c], idx[index][0], idx[index][1]);
  double nref = (L * (L - 1)) / 2;
  model_sp += c / nref;
  cout << "Sparsity after decimation is " << scientific << setprecision(3) << model_sp << endl;
  return 0;
}

int Model::activate_compwise(int c, int iter, vector<vector<MYFLOAT>> &sm)
{
  int i, j, a, b, index, m = 0;
  double smalln = min(1e-30, params->pseudocount * 0.03);
  double mindl = 1e50;
  double pc = 0.01;
  printf("Activating %d couplings\n", c);
  for (int k = 0; k < int(tmp_idx.size()); k++)
  {
    tmp_idx[k] = k;
    i = idx[k][0];
    j = idx[k][1];
    a = idx[k][2];
    b = idx[k][3];
    if (decJ[i * q + a][j * q + b] < 0.5)
    {
      m += 1;
      double f_d = (1 - pc) * sm[i * q + a][j * q + b] + pc / (q * q);
      double f_m = (1 - pc) * mstat->sm_s[i * q + a][j * q + b] + pc / (q * q);
      sorted_struct[k] = f_d * log(f_d / f_m) + (1 - f_d) * log((1 - f_d) / (1 - f_m));
      sorted_struct[k] += rand01() * smalln;
      mindl = min(mindl, sorted_struct[k]);
      // fprintf(fileout, "J %i %i %i %i %.2e %f %f\n", i, j, a, b, sorted_struct[k], J[i*q+a][j*q+b], sm_s[i*q+a][j*q+b]);
    }
    else
    {
      double f = rand01();
      sorted_struct[k] = -int(tmp_idx.size()) + f; // to be optimized: elements should be removed instead of putting large numbers
    }
  }
  cout << "Inactive couplings before decimation " << m << endl;
  quicksort(sorted_struct, tmp_idx, 0, int(tmp_idx.size()) - 1);
  for (int k = 0; k < c; k++)
  {
    index = tmp_idx[int(tmp_idx.size()) - 1 - k];
    i = idx[index][0];
    j = idx[index][1];
    a = idx[index][2];
    b = idx[index][3];
    if (decJ[i * q + a][j * q + b] > 0.5)
    {
      printf("Error: coupling %d %d %d %d was already present\n", i, j, a, b);
      exit(1.);
    }
    J[i * q + a][j * q + b] = 0;
    J[j * q + b][i * q + a] = 0;
    decJ[i * q + a][j * q + b] = 1;
    decJ[j * q + b][i * q + a] = 1;
  }
  index = tmp_idx[int(tmp_idx.size()) - 1 - c];
  fprintf(stdout, "Smallest DeltaL associated with added couplings is %.2e (i: %i j: %i a: %i b: %i) \n", sorted_struct[int(tmp_idx.size()) - 1 - c], idx[index][0], idx[index][1], idx[index][2], idx[index][3]);
  double nref = (L * (L - 1) * q * q) / 2;
  model_sp -= c / nref;
  cout << "Sparsity after decimation is " << scientific << setprecision(3) << model_sp << endl;
  return 0;
}
