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
#include <fstream>
#include <cctype>
#include "BMaux.h"

using namespace std;
typedef double MYFLOAT;

double rand01()
{
	return (double)rand() / (double)((unsigned)RAND_MAX + 1);
}

double randrange(int xmin, int xmax)
{
	return xmin + (xmax - xmin) * rand() / (double)RAND_MAX;
}

double max(double x, double y)
{
	return (x > y) ? x : y;
}

double min(double x, double y)
{
	return (x < y) ? x : y;
}

void permute(int *vector, int s)
{
	int i, j, aux;
	for (i = 0; i < s; i++)
	{
		j = i + rand() % (s - i);
		aux = vector[i];
		vector[i] = vector[j];
		vector[j] = aux;
	}
}

int overlap(vector<unsigned char> &x1, vector<unsigned char> &x2)
{
	int q = 0;
	for (int i = 0; i < int(x1.size()); i++)
	{
		if (x1[i] == x2[i])
			q++;
	}
	return q;
}

int overlap_two(vector<int> &x1, vector<int> &x2)
{
	int Q = 0;
	for (int i = 0; i < int(x1.size()); i++)
	{
		if (x1[i] == x2[i])
		{
			Q++;
			for (int j = i + 1; j < int(x1.size()); j++)
			{
				Q += int(x1[j] == x2[j]);
			}
		}
	}
	return Q;
}

int convert_char_amino(char a)
{
	int i;
	switch (a)
	{
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
		i = 18;
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
		cerr << a << "not recognized" << endl;
		return (EXIT_FAILURE);
	}
	return (unsigned char)i;
}

int convert_char_nbase(char a)
{
	int i;
	switch (a)
	{
	case '-':
		i = 0;
		break;
	case 'A':
		i = 1;
		break;
	case 'C':
		i = 2;
		break;
	case 'G':
		i = 3;
		break;
	case 'U':
		i = 4;
		break;
	case 'T':
		i = 4;
		break;
	default:
		cerr << a << "not recognized" << endl;
		i = 0;
		break;
	}
	return (unsigned char)i;
}

int convert_char_epi(char a)
{
	int i;
	switch (a)
	{
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
		cerr << a << "not recognized, assuming '-'" << endl;
		i = 0;
		break;
		// return(EXIT_FAILURE);
	}
	return (unsigned char)i;
}

int convert_char_ising(char a)
{
	int i;
	switch (a)
	{
	case '0':
		i = 0;
		break;
	case 'A':
		i = 0;
		break;
	case 'P':
		i = 1;
		break;
	case '1':
		i = 1;
		break;
	case 'u':
		i = 1;
		break;
	case 'd':
		i = 0;
		break;
	default:
		cerr << a << "not recognized" << endl;
		return (EXIT_FAILURE);
	}
	return (unsigned char)i;
}

int convert_char_potts(char a)
{
	return (int)(a - '0');
}

int convert_char_adhoc(char a, char* abc)
{
	std::string abcd;
	abcd = abc;
	return (int)abcd.find(a);
}

int convert_char(char ch, char *ctype)
{
	int ris;
	std::string abc;
	abc = ctype;
	if (!strcmp(ctype, "a"))
		ris = convert_char_amino(ch);
	else if (!strcmp(ctype, "n"))
		ris = convert_char_nbase(ch);
	else if (!strcmp(ctype, "i"))
		ris = convert_char_ising(ch);
	else if (!strcmp(ctype, "e"))
		ris = convert_char_epi(ch);
	else if (isdigit(*ctype))
		ris = convert_char_potts(ch);
	else if (abc.length() > 1)
		ris = convert_char_adhoc(ch, ctype);
	else
	{
		cerr << "Error in alphabet specification" << endl;
		exit(EXIT_FAILURE);
	}
	return ris;
}

int print_alphabet(char *ctype)
{

	int q = 0;
	std::string abc;
	abc = ctype;
	if (!strcmp(ctype, "a"))
	{
		cout << "Using alphabet: -ACDEFGHIKLMNPQRSTVWY" << endl;
		q = 21;
	}
	else if (!strcmp(ctype, "n"))
	{
		cout << "Using alphabet: -ACGU" << endl;
		q = 5;
	}
	else if (!strcmp(ctype, "i"))
	{
		cout << "Using alphabet: {-1,1} spins. Input: AP (absent/present) or ud (down/up) or binary {0,1}" << endl;
		q = 1;
	}
	else if (!strcmp(ctype, "e"))
	{
		cout << "Using alphabet: -AF5TtGEZhBbeRrq \n"
			 << endl;
		q = 16;
	}
	else if (isdigit(*ctype) && sscanf(ctype, "%d", &q) == 1)
	{
		cout << "Using alphabet: Potts with q = " << q << endl;
	}
	else if (abc.length() > 1)
	{
		cout << "Using alphabet: " << abc << endl;
		q = abc.length();
	}
	else
	{
		cerr << "Use 'a' for amino-acids, 'n' for nitrogenous bases or 'i' for Ising-like variables" << endl;
		return EXIT_FAILURE;
	}
	return q;
}

vector<char> alphabet(char *ctype)
{

	vector<char> ris;
	std::string abc;
	abc = ctype;
	int q;
	if (!strcmp(ctype, "a"))
	{
		char app[] = "-ACDEFGHIKLMNPQRSTVWY\n";
		ris = vector<char>(app, app + sizeof(app) / sizeof(char));
	}
	else if (!strcmp(ctype, "n"))
	{
		char app[] = "-ACGU";
		ris = vector<char>(app, app + sizeof(app) / sizeof(char));
	}
	else if (!strcmp(ctype, "i"))
	{
		char app[] = "01";
		ris = vector<char>(app, app + sizeof(app) / sizeof(char));
	}
	else if (!strcmp(ctype, "e"))
	{
		char app[] = "-AF5TtGEZhBbeRrq";
		ris = vector<char>(app, app + sizeof(app) / sizeof(char));
	}
	else if (isdigit(*ctype) && sscanf(ctype, "%d", &q) == 1)
	{
		ris.clear();
		for (int i = 0; i < q; i++)
			ris.push_back(i + '0');
	}
	else if (abc.length() > 1)
	{
		ris.resize(abc.length());
		memcpy(&ris[0], ctype, abc.length());
	}
	else
	{
		cerr << "Error in alphabet specification" << endl;
		cerr << "Use 'a' for amino-acids or 'n' for nitrogenous bases" << endl;
		cerr << "Use 'i' for ising variables" << endl;
		exit(EXIT_FAILURE);
	}
	return ris;
}
