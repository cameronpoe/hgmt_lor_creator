#ifndef helper_functions_h
#define helper_functions_h

#include <stdio.h>

// gives a random number following a guassian distribution
double gaussian(double sd, int num_additions);
char **get_args(int argc, char **argv);
char **get_flags(int argc, char **argv);
int num_flags(int argc, char **argv);
int num_args(int argc, char **argv);
void printm(char message[], int counter, int mod);
int factorial(int num);
typedef struct perm_ {
  int *perm;
  int *index;
  int *places;
  int *parity;
  int length;
} perm;
perm *first_perm(int n);
void free_perm(perm *permutation);
void increment_perm(perm *permutation);
#endif
