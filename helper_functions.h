#ifndef helper_functions_h
#define helper_functions_h

#include <stdio.h>

// gives a random number following a guassian distribution
double gaussian(double sd, int num_additions);
char **get_args(int argc, char **argv);
char **get_flags(int argc, char **argv);
int num_flags(int argc, char **argv);
int num_args(int argc, char **argv);
#endif
