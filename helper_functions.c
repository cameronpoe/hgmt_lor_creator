#include "helper_functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// gives a random number following a guassian distribution
double gaussian(double sd, int num_additions) {
  double rand = 0.0;
  for (int i = 0; i < num_additions; i++) {
    rand += drand48();
  }
  rand -= (double)num_additions / 2.0;
  rand *= sqrt(12.0) * sd;
  rand /= sqrt((double)num_additions);
  return rand;
}
int num_args(int argc, char **argv) {
  int numargs = 0;
  for (int i = 1; i < argc; i++) {
    if ('-' != argv[i][0]) {
      numargs++;
    }
  }
  return numargs;
}
int num_flags(int argc, char **argv) {
  int numflags = 0;
  for (int i = 1; i < argc; i++) {
    if ('-' == argv[i][0]) {
      numflags++;
    }
  }
  return numflags;
}
char **get_args(int argc, char **argv) {
  int numargs = num_args(argc, argv);
  char **args = malloc(sizeof(char *) * numargs);
  int j = 0;
  for (int i = 1; i < argc; i++) {
    if ('-' != argv[i][0]) {
      int len = strlen(argv[i]) + 1;
      args[j] = malloc(sizeof(char) * len);
      strcpy(args[j], argv[i]);
      j++;
    }
  }
  return args;
}
char **get_flags(int argc, char **argv) {
  int numflags = num_flags(argc, argv);
  char **flags = malloc(sizeof(char *) * numflags);
  int j = 0;
  for (int i = 1; i < argc; i++) {
    if ('-' == argv[i][0]) {
      int len = strlen(argv[i]) + 1;
      flags[j] = malloc(sizeof(char) * len);
      strcpy(flags[j], argv[i]);
      j++;
    }
  }
  return flags;
}
