
// including standard files
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// including custom files
#include "compton_chain_ordering.h"
#include "helper_functions.h"
#include "hgmt_structs.h"
#include "llist.h"
#include "vector_ops.h"
uint num_datas;
FILE *debug_out;
void print_data(double *data, FILE *output) {
  fwrite(data, sizeof(double), 1, output);
}
double read_double(FILE *input) {
  double num;
  int worked = fread(&num, sizeof(double), 1, input);

  if (worked != 1) {
    return -1;
  }
  // make a new event to be passed out
  return num;
}
void error_debug(FILE *input) {
  histogram *hist = new_histogram(0.0, 100.0, 40);

  double num = read_double(input);
  double tot = 0;
  while (num != -1) {
    num_datas++;
    tot += num;
    add_to_histogram(num, hist);
    num = read_double(input);
  }
  printf("number of lors: %i\n", num_datas);
  printf("average error: %lf\n", (double)tot / num_datas);
  print_histogram(hist);
}
int main(int argc, char **argv) {
  debug_out = fopen("debug_out.data", "wb");
  if (debug_out == NULL) {
    printf("Unable to open output file for writing\n");
    return 1;
  }
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  // defines the help function and how to call it (by using -h or --help)
  for (int i = 0; i < num_flags(argc, argv); i++) {
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./hgmt_debug -task debug_file_loc arg2\n");
      printf("-h: print this help\n");
      printf("-e: run lor error debug\n");
      exit(0);
    }
    if (strcmp(flags[i], "-e") == 0) {
      printf("running with lor error debug\n");
      FILE *input = fopen(args[0], "rb");
      error_debug(input);
    }
  }
}
