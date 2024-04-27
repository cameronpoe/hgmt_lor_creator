
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
double tot_error;
uint num_lors;
lor *read_lor(FILE *input) {
  if (input == NULL) {
    return NULL;
  }
  lor *new = (lor *)malloc(sizeof(lor));
  int worked = fread(new, sizeof(lor), 1, input);

  if (worked != 1) {
    return NULL;
  }
  // make a new event to be passed out
  return new;
}
void error_debug(FILE *perfects, FILE *fuzzies) {
  histogram *hist = new_histogram(0.0, 100.0, 40);
  lor *perfect = read_lor(perfects);
  lor *fuzzy = read_lor(fuzzies);
  while (fuzzy != NULL && perfect != NULL) {
    num_lors++;
    double dist = vec_dist(perfect->center, fuzzy->center);
    tot_error += dist;
    add_to_histogram(dist, hist);
    free(perfect);
    free(fuzzy);
    perfect = read_lor(perfects);
    fuzzy = read_lor(fuzzies);
  }
  printf("number of lors: %i\n", num_lors);
  printf("average error: %lf\n", (double)tot_error / num_lors);
  print_histogram(hist);
}
int main(int argc, char **argv) {
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
      FILE *perfects = fopen(args[0], "rb");
      FILE *wrongs = fopen(args[1], "rb");
      error_debug(perfects, wrongs);
    }
  }
}
