
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
event *read_event(FILE *source) {

  uint event_id;
  double energy_deposit;
  float x;
  float y;
  float z;
  // float mx;
  // float my;
  // float mz;
  float tof;
  int particle_type;
  int track_id;
  int worked = 0;

  worked += fread(&event_id, sizeof(uint), 1, source);
  worked += fread(&energy_deposit, sizeof(double), 1, source);
  worked += fread(&x, sizeof(float), 1, source);
  worked += fread(&y, sizeof(float), 1, source);
  worked += fread(&z, sizeof(float), 1, source);
  // worked += fread(&mx, sizeof(float), 1, source);
  // worked += fread(&my, sizeof(float), 1, source);
  // worked += fread(&mz, sizeof(float), 1, source);
  worked += fread(&tof, sizeof(float), 1, source);
  worked += fread(&particle_type, sizeof(int), 1, source);
  worked += fread(&track_id, sizeof(int), 1, source);

  if (worked != 8) {
    return NULL;
  }

  // make a new event to be passed out
  event *new_event = (event *)malloc(sizeof(event));
  new_event->event_id = event_id;
  new_event->energy_deposit = energy_deposit;
  new_event->location = three_vec((double)x, (double)y, (double)z);
  // new_event->momentum = three_vec((double)mx, (double)my, (double)mz);
  new_event->tof = (double)tof;
  new_event->particle_type = particle_type;
  new_event->track_id = track_id;
  return new_event;
}
event *read_gamma(FILE *source) {
  event *new_event = read_event(source);
  while (new_event != NULL && new_event->particle_type != 22) {
    free(new_event);
    new_event = read_event(source);
  }
  // printf("%i %i\n", new_event->track_id, new_event->event_id);
  return new_event;
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
event *read_history(int event_id, FILE *source) {
  event *new_event = read_gamma(source);
  while (new_event->event_id == event_id && new_event != NULL) {
    free(new_event);
    new_event = read_gamma(source);
  }
  return new_event;
}
void phsp_diagnostics(FILE *source) {
  event *new_event = read_gamma(source);
  int num_events = 0;
  while (new_event != NULL) {
    free(new_event);
    event_id = read_history(event_id, source);
    num_events++;
  }
  printf("number of events: %i\n", num_events);
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
      printf("-hi: run with histogram\n");
      printf("-p: run diagnostics on phsp file\n");
      exit(0);
    }
    if (strcmp(flags[i], "-hi") == 0) {
      printf("running histogram\n");
      FILE *input = fopen(args[0], "rb");
      error_debug(input);
    }
    if (strcmp(flags[i], "-p") == 0) {
      printf("diagnostics on phsp file\n");
      FILE *input = fopen(args[0], "rb");
      phsp_diagnostics(input);
    }
  }
}
