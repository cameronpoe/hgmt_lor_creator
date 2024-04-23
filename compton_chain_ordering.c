
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
#include "vector_ops.h"
#define SPD_LGHT 29.9792458 // cm/ns
#define SPC_UNC 0.1         // cm
#define TIME_UNC 0.042463   // ns, sigma (0.100 ns FWHM)
#define DETECTOR_THICKNESS 2.54
#define DETECTOR_SEGMENTATION 1
double time_FOM(photon_path path, int *order) {
  double var_rad = DETECTOR_THICKNESS * DETECTOR_THICKNESS / 48.0;
  double var_tan = SPC_UNC * SPC_UNC;
  double FOM = 0;
  for (int i = 0; i < path.num_hits - 1; i++) {
    hit *hit1 = &path.hits[order[i]];
    hit *hit2 = &path.hits[order[i + 1]];
    // calculating sigma values
    double rad1 = radial_dist(hit1->location);
    double rad2 = radial_dist(hit2->location);
    double d = vec_dist(hit1->location, hit2->location);
    double delta_x = hit1->location.x - hit2->location.x;
    double delta_y = hit1->location.y - hit2->location.y;
    double delta_z = hit1->location.z - hit2->location.z;
    double cos_theta1 = hit1->location.x / rad1;
    double sin_theta1 = hit1->location.y / rad1;
    double cos_theta2 = hit2->location.x / rad2;
    double sin_theta2 = hit1->location.y / rad2;
    double var_delta_t = TIME_UNC * TIME_UNC * 2;
    double var_delta_x =
        cos_theta1 * cos_theta1 * var_rad + sin_theta1 * sin_theta1 * var_tan +
        cos_theta2 * cos_theta2 * var_rad + sin_theta2 * sin_theta2 * var_tan;
    double var_delta_y =
        sin_theta1 * sin_theta1 * var_rad + cos_theta1 * cos_theta1 * var_tan +
        sin_theta2 * sin_theta2 * var_rad + cos_theta2 * cos_theta2 * var_tan;
    double var_delta_z = 2 * SPC_UNC * SPC_UNC;
    double var_d = delta_x * delta_x * var_delta_x / (d * d) +
                   delta_y * delta_y * var_delta_y / (d * d) +
                   delta_z * delta_z * var_delta_z / (d * d);
    double var_dt = var_d / (SPD_LGHT * SPD_LGHT);
    // calculating FOM
    double delta_t = fabs(hit1->tof - hit2->tof);
    FOM += fabs(delta_t - d / SPD_LGHT) / sqrt(var_delta_t + var_dt);
  }
  return FOM;
}
int *best_order(photon_path path, double (*FOM)(photon_path, int *)) {
  int num_hits = path.num_hits;
  int factor = factorial(num_hits);
  int *the_best_order = (int *)calloc(sizeof(int), num_hits);
  double bestFOM = 0;
  perm *permutation = first_perm(num_hits);
  // print_perm(permutation);
  for (int i = 0; i < factor; i++) {
    double figure_of_merit = (*FOM)(path, permutation->perm);
    if (figure_of_merit > bestFOM) {
      bestFOM = figure_of_merit;
      memcpy(the_best_order, permutation->perm, sizeof(int) * num_hits);
    }
    if (i != factor - 1) {
      increment_perm(permutation);
    }
  }
  free_perm(permutation);
  return the_best_order;
}
hit *initial_by_best_order(photon_path path,
                           double (*FOM)(photon_path, int *)) {
  int *the_best_order = best_order(path, FOM);
  int index;
  if (path.hits[the_best_order[0]].tof <
      path.hits[the_best_order[path.num_hits - 1]].tof) {
    index = the_best_order[0];
  } else {
    index = the_best_order[path.num_hits - 1];
  }
  free(the_best_order);
  return &path.hits[index];
}
hit *initial_by_best_time(photon_path path) {
  double best_tof = path.hits[0].tof;
  hit *initial = &path.hits[0];
  for (int i = 1; i < path.num_hits; i++) {
    if (path.hits[i].tof < best_tof) {
      best_tof = path.hits[i].tof;
      initial = &path.hits[i];
    }
  }
  return initial;
}
