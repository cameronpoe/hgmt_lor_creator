
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
#define TIME_SD_TOLERANCE 3
#define OPTIMIZE 1
double variance_dist(vec3d loc1, vec3d loc2, double d) {
  double delta_x = loc1.x - loc2.x;
  double delta_y = loc1.y - loc2.y;
  double delta_z = loc1.z - loc2.z;
  double rad1 = radial_dist(loc1);
  double rad2 = radial_dist(loc2);
  double var_rad;
  if (DETECTOR_SEGMENTATION) {
    double var_rad = DETECTOR_THICKNESS * DETECTOR_THICKNESS / 48.0;
  } else {
    var_rad = SPC_UNC;
  }
  double var_tan = SPC_UNC * SPC_UNC;
  double cos_theta1 = loc1.x / rad1;
  double sin_theta1 = loc1.y / rad1;
  double cos_theta2 = loc2.x / rad2;
  double sin_theta2 = loc1.y / rad2;
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
  return var_d;
}
double time_FOM_cum(hit *hits, int *order, int num_hits) {
  double FOM = 0;
  double var_t = TIME_UNC * TIME_UNC;
  double var_dt = 0;
  double dt = hits[order[0]].tof;
  for (int i = 0; i < num_hits - 1; i++) {
    hit *hit1 = &hits[order[i]];
    hit *hit2 = &hits[order[i + 1]];
    // calculating sigma values
    double d = vec_dist(hit1->location, hit2->location);
    dt += d / SPD_LGHT;
    double var_d = variance_dist(hit1->location, hit2->location, d);
    var_dt += var_d / (SPD_LGHT * SPD_LGHT);
    // calculating FOM
    FOM += fabs(hit2->tof - dt) / sqrt(var_t + var_dt);
  }
  return FOM;
}
double time_FOM(hit *hits, int *order, int num_hits) {
  double FOM = 0;
  for (int i = 0; i < num_hits - 1; i++) {
    hit *hit1 = &hits[order[i]];
    hit *hit2 = &hits[order[i + 1]];
    // calculating sigma values
    double d = vec_dist(hit1->location, hit2->location);
    double var_d = variance_dist(hit1->location, hit2->location, d);
    double var_delta_t = TIME_UNC * TIME_UNC * 2;
    double var_dt = var_d / (SPD_LGHT * SPD_LGHT);
    // calculating FOM
    double delta_t = hit2->tof - hit1->tof;
    FOM += fabs(delta_t - d / SPD_LGHT) / sqrt(var_delta_t + var_dt);
  }
  return FOM;
}
int compare_hits(const void *hit1, const void *hit2) {
  return ((hit *)hit1)->tof > ((hit *)hit2)->tof;
}
hit *initial_by_best_order(photon_path *path,
                           double (*FOM)(hit *hits, int *order, int num_hits)) {
  int num_hits = path->num_hits;
  hit *hits = path->hits;
  if (OPTIMIZE) {
    qsort(hits, num_hits, sizeof(hit), compare_hits);
    for (int i = 0; i < num_hits; i++) {
      if (hits[i].tof - hits[0].tof > TIME_SD_TOLERANCE * TIME_UNC) {
        num_hits = i + 1;
        break;
      }
    }
  }
  int factor = factorial(num_hits);
  perm *permutation = first_perm(num_hits);
  hit *the_best_hit = hits;
  double bestFOM = (*FOM)(hits, permutation->perm, num_hits);
  for (int i = 0; i < factor - 1; i++) {
    increment_perm(permutation);
    double figure_of_merit = (*FOM)(hits, permutation->perm, num_hits);
    if (figure_of_merit < bestFOM) {
      bestFOM = figure_of_merit;
      the_best_hit = &hits[permutation->perm[0]];
    } else if (figure_of_merit == bestFOM &&
               the_best_hit->tof > hits[permutation->perm[0]].tof) {
      the_best_hit = &hits[permutation->perm[0]];
    }
  }
  free_perm(permutation);
  return the_best_hit;
}
hit *initial_by_best_time(photon_path *path) {
  double best_tof = path->hits[0].tof;
  hit *initial = &path->hits[0];
  for (int i = 1; i < path->num_hits; i++) {
    if (path->hits[i].tof < best_tof) {
      best_tof = path->hits[i].tof;
      initial = &path->hits[i];
    }
  }
  return initial;
}
