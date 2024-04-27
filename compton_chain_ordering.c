
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
double variance_dist(vec3d loc1, vec3d loc2, double d) {
  double var_rad = DETECTOR_THICKNESS * DETECTOR_THICKNESS / 48.0;
  double var_tan = SPC_UNC * SPC_UNC;
  double rad1 = radial_dist(loc1);
  double rad2 = radial_dist(loc2);
  double delta_x = loc1.x - loc2.x;
  double delta_y = loc1.y - loc2.y;
  double delta_z = loc1.z - loc2.z;
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
double time_FOM_cum(photon_path *path, int *order) {
  double FOM = 0;
  double var_t = TIME_UNC * TIME_UNC;
  double var_dt = 0;
  double dt = 0;
  double t_0 = path->hits[order[0]].tof;
  for (int i = 0; i < path->num_hits - 1; i++) {
    hit *hit1 = &path->hits[order[i]];
    hit *hit2 = &path->hits[order[i + 1]];
    // calculating sigma values
    double d = vec_dist(hit1->location, hit2->location);
    dt += d / SPD_LGHT;
    double var_d = variance_dist(hit1->location, hit2->location, d);
    var_dt += var_d / (SPD_LGHT * SPD_LGHT);
    double rel_t = hit2->tof - t_0;
    // calculating FOM
    FOM += fabs(rel_t - dt) / sqrt(var_t + var_dt);
  }
  return FOM;
}
double time_FOM(photon_path *path, int *order) {
  double FOM = 0;
  for (int i = 0; i < path->num_hits - 1; i++) {
    hit *hit1 = &path->hits[order[i]];
    hit *hit2 = &path->hits[order[i + 1]];
    // calculating sigma values
    double d = vec_dist(hit1->location, hit2->location);
    double var_d = variance_dist(hit1->location, hit2->location, d);
    double var_delta_t = TIME_UNC * TIME_UNC * 2;
    double var_dt = var_d / (SPD_LGHT * SPD_LGHT);
    // calculating FOM
    double delta_t = fabs(hit1->tof - hit2->tof);
    FOM += fabs(delta_t - d / SPD_LGHT) / sqrt(var_delta_t + var_dt);
  }
  return FOM;
}
hit *initial_by_best_order(photon_path *path,
                           double (*FOM)(photon_path *, int *)) {
  int num_hits = path->num_hits;
  int factor = factorial(num_hits);
  perm *permutation = first_perm(num_hits);
  hit *the_best_hit = path->hits;
  double bestFOM = (*FOM)(path, permutation->perm);
  for (int i = 0; i < factor - 1; i++) {
    increment_perm(permutation);
    double figure_of_merit = (*FOM)(path, permutation->perm);
    if (figure_of_merit < bestFOM) {
      bestFOM = figure_of_merit;
      the_best_hit = &path->hits[permutation->perm[0]];
    } else if (figure_of_merit == bestFOM &&
               the_best_hit->tof > path->hits[permutation->perm[0]].tof) {
      the_best_hit = &path->hits[permutation->perm[0]];
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
hit *initial_by_weighted(photon_path *path, double (*FOM)(photon_path *, int *),
                         double time_weight) {
  int num_hits = path->num_hits;
  int factor = factorial(num_hits);
  perm *permutation = first_perm(num_hits);
  double best_score = path->hits[0].tof * time_weight +
                      (1 - time_weight) * (*FOM)(path, permutation->perm);
  hit *best_hit = path->hits;
  for (int i = 0; i < factor - 1; i++) {
    increment_perm(permutation);
    double score = path->hits[permutation->perm[0]].tof * time_weight +
                   (1 - time_weight) * (*FOM)(path, permutation->perm);
    if (score < best_score) {
      best_score = score;
      best_hit = &path->hits[permutation->perm[0]];
    }
  }
  free_perm(permutation);
  return best_hit;
}
