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

// global variables
#define PI 3.141592653589
#define COLS 105
#define SPD_LGHT 29.9792458 // cm/ns
#define UNCERT_REP 30
#define SPC_UNC 0.1
#define TIME_UNC 0.042463 // ns, sigma (0.100 ns FWHM)
#define DETECTOR_THICKNESS 2.54
#define DETECTOR_SEGMENTATION 1
uint num_scatters = 0;
uint scatters_detected = 0;
uint annihilations_occured = 0;
uint lors_created = 0;
uint first_correct = 0;
uint first_guessed = 0;
event *first_event;
double eff_by_energy[COLS];
double E_max = 520.0;
double E_min = 0.0;

double linear_interpolation(double nums[COLS], double min, double max,
                            double value) {
  double i = (COLS - 1) * (value - min) / (max - min);
  int i_l = (int)i;
  int i_r = i_l + 1;
  double i_space = i - i_l;
  // printf("%lf \n", nums[i_l] * i_space + nums[i_r] * (1.0 - i_space));
  return nums[i_l] * i_space + nums[i_r] * (1.0 - i_space);
}
int read_eff(FILE *source) {

  // exits with code 1 if the source pointer isn't pointing to the file
  if (source == NULL) {
    return 1;
  }

  // loops through all the entries in a row
  for (int i = 0; i < COLS; i++) {
    int worked = fscanf(source, "%lf,", &eff_by_energy[i]);
  }
  return 0;
}
event *read_event(FILE *source) {

  uint event_id;
  double energy_deposit;
  float x;
  float y;
  float z;
  float mx;
  float my;
  float mz;
  float tof;
  int particle_type;
  int track_id;
  int worked = 0;

  worked += fread(&event_id, sizeof(uint), 1, source);
  worked += fread(&energy_deposit, sizeof(double), 1, source);
  worked += fread(&x, sizeof(float), 1, source);
  worked += fread(&y, sizeof(float), 1, source);
  worked += fread(&z, sizeof(float), 1, source);
  worked += fread(&mx, sizeof(float), 1, source);
  worked += fread(&my, sizeof(float), 1, source);
  worked += fread(&mz, sizeof(float), 1, source);
  worked += fread(&tof, sizeof(float), 1, source);
  worked += fread(&particle_type, sizeof(int), 1, source);
  worked += fread(&track_id, sizeof(int), 1, source);

  if (worked != 11) {
    return NULL;
  }

  // make a new event to be passed out
  event *new_event = (event *)malloc(sizeof(event));
  new_event->event_id = event_id;
  new_event->energy_deposit = energy_deposit;
  new_event->location = three_vec((double)x, (double)y, (double)z);
  new_event->momentum = three_vec((double)mx, (double)my, (double)mz);
  new_event->tof = (double)tof;
  new_event->particle_type = particle_type;
  new_event->track_id = track_id;

  if (new_event->particle_type == 22) {
    num_scatters++;
    printm("Number of recorded gammas read: ", num_scatters, 100000);
  }
  return new_event;
}

uint is_scatter_detected(event *single_event) {
  if (single_event->particle_type == 22 &&
      drand48() < linear_interpolation(eff_by_energy, E_min, E_max,
                                       single_event->energy_deposit)) {
    scatters_detected++;
    return 1;
  }
  return 0;
}
hit *event_to_hit(event *single_event) {
  vec_mag(three_vec(single_event->location.x, single_event->location.y, 0.0));
  vec3d z_hat = three_vec(0.0, 0.0, 1.0);
  vec3d circ_hat = vec_norm(vec_cross(z_hat, single_event->location));
  vec3d offset = vec_add(vec_scale(z_hat, gaussian(SPC_UNC, 30)),
                         vec_scale(circ_hat, gaussian(SPC_UNC, 30)));
  hit *new_hit = (hit *)malloc(sizeof(hit));
  new_hit->location = single_event->location;
  new_hit->tof = single_event->tof + gaussian(TIME_UNC, 30);
  double rad_dist = radial_dist(new_hit->location);
  if (DETECTOR_SEGMENTATION) {
    // we move the radial component to one of the planes of the detector
    new_hit->location = radial_scale(
        new_hit->location,
        DETECTOR_THICKNESS * ((double)((int)(rad_dist / DETECTOR_THICKNESS))) /
                rad_dist +
            DETECTOR_THICKNESS / (2 * rad_dist));
  } else {
    vec3d r_hat =
        vec_scale(three_vec(new_hit->location.x, new_hit->location.y, 0),
                  1.0 / radial_dist(new_hit->location));
    offset = vec_add(offset, vec_scale(r_hat, gaussian(SPC_UNC, 30)));
  }
  new_hit->location = vec_add(new_hit->location, offset);
  return new_hit;
}

annihilation *read_annihilation(FILE *source) {
  // reading all the scatters within the time window
  event *first_event2 = NULL;
  llist *path1 = NULL;
  llist *path2 = NULL;
  int path1_count = 0;
  int path2_count = 0;
  int need_to_free_e1 = 1;
  if (first_event == NULL) {
    return NULL;
  }
  // organizing the scatters that are detected
  if (is_scatter_detected(first_event)) {
    path1 = add_to_top(path1, first_event);
    path1_count++;
    need_to_free_e1 = 0;
  }
  while (1) {
    event *new_event = read_event(source);
    if (new_event == NULL) {
      free(new_event);
      break;
    }
    if (first_event->event_id != new_event->event_id) {
      first_event2 = new_event;
      break;
    }
    if (is_scatter_detected(new_event)) {
      if (new_event->track_id == first_event->track_id) {
        path1 = add_to_top(path1, new_event);
        path1_count++;
      } else {
        path2 = add_to_top(path2, new_event);
        path2_count++;
      }
    } else {
      free(new_event);
    }
  }
  // generating the annihilation
  annihilation *new_annihilation = (annihilation *)malloc(sizeof(annihilation));
  new_annihilation->photon1_path.hits =
      (hit *)malloc(sizeof(hit) * path1_count);
  new_annihilation->photon2_path.hits =
      (hit *)malloc(sizeof(hit) * path2_count);
  new_annihilation->photon1_path.num_hits = path1_count;
  new_annihilation->photon2_path.num_hits = path2_count;
  // filling in all the values
  for (int i = path1_count - 1; i >= 0; i--) {
    hit *detector_hit = event_to_hit(path1->data);
    new_annihilation->photon1_path.hits[i] = *detector_hit;
    free(detector_hit);
    if (i != path1_count - 1) {
      path1 = path1->down;
    }
  }
  for (int i = path2_count - 1; i >= 0; i--) {
    hit *detector_hit = event_to_hit(path2->data);
    new_annihilation->photon2_path.hits[i] = *detector_hit;
    free(detector_hit);
    if (i != path2_count - 1) {
      path2 = path2->down;
    }
  }
  // freeing all the redundant data
  if (need_to_free_e1) {
    free(first_event);
  }
  first_event = first_event2;
  wipe_list(path1);
  wipe_list(path2);
  // returning
  return new_annihilation;
}
prim_lor *create_prim_lor(annihilation *new_annihilation) {
  // hit *hit1 = initial_by_best_order(new_annihilation->photon1_path,
  // time_FOM); hit *hit2 =
  // initial_by_best_order(new_annihilation->photon2_path, time_FOM);
  hit *hit1 = initial_by_best_time(new_annihilation->photon1_path);
  hit *hit2 = initial_by_best_time(new_annihilation->photon2_path);
  if (hit1 == new_annihilation->photon1_path.hits) {
    first_correct++;
  }
  if (hit2 == new_annihilation->photon2_path.hits) {
    first_correct++;
  }
  first_guessed += 2;
  prim_lor *new_prim_lor = (prim_lor *)malloc(sizeof(prim_lor));
  new_prim_lor->first_loc = hit1->location;
  new_prim_lor->second_loc = hit2->location;
  new_prim_lor->first_time = hit1->tof;
  new_prim_lor->second_time = hit2->tof;
  return new_prim_lor;
}
lor *create_lor(prim_lor *primitive_lor) {

  vec3d a = primitive_lor->first_loc;
  // printf("create_lor: a: \n");
  // vec_print(a, stdout);
  // printf("\n");
  // printf("b: \n");
  vec3d b = primitive_lor->second_loc;
  // vec_print(b, stdout);
  // printf("\n");
  vec3d c = vec_sub(a, b);
  vec3d geometric_center = vec_add(b, vec_scale(c, 0.5));
  // printf("geometic center: \n");
  // vec_print(geometric_center, stdout);
  // printf("\n");
  vec3d c_hat = vec_norm(c);
  double delta_t = -(primitive_lor->first_time - primitive_lor->second_time);
  vec3d displacement_from_center = vec_scale(c_hat, SPD_LGHT * delta_t * 0.5);
  vec3d annihilation_loc = vec_add(geometric_center, displacement_from_center);

  double transverse_uncert = sqrt(2 * SPC_UNC * SPC_UNC);
  double longtidudinal_uncert = sqrt(DETECTOR_THICKNESS / 24);

  lor *new = (lor *)malloc(sizeof(lor));
  new->center = annihilation_loc;
  new->dir = c_hat;
  new->long_uncert = longtidudinal_uncert;
  new->transverse_uncert = transverse_uncert;

  return new;
}
void print_lor(lor *new_lor, FILE *output) {
  fwrite(&new_lor->center.x, sizeof(double), 1, output);
  fwrite(&new_lor->center.y, sizeof(double), 1, output);
  fwrite(&new_lor->center.z, sizeof(double), 1, output);
  fwrite(&new_lor->dir.x, sizeof(double), 1, output);
  fwrite(&new_lor->dir.y, sizeof(double), 1, output);
  fwrite(&new_lor->dir.z, sizeof(double), 1, output);
  fwrite(&new_lor->long_uncert, sizeof(double), 1, output);
  fwrite(&new_lor->transverse_uncert, sizeof(double), 1, output);
}

int main(int argc, char **argv) {
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  // defines the help function and how to call it (by using -h or --help)
  for (int i = 0; i < num_flags(argc, argv); i++) {
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./hgmt_lor_creator [TOPAS_file_location.phsp] "
             "[efficiency_table_location.csv] [LOR_output_location]\n");
      printf("-h: print this help\n");
      exit(0);
    }
  }
  // checks to make sure you have correct number of args
  if (num_args(argc, argv) != 3) {
    printf("Incorrect number of arguments, three arguments required.\n");
    printf("Use the -h command to get options.\n\n");
    exit(1);
  }

  // reads in efficiency table into 2D array called eff_by_ang
  printf("HGMT LOR Creator\n\nLoading in '%s' as efficiencies table...\n",
         args[1]);
  FILE *eff_table_file = fopen(args[1], "r");
  int eff_file_read = read_eff(eff_table_file);
  fclose(eff_table_file);
  printf("Done!\n\n");

  // opens up a .lor file to output each LOR into
  char *lor_file_loc = (char *)malloc(sizeof(char) * (strlen(args[2]) + 10));
  strcpy(lor_file_loc, args[2]);
  lor_file_loc = strcat(lor_file_loc, ".lor");
  FILE *lor_output = fopen(lor_file_loc, "wb");
  if (lor_output == NULL) {
    printf("Unable to open output file for writing\n");
    return 1;
  }

  FILE *phsp_file = fopen(args[0], "rb");

  printf("Constructing the hits...\n");
  first_event = read_event(phsp_file);
  annihilation *new_annihilation = read_annihilation(phsp_file);
  while (new_annihilation != NULL) {
    annihilations_occured++;
    if (new_annihilation->photon1_path.num_hits != 0 &&
        new_annihilation->photon2_path.num_hits != 0) {
      prim_lor *primitive_lor = create_prim_lor(new_annihilation);
      if (primitive_lor != NULL) {
        lor *new_lor = create_lor(primitive_lor);
        print_lor(new_lor, lor_output);
        lors_created++;
        free(new_lor);
      }
      free(primitive_lor);
    }
    free(new_annihilation->photon1_path.hits);
    free(new_annihilation->photon2_path.hits);
    free(new_annihilation);
    new_annihilation = read_annihilation(phsp_file);
  }
  printf("Done!\n\n");

  printf("Scatters that occurred: %u\n", num_scatters);
  printf("Scatters detected: %u\n", scatters_detected);
  printf("Lors created: %u\n", lors_created);
  printf("Annihilations ocurred: %u\n", annihilations_occured);
  printf("Lor creation efficiency: %lf\n",
         (double)lors_created / annihilations_occured);
  printf("First detected scatter correctness: %lf\n",
         (double)first_correct / first_guessed);
  return 0;
}
