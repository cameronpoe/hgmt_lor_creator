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
uint num_scatters = 0;
uint scatters_detected = 0;
uint annihilations_occured = 0;
uint lors_created = 0;
uint first_correct = 0;
uint events_occurred = 0;
uint paths_created = 0;
uint first_detected = 0;
uint first_guessed = 0;
uint first_scatter_correct = 0;
uint first_scatter_in_detector = 0;
uint error_debug = 0;
uint writing_to_lor = 1;
event *first_event;
double eff_by_energy[COLS];
double E_max = 520.0;
double E_min = 0.0;
FILE *debug;
double detector_locations[12] = {
    45,    47.54, 50.08, 52.62, 55.16, 57.7, 60.24,
    62.78, 65.32, 67.86, 70.4,  72.94}; // inner radii of detectors, MUST BE
//  SORTED
// double detector_locations[6] = {45.0, 50, 55, 60, 65, 70};
void print_lor(lor *new_lor, FILE *output) {
  fwrite(new_lor, sizeof(lor), 1, output);
}

void print_double(double numb, FILE *output) {
  fwrite(&numb, sizeof(double), 1, output);
}
prim_lor *create_prim_lor(annihilation *new_annihilation) {
  hit *hit1 = initial_by_best_order(new_annihilation->photon1_path, time_FOM);
  hit *hit2 = initial_by_best_order(new_annihilation->photon2_path, time_FOM);
  // hit *hit1 = initial_by_best_time(new_annihilation->photon1_path);
  // hit *hit2 = initial_by_best_time(new_annihilation->photon2_path);
  if (hit1->first) {
    first_correct++;
    if (new_annihilation->photon1_path->has_first) {
      first_scatter_correct++;
    }
  } else if (error_debug == 2) {
    print_double(vec_dist(hit1->location,
                          new_annihilation->photon1_path->hits[0].location),
                 debug);
  }
  if (hit2->first) {
    first_correct++;
    if (new_annihilation->photon2_path->has_first) {
      first_scatter_correct++;
    }
  } else if (error_debug == 2) {
    print_double(vec_dist(hit2->location,
                          new_annihilation->photon2_path->hits[0].location),
                 debug);
  }
  if (error_debug == 1) {
    print_double(
        vec_dist(new_annihilation->photon1_path->debug_path[0].location,
                 hit1->location),
        debug);
    print_double(
        vec_dist(new_annihilation->photon2_path->debug_path[0].location,
                 hit2->location),
        debug);
  }
  first_guessed += 2;
  prim_lor *new_prim_lor = (prim_lor *)malloc(sizeof(prim_lor));
  new_prim_lor->hit1 = hit1;
  new_prim_lor->hit2 = hit2;
  return new_prim_lor;
}
lor *create_lor(prim_lor *primitive_lor) {

  vec3d a = primitive_lor->hit1->location;
  // printf("create_lor: a: \n");
  // vec_print(a, stdout);
  // printf("\n");
  // printf("b: \n");
  vec3d b = primitive_lor->hit2->location;
  // vec_print(b, stdout);
  // printf("\n");
  vec3d c = vec_sub(a, b);
  vec3d geometric_center = vec_add(b, vec_scale(c, 0.5));
  // printf("geometic center: \n");
  // vec_print(geometric_center, stdout);
  // printf("\n");
  vec3d c_hat = vec_norm(c);
  double delta_t = -(primitive_lor->hit1->tof - primitive_lor->hit2->tof);
  vec3d displacement_from_center = vec_scale(c_hat, SPD_LGHT * delta_t * 0.5);
  vec3d annihilation_loc = vec_add(geometric_center, displacement_from_center);

  double transverse_uncert = sqrt(2 * SPC_UNC * SPC_UNC);
  double longtidudinal_uncert = sqrt(TIME_UNC * TIME_UNC * 2);

  lor *new = (lor *)malloc(sizeof(lor));
  new->center = annihilation_loc;
  new->dir = c_hat;
  new->long_uncert = longtidudinal_uncert;
  new->transverse_uncert = transverse_uncert;

  return new;
}
double linear_interpolation(double nums[COLS], double min, double max,
                            double value) {
  double i = (COLS - 1) * (value - min) / (max - min);
  int i_l = (int)i;
  int i_r = i_l + 1;
  double i_space = i - i_l;
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
// gets the detector an event happened in. return -1 if it didn't happen in
// detector
int get_detector(vec3d location) {
  double rad_dist = radial_dist(location);
  for (int i = 0; i < sizeof(detector_locations) / sizeof(double); i++) {
    if (rad_dist > detector_locations[i] &&
        rad_dist < detector_locations[i] + DETECTOR_THICKNESS) {
      return i;
    }
  }
  return -1;
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
  int parent_id;
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
  worked += fread(&parent_id, sizeof(int), 1, source);
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
  new_event->parent_id = parent_id;
  new_event->track_id = track_id;
  new_event->detector_id = get_detector(new_event->location);
  num_scatters++;
  return new_event;
}

uint is_scatter_detected(event *single_event) {
  if (single_event->detector_id != -1 &&
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
    // we move the radial component to the midpoint of the detector which it hit
    new_hit->location = radial_scale(
        new_hit->location, (detector_locations[single_event->detector_id] +
                            DETECTOR_THICKNESS / 2) /
                               rad_dist);
  } else {
    vec3d r_hat =
        vec_scale(three_vec(new_hit->location.x, new_hit->location.y, 0),
                  1.0 / radial_dist(new_hit->location));
    offset = vec_add(offset, vec_scale(r_hat, gaussian(SPC_UNC, 30)));
  }
  new_hit->location = vec_add(new_hit->location, offset);
  return new_hit;
}
int compare_hits(const void *hit1, const void *hit2) {
  return ((hit *)hit1)->tof > ((hit *)hit2)->tof;
}
photon_path *read_photon_path(FILE *source) {
  // reading all the scatters within the time window
  // scatters are given backwards so it is difficult to figure stuff out
  llist *path_perfect = NULL; // error debug stuff
  photon_path *photon = (photon_path *)malloc(sizeof(photon_path));
  photon->num_events = 0;
  photon->num_hits = 0;
  event *new_event = first_event;
  // constructing list of events
  while (new_event != NULL && first_event->event_id == new_event->event_id &&
         first_event->parent_id == new_event->parent_id) {
    path_perfect = add_to_top(path_perfect, new_event);
    photon->num_events++;
    new_event = read_event(source);
  }
  photon->debug_path = (event *)malloc(sizeof(event) * photon->num_events);
  for (int i = 0; i < photon->num_events; i++) {
    photon->debug_path[i] = *path_perfect->data;
    if (i != photon->num_events - 1) {
      path_perfect = path_perfect->down;
    }
  }
  // determining which hits are detected
  int *detected = (int *)calloc(photon->num_events, sizeof(int));
  for (int i = 0; i < photon->num_events; i++) {
    if (is_scatter_detected(&photon->debug_path[i])) {
      detected[i] = 1;
      photon->num_hits++;
    }
    if (i == 0) {
      photon->has_first = *detected;
      if (photon->debug_path->detector_id != -1) {
        first_scatter_in_detector++;
      }
    }
  }
  // constructing list of hits
  photon->hits = (hit *)malloc(sizeof(hit) * photon->num_hits);
  int j = 0;
  for (int i = 0; i < photon->num_events; i++) {
    if (detected[i]) {
      hit *detector_hit = event_to_hit(&photon->debug_path[i]);
      photon->hits[j] = *detector_hit;
      photon->hits[j].first = j == 0;
      free(detector_hit);
      j++;
      if (error_debug == 9)
        print_double((double)photon->debug_path[i].detector_id, debug);
    }
  }
  free(detected);
  // freeing all the redundant data
  wipe_list(path_perfect);
  first_event = new_event;
  // sorting by detected time
  if (photon->num_hits > 0)
    qsort(photon->hits, photon->num_hits, sizeof(hit), compare_hits);
  // error debug stuff
  events_occurred++;
  if (photon->num_hits > 0) {
    paths_created++;
  }
  if (photon->has_first) {
    first_detected++;
  }
  if (error_debug == 3) {
    if (photon->num_events >= 2) {
      print_double(vec_dist(photon->debug_path[0].location,
                            photon->debug_path[1].location),
                   debug);
    }
  }
  if (error_debug == 4) {
    if (photon->num_events >= 2) {
      print_double(photon->debug_path[1].tof - photon->debug_path[0].tof,
                   debug);
    }
  }
  if (error_debug == 5) {
    if (photon->num_hits >= 2) {
      print_double(vec_dist(photon->hits[0].location, photon->hits[1].location),
                   debug);
    }
  }
  if (error_debug == 6) {
    if (photon->num_hits >= 2) {
      print_double(photon->hits[1].tof - photon->hits[0].tof, debug);
    }
  }
  if (error_debug == 7) {
    print_double(photon->debug_path[0].energy_deposit, debug);
  }
  if (error_debug == 8) {
    print_double(photon->num_hits, debug);
  }
  // returning
  return photon;
}
void skip_photon_path(FILE *source) {
  event *new_event = read_event(source);
  while (new_event != NULL && new_event->parent_id == first_event->parent_id &&
         new_event->event_id == first_event->event_id) {
    free(new_event);
    new_event = read_event(source);
  }
  free(first_event);
  first_event = new_event;
  return;
}
annihilation *read_annihilation(FILE *source) {
  // this is so complicated because we only get photon paths with trackid 2 or 3
  // so we must filter out all the other junk
  if (first_event == NULL) {
    return NULL;
  }
  annihilation *photon_pair = (annihilation *)malloc(sizeof(annihilation));
  int event_id = first_event->event_id;
  // printf("%i\n", event_id);
  photon_pair->photon1_path = read_photon_path(source);
  if (first_event->event_id == event_id) {
    photon_pair->photon2_path = read_photon_path(source);
  } else {
    photon_pair->photon2_path = NULL;
  }
  while (first_event != NULL && first_event->event_id == event_id) {
    // skip_photon_path(source);
  }
  // printf("\n");
  return photon_pair;
}
void free_annihilation(annihilation *new_annihilation) {
  if (new_annihilation->photon1_path != NULL) {
    free(new_annihilation->photon1_path->hits);
    free(new_annihilation->photon1_path);
    if (error_debug) {
      free(new_annihilation->photon1_path->debug_path);
    }
  }
  if (new_annihilation->photon2_path != NULL) {
    free(new_annihilation->photon2_path->hits);
    free(new_annihilation->photon2_path);
    if (error_debug) {
      free(new_annihilation->photon2_path->debug_path);
    }
  }
  free(new_annihilation);
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
      printf("-e1: run with compton error debug\n");
      printf("-e2: run with compton ordering error debug\n");
      exit(0);
    } else if (strcmp(flags[i], "-e1") == 0) {
      printf("running with compton error debug\n");
      error_debug = 1;
    } else if (strcmp(flags[i], "-e2") == 0) {
      printf("running with compton ordering error debug\n");
      error_debug = 2;
    } else if (strcmp(flags[i], "-e3") == 0) {
      printf("running with scatter distance\n");
      error_debug = 3;
    } else if (strcmp(flags[i], "-e4") == 0) {
      printf("running with scatter time difference\n");
      error_debug = 4;
    } else if (strcmp(flags[i], "-e5") == 0) {
      printf("running with hit distance\n");
      error_debug = 5;
    } else if (strcmp(flags[i], "-e6") == 0) {
      printf("running with hit time distance\n");
      error_debug = 6;
    } else if (strcmp(flags[i], "-e7") == 0) {
      printf("running with energy distribution\n");
      error_debug = 7;
    } else if (strcmp(flags[i], "-e8") == 0) {
      printf("running with chain length debug\n");
      error_debug = 8;
    } else if (strcmp(flags[i], "-e9") == 0) {
      printf("running with detector activity debug\n");
      error_debug = 9;
    } else if (strcmp(flags[i], "-d") == 0) {
      printf("running in debug mode, won't write to a lor file\n");
      writing_to_lor = 0;
    }
  }
  // checks to make sure you have correct number of args
  if (num_args(argc, argv) != 3 - !writing_to_lor) {
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
  FILE *lor_output = NULL;
  if (writing_to_lor) {
    char *lor_file_loc = (char *)malloc(sizeof(char) * (strlen(args[2]) + 10));
    strcpy(lor_file_loc, args[2]);
    lor_file_loc = strcat(lor_file_loc, ".lor");
    lor_output = fopen(lor_file_loc, "wb");
    if (lor_output == NULL) {
      printf("Unable to open output file for writing\n");
      return 1;
    }
  }
  if (error_debug) {
    debug = fopen("debug.data", "wb");
    if (debug == NULL) {
      printf("Unable to open debug file for writing\n");
      return 1;
    }
  }
  FILE *phsp_file = fopen(args[0], "rb");

  printf("Constructing the hits...\n");
  first_event = read_event(phsp_file);
  annihilation *new_annihilation = read_annihilation(phsp_file);
  while (new_annihilation != NULL) {
    annihilations_occured++;
    if (new_annihilation->photon2_path != NULL &&
        new_annihilation->photon1_path->num_hits != 0 &&
        new_annihilation->photon2_path->num_hits != 0) {
      // printf("%i \n", new_annihilation->photon1_path.num_hits);
      // printf("%i \n\n", new_annihilation->photon2_path.num_hits);
      prim_lor *primitive_lor = create_prim_lor(new_annihilation);
      if (primitive_lor != NULL) {
        lors_created++;
        if (writing_to_lor) {
          lor *new_lor = create_lor(primitive_lor);
          print_lor(new_lor, lor_output);
          free(new_lor);
        }
      }
      free(primitive_lor);
    }
    free_annihilation(new_annihilation);
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
  printf("Proportion of first scatters detected: %lf\n",
         (double)first_detected / events_occurred);
  printf("Proportion of first scatters occurred in a detector: %lf\n",
         (double)first_scatter_in_detector / events_occurred);
  printf("Proportion of first scatters detected in chains with atleast 1 "
         "scatter detected: %lf\n",
         (double)first_detected / paths_created);
  printf("First scatter correctness: %lf\n",
         (double)first_scatter_correct / first_guessed);
  return 0;
}
