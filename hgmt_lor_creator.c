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

// params
bool writing_to_lor = true;
uint vis_events = 0;
double detector_locations[6] = {45.0, 50, 55, 60, 65, 70}; // MUST BE SORTED

// global variables
#define NUM_CUTS 5
// cuts are {occured, interacted with something, wasn't inpatient, detected,
// first or second detected}
char *cut_descriptions[] = {"occured", "interected with something",
                            "wasn't inpatient", "detected",
                            "first or second detected"};
uint cuts[NUM_CUTS] = {0};
// dual cuts are the same but require both to happen in an annihilation
uint dual_cuts[NUM_CUTS] = {0};
uint num_scatters = 0;
uint num_hits = 0;
event *first_event;
double eff_by_energy[COLS];
double E_max = 520.0;
double E_min = 0.0;
bool debug_options[NUM_DEBUG_OPTIONS];
FILE *debug[NUM_DEBUG_OPTIONS];
FILE *visualization;
void print_lor(lor *new_lor, FILE *output) {
  fwrite(new_lor, sizeof(lor), 1, output);
}

void print_double(double numb, FILE *output) {
  fwrite(&numb, sizeof(double), 1, output);
}
prim_lor *create_prim_lor(annihilation *new_annihilation) {
  // hit *hit1 = initial_by_best_order(new_annihilation->photon1_path,
  // time_FOM); hit *hit2 =
  // initial_by_best_order(new_annihilation->photon2_path, time_FOM);
  hit *hit1 = initial_by_least_radial(new_annihilation->photon1_path);
  hit *hit2 = initial_by_least_radial(new_annihilation->photon2_path);

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
void read_eff(FILE *source) {
  // loops through all the entries in a row
  for (int i = 0; i < COLS; i++) {
    fscanf(source, "%lf,", &eff_by_energy[i]);
  }
}
// gets the detector an event happened in. return -1 if it didn't happen in a
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
  return new_event;
}

uint is_detected(event *single_event) {
  if (single_event->detector_id != -1 &&
      drand48() < linear_interpolation(eff_by_energy, E_min, E_max,
                                       single_event->energy_deposit)) {
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
  new_hit->source = single_event;
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
    offset = vec_add(offset, vec_scale(r_hat, gaussian(RAD_UNC, 30)));
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
  photon->events = (event *)malloc(sizeof(event) * photon->num_events);
  for (int i = 0; i < photon->num_events; i++) {
    photon->events[i] = *path_perfect->data;
    if (i != photon->num_events - 1) {
      path_perfect = path_perfect->down;
    }
  }
  // determining which hits are detected
  int *detected = (int *)calloc(photon->num_events, sizeof(int));
  for (int i = 0; i < photon->num_events; i++) {
    if (is_detected(&photon->events[i])) {
      detected[i] = 1;
      photon->num_hits++;
    }
  }
  // constructing list of hits
  photon->hits = (hit *)malloc(sizeof(hit) * photon->num_hits);
  for (int i = 0, j = 0; i < photon->num_events; i++) {
    if (detected[i]) {
      hit *detector_hit = event_to_hit(&photon->events[i]);
      photon->hits[j] = *detector_hit;
      free(detector_hit);
      j++;
    }
  }
  free(detected);
  // freeing all the redundant data
  wipe_list(path_perfect);
  first_event = new_event;
  // sorting by detected time
  if (photon->num_hits > 0)
    qsort(photon->hits, photon->num_hits, sizeof(hit), compare_hits);
  // returning
  return photon;
}
annihilation *read_annihilation(FILE *source) {
  // this is so complicated because we only get photon paths with trackid 2 or 3
  // so we must filter out all the other junk
  if (first_event == NULL) {
    return NULL;
  }
  annihilation *new_annihilation = (annihilation *)malloc(sizeof(annihilation));
  int event_id = first_event->event_id;
  // printf("%i\n", event_id);
  new_annihilation->photon1_path = read_photon_path(source);
  if (first_event->event_id == event_id) {
    new_annihilation->photon2_path = read_photon_path(source);
    vec3d a = new_annihilation->photon1_path->events->location;
    vec3d b = new_annihilation->photon2_path->events->location;
    vec3d c = vec_sub(a, b);
    vec3d center = vec_add(b, vec_scale(c, 0.5));
    vec3d c_hat = vec_norm(c);
    double delta_t = -(new_annihilation->photon1_path->events->tof -
                       new_annihilation->photon2_path->events->tof);
    vec3d displacement_from_center = vec_scale(c_hat, SPD_LGHT * delta_t * 0.5);
    vec3d annihilation_loc = vec_add(center, displacement_from_center);
    new_annihilation->center = annihilation_loc;
  } else
    new_annihilation->photon2_path = NULL;
  return new_annihilation;
}
void free_annihilation(annihilation *new_annihilation) {
  if (new_annihilation->photon1_path != NULL) {
    free(new_annihilation->photon1_path->hits);
    free(new_annihilation->photon1_path);
    free(new_annihilation->photon1_path->events);
  }
  if (new_annihilation->photon2_path != NULL) {
    free(new_annihilation->photon2_path->hits);
    free(new_annihilation->photon2_path);
    free(new_annihilation->photon2_path->events);
  }
  free(new_annihilation);
}
void print_path(photon_path *path) {
  bool *detected = calloc(path->num_events, sizeof(bool));
  for (int i = 0; i < path->num_hits; i++) {
    int index = path->hits[i].source - path->events;
    detected[index] = 1;
  }
  for (int i = 0; i < path->num_events; i++) {
    // format: x,y,z, energy deposit, detected
    fprintf(visualization, "%lf %lf %lf %lf %d \n", path->events[i].location.x,
            path->events[i].location.y, path->events[i].location.z,
            path->events[i].energy_deposit, detected[i] ? 1 : 0);
  }
}
void print_annihilation(annihilation *new_annihilation) {
  fprintf(visualization, "%lf %lf %lf \n\n", new_annihilation->center.x,
          new_annihilation->center.y, new_annihilation->center.z);
  print_path(new_annihilation->photon1_path);
  fprintf(visualization, "\n");
  print_path(new_annihilation->photon2_path);
  fprintf(visualization, "\n\n");
}
// provide debug statistics
int debug_path(photon_path *path) {
  if (path == NULL) {
    cuts[0]++;
    return 0;
  }
  bool *detected = calloc(path->num_events, sizeof(bool));
  for (int i = 0; i < path->num_hits; i++) {
    int index = path->hits[i].source - path->events;
    detected[index] = 1;
  }
  // getting all the important statistics
  num_scatters += path->num_events;
  num_hits += path->num_hits;
  if (debug_options[0]) {
    for (int j = 0; j < path->num_events; j++)
      print_double(path->events[j].detector_id, debug[0]);
  }
  if (debug_options[1])
    print_double(path->events->detector_id, debug[1]);
  // figure out which cut the photon got to, format is: if (not cut n) cut=n-1
  int cut;
  if (path->events->detector_id == -1)
    cut = 1;
  else if (path->num_hits == 0)
    cut = 2;
  else if (path->hits->source != path->events)
    cut = 3;
  else
    cut = 4;
  cuts[cut]++;
  free(detected);
  return cut;
}
void debug_annihilation(annihilation *new_annihilation) {
  // fprintf(visualization, "%i\n", num_scatters);
  if (vis_events > 0 && new_annihilation->photon1_path != NULL &&
      new_annihilation->photon2_path != NULL) {
    print_annihilation(new_annihilation);
    vis_events--;
  }
  int cut1 = debug_path(new_annihilation->photon1_path);
  int cut2 = debug_path(new_annihilation->photon2_path);
  dual_cuts[MIN(cut1, cut2)]++;
}
int main(int argc, char **argv) {
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  // handling all flags and arguments
  for (int i = 0; i < num_flags(argc, argv); i++) {
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./hgmt_lor_creator [TOPAS_file_location.phsp] "
             "[efficiency_table_location.csv] [LOR_output_location]\n");
      printf("-h: print this help\n");
      printf("-d: run in debug mode, do not write to lor file\n");
      printf("-v#: visualize # events\n");
      printf("-e#: run with debug option #\n");
      printf("\t0: histogram of detector vs number of scatters\n");
      printf("\t1: histogram of detector vs number of first scatters\n");
      exit(0);
    } else if (strcmp(flags[i], "-d") == 0) {
      printf("running in debug mode, won't write to a lor file\n");
      writing_to_lor = false;
    } else if (strncmp(flags[i], "-e", 2) == 0) {
      uint debug_option;
      sscanf(flags[i], "-e%u", &debug_option);
      debug_options[debug_option] = true;
      printf("running with debug option: %u\n", debug_option); // Output: 42
    } else if (strncmp(flags[i], "-v", 2) == 0) {
      sscanf(flags[i], "-v%u", &vis_events);
      printf("outputting data to visualize %u events\n", vis_events);
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
  read_eff(eff_table_file);
  fclose(eff_table_file);

  // opens up a .lor file to output each LOR into
  FILE *lor_output = NULL;
  if (writing_to_lor) {
    printf("Unable to open output file for writing1\n");
    char *lor_file_loc;
    asprintf(&lor_file_loc, "data/%s.lor", args[2]);
    lor_output = fopen(lor_file_loc, "wb");
    free(lor_file_loc);
  }
  for (int i = 0; i < NUM_DEBUG_OPTIONS; i++) {
    if (debug_options[i]) {
      char *filename;
      asprintf(&filename, "data/debug%d.data", i);
      debug[i] = fopen(filename, "wb");
      free(filename);
    }
  }
  if (vis_events)
    visualization = fopen("data/visualization.data", "w");
  FILE *phsp_file = fopen(args[0], "rb");
  printf("Loading in '%s' as the phsp file...\n", args[0]);

  printf("Constructing the hits...\n\n");
  first_event = read_event(phsp_file);
  annihilation *new_annihilation = read_annihilation(phsp_file);
  while (new_annihilation != NULL) {
    if (new_annihilation->photon2_path != NULL &&
        new_annihilation->photon1_path->num_hits != 0 &&
        new_annihilation->photon2_path->num_hits != 0) {
      // printf("%i \n", new_annihilation->photon1_path.num_hits);
      // printf("%i \n\n", new_annihilation->photon2_path.num_hits);
      prim_lor *primitive_lor = create_prim_lor(new_annihilation);
      if (primitive_lor != NULL) {
        if (writing_to_lor) {
          lor *new_lor = create_lor(primitive_lor);
          print_lor(new_lor, lor_output);
          free(new_lor);
        }
      }
      free(primitive_lor);
    }
    debug_annihilation(new_annihilation);
    free_annihilation(new_annihilation);
    new_annihilation = read_annihilation(phsp_file);
  }
  // fixing cuts formating to be cumulative
  for (int i = NUM_CUTS - 2; i >= 0; i--) {
    cuts[i] += cuts[i + 1];
    dual_cuts[i] += dual_cuts[i + 1];
  }
  printf("total annihilations: %u\n", dual_cuts[0]);
  printf("total scatters: %u\n", num_scatters);
  printf("scatters detected: %u\n", num_hits);
  printf(
      "(DUAL)CUT 'N': 'num' 'percent passing' 'cumulative percent passing'\n");
  printf("\n");
  for (int i = 1; i < NUM_CUTS; i++) {
    printf("%u: %s\n", i, cut_descriptions[i]);
  }
  printf("\n");
  for (int i = 1; i < NUM_CUTS; i++) {
    printf("CUT %u: %u %lf %lf\n", i, cuts[i], (double)cuts[i] / cuts[i - 1],
           (double)cuts[i] / cuts[0]);
  }
  printf("\n");
  for (int i = 1; i < NUM_CUTS; i++) {
    printf("DUALCUT %u: %u %lf %lf\n", i, dual_cuts[i],
           (double)dual_cuts[i] / dual_cuts[i - 1],
           (double)dual_cuts[i] / dual_cuts[0]);
  }
  return 0;
}
