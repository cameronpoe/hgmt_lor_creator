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
uint counter = 0;
double detector_positions[12] = {45, 50, 55, 60, 65, 70,
                                 75, 80, 85, 90, 95, 100}; // MUST BE SORTED

// global variables
#define NUM_CUTS 5
#define NUM_DEBUG_OPTIONS 5
// cuts are {occured, interacted with something, wasn't inpatient, detected,
// first or second detected}
char *cut_descriptions[] = {"occured", "interected with something", "detected",
                            "wasn't inpatient", "first scatter detected"};
uint cuts[NUM_CUTS] = {0};
int array[2][NUM_DEBUG_OPTIONS];
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
void print_int(int numb, FILE *output) {
  fwrite(&numb, sizeof(int), 1, output);
}
prim_lor *create_prim_lor(annihilation *new_annihilation) {
  // hit *hit1 = initial_by_best_order(new_annihilation->photon1_path,
  // time_FOM); hit *hit2 =
  // initial_by_best_order(new_annihilation->photon2_path, time_FOM);
  hit hit1 =
      new_annihilation
          ->hits[0]; // initial_by_best_time(new_annihilation->photon1_path);
  hit hit2 =
      new_annihilation
          ->hits[1]; // initial_by_best_time(new_annihilation->photon2_path);

  prim_lor *new_prim_lor = (prim_lor *)malloc(sizeof(prim_lor));
  new_prim_lor->hit1 = hit1;
  new_prim_lor->hit2 = hit2;
  return new_prim_lor;
}

lor *create_lor(prim_lor *primitive_lor) {

  vec3d a = primitive_lor->hit1.position;
  vec3d b = primitive_lor->hit2.position;
  vec3d c = vec_sub(a, b);
  vec3d geometric_center = vec_add(b, vec_scale(c, 0.5));
  vec3d c_hat = vec_norm(c);
  double delta_t = -(primitive_lor->hit1.tof - primitive_lor->hit2.tof);
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
double impact_parameter(vec3d loc1, vec3d loc2, double tof1, double tof2,
                        vec3d true_center) {
  vec3d c = vec_sub(loc1, loc2);
  vec3d geometric_center = vec_add(loc2, vec_scale(c, 0.5));
  vec3d c_hat = vec_norm(c);
  double delta_t = tof2 - tof1;
  vec3d displacement_from_center = vec_scale(c_hat, SPD_LGHT * delta_t * 0.5);
  vec3d estimated_loc = vec_add(geometric_center, displacement_from_center);
  return vec_mag(vec_rejection(vec_sub(estimated_loc, true_center), c));
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
int get_detector(vec3d position) {
  double rad_dist = radial_dist(position);
  for (int i = 0; i < sizeof(detector_positions) / sizeof(double); i++) {
    if (rad_dist > detector_positions[i] &&
        rad_dist < detector_positions[i] + DETECTOR_THICKNESS) {
      return i;
    }
  }
  return -1;
}
bool is_detected(event *single_event) {
  if (single_event->detector_id != -1 &&
      drand48() < linear_interpolation(eff_by_energy, E_min, E_max,
                                       single_event->energy)) {
    return true;
  }
  return false;
}
hit *event_to_hit(event *single_event) {
  vec_mag(three_vec(single_event->position.x, single_event->position.y, 0.0));
  vec3d z_hat = three_vec(0.0, 0.0, 1.0);
  vec3d circ_hat = vec_norm(vec_cross(z_hat, single_event->position));
  vec3d offset = vec_add(vec_scale(z_hat, gaussian(SPC_UNC, 30)),
                         vec_scale(circ_hat, gaussian(SPC_UNC, 30)));
  hit *new_hit = (hit *)malloc(sizeof(hit));
  new_hit->source = single_event;
  new_hit->position = single_event->position;
  new_hit->tof = single_event->tof + gaussian(TIME_UNC, 30);
  double rad_dist = radial_dist(new_hit->position);
  if (DETECTOR_SEGMENTATION) {
    // we move the radial component to the midpoint of the detector which it hit
    new_hit->position = radial_scale(
        new_hit->position, (detector_positions[single_event->detector_id] +
                            DETECTOR_THICKNESS / 2) /
                               rad_dist);
  } else {
    vec3d r_hat =
        vec_scale(three_vec(new_hit->position.x, new_hit->position.y, 0),
                  1.0 / radial_dist(new_hit->position));
    offset = vec_add(offset, vec_scale(r_hat, gaussian(RAD_UNC, 30)));
  }
  new_hit->position = vec_add(new_hit->position, offset);
  return new_hit;
}
vec3d read_vec3d(FILE *source) {
  vec3d vec;
  fread(&vec.x, sizeof(double), 1, source);
  fread(&vec.y, sizeof(double), 1, source);
  fread(&vec.z, sizeof(double), 1, source);
  return vec;
}
bool read_annihilation(annihilation *annihilation_pointer, FILE *source) {
  if (!fread(&annihilation_pointer->time, sizeof(double), 1, source))
    return false;
  annihilation_pointer->origin = read_vec3d(source);
  annihilation_pointer->center = read_vec3d(source);
  uint num_events;
  uint num_primary1 = 0;
  uint num_primary2 = 0;
  fread(&num_events, sizeof(uint), 1, source);
  annihilation_pointer->events = (event *)malloc(sizeof(event) * num_events);
  uint num_hits = 0;
  uint num_primary1_hits = 0;
  uint num_primary2_hits = 0;
  for (int i = 0; i < num_events; i++) {
    event *next_event = &annihilation_pointer->events[i];
    fread(&next_event->tof, sizeof(double), 1, source);
    fread(&next_event->energy, sizeof(double), 1, source);
    next_event->position = read_vec3d(source);
    next_event->direction = read_vec3d(source);
    fread(&next_event->primary, sizeof(uint), 1, source);
    next_event->detector_id = get_detector(next_event->position);
    if (next_event->primary == 1) {
      next_event->number = num_primary1;
      num_primary1++;
    } else if (next_event->primary == 2) {
      next_event->number = num_primary2;
      num_primary2++;
    }
    if (is_detected(next_event)) {
      num_hits++;
      next_event->detected = true;
      if (next_event->primary == 1)
        num_primary1_hits++;
      else if (next_event->primary == 2)
        num_primary2_hits++;
    } else
      next_event->detected = false;
  }
  annihilation_pointer->hits = (hit *)malloc(sizeof(hit) * num_hits);
  annihilation_pointer->photon1_path.events =
      (event **)malloc(sizeof(event *) * num_primary1);
  annihilation_pointer->photon2_path.events =
      (event **)malloc(sizeof(event *) * num_primary2);
  annihilation_pointer->photon1_path.hits =
      (hit **)malloc(sizeof(hit *) * num_primary1_hits);
  annihilation_pointer->photon2_path.hits =
      (hit **)malloc(sizeof(hit *) * num_primary2_hits);
  num_primary1 = 0;
  num_primary2 = 0;
  num_hits = 0;
  num_primary1_hits = 0;
  num_primary2_hits = 0;
  for (int i = 0; i < num_events; i++) {
    event *current_event = &annihilation_pointer->events[i];
    if (annihilation_pointer->events[i].primary == 1) {
      annihilation_pointer->photon1_path.events[num_primary1] = current_event;
      num_primary1++;
    } else if (annihilation_pointer->events[i].primary == 2) {
      annihilation_pointer->photon2_path.events[num_primary2] = current_event;
      num_primary2++;
    }
    if (current_event->detected) {
      hit *detector_hit = event_to_hit(current_event);
      annihilation_pointer->hits[num_hits] = *detector_hit;
      free(detector_hit);
      if (current_event->primary == 1) {
        annihilation_pointer->photon1_path.hits[num_primary1_hits] =
            &annihilation_pointer->hits[num_hits];
        num_primary1_hits++;
      } else if (current_event->primary == 2) {
        annihilation_pointer->photon2_path.hits[num_primary2_hits] =
            &annihilation_pointer->hits[num_hits];
        num_primary2_hits++;
      }
      num_hits++;
    }
  }
  annihilation_pointer->num_events = num_events;
  annihilation_pointer->num_hits = num_hits;
  annihilation_pointer->photon1_path.num_events = num_primary1;
  annihilation_pointer->photon2_path.num_events = num_primary2;
  annihilation_pointer->photon1_path.num_hits = num_primary1_hits;
  annihilation_pointer->photon2_path.num_hits = num_primary2_hits;
  return true;
}

int compare_hits(const void *hit1, const void *hit2) {
  return ((hit *)hit1)->tof > ((hit *)hit2)->tof;
}
void free_annihilation(annihilation *annihilation_pointer) {
  free(annihilation_pointer->events);
  free(annihilation_pointer->hits);
  free(annihilation_pointer->photon1_path.events);
  free(annihilation_pointer->photon1_path.hits);
  free(annihilation_pointer->photon2_path.events);
  free(annihilation_pointer->photon2_path.hits);
}
void print_path(photon_path *path) {
  for (int i = 0; i < path->num_events; i++)
    // format: x,y,z, energy deposit, detected
    fprintf(visualization, "%lf %lf %lf %lf %d \n", path->events[i]->position.x,
            path->events[i]->position.y, path->events[i]->position.z,
            path->events[i]->energy, path->events[i]->detected ? 1 : 0);
}
void print_annihilation(annihilation *new_annihilation) {
  fprintf(visualization, "%lf %lf %lf \n\n", new_annihilation->center.x,
          new_annihilation->center.y, new_annihilation->center.z);
  print_path(&new_annihilation->photon1_path);
  fprintf(visualization, "\n");
  print_path(&new_annihilation->photon2_path);
  fprintf(visualization, "\n\n");
}
// provide debug statistics
int debug_path(photon_path *path) {
  if (path->num_events == 0) {
    cuts[0]++;
    return 0;
  }
  // getting all the important statistics
  if (debug_options[1])
    print_double(path->events[0]->detector_id, debug[1]);
  // figure out which cut the photon got to, format is: if (not cut n) cut=n-1
  int cut;
  if (path->num_hits == 0)
    cut = 1;
  else if (path->events[0]->detector_id == -1)
    cut = 2;
  else if (path->hits[0]->source != path->events[0])
    cut = 3;
  else
    cut = 4;
  cuts[cut]++;
  return cut;
}
int debug_annihilation(annihilation *new_annihilation) {
  num_scatters += new_annihilation->num_events;
  num_hits += new_annihilation->num_hits;
  if (debug_options[0])
    for (int j = 0; j < new_annihilation->num_events; j++)
      print_double(new_annihilation->events[j].detector_id, debug[0]);

  // fprintf(visualization, "%i\n", num_scatters);
  if (vis_events > 0) {
    print_annihilation(new_annihilation);
    vis_events--;
  }
  int cut1 = debug_path(&new_annihilation->photon1_path);
  int cut2 = debug_path(&new_annihilation->photon2_path);
  int cut = MIN(cut1, cut2);
  dual_cuts[cut]++;

  if (debug_options[4] && cut >= 1) {
    array[0][cut1]++;
    array[1][cut2]++;
    for (int i = 0; i < MIN(new_annihilation->photon1_path.num_events, 4); i++)
      for (int j = 0; j < MIN(new_annihilation->photon2_path.num_events, 4);
           j++) {
        vec3d true_center = new_annihilation->center;
        vec3d loc1 = new_annihilation->photon1_path.events[i]->position;
        vec3d loc2 = new_annihilation->photon2_path.events[j]->position;
        double tof1 = new_annihilation->photon1_path.events[i]->tof;
        double tof2 = new_annihilation->photon2_path.events[j]->tof;
        print_int(i + 1, debug[4]);
        print_int(j + 1, debug[4]);
        print_double(impact_parameter(loc1, loc2, tof1, tof2, true_center),
                     debug[4]);
      }
  }
  return cut;
}
void debug_lor(lor *new_lor, vec3d truecenter) {
  if (debug_options[2]) {
    print_double(vec_mag(vec_rejection(vec_sub(new_lor->center, truecenter),
                                       new_lor->dir)),
                 debug[2]);
  }
  if (debug_options[3]) {
    print_double(vec_mag(vec_projection(vec_sub(new_lor->center, truecenter),
                                        new_lor->dir)),
                 debug[3]);
  }
}
int main(int argc, char **argv) {
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  // handling all flags and arguments
  for (int i = 0; i < num_flags(argc, argv); i++) {
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./hgmt_lor_creator [TOPAS_file_position.phsp] "
             "[efficiency_table_position.csv] [output_directory]\n");
      printf("-h: print this help\n");
      printf("-d: run in debug mode, do not write to lor file\n");
      printf("-v#: visualize # events\n");
      printf("-e#: run with debug option #\n");
      printf("\t0: histogram of detector vs number of scatters\n");
      printf("\t1: histogram of detector vs number of first scatters\n");
      printf("\t2: lor reconstruction error to real center (transverse)\n");
      printf("\t3: lor reconstruction error to real center (longitudinal)\n");
      printf("\t4: Henry Plot (errors truth study)\n");
      exit(0);
    } else if (strcmp(flags[i], "-d") == 0) {
      printf("running in debug mode, won't write to a lor file\n");
      writing_to_lor = false;
    } else if (strncmp(flags[i], "-e", 2) == 0) {
      uint debug_option;
      sscanf(flags[i], "-e%u", &debug_option);
      debug_options[debug_option] = true;
    } else if (strncmp(flags[i], "-v", 2) == 0) {
      sscanf(flags[i], "-v%u", &vis_events);
      printf("outputting data to visualize %u events\n", vis_events);
    }
  }

  // checks to make sure you have correct number of args
  if (num_args(argc, argv) != 3) {
    printf("Incorrect number of arguments, three arguments required.\n");
    printf("Use the -h command to get options.\n\n");
    exit(1);
  }
  // opens files for debug output
  printf("running with debug options:"); // Output: 42
  for (int i = 0; i < NUM_DEBUG_OPTIONS; i++)
    if (debug_options[i]) {
      printf(" %i", i);
      char *filename;
      asprintf(&filename, "%sdebug%d.data", args[2], i);
      debug[i] = fopen(filename, "wb");
      free(filename);
    }
  printf("\n");
  // reads in efficiency table into 2D array called eff_by_ang
  printf("HGMT LOR Creator\n\nLoading in '%s' as efficiencies table...\n",
         args[1]);
  FILE *eff_table_file = fopen(args[1], "r");
  read_eff(eff_table_file);
  fclose(eff_table_file);

  // opens up a .lor file to output each LOR into
  FILE *lor_output = NULL;
  if (writing_to_lor) {
    printf("Unable to open output file for writing\n");
    char *lor_file_loc;
    asprintf(&lor_file_loc, "%sHGMTDerenzo.lor", args[2]);
    lor_output = fopen(lor_file_loc, "wb");
    free(lor_file_loc);
  }
  if (vis_events) {
    char *filename;
    asprintf(&filename, "%svisualization.data", args[2]);
    visualization = fopen(filename, "w");
    free(filename);
  }
  FILE *phsp_file = fopen(args[0], "rb");
  printf("Loading in '%s' as the phsp file...\n", args[0]);

  printf("Constructing the hits...\n\n");
  annihilation new_annihilation;
  bool worked = read_annihilation(&new_annihilation, phsp_file);
  while (worked) {
    if (debug_annihilation(&new_annihilation) >= 2) {
      prim_lor *primitive_lor = create_prim_lor(&new_annihilation);
      lor *new_lor = create_lor(primitive_lor);
      if (writing_to_lor)
        print_lor(new_lor, lor_output);
      debug_lor(new_lor, new_annihilation.center);
      free(primitive_lor);
      free(new_lor);
    }
    free_annihilation(&new_annihilation);
    worked = read_annihilation(&new_annihilation, phsp_file);
  }
  // fixing cuts formating to be cumulative
  for (int i = NUM_CUTS - 2; i >= 0; i--) {
    cuts[i] += cuts[i + 1];
    dual_cuts[i] += dual_cuts[i + 1];
  }
  printf("total annihilations: %u\n", dual_cuts[0]);
  printf("total scatters: %u\n", num_scatters);
  printf("total hits: %u\n\n", num_hits);
  printf(
      "(DUAL)CUT 'N': 'num' 'percent passing' 'cumulative percent passing'\n");
  for (int i = 1; i < NUM_CUTS; i++)
    printf("%u: %s\n", i, cut_descriptions[i]);
  printf("\n");
  for (int i = 1; i < NUM_CUTS; i++)
    printf("CUT %u: %u %lf %lf\n", i, cuts[i], (double)cuts[i] / cuts[i - 1],
           (double)cuts[i] / cuts[0]);
  printf("\n");
  for (int i = 1; i < NUM_CUTS; i++)
    printf("DUALCUT %u: %u %lf %lf\n", i, dual_cuts[i],
           (double)dual_cuts[i] / dual_cuts[i - 1],
           (double)dual_cuts[i] / dual_cuts[0]);
  printf("\nThe 2D array is:\n");
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < NUM_DEBUG_OPTIONS; j++) {
      printf("%u ", array[i][j]);
    }
    printf("\n"); // Move to the next row
  }
  printf("%f \n", (double)array[0][2] / (array[0][2] + array[0][1]));
  printf("%f \n", (double)array[1][2] / (array[1][2] + array[1][1]));

  // closing stuff out

  for (int i = 0; i < NUM_DEBUG_OPTIONS; i++)
    if (debug_options[i])
      fclose(debug[i]);
  if (visualization != NULL)
    fclose(visualization);
  return 0;
}
