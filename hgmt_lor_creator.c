// including standard files
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// including custom files
#include "helper_functions.h"
#include "hgmt_event.h"
#include "vector_ops.h"

// global variables
#define ROWS 31
#define COLS 31
#define ROWS_DIV 3 // deg/row
#define COLS_DIV 3 // deg/col
#define PI 3.141592653589
#define SPD_LGHT 29.9792458 // cm/ns
#define UNCERT_REP 30
#define SPC_UNC 0.1       // cm
#define TIME_UNC 0.042463 // ns, sigma (0.100 ns FWHM)
#define DIAGNOSTIC_TOTAL_GAMMAS_DETECTED 1
#define MAX_THREADS 7
double eff_by_ang[ROWS][COLS];
uint gammas_entered_scanner = 0;
uint binary = 0;
uint gamma_pairs_entered_scanner = 0;
uint diagnostic_single_gammas_detected = 0;
uint gamma_pairs_detected = 0;
uint gammas_detected_ips = 0;
uint gamma_pairs_detected_atleast_one_ips = 0;
uint gamma_pairs_detected_both_ips = 0;

double bilinear_interp(double arr[][COLS], int num_rows, int num_cols,
                       double i_full, double j_full) {

  int i = (int)i_full;
  int j = (int)j_full;

  int i_l = i;
  int i_r = i + 1;
  int j_l = j;
  int j_r = j + 1;

  double di_l = (double)i_l;
  double di_r = (double)i_r;
  double dj_l = (double)j_l;
  double dj_r = (double)j_r;

  double prefactor;

  double interpolated_val;

  // checks if at corners
  if ((i == 0 && j == 0) || (i == num_rows - 1 && j == 0) ||
      (i == 0 && j == num_cols - 1) ||
      (i == num_rows - 1 && j == num_cols - 1)) {
    interpolated_val = arr[i][j];
    // checks if on top  or bottom side
  } else if (i == 0 || i == num_rows - 1) {
    prefactor = 1.0 / (dj_r - dj_l);
    interpolated_val = prefactor * (arr[i][j_l] * (dj_r - j_full) +
                                    arr[i][j_r] * (j_full - dj_l));
    // checks if on left or right side
  } else if (j == 0 || j == num_cols - 1) {
    prefactor = 1.0 / (di_r - di_l);
    interpolated_val = prefactor * (arr[i_l][j] * (di_r - i_full) +
                                    arr[i_r][j] * (i_full - di_l));
    // handles point being in the center of the rectangle
  } else {
    prefactor = 1.0 / (di_r - di_l) * (dj_r - dj_l);
    double term_1 = (di_r - i_full) * (arr[i_l][j_l] * (dj_r - j_full) +
                                       arr[i_l][j_r] * (j_full - dj_l));
    double term_2 = (i_full - di_l) * (arr[i_r][j_l] * (dj_r - j_full) +
                                       arr[i_r][j_r] * (j_full - dj_l));
    interpolated_val = prefactor * (term_1 + term_2);
  }

  return interpolated_val;
}

void test_bl_interp() {

  // results of tests with a table are to the right of the printf statements
  printf("Checking corners:\n");
  printf("i=0.1, j=0.1 -> %lf\n",
         bilinear_interp(eff_by_ang, ROWS, COLS, 0.1, 0.1)); // correct (exact)
  printf("i=30.1, j=0.1 -> %lf\n",
         bilinear_interp(eff_by_ang, ROWS, COLS, 30.1, 0.1)); // correct (exact)
  printf("i=0.1, j=30.1 -> %lf\n",
         bilinear_interp(eff_by_ang, ROWS, COLS, 0.1, 30.1)); // correct (exact)
  printf(
      "i=30.1, j=30.1 -> %lf\n",
      bilinear_interp(eff_by_ang, ROWS, COLS, 30.1, 30.1)); // correct (exact)
  printf("Checking sides:\n");
  printf("i=0.1, j=5.48 -> %lf\n", bilinear_interp(eff_by_ang, ROWS, COLS, 0.1,
                                                   5.48)); // correct (in range)
  printf("i=2.241, j=0.1 -> %lf\n",
         bilinear_interp(eff_by_ang, ROWS, COLS, 2.241,
                         0.1)); // correct (in range)
  printf("i=30.1, j=5.48 -> %lf\n",
         bilinear_interp(eff_by_ang, ROWS, COLS, 30.1,
                         5.48)); // correct (in range)
  printf("i=2.241, j=30.1 -> %lf\n",
         bilinear_interp(eff_by_ang, ROWS, COLS, 2.241,
                         30.1)); // correct (in range)
  printf("Checking center:\n");
  printf("i=2.241, j=5.48 -> %lf\n",
         bilinear_interp(eff_by_ang, ROWS, COLS, 2.241,
                         5.48)); // correct (in range)

  return;
}

int read_eff(FILE *source) {

  // exits with code 1 if the source pointer isn't pointing to the file
  if (source == NULL) {
    return 1;
  }

  // loops through all the entries in a row
  for (int i = 0; i < ROWS; i++) {
    int j = 0;
    fscanf(source, "%lf", &eff_by_ang[i][j]);
    j++;
    for (; j < COLS; j++) {
      fscanf(source, ",%lf", &eff_by_ang[i][j]);
    }
  }

  return 0;
}

void test_eff_array(float eff_table[ROWS][COLS]) {

  for (int i = 0; i < ROWS; i++) {
    int j = 0;
    printf("%f", eff_table[i][j]);
    j++;
    for (; j < COLS; j++) {
      printf(", %f", eff_table[i][j]);
    }
    printf("\n");
  }

  return;
}

event *read_event(FILE *source) {

  double x;          // 1: Position X [cm]
  double y;          // 2: Position Y [cm]
  double z;          // 3: Position Z [cm]
  double dir_cos_x;  // 4: Direction Cosine X
  double dir_cos_y;  // 5: Direction Cosine Y
  double energy;     // 6: Energy [MeV]
  double weight;     // 7: Weight
  int particle;      // 8: Particle Type (in PDG Format)
  int dir_cos_neg;   // 9: Flag to tell if Third Direction Cosine is Negative (1
                     // means true)
  int first_in_hist; // 10: Flag to tell if this is the First Scored Particle
                     // from this History (1 means true)
  double tof;        // 11: Time of Flight [ns]
  int run_id;
  int event_id;
  int track_id;

  int worked = 0;
  if (binary) {
    worked += fread(&x, sizeof(double), 1, source);
    worked += fread(&y, sizeof(double), 1, source);
    worked += fread(&z, sizeof(double), 1, source);
    worked += fread(&dir_cos_x, sizeof(double), 1, source);
    worked += fread(&dir_cos_y, sizeof(double), 1, source);
    worked += fread(&energy, sizeof(double), 1, source);
    worked += fread(&weight, sizeof(double), 1, source);
    worked += fread(&particle, sizeof(int), 1, source);
    worked += fread(&dir_cos_neg, sizeof(uint), 1, source);
    worked += fread(&first_in_hist, sizeof(uint), 1, source);
    worked += fread(&tof, sizeof(double), 1, source);
    worked += fread(&run_id, sizeof(uint), 1, source);
    worked += fread(&event_id, sizeof(uint), 1, source);
    worked += fread(&track_id, sizeof(uint), 1, source);
  } else {
    // reads in position as long float
    worked = fscanf(source, "%lf", &x);
    worked = fscanf(source, "%lf", &y);
    worked = fscanf(source, "%lf", &z);
    // reads in cosines as long floats
    worked = fscanf(source, "%lf", &dir_cos_x);
    worked = fscanf(source, "%lf", &dir_cos_y);
    // adds energy and weight as long floats
    worked = fscanf(source, "%lf", &energy);
    worked = fscanf(source, "%lf", &weight);
    // adds particle ID as an int (all PIDs are ints per PDG)
    worked = fscanf(source, "%i", &particle);
    // adds flags as uints (since they're 0 or 1)
    worked = fscanf(source, "%u", &dir_cos_neg);
    worked = fscanf(source, "%u", &first_in_hist);
    // adds tof as long float
    worked = fscanf(source, "%lf", &tof);
    // adds identifying info as uints
    worked = fscanf(source, "%u", &run_id);
    worked = fscanf(source, "%u", &event_id);
    worked = fscanf(source, "%u", &track_id);
  }

  if (worked == EOF) {
    return NULL;
  }

  // make a new event to be passed out
  event *new_event = (event *)malloc(sizeof(event));
  if (new_event == NULL) {
    return NULL;
  }

  new_event->location = three_vec(x, y, z);
  new_event->dir_cos_x = dir_cos_x;
  new_event->dir_cos_y = dir_cos_y;
  new_event->energy = energy;
  new_event->weight = weight;
  new_event->particle = particle;
  new_event->dir_cos_neg = dir_cos_neg;
  new_event->first_in_hist = first_in_hist;
  new_event->tof = tof;
  new_event->run_id = run_id;
  new_event->event_id = event_id;
  new_event->track_id = track_id;

  if (new_event->particle == 22) {
    gammas_entered_scanner++;
  }

  return new_event;
}

event *get_event_pair(FILE *source, event *(*f)(FILE *)) {

  int found_pair = 1;

  event *first_event = f(source);
  event *second_event = f(source);

  if (first_event == NULL || second_event == NULL) {
    return NULL;
  }

  while (found_pair) {
    // are these from the same history and run
    if ((first_event->run_id != second_event->run_id) ||
        (first_event->event_id != second_event->event_id) ||
        (first_event->particle != 22) || (second_event->particle != 22)) {
      free(first_event);
      first_event = second_event;
      second_event = f(source);
      if (second_event == NULL) {
        return NULL;
      }
    } else {
      // they are!
      found_pair = 0;
    }
  }

  event *event_pair = (event *)malloc(2 * sizeof(event));

  // orders the events such that event_pair[0] arrived first and event_pair[1]
  // arrived second
  if (first_event->tof <= second_event->tof) {
    event_pair[0] = *first_event;
    event_pair[1] = *second_event;
  } else {
    event_pair[0] = *second_event;
    event_pair[1] = *first_event;
  }

  free(first_event);
  free(second_event);

  // increments number of pairs that entered the detector
  gamma_pairs_entered_scanner++;

  return event_pair;
}

void print_event_pair(event *pair) {

  for (int i = 0; i < 2; i++) {
    event single_event = pair[i];
    printf("(x,y,z): (%lf, %lf, %lf) cm, TOF: %lf, Run ID: %u, History ID: %u, "
           "Particle in History ID: %u\n",
           single_event.location.x, single_event.location.y,
           single_event.location.z, single_event.tof, single_event.run_id,
           single_event.event_id, single_event.track_id);
  }

  return;
}

int is_gamma_detected(event *single_event) {

  int is_detected = 0;

  double alpha = single_event->dir_cos_x;
  double beta = single_event->dir_cos_y;
  int sign_of_gamma;
  double rad_to_deg = 180.0 / PI;

  double phi;
  double theta;

  // converts the flag for the sign of the third direction cosine into a factor
  // of 1 or -1
  if (single_event->dir_cos_neg) {
    sign_of_gamma = 1;
  } else {
    sign_of_gamma = -1;
  }

  // dir is the momentum unit vector of the gamma; has to check that the square
  // root is not imaginary (happens with floats)
  double operand = 1.0 - (alpha * alpha) + (beta * beta);
  if (operand < 0.0) {
    operand = 0.0;
  }
  vec3d dir = three_vec(alpha, beta, -sqrt(operand) * (float)sign_of_gamma);

  // constructing the normalized normal vector to the detector
  vec3d normal = single_event->location;
  normal.z = 0;
  normal = vec_norm(normal);

  // uses vector ops to find the angle between various vectors
  theta = vec_angle(normal, dir) * rad_to_deg;
  phi = vec_angle(three_vec(0, 0, 1), vec_rejection(dir, normal)) * rad_to_deg;

  // converts theta and phi into their respsective indices on the table
  int theta_ind = (int)(theta / COLS_DIV);
  if (phi > 90.0) {
    phi = phi - 90.0;
  }
  int phi_ind = (int)(phi / ROWS_DIV);

  // printf("theta index: %i, phi index: %i\n", theta_ind, phi_ind);
  double rand_num = rand() / ((double)RAND_MAX);

  double efficiency =
      bilinear_interp(eff_by_ang, ROWS, COLS, phi_ind, theta_ind);

  if (rand_num <= efficiency) {
    is_detected = 1;
  }

  return is_detected;
}

prim_lor *create_primitive_lor(event *event_pair) {

  // srand48(4);

  // checks to see if in-patient scatter occurs for the gamma pair
  double energy_1 = event_pair[0].energy;
  double energy_2 = event_pair[1].energy;
  energy_1 = (1000.0 * energy_1 + 0.0015);
  energy_2 = (1000.0 * energy_2 + 0.0015);
  int energy_1_int = (int)energy_1;
  int energy_2_int = (int)energy_2;
  if (energy_1_int != 511 || energy_2_int != 511) {
    gamma_pairs_detected_atleast_one_ips++;
  }
  if (energy_1_int != 511 && energy_2_int != 511) {
    gamma_pairs_detected_both_ips++;
  }

  vec3d z_hat = three_vec(0.0, 0.0, 1.0);

  prim_lor *new_prim_lor = (prim_lor *)malloc(sizeof(prim_lor));

  for (int i = 0; i < 2; i++) {

    // loc_i is vector from origin to event location on the detector
    event gamma_i = event_pair[i];
    vec3d loc_i = gamma_i.location;

    // circ_i_hat is the normalized vector that is tangent to a circular cross
    // section of the detector btw, this is where the assumption that the
    // detector is pointed along and is cocentric with the z-axis
    vec3d circ_i = vec_cross(loc_i, z_hat);
    vec3d circ_i_hat = vec_norm(circ_i);

    // creates three separate random numbers that are gaussian distributed about
    // zero
    double i_z_rand = gaussian(SPC_UNC, UNCERT_REP);
    double i_circ_rand = gaussian(SPC_UNC, UNCERT_REP);
    double i_time_rand = gaussian(TIME_UNC, UNCERT_REP);
    // creates vectors in the two distinct directions tangential to the detector
    // with magnitudes a random fraction of the spatial uncertainty
    vec3d i_z_rand_vec = vec_scaler(z_hat, i_z_rand);
    // printf("z randomness: ");
    // vec_print(i_z_rand_vec, stdout);
    // printf("\n");
    vec3d i_circ_rand_vec = vec_scaler(circ_i_hat, i_circ_rand);
    // any uncertainty applied in the r_hat direction?

    // creates the new (blurred) location of the gamma interaction, and inserts
    // it into new_prim_lor
    vec3d new_i_loc = vec_add(vec_add(loc_i, i_z_rand_vec), i_circ_rand_vec);
    if (i == 0) {
      new_prim_lor->first_loc = new_i_loc;
      new_prim_lor->first_energy = gamma_i.energy;
      new_prim_lor->first_time = gamma_i.tof + i_time_rand;
    } else {
      new_prim_lor->second_loc = new_i_loc;
      new_prim_lor->second_energy = gamma_i.energy;
      new_prim_lor->second_time = gamma_i.tof + i_time_rand;
    }
  }

  return new_prim_lor;
}

void test_create_primitive_lor(event *event_pair, prim_lor *primitive_lor) {

  event *gamma_1 = &event_pair[0];
  event *gamma_2 = &event_pair[1];
  vec3d loc_1 = gamma_1->location;
  vec3d loc_2 = gamma_2->location;
  vec3d loc_1_new = primitive_lor->first_loc;
  vec3d loc_2_new = primitive_lor->second_loc;

  printf("first gamma: \n");
  printf("(x,y,z) = (%lf, %lf, %lf) cm; TOF = %lf ns\n", loc_1.x, loc_1.y,
         loc_1.z, gamma_1->tof);
  printf("(x',y',z') = (%lf, %lf, %lf) cm; TOF' = %lf ns\n\n", loc_1_new.x,
         loc_1_new.y, loc_1_new.z, primitive_lor->first_time);

  printf("second gamma: \n");
  printf("(x,y,z) = (%lf, %lf, %lf) cm; TOF = %lf ns\n", loc_2.x, loc_2.y,
         loc_2.z, gamma_2->tof);
  printf("(x',y',z') = (%lf, %lf, %lf) cm; TOF' = %lf ns\n\n", loc_2_new.x,
         loc_2_new.y, loc_2_new.z, primitive_lor->second_time);

  return;
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
  vec3d geometric_center = vec_add(b, vec_scaler(c, 0.5));
  // printf("geometic center: \n");
  // vec_print(geometric_center, stdout);
  // printf("\n");
  vec3d c_hat = vec_norm(c);
  double delta_t = -(primitive_lor->first_time - primitive_lor->second_time);
  vec3d displacement_from_center = vec_scaler(c_hat, SPD_LGHT * delta_t * 0.5);
  vec3d annihilation_loc = vec_add(geometric_center, displacement_from_center);

  double transverse_uncert = sqrt(SPC_UNC * SPC_UNC + SPC_UNC * SPC_UNC);
  double longtidudinal_uncert =
      sqrt((SPD_LGHT * TIME_UNC) * (SPD_LGHT * TIME_UNC) +
           (SPD_LGHT * TIME_UNC) * (SPD_LGHT * TIME_UNC));

  lor *new = (lor *)malloc(sizeof(lor));
  new->center = annihilation_loc;
  new->dir = c_hat;
  new->long_uncert = longtidudinal_uncert;
  new->transverse_uncert = transverse_uncert;

  return new;
}

void test_create_lor(lor *new_lor) {

  vec3d center = new_lor->center;
  printf("LOR center (x,y,z) = (%lf, %lf, %lf) cm\n", center.x, center.y,
         center.z);

  return;
}

void print_lor(FILE *output, lor *lor) {
  fprintf(output, "%i, ", gamma_pairs_detected);
  fprintf(output, "%f, %f, %f,", lor->center.x, lor->center.y, lor->center.z);
  fprintf(output, " %f, %f, %f,", lor->dir.x, lor->dir.y, lor->dir.z);
  fprintf(output, " %f, %f\n", lor->long_uncert, lor->transverse_uncert);
}

int main(int argc, char **argv) {
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  // defines the help function and how to call it (by using -h or --help)
  for (int i = 0; i < num_flags(argc, argv); i++) {
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./hgm_lor_creator [TOPAS_file_location.phsp] "
             "[efficiency_table_location.csv] [LOR_output_location]\n");
      printf("-h: print this help\n");
      printf("-b: read phsp in binary format\n");
      exit(0);
    }
    if (strcmp(flags[i], "-b") == 0) {
      binary = 1;
    }
  }
  // checks to make sure you have correct number of args
  if (num_args(argc, argv) != 3) {
    printf("Incorrect number of arguments, three arguments required.\n");
    printf("Use the -h or --help command to get options.\n\n");
    exit(1);
  }

  // reads in efficiency table into 2D array called eff_by_ang
  printf("HGM LOR Creator\n\nLoading in '%s' as efficiencies table...\n",
         args[1]);
  FILE *eff_table_file = fopen(args[1], "r");
  int eff_file_read = read_eff(eff_table_file);
  fclose(eff_table_file);
  printf("Done!\n\n");

  // opens up a .lor file to output each LOR into
  char *lor_file_loc = (char *)malloc(sizeof(char) * (strlen(args[2]) + 10));
  strcpy(lor_file_loc, args[2]);
  lor_file_loc = strcat(lor_file_loc, ".lor");
  FILE *lor_output = fopen(lor_file_loc, "w");
  if (lor_output == NULL) {
    printf("Unable to open output file for writing\n");
    return 1;
  }

  FILE *phsp_file = fopen(args[0], "r");

  printf("Constructing the LORs...\n");

  event *pair_of_events = get_event_pair(phsp_file, read_event);

  while (pair_of_events != NULL) {

    // prints a status update every 1,000,000 gammas we read through
    // (effectively the line number the code is on)
    if ((gammas_entered_scanner / 1000000) * 1000000 ==
        gammas_entered_scanner) {
      printf("Number of recorded gammas read: %u\n", gammas_entered_scanner);
    }

    // works out if the pair of gammas that enter the detector are actually
    // detected
    int pair_detected = 0;
    for (int i = 0; i < 2; i++) {
      event *single_event = &pair_of_events[i];
      int detected = is_gamma_detected(single_event);
      // printf("Was this gamma detected? %u\n", detected);
      if (!detected) {

        if (DIAGNOSTIC_TOTAL_GAMMAS_DETECTED && i == 0) {
          int second_detected = is_gamma_detected(&pair_of_events[1]);
          if (second_detected) {
            diagnostic_single_gammas_detected++;
          }
        }

        break;
      }
      pair_detected++;
    }
    pair_detected = pair_detected / 2; // interesting question: is 1/2 = 0 or 1?

    // if the pair was detected, begins to construct the LOR
    if (pair_detected) {
      // increments metadata about gamma pairs that were detected
      // printf("\nReadin %i:\nGamma 1 location: ", pair_of_events[0].event_id);
      // vec_print(pair_of_events[0].location, stdout);
      // printf("\nGamma 2 location: ");
      // vec_print(pair_of_events[1].location, stdout);
      // printf("\n");
      gamma_pairs_detected++;

      prim_lor *primitive_lor = create_primitive_lor(pair_of_events);
      // printf("\nCreated primitive lor:\nGamma 1 location: ");
      // vec_print(primitive_lor->first_loc, stdout);
      // printf("\nGamma 2 location: ");
      // vec_print(primitive_lor->second_loc, stdout);
      // printf("\n");
      lor *new_lor = create_lor(primitive_lor);

      print_lor(lor_output, new_lor);

      free(primitive_lor);
      free(new_lor);
    }

    // frees memory pointed to by pair_of_events and proceeds on to get next
    // pair of events
    free(pair_of_events);
    pair_of_events = get_event_pair(phsp_file, read_event);
  }

  printf("Done!\n\n");

  // prints metadata
  printf("Gammas that entered the detector: %u\nGamma pairs that entered the "
         "detector: %u\nGamma pairs that were detected by the detector: %u\n",
         gammas_entered_scanner, gamma_pairs_entered_scanner,
         gamma_pairs_detected);
  uint one_ips_scatter =
      gamma_pairs_detected_atleast_one_ips - gamma_pairs_detected_both_ips;
  double percent_one_ips =
      100.0 * (double)one_ips_scatter / ((double)(gamma_pairs_detected));
  printf("Detected gamma pairs with only one gamma doing in-patient "
         "scattering: %i (%.2f%%)\n",
         one_ips_scatter, percent_one_ips);
  double percent_both_ips = 100.0 * (double)gamma_pairs_detected_both_ips /
                            ((double)(gamma_pairs_detected));
  printf("Detected gamma pairs with both gammas doing in-patient scattering: "
         "%i (%.2f%%)\n",
         gamma_pairs_detected_both_ips, percent_both_ips);
  double percent_atleastone_ips = 100.0 *
                                  (double)gamma_pairs_detected_atleast_one_ips /
                                  ((double)(gamma_pairs_detected));
  printf("Total number of detected gamma pairs with an in-patient scatter: %i "
         "(%.2f%%)\n",
         gamma_pairs_detected_atleast_one_ips, percent_atleastone_ips);
  if (DIAGNOSTIC_TOTAL_GAMMAS_DETECTED) {
    uint total_gammas_detected =
        2 * gamma_pairs_detected + diagnostic_single_gammas_detected;
    printf("Total gammas detected (including pairs and single gammas "
           "detected): %u\n",
           total_gammas_detected);
    double total_eff = 100.0 * (double)total_gammas_detected /
                       ((double)gammas_entered_scanner);
    printf("Total efficiency: (%.2f%%)\n", total_eff);
  }

  return 0;
}
