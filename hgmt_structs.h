#ifndef hgmt_structs_h
#define hgmt_structs_h

#include "vector_ops.h"
#include <stdio.h>

#define PI 3.141592653589
#define COLS 105
#define SPD_LGHT 29.9792458 // cm/ns
#define UNCERT_REP 30
#define SPC_UNC 0.1
#define TIME_UNC 0.042463 // ns, sigma (0.100 ns FWHM)
#define DETECTOR_THICKNESS 2.54
#define DETECTOR_SEGMENTATION 1
typedef unsigned int uint;

typedef struct event_ {
  uint event_id;
  double energy_deposit;
  vec3d location;
  vec3d momentum;
  double tof;
  int particle_type;
  int track_id;
} event;
typedef struct hit_ {
  vec3d location;
  double tof;
} hit;

typedef struct prim_lor_ {
  hit *hit1;
  hit *hit2;
} prim_lor;

typedef struct _lor {
  vec3d center;
  vec3d dir;
  double long_uncert;
  double transverse_uncert;
} lor;
typedef struct _photon_path {
  hit *hits;
  int num_hits;
  int has_first;
} photon_path;
typedef struct _annihilation {
  photon_path *photon1_path;
  photon_path *photon2_path;
} annihilation;
#endif
