#ifndef hgmt_structs_h
#define hgmt_structs_h

#include "vector_ops.h"
#include <stdio.h>

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
  vec3d first_loc;
  vec3d second_loc;
  double first_time;
  double second_time;
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
  photon_path photon1_path;
  photon_path photon2_path;
} annihilation;
#endif
