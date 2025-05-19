#ifndef hgmt_structs_h
#define hgmt_structs_h

#include "vector_ops.h"
#include <stdio.h>

#define PI 3.141592653589
#define COLS 105
#define SPD_LGHT 29.9792458 // cm/ns
#define UNCERT_REP 30
#define SPC_UNC 0     // 0.1 // cm
#define RAD_UNC 0     // 0.5
#define TIME_UNC 0.05 // 0.1 // 0.042463 // ns, sigma (0.100 ns FWHM)
#define DETECTOR_THICKNESS 2.54
#define DETECTOR_SEGMENTATION 0
typedef unsigned int uint;

typedef struct event_ {
  double tof;
  double energy;
  vec3d position;
  vec3d direction;
  int detector_id;
  uint primary; // created by a primary scattering (1=gamma_1 and 2=gamma_2)
  uint number;  // 0=first scatter, 1 = second, etc, only valid if above nonzero
  bool detected;
} event;
typedef struct hit_ {
  vec3d position;
  double tof;
  event *source;
} hit;

typedef struct prim_lor_ {
  hit hit1;
  hit hit2;
} prim_lor;

typedef struct _lor {
  vec3d center;
  vec3d dir;
  double long_uncert;
  double transverse_uncert;
} lor;
typedef struct _photon_path {
  int num_events;
  event **events;
  int num_hits;
  hit **hits;
} photon_path;

typedef struct _annihilation {
  vec3d origin;
  vec3d center;
  double time;
  uint num_events;
  event *events;
  uint num_hits;
  hit *hits;
  photon_path photon1_path;
  photon_path photon2_path;
} annihilation;
#endif
